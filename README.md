# shuffle

Anonymise genomic VCF files via simulated meiotic recombination. Each
synthetic output individual is a diploid mosaic of genotype segments drawn
from different real donors. Local linkage disequilibrium is preserved within
segments; no original individual's genome appears in the output intact.

---

## Table of contents

1. [How it works](#how-it-works)
   - [Continuous mode (whole-chromosome)](#continuous-mode-whole-chromosome)
   - [Region-sampling mode (targeted panels)](#region-sampling-mode-targeted-panels)
   - [Minimum-donors constraint](#minimum-donors-constraint)
   - [Sex chromosome donor filtering](#sex-chromosome-donor-filtering)
   - [Chromosome name normalisation](#chromosome-name-normalisation)
2. [Installation](#installation)
3. [Genetic map](#genetic-map)
4. [How-to guides](#how-to-guides)
   - [Shuffle a whole-chromosome cohort](#shuffle-a-whole-chromosome-cohort)
   - [Shuffle a targeted gene panel](#shuffle-a-targeted-gene-panel)
   - [Guarantee mixing with --min-donors](#guarantee-mixing-with---min-donors)
   - [Shuffle sex chromosomes](#shuffle-sex-chromosomes)
   - [Process multiple chromosomes](#process-multiple-chromosomes)
   - [Validate output](#validate-output)
   - [Programmatic use](#programmatic-use)
5. [CLI reference](#cli-reference)
   - [shuffle](#shuffle)
   - [validate](#validate)
6. [Inputs and assumptions](#inputs-and-assumptions)
   - [Normalisation (required)](#normalisation-required)
   - [One VCF per sample](#one-vcf-per-sample)
   - [Phasing](#phasing)
   - [Site concordance across VCFs](#site-concordance-across-vcfs)
   - [Per-sample variant-only VCFs (panel and WES data)](#per-sample-variant-only-vcfs-panel-and-wes-data)
   - [File format and indexing](#file-format-and-indexing)
   - [Cohort size](#cohort-size)
   - [Ploidy](#ploidy)
7. [Module and API reference](#module-and-api-reference)
   - [config.py — ShufflerConfig](#configpy--shufflerconfig)
   - [io/genetic_map.py — GeneticMap](#iogenetic_mapy--geneticmap)
   - [io/vcf_reader.py — PerSampleVCFReader](#iovcf_readery--persamplevcfreader)
   - [io/sex_map.py — sex map utilities](#iosex_mapy--sex-map-utilities)
   - [io/vcf_writer.py — SyntheticVCFWriter](#iovcf_writery--syntheticvcfwriter)
   - [core/genotype_pool.py — GenotypePool](#coregenotype_poolpy--genotypepool)
   - [core/recombination.py](#corerecombinationpy)
   - [core/mosaic_builder.py](#coremosaic_builderpy)
   - [validate.py](#validatepy)
   - [cli.py](#clipy)
8. [Test suite](#test-suite)
   - [Running tests](#running-tests)
   - [Unit tests](#unit-tests)
   - [Tier 1 — recombination model fidelity](#tier-1--recombination-model-fidelity)
   - [Tier 2 — privacy and biological plausibility](#tier-2--privacy-and-biological-plausibility)
   - [Tier 3 — PCA and LD (require external tools)](#tier-3--pca-and-ld-require-external-tools)
   - [Patient end-to-end integration test](#patient-end-to-end-integration-test)
   - [Full metrics summary](#full-metrics-summary)
9. [Known limitations](#known-limitations)
10. [Interpretation guide](#interpretation-guide)
11. [Further details](#further-details)

---

## How it works

### Continuous mode (whole-chromosome)

For each synthetic output individual, v-shuffler simulates meiotic crossovers
along the chromosome using a **Poisson process**:

- **Rate:** λ = `total_cM / 100` (one expected crossover per Morgan per
  meiosis). Overridable via `ShufflerConfig.n_crossovers_lambda`.
- **Positions:** drawn uniformly in cM space, i.e.
  `uniform(start_cM, end_cM)`. The genetic map encodes the spatially varying
  recombination rate, so breakpoints land more densely in hotspot regions
  when converted back to physical coordinates.

The chromosome is divided into non-overlapping segments at the crossover
positions. Each segment is assigned a randomly chosen donor (consecutive
segments always come from different donors). Genotypes are then streamed in
memory-efficient chunks; each synthetic individual's dosage at every variant
is copied from the donor assigned to that segment.

Continuous mode is intended for **whole-chromosome or whole-genome** data,
where the captured fraction of any chromosome is large enough that the
Poisson process produces meaningful crossover counts. Enable it explicitly
with `--no-region-sampling`.

### Region-sampling mode (targeted panels)

> **This is the default mode.** Use `--no-region-sampling` to switch to
> continuous mode.

Targeted NGS panels (exomes, gene panels, custom captures) only cover a small
scattered fraction of each chromosome. The effective λ across the captured
regions is near zero, so the Poisson crossover process almost never fires.
Most synthetics end up being copies of a single donor — providing no
anonymisation.

Region-sampling solves this by treating each **captured genomic region** as
an independent mixing unit:

1. **First pass** — `iter_positions()` scans the first input VCF to collect
   all variant positions (one lightweight pass, no genotype loading).
2. **Region detection** — `detect_regions()` groups positions into contiguous
   clusters using a configurable gap threshold (default: 10 kb). Each cluster
   corresponds crudely to one captured gene.
3. **Plan generation** — `build_region_segment_plan()` assigns one randomly
   chosen donor to each region independently. A without-replacement phase
   ensures that the first `min_donors` regions each get a distinct donor;
   the remaining regions are sampled freely subject to the adjacency
   constraint (consecutive regions never share the same donor).
4. **Mosaic assembly** — variants are streamed in chunks and each synthetic
   individual's dosage is copied from whichever donor the segment plan assigns
   to that variant's cM position.

With 100 captured regions and 200 donors, each synthetic draws ~63 distinct
donors — thorough mixing regardless of chromosome length. The primary-donor
attack success rate (see below) drops from ~89% (continuous mode, realistic λ ≈ 0.5) to
~5% (region mode, 100 regions).

### Minimum-donors constraint

Both modes accept `--min-donors N` (default: 1, meaning unconstrained).
When set, every synthetic individual is guaranteed to have at least N distinct
donors.

- **Region mode:** the first `effective_min = min(N, n_donors, n_regions)`
  regions are filled from a shrinking without-replacement pool; the rest are
  sampled freely.
- **Continuous mode:** extra breakpoints are injected until the plan contains
  at least N distinct donors. (The adjacency-only constraint can produce
  cycles like A→B→A, so a retry loop adds breakpoints one at a time until
  the count is satisfied.)

In both cases `effective_min` is capped at `min(N, n_donors, n_regions_or_segments)`
so the constraint always resolves without error.

### Sex chromosome donor filtering

Pass `--sex-file` to ensure the correct donor subset is used for sex
chromosomes:

- **chrX** → only female donors (diploid throughout, dosage 0/1/2 valid)
- **chrY** → only male donors
- **All other chromosomes** → full donor pool, sex file ignored

Without `--sex-file`, a warning is logged and the full pool is used; for
autosomes this is always correct. For chrX in a mixed-sex cohort it would
mix diploid female calls with haploid-or-homozygous male calls, producing
biologically implausible output.

The sex file is a plain two-column whitespace-delimited text file:

```
# path                           sex
/data/samples/sample001.vcf.gz   F
/data/samples/sample002.vcf.gz   M
sample003.vcf.gz                 female   # basename match works too
sample004.vcf.gz                 2        # PLINK convention (1=male, 2=female)
```

Lines starting with `#` and blank lines are ignored. An optional header row
(`path sex`) is skipped automatically. Accepted sex labels (case-insensitive):
`F` / `female` / `2` for female; `M` / `male` / `1` for male.

### Chromosome name normalisation

Different bioinformatics pipelines name chromosomes differently: some write
`chr22`, others write `22`. v-shuffler accepts either form for `--chromosome`
and automatically translates to whichever convention the input VCFs actually
use.

At the very start of a run, before the genetic map is loaded or any VCF
records are read, `resolve_chromosome_name()` opens the first input VCF and
inspects its sequence names (from `##contig` headers or the tabix index).
It returns the matching form — adding or stripping the `chr` prefix as
needed. If the VCF provides no sequence-name information (plain unindexed
files without `##contig` headers) the supplied name is used unchanged, and
the reader's existing fallback filter handles both forms automatically.

The resolved name flows through consistently to the genetic map lookup, tabix
queries, the writer's provenance header, and the output filename in
`multi_sample` mode. A log message is emitted whenever a translation occurs:

```
INFO  Chromosome name normalised from 'chr22' to '22' to match VCF convention.
```

---

## Installation

```bash
pip install -e .
# development (includes pytest, scipy, hypothesis, etc.)
pip install -e ".[dev]"
```

**Runtime dependencies:** `cyvcf2 ≥ 0.30`, `numpy ≥ 1.24`, `click ≥ 8.1`,
`pysam ≥ 0.21`, `scipy ≥ 1.10`, `tqdm ≥ 4.65`.

**External tools** (must be on `PATH` for compression and indexing):
`bgzip`, `tabix`. If absent, output VCFs are written uncompressed and
unindexed with a warning; the shuffle itself still completes.

---

## Genetic map

v-shuffler accepts two genetic map formats:

| Format | Columns | Header example |
|--------|---------|----------------|
| **SHAPEIT5 / Oxford** | `pos chr cM` | `pos chr cM` |
| **HapMap** | `Chromosome Position(bp) Rate(cM/Mb) Map(cM)` | `Chromosome Position(bp) Rate(cM/Mb) Map(cM)` |

Both plain text and `.gz` files are supported. The chromosome column may use
`chr22` or `22` spellings interchangeably.

**Recommended source — SHAPEIT5 GRCh38 maps:**

```
https://github.com/odelaneau/shapeit5/tree/main/maps/b38
```

Download the per-chromosome file (e.g. `chr22.b38.gmap.gz`) and pass it to
`--genetic-map`. In region-sampling mode the map is used only to convert bp
positions to cM for the segment plan; map accuracy is less critical than in
continuous mode.

---

## How-to guides

### Shuffle a whole-chromosome cohort

Use `--no-region-sampling` for whole-chromosome WGS or WES data where there
are enough variants per chromosome for the Poisson crossover model to produce
meaningful mixing.

```bash
v-shuffler shuffle \
    --input "data/per_sample/*.vcf.gz" \
    --output-dir shuffled/chr22/ \
    --genetic-map chr22.b38.gmap.gz \
    --chromosome chr22 \
    --n-samples 1000 \
    --seed 42 \
    --no-region-sampling
```

Or use a filelist instead of a glob:

```bash
ls data/per_sample/*.vcf.gz > samples.txt
v-shuffler shuffle \
    --input @samples.txt \
    --output-dir shuffled/chr22/ \
    --genetic-map chr22.b38.gmap.gz \
    --chromosome chr22 \
    --no-region-sampling
```

### Shuffle a targeted gene panel

Region-sampling is the default, so no extra flags are needed for panel data.
The tool auto-detects captured regions from the variant positions.

```bash
v-shuffler shuffle \
    --input "panel_data/*.vcf.gz" \
    --output-dir shuffled/ \
    --genetic-map chr22.b38.gmap.gz \
    --chromosome chr22 \
    --n-samples 500 \
    --seed 42
```

If your panel has unusually large intra-gene gaps you may need to increase
`--region-gap` beyond the default 10 kb:

```bash
v-shuffler shuffle \
    --input "panel_data/*.vcf.gz" \
    --output-dir shuffled/ \
    --genetic-map chr22.b38.gmap.gz \
    --chromosome chr22 \
    --region-gap 50000    # 50 kb threshold
```

Conversely, if variants within a single gene are not being grouped together,
decrease `--region-gap`.

### Guarantee mixing with --min-donors

Use `--min-donors` to ensure every synthetic draws from at least N distinct
individuals. This is useful when the donor pool is small or when you want a
hard lower bound on diversity regardless of the number of detected regions.

```bash
# Require at least 5 distinct donors per synthetic
v-shuffler shuffle \
    --input "data/*.vcf.gz" \
    --output-dir shuffled/ \
    --genetic-map chr22.b38.gmap.gz \
    --chromosome chr22 \
    --min-donors 5

# Combine with continuous mode for whole-chromosome data
v-shuffler shuffle \
    --input "data/*.vcf.gz" \
    --output-dir shuffled/ \
    --genetic-map chr22.b38.gmap.gz \
    --chromosome chr22 \
    --no-region-sampling \
    --min-donors 3
```

If `--min-donors` exceeds the number of available donors or detected regions
it is silently capped — the run always completes without error.

### Shuffle sex chromosomes

Use `--sex-file` to route chrX to female donors and chrY to male donors.
The simplest approach for most pipelines is to declare all synthetics as
female and shuffle chrX from female donors only, skipping chrY.

**Step 1 — Create a sex file.**
The file maps each donor VCF to its sex. If you have a PLINK `.fam` file
(columns: FID IID PAT MAT SEX PHENO, sex: 1=male, 2=female):

```bash
# Extract IID and sex, then match to VCF filenames
awk '{print $2".vcf.gz", $5}' cohort.fam > sex.txt
```

Or write it manually:

```
/data/donors/sample001.vcf.gz  F
/data/donors/sample002.vcf.gz  M
/data/donors/sample003.vcf.gz  F
```

**Step 2 — Run chrX with the sex file.**
Only female donors are used; male donors are silently excluded.

```bash
v-shuffler shuffle \
    --input "data/donors/*.vcf.gz" \
    --output-dir shuffled/chrX/ \
    --genetic-map chrX.b38.gmap.gz \
    --chromosome chrX \
    --n-samples 500 \
    --seed 42 \
    --sex-file sex.txt
```

**Step 3 — Skip chrY** (or shuffle it with male donors if needed):

```bash
# Male donors only for chrY
v-shuffler shuffle \
    --input "data/donors/*.vcf.gz" \
    --output-dir shuffled/chrY/ \
    --genetic-map chrY.b38.gmap.gz \
    --chromosome chrY \
    --n-samples 500 \
    --seed 42 \
    --sex-file sex.txt
```

**Note on donor pool size:** using only female donors for chrX halves (or
more) the pool. With region-sampling mode this is well-tolerated, but ensure
you have at least ~50 female donors for meaningful anonymisation.

### Process multiple chromosomes

Run v-shuffler once per chromosome. Use the same `--seed` for reproducibility
and `bcftools concat` to merge into per-synthetic genome-wide VCFs.

```bash
for chr in $(seq 1 22); do
    v-shuffler shuffle \
        --input @samples.txt \
        --output-dir "shuffled/chr${chr}/" \
        --genetic-map "maps/chr${chr}.b38.gmap.gz" \
        --chromosome "chr${chr}" \
        --n-samples 1000 \
        --seed 42
done

# Combine per-chromosome files into one file per synthetic individual
for i in $(seq 0 999); do
    bcftools concat shuffled/chr*/synthetic_${i}.vcf.gz \
        -Oz -o combined/synthetic_${i}.vcf.gz
    tabix -p vcf combined/synthetic_${i}.vcf.gz
done
```

For a multi-sample output (one VCF with all synthetics per chromosome):

```bash
v-shuffler shuffle \
    --input @samples.txt \
    --output-dir shuffled/chr22/ \
    --genetic-map chr22.b38.gmap.gz \
    --chromosome chr22 \
    --output-mode multi_sample
# produces shuffled/chr22/synthetic_chr22.vcf.gz
```

### Validate output

Run the built-in validator to check allele frequency preservation and confirm
no synthetic is a near-identical copy of any donor:

```bash
# Merge all donor VCFs into a single reference file first
bcftools merge -Oz -o merged_donors.vcf.gz data/per_sample/*.vcf.gz
tabix -p vcf merged_donors.vcf.gz

v-shuffler validate \
    --input "shuffled/synthetic_*.vcf.gz" \
    --reference-vcf merged_donors.vcf.gz \
    --chromosome chr22
```

**For panel or WES data with many rare variants:** the AF Pearson r ≥ 0.99
threshold applies to common variants only. When the site set is dominated by
rare private variants (MAF < 5%), the AF correlation across all sites will be
substantially lower than 0.99 even when the tool is working correctly — this
is a statistical property of estimating rare-variant frequencies from a small
synthetic cohort, not an algorithmic failure. Restrict the reference VCF to
MAF > 5% sites before validation:

```bash
# Filter reference to MAF >5% for AF validation
bcftools view -q 0.05:minor merged_donors.vcf.gz \
    -Oz -o merged_donors_common.vcf.gz
tabix -p vcf merged_donors_common.vcf.gz

v-shuffler validate \
    --input "shuffled/synthetic_*.vcf.gz" \
    --reference-vcf merged_donors_common.vcf.gz \
    --chromosome chr22
```

For a more thorough quantitative assessment, run the empirical test suite
(see [Test suite](#test-suite)).

### Programmatic use

The pipeline can be driven entirely from Python without the CLI:

```python
import numpy as np
from pathlib import Path
from v_shuffler.config import ShufflerConfig
from v_shuffler.cli import _run_shuffle

config = ShufflerConfig(
    input_vcfs=list(Path("data/").glob("*.vcf.gz")),
    output_dir=Path("shuffled/"),
    genetic_map=Path("chr22.b38.gmap.gz"),
    chromosome="chr22",
    n_output_samples=500,
    seed=42,
    region_sampling=True,           # default; False for whole-chromosome mode
    region_gap_bp=10_000,           # bp gap that starts a new region
    min_donors_per_synthetic=5,     # minimum distinct donors per synthetic
    sex_file=Path("sex.txt"),       # optional; filters donors for chrX/chrY
)
_run_shuffle(config)
```

For lower-level access, assemble the pipeline components directly:

```python
import numpy as np
from v_shuffler.io.genetic_map import GeneticMap
from v_shuffler.io.vcf_reader import PerSampleVCFReader
from v_shuffler.io.vcf_writer import SyntheticVCFWriter
from v_shuffler.core.recombination import (
    detect_regions,
    generate_all_region_plans,
)
from v_shuffler.core.mosaic_builder import build_synthetic_genotypes

gmap = GeneticMap("chr22.b38.gmap.gz", "chr22")
vcf_paths = list(Path("data/").glob("*.vcf.gz"))

reader = PerSampleVCFReader(vcf_paths, "chr22", gmap)

# First pass: detect captured regions
all_positions = reader.iter_positions()
regions_bp = detect_regions(all_positions, gap_threshold_bp=10_000)

# Convert to cM and generate plans
bp_flat = np.array([[r[0], r[1]] for r in regions_bp], dtype=np.int64).ravel()
cm_flat = gmap.bp_to_cm(bp_flat).reshape(-1, 2)
regions_cm = [(float(row[0]), float(row[1])) for row in cm_flat]

rng = np.random.default_rng(42)
plans = generate_all_region_plans(
    n_output_samples=100,
    regions_cm=regions_cm,
    n_pool_samples=len(vcf_paths),
    rng=rng,
    min_donors=5,
)

# Second pass: stream variants and write output
writer = SyntheticVCFWriter(
    output_dir=Path("shuffled/"),
    sample_names=[f"synthetic_{i}" for i in range(100)],
    template_vcf_path=vcf_paths[0],
    output_mode="per_sample",
    seed=42,
    chromosome="chr22",
    version="0.1.0",
)
for pool in reader.iter_chunks():
    writer.write_chunk(pool, build_synthetic_genotypes(pool, plans))
writer.finalize()
```

---

## CLI reference

### shuffle

```
v-shuffler shuffle [OPTIONS]
```

Shuffles genotypes between donor VCFs to produce anonymous synthetic
individuals.

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--input / -i` | required | Glob pattern or `@filelist.txt` of per-sample VCFs |
| `--output-dir / -o` | required | Output directory (created if absent) |
| `--genetic-map / -m` | required | SHAPEIT5 or HapMap format map file |
| `--chromosome / -c` | required | Chromosome to process. Both `chr22` and `22` forms are accepted; the tool normalises to match the VCF automatically. |
| `--n-samples / -n` | = n inputs | Number of synthetic individuals to produce |
| `--seed / -s` | None | Random seed; omit for a non-reproducible run |
| `--output-mode` | `per_sample` | `per_sample` (one VCF per synthetic) or `multi_sample` |
| `--chunk-size` | 50 000 | Variants loaded per memory chunk |
| `--max-missing` | 0.05 | Skip variants with > this fraction of missing calls |
| `--no-region-sampling` | off | Disable region mode; use classic Poisson crossover model |
| `--region-gap` | 10 000 | bp gap between variants that starts a new captured region |
| `--min-donors` | 1 | Minimum distinct donors per synthetic (1 = unconstrained) |
| `--sex-file` | None | Two-column donor-sex file; restricts chrX to female donors and chrY to male donors |
| `--threads` | 4 | Reserved for future parallelism; not yet used |
| `--verbose` | False | Enable DEBUG-level logging |

**Output (per_sample mode):** `{output-dir}/synthetic_0.vcf.gz`,
`synthetic_1.vcf.gz`, … plus corresponding `.tbi` index files.

**Output (multi_sample mode):** `{output-dir}/synthetic_{chromosome}.vcf.gz`
and `.tbi`.

---

### validate

```
v-shuffler validate [OPTIONS]
```

Sanity-checks synthetic output against the original donor VCFs.

```bash
v-shuffler validate \
    --input "shuffled/synthetic_*.vcf.gz" \
    --reference-vcf merged_donors.vcf.gz \
    --chromosome chr22
```

**Options:**

| Option | Description |
|--------|-------------|
| `--input / -i` | Glob or `@filelist.txt` of synthetic VCFs |
| `--reference-vcf / -r` | Merged multi-sample donor VCF |
| `--chromosome / -c` | Chromosome to validate |
| `--verbose` | Enable DEBUG logging |

**Checks performed:**

1. **Allele frequency correlation** — Pearson r between donor AF and
   synthetic AF across all variants (threshold: r ≥ 0.99).
2. **Identity leak** — pairwise concordance between every synthetic
   individual and every donor (threshold: no pair > 0.99).

For a more complete quantitative assessment use the empirical test suite
described in [Test suite](#test-suite).

---

## Inputs and assumptions

### Normalisation (required)

v-shuffler represents every genotype as a single integer dosage (0/1/2 for
hom-ref, het, hom-alt). This is unambiguous for biallelic sites, but
**multi-allelic records collapse to the wrong representation**: a donor with
`0/2` (one copy of ALT allele 2) and a donor with `0/1` (one copy of ALT
allele 1) both encode to dosage 1, and are written back as `0/1` — the
wrong allele for the `0/2` donor.

**Input VCFs must be normalised before running v-shuffler.** The required
preprocessing step is:

```bash
bcftools norm \
    -m -any \        # split all multiallelic records into separate biallelic records
    --keep-sum AD \  # sum allelic depth across split records rather than discarding
    input.vcf.gz \
    -Oz -o normalised.vcf.gz
tabix -p vcf normalised.vcf.gz
```

To also left-align indels (strongly recommended whenever a reference FASTA is
available, as it ensures equivalent indel representations are merged to
identical records):

```bash
bcftools norm \
    -m -any \
    --keep-sum AD \
    -f /path/to/reference.fa \   # left-aligns indels against the reference
    input.vcf.gz \
    -Oz -o normalised.vcf.gz
tabix -p vcf normalised.vcf.gz
```

For a cohort of per-sample VCFs:

```bash
for f in data/per_sample/*.vcf.gz; do
    base=$(basename "$f" .vcf.gz)
    bcftools norm -m -any --keep-sum AD -f ref.fa "$f" \
        -Oz -o "normalised/${base}.vcf.gz"
    tabix -p vcf "normalised/${base}.vcf.gz"
done
```

After normalisation every ALT column contains exactly one allele, dosage 1
unambiguously means one copy of that allele, and v-shuffler's dosage model is
fully correct.

### One VCF per sample

v-shuffler reads one file per donor, each containing exactly one sample
column. Multi-sample VCFs must be split first:

```bash
bcftools +split -Oz merged.vcf.gz -o per_sample/
```

Typically 100–2 000 donor files are used.

### Phasing

VCFs need **not** be phased. Unphased `0/1` and phased `0|1` calls are both
treated as dosage 1 — the phase separator is ignored on read. The shuffled
output is always written as unphased (`/` separator).

### Site concordance across VCFs

All per-sample VCFs must cover **exactly the same set of variants in the same
order** (identical CHROM/POS/REF/ALT tuples). Site consistency is validated
across files on the first chunk and a clear error is raised on mismatch.

If samples were genotyped at different site sets, intersect them before
running:

```bash
bcftools isec -p isec_out/ -Oz \
    normalised/sample0.vcf.gz \
    normalised/sample1.vcf.gz \
    ... \
    normalised/sampleN.vcf.gz
```

### Per-sample variant-only VCFs (panel and WES data)

Per-sample HaplotypeCaller (or equivalent) VCFs contain only the sites where
that individual had a variant call — homozygous-reference sites are absent.
This means different samples have different site sets, causing the site
consistency check to fail immediately.

The correct preprocessing for this data type is to apply a depth filter,
merge with reference-fill, re-normalise, and split back:

**Step 1 — Filter low-coverage variants.**
Calls with fewer than 20 reads are unreliable and should be removed before
any downstream use:

```bash
mkdir -p normalised_dp20/
for f in per_sample/*.vcf.gz; do
    base=$(basename "$f")
    bcftools view -i 'FORMAT/DP>=20' "$f" \
        -Oz -o "normalised_dp20/${base}"
    tabix -p vcf "normalised_dp20/${base}"
done
```

**Step 2 — Merge with reference-fill.**
`bcftools merge --missing-to-ref` fills absent sites as `0/0`, consolidating
all per-sample site sets into one concordant multi-sample VCF. This is valid
when all samples were generated from the **same capture panel** on the same
sequencing platform — the assumption that "absent variant = homozygous
reference" holds for well-covered panel sites. Sites where a sample had
genuinely insufficient coverage were removed in Step 1.

```bash
ls normalised_dp20/*.vcf.gz > sample_list.txt

bcftools merge \
    --missing-to-ref \
    --file-list sample_list.txt \
    -r 22 \                    # restrict to chromosome of interest
    -Oz -o merged_chr22.vcf.gz
tabix -p vcf merged_chr22.vcf.gz
```

**Step 3 — Re-normalise the merged VCF.**
Merging can create multiallelic records at positions where different samples
carry different ALT alleles. Split these back to biallelic before shuffling:

```bash
bcftools norm \
    -m -any \
    --keep-sum AD \
    merged_chr22.vcf.gz \
    -Oz -o merged_chr22_norm.vcf.gz
tabix -p vcf merged_chr22_norm.vcf.gz
```

**Step 4 — Split back to per-sample VCFs.**
v-shuffler requires one file per donor. The split files now share an
identical site set across all samples:

```bash
mkdir -p per_sample_merged/
bcftools +split -Oz merged_chr22_norm.vcf.gz -o per_sample_merged/
for f in per_sample_merged/*.vcf.gz; do tabix -p vcf "$f"; done
```

The resulting per-sample VCFs are suitable as direct input to v-shuffler.
See [AF correlation on rare-variant-dominated data](#af-correlation-on-rare-variant-dominated-data)
for an important limitation that applies when validating the output.

### File format and indexing

VCFs should be **bgzipped and tabix-indexed** (`.vcf.gz` + `.vcf.gz.tbi` or
`.csi`). Plain `.vcf` files work but require full-file iteration rather than a
tabix region query, which is substantially slower for large datasets.

Both `chr22` and `22` chromosome naming conventions are accepted; the tool
detects which form the VCFs use and normalises automatically.

### Cohort size

A minimum of **100 donors** is recommended. With fewer, one individual tends
to dominate each synthetic individual and re-identification risk is high. The
re-identification attack metrics (P2, P4) worsen non-linearly as pool size
falls below ~50.

### Ploidy

v-shuffler assumes **diploid autosomes** throughout. Dosage semantics (0/1/2)
are defined for diploid calls. For sex chromosomes, see
[Sex chromosome donor filtering](#sex-chromosome-donor-filtering).

---

## Module and API reference

### `config.py` — `ShufflerConfig`

All parameters for a single shuffle run. Used internally by the CLI and
available for programmatic use.

```python
from v_shuffler.config import ShufflerConfig
from pathlib import Path

config = ShufflerConfig(
    input_vcfs=[Path("sample0.vcf.gz"), Path("sample1.vcf.gz")],
    output_dir=Path("out/"),
    genetic_map=Path("chr22.b38.gmap.gz"),
    chromosome="chr22",
    n_output_samples=100,
    seed=42,
    # Recombination
    n_crossovers_lambda=None,       # None = total_cM / 100
    # Region sampling
    region_sampling=True,           # False → continuous mode
    region_gap_bp=10_000,           # bp gap that starts a new region
    min_donors_per_synthetic=1,     # 1 = no constraint
    # Sex chromosome filtering
    sex_file=None,                  # Path to donor-sex file; None = use all donors
    # Processing
    max_missing_rate=0.05,
    chunk_size_variants=50_000,
    output_mode="per_sample",
)
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_vcfs` | `list[Path]` | required | Per-sample donor VCF paths |
| `output_dir` | `Path` | required | Output directory (auto-created) |
| `genetic_map` | `Path` | required | Map file path |
| `chromosome` | `str` | required | Chromosome name |
| `n_output_samples` | `int` | required | Number of synthetics to produce |
| `seed` | `int \| None` | `None` | RNG seed |
| `n_crossovers_lambda` | `float \| None` | `None` | Override Poisson λ (continuous mode) |
| `min_segment_cM` | `float` | `0.5` | Minimum segment length (reserved) |
| `region_sampling` | `bool` | `True` | Enable region-based sampling |
| `region_gap_bp` | `int` | `10_000` | bp gap that separates two regions |
| `min_donors_per_synthetic` | `int` | `1` | Minimum distinct donors per synthetic |
| `sex_file` | `Path \| None` | `None` | Donor-sex file; restricts chrX/chrY to the appropriate sex |
| `max_missing_rate` | `float` | `0.05` | Per-variant missing-call filter |
| `chunk_size_variants` | `int` | `50_000` | Streaming chunk size |
| `n_threads` | `int` | `4` | Reserved for future use |
| `output_mode` | `str` | `"per_sample"` | `"per_sample"` or `"multi_sample"` |

---

### `io/genetic_map.py` — `GeneticMap`

Parses a SHAPEIT5 or HapMap genetic map file for one chromosome.

```python
from v_shuffler.io.genetic_map import GeneticMap

gmap = GeneticMap("chr22.b38.gmap.gz", "chr22")

gmap.total_length_cm   # float — total cM span
gmap.start_cm          # float — cM at first map position
gmap.end_cm            # float — cM at last map position

cm_positions = gmap.bp_to_cm(np.array([1_000_000, 5_000_000, 20_000_000]))
# → np.ndarray of float64, same shape as input
# Positions outside the map range are clamped to the boundary value.
```

`bp_to_cm` uses `numpy.interp`, which clamps at boundaries. Positions before
the first map entry receive `start_cm`; positions after the last entry
receive `end_cm`.

---

### `io/vcf_reader.py` — `PerSampleVCFReader` and `resolve_chromosome_name`

#### `resolve_chromosome_name`

```python
from v_shuffler.io.vcf_reader import resolve_chromosome_name

chrom = resolve_chromosome_name(Path("donor0.vcf.gz"), "chr22")
# → "chr22" if the VCF uses chr-prefixed names
# → "22"    if the VCF uses bare names
# → "chr22" unchanged if the VCF has no ##contig headers and no index
```

Inspects the sequence names reported by the first VCF's `##contig` headers
or tabix index — whichever is available. Tries the supplied name, then the
prefixed form, then the bare form; returns the first match. Falls back to the
supplied name if no sequence information is present. Called automatically at
the start of `_run_shuffle`; importable for programmatic pipelines.

#### `PerSampleVCFReader`

Streams one VCF file per donor simultaneously, yielding `GenotypePool`
chunks of shape `(chunk_size, n_donors)`.

```python
from v_shuffler.io.vcf_reader import PerSampleVCFReader

reader = PerSampleVCFReader(
    vcf_paths=[Path("donor0.vcf.gz"), Path("donor1.vcf.gz"), ...],
    chromosome="chr22",
    genetic_map=gmap,
    chunk_size=50_000,
    max_missing_rate=0.05,
)

# Lightweight first pass: collect positions for region detection
positions = reader.iter_positions()  # np.ndarray, int64

# Main pass: stream genotype chunks
for pool in reader.iter_chunks():
    # pool.dosages  : np.ndarray (V, N), uint8  — 0/1/2 or MISSING=255
    # pool.positions: np.ndarray (V,),   int64  — 1-based bp positions
    # pool.cm_pos   : np.ndarray (V,),   float64
    # pool.variant_info: list[VariantInfo], length V
    process(pool)
```

**`iter_positions()`** opens only the first VCF file, applies no
missing-rate filter, and returns a plain position array. It is used by the
CLI in region-sampling mode to detect captured regions before loading any
genotype data.

**`iter_chunks()`** notes:
- All VCF files must have identical variants in the same order. Site
  consistency (CHROM/POS/REF/ALT) is validated on the first chunk.
- Multi-allelic sites are encoded as total alt-allele count (dosage sum).
- Variants where the missing-call fraction exceeds `max_missing_rate` are
  silently skipped.
- Indexed (`.tbi`/`.csi`) files use a tabix region query; plain VCFs fall
  back to full-file iteration with a chromosome filter.

---

### `io/sex_map.py` — sex map utilities

Parses donor-sex files and filters VCF paths for sex chromosome runs.

```python
from v_shuffler.io.sex_map import (
    load_sex_map,
    filter_vcfs_by_sex,
    sex_filter_for_chromosome,
)
from pathlib import Path

vcf_paths = list(Path("data/").glob("*.vcf.gz"))

# Load the sex file → {Path: 'F'|'M'}
sex_map = load_sex_map(Path("sex.txt"), vcf_paths)

# Keep only female donors for chrX
female_vcfs = filter_vcfs_by_sex(vcf_paths, sex_map, keep_sex="F")

# Determine which sex to keep for a given chromosome
sex_filter_for_chromosome("chrX")   # → 'F'
sex_filter_for_chromosome("chrY")   # → 'M'
sex_filter_for_chromosome("chr22")  # → None  (use all donors)
```

**`parse_sex_label(raw)`** — normalises any accepted label to `'F'` or
`'M'`. Accepted values (case-insensitive): `F` / `female` / `2` for female;
`M` / `male` / `1` for male. Raises `ValueError` for anything else.

**`load_sex_map(sex_file, vcf_paths)`** — parses the two-column file. Tries
exact path match first, then basename fallback. Warns for any VCF paths not
found in the file. Gracefully skips comment lines (`#`), blank lines, and a
first-line header row.

**`filter_vcfs_by_sex(vcf_paths, sex_map, keep_sex)`** — returns only the
paths whose sex in `sex_map` equals `keep_sex`, preserving the original order.

**`sex_filter_for_chromosome(chromosome)`** — returns `'F'`, `'M'`, or
`None` based on the chromosome name. Handles `chrX`/`X`/`CHRX` variants for
X and `chrY`/`Y`/`CHRY` variants for Y.

---

### `io/vcf_writer.py` — `SyntheticVCFWriter`

Writes synthetic genotype data to one or more VCF files.

```python
from v_shuffler.io.vcf_writer import SyntheticVCFWriter

writer = SyntheticVCFWriter(
    output_dir=Path("out/"),
    sample_names=["synthetic_0", "synthetic_1"],
    template_vcf_path=Path("donor0.vcf.gz"),
    output_mode="per_sample",
    seed=42,
    chromosome="chr22",
    version="0.1.0",
)

for pool in reader.iter_chunks():
    dosages = build_synthetic_genotypes(pool, plans)  # (V, S) uint8
    writer.write_chunk(pool, dosages)

final_paths = writer.finalize()  # bgzip + tabix; returns list of .vcf.gz paths
```

**Dosage encoding:** `0 → "0/0"`, `1 → "0/1"`, `2 → "1/1"`, `255 → "./."`.
All output genotypes are **unphased** (`/` separator).

**Header handling:** copies all header lines from the template VCF, strips
`##sample=` metadata, replaces sample column names, and prepends a
`##v-shuffler=` provenance line (tool version, chromosome, seed, UTC
timestamp).

---

### `core/genotype_pool.py` — `GenotypePool`

Container for one chunk of genotype data across all donors.

```python
from v_shuffler.core.genotype_pool import GenotypePool, VariantInfo, MISSING

# MISSING = 255  (uint8 sentinel for missing genotype calls)

pool = GenotypePool(
    dosages=np.array(..., dtype=np.uint8),  # shape (V, N)
    positions=np.array(..., dtype=np.int64),
    cm_pos=np.array(..., dtype=np.float64),
    variant_info=[VariantInfo(chrom, pos, ref, alts, id, qual, filters, cm_pos), ...],
)

pool.n_variants  # int
pool.n_samples   # int
```

`VariantInfo` stores per-site metadata: `chrom`, `pos` (1-based), `ref`,
`alts` (list), `id`, `qual`, `filters`, `cm_pos`.

---

### `core/recombination.py`

Implements both the Poisson crossover model (continuous mode) and the
region-based sampling model (panel mode). All plan functions return lists of
`Segment(cm_start, cm_end, sample_idx)` frozen dataclasses.

#### `Segment`

```python
from v_shuffler.core.recombination import Segment

# Segment(cm_start: float, cm_end: float, sample_idx: int)
# Frozen dataclass. sample_idx is the column index in the GenotypePool matrix.
```

#### Continuous-mode functions

**`simulate_crossover_breakpoints`** — draws Poisson(λ) crossover positions
uniformly in cM space:

```python
from v_shuffler.core.recombination import simulate_crossover_breakpoints

rng = np.random.default_rng(42)
breakpoints = simulate_crossover_breakpoints(gmap, rng, lambda_override=None)
# → np.ndarray, float64, shape (k,) — sorted cM positions; may be empty
```

**`build_segment_plan`** — converts breakpoints into a segment plan:

```python
from v_shuffler.core.recombination import build_segment_plan

plan = build_segment_plan(breakpoints, gmap, n_samples=200, rng=rng)
# → list[Segment]  (covers [gmap.start_cm, gmap.end_cm] completely)
# Consecutive segments always come from different donors (n_samples > 1).
```

**`generate_all_segment_plans`** — generates plans for all synthetics
upfront:

```python
from v_shuffler.core.recombination import generate_all_segment_plans

plans = generate_all_segment_plans(
    n_output_samples=1000,
    genetic_map=gmap,
    n_pool_samples=200,
    rng=rng,
    lambda_override=None,  # None = total_cM / 100
    min_donors=1,          # 1 = unconstrained
)
# → list[list[Segment]], length n_output_samples
```

When `min_donors > 1`, extra breakpoints are injected until each plan
contains at least `min(min_donors, n_pool_samples)` distinct donors.

#### Region-sampling functions

**`detect_regions`** — groups sorted variant positions into captured regions:

```python
from v_shuffler.core.recombination import detect_regions

positions = np.array([1000, 1100, 1200, 600_000, 600_100], dtype=np.int64)
regions = detect_regions(positions, gap_threshold_bp=10_000)
# → [(1000, 1200), (600000, 600100)]
```

A new region starts whenever `positions[i] - positions[i-1] > gap_threshold_bp`.
Returns `[]` for empty input; `[(pos, pos)]` for a single position.

**`build_region_segment_plan`** — assigns one donor per region:

```python
from v_shuffler.core.recombination import build_region_segment_plan

regions_cm = [(0.0, 0.004), (1.0, 1.004), (2.0, 2.004)]
plan = build_region_segment_plan(
    regions_cm,
    n_samples=200,
    rng=rng,
    min_donors=3,   # first 3 regions get distinct donors
)
# → list[Segment], one per region
```

The first `effective_min = min(min_donors, n_samples, n_regions)` regions
are assigned donors without replacement (adjacency constraint relaxed only
when necessary); remaining regions use free sampling with adjacency
constraint only.

**`generate_all_region_plans`** — generates region-based plans for all
synthetics:

```python
from v_shuffler.core.recombination import generate_all_region_plans

plans = generate_all_region_plans(
    n_output_samples=500,
    regions_cm=regions_cm,
    n_pool_samples=200,
    rng=rng,
    min_donors=1,
)
# → list[list[Segment]], length n_output_samples
```

---

### `core/mosaic_builder.py`

Applies segment plans to a pool chunk to produce synthetic dosages.

**`apply_segment_plan`** — copies dosages from each segment's assigned donor:

```python
from v_shuffler.core.mosaic_builder import apply_segment_plan

dosages = apply_segment_plan(pool, plan)
# → np.ndarray, shape (pool.n_variants,), uint8
# Variants not covered by any segment receive MISSING = 255.
```

Matching condition: `cm_pos >= seg.cm_start AND cm_pos <= seg.cm_end`.
Single-variant regions (`cm_start == cm_end`) match correctly by equality.

**`build_synthetic_genotypes`** — vectorised over all plans:

```python
from v_shuffler.core.mosaic_builder import build_synthetic_genotypes

synth = build_synthetic_genotypes(pool, plans)
# → np.ndarray, shape (pool.n_variants, n_output_samples), uint8
```

---

### `validate.py`

Programmatic validation utilities used by the `validate` CLI subcommand.

```python
from v_shuffler.validate import run_validate
from pathlib import Path

run_validate(
    synth_paths=[Path("synthetic_0.vcf.gz"), ...],
    reference_vcf=Path("merged_donors.vcf.gz"),
    chromosome="chr22",
)
```

**Internal helpers** (importable but considered private):

| Function | Returns |
|----------|---------|
| `_read_reference(vcf_path, chrom)` | `(ref_afs, positions, sample_dosages)` — reads a multi-sample VCF in one pass |
| `_read_dosages(vcf_path, chrom, target_positions)` | `(positions, dosages)` — reads a single-sample VCF at specified positions |

---

### `cli.py`

Entry point for the `v-shuffler` command. Exposes two Click subcommands:
`shuffle` and `validate`.

**`_resolve_inputs(input_spec: str) → list[Path]`** — resolves the `--input`
argument. Accepts a glob pattern (`"data/*.vcf.gz"`), a `@filelist.txt`
reference, or a single file path. Raises `click.BadParameter` if no files
match or any listed path does not exist. This function is importable and
reused in the test suite.

**`_run_shuffle(config: ShufflerConfig) → None`** — the main pipeline
function. First resolves the chromosome name against the VCF convention via
`resolve_chromosome_name`, then optionally sex-filters the donor pool, then
constructs the reader (needed for the `iter_positions()` call in region mode),
generates plans, and streams variant chunks through to the writer. Importable
for programmatic use.

---

## Test suite

### Running tests

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run everything (unit + Tier 1 + Tier 2; Tier 3 auto-skipped without tools)
pytest -v

# Run only a specific tier
pytest tests/test_empirical_tier1.py -v
pytest tests/test_empirical_tier2.py -v

# Run with output printed (useful for inspecting informational metrics)
pytest tests/test_empirical_tier2.py -v -s
```

**Current status:** 145 tests pass, 6 skip (2 require plink2/bcftools; 4
require patient VCF env vars). Run time ≈ 18 s.

---

### Unit tests

| File | What it tests |
|------|---------------|
| `tests/test_recombination.py` | `simulate_crossover_breakpoints`, `build_segment_plan`, `generate_all_segment_plans`, `detect_regions`, `build_region_segment_plan`, `generate_all_region_plans` — edge cases, determinism, adjacency constraint, min_donors |
| `tests/test_genetic_map.py` | `GeneticMap` — SHAPEIT5 and HapMap format parsing, `bp_to_cm` interpolation, boundary clamping, error handling |
| `tests/test_mosaic_builder.py` | `apply_segment_plan`, `build_synthetic_genotypes` — correct segment assignment, MISSING fill, chunk boundaries |
| `tests/test_vcf_io.py` | `PerSampleVCFReader`, `SyntheticVCFWriter`, `resolve_chromosome_name` — round-trip VCF read/write, missing-rate filter, site consistency check, chr-prefix normalisation (four unit tests + two integration tests for cross-convention reading) |
| `tests/test_cli.py` | Full CLI via `click.testing.CliRunner` — shuffle run, determinism, multi-sample mode, unphased output, error handling; `TestChromosomeNormalisation` (bare VCF + chr-prefix flag, bare flag, normalisation log) |
| `tests/test_sex_filter.py` | `parse_sex_label`, `sex_filter_for_chromosome`, `load_sex_map`, `filter_vcfs_by_sex` — all label variants, path/basename matching, header/comment handling, error cases; CLI integration for autosome passthrough, chrX female-only filtering, no-sex-file warning, empty-pool error |

---

### Tier 1 — recombination model fidelity

**File:** `tests/test_empirical_tier1.py`
**Requirements:** none (numpy + scipy only)
**Run time:** ≈ 30 s (10 000 Monte Carlo repetitions for R1–R3; 1 000 for R5–R6)

Verifies the statistical properties of both the Poisson crossover model and
the region-sampling model.

#### R1–R3: Poisson crossover model

| Test ID | Test name | What it checks | Pass threshold |
|---------|-----------|----------------|---------------|
| R1 | `test_r1_mean` | `\|empirical mean − λ\| / λ` | < 1% |
| R1 | `test_r1_dispersion` | Poisson dispersion (var/mean) | 0.95 – 1.05 |
| R1 | `test_r1_ks_vs_poisson` | Discrete KS vs Poisson(λ) CDF | < 0.02 |
| R2 | `test_r2_cM_uniformity` | KS vs Uniform(0,1) in cM space | < 0.02 |
| R2 | `test_r2_mean_position` | Mean normalised cM position | 0.49 – 0.51 |
| R3 | `test_r3_mean_segment_length` | `\|mean length − L/(λ+1)\| / expected` | < 2% |

**R2 note:** breakpoints are drawn uniformly in *cM* space, not bp space. The
test confirms cM-space uniformity. Physical positions are deliberately
non-uniform (hotspot clustering) — this is correct behaviour captured by the
genetic map.

#### R4–R6: Region-sampling model

| Test ID | Test name | What it checks | Pass threshold |
|---------|-----------|----------------|---------------|
| R4 | `test_r4_region_detection_count` | `detect_regions` returns correct count | exact |
| R4 | `test_r4_region_boundaries` | Each region spans correct bp range | exact |
| R5 | `test_r5_cross_region_independence` | Mean cross-region co-assignment ≈ 1/N_DONORS | within 2% |
| R6 | `test_r6_min_donors_always_satisfied` | 100% of plans meet `min_donors=10` | 0 failures |
| R6 | `test_r6_min_donors_capped_by_n_regions` | Plans with `min_donors=200` capped at `n_regions` | 100% exact |

---

### Tier 2 — privacy and biological plausibility

**File:** `tests/test_empirical_tier2.py`
**Requirements:** none (fully in-memory synthetic data)
**Run time:** ≈ 15–20 s

Uses an in-memory panel-like fixture: 200 donors, 2 000 variants arranged in
100 gene regions (20 variants each, 100 bp intra-region spacing, 500 kb
inter-region gaps), 100 synthetics, 20 held-out individuals. Region-sampling
mode is the default; plans are generated via `generate_all_region_plans`.

#### Privacy tests

| Test ID | Test name | What it checks | Pass threshold |
|---------|-----------|----------------|---------------|
| P1 | `test_p1_max_concordance` | Max pairwise concordance (synth × donor) | < 0.99 |
| P1 | `test_p1_99th_percentile_concordance` | 99th pct of per-synth max concordance | < 0.85 |
| P1 | `test_p1_mean_concordance_vs_baseline` | Mean synth–donor vs donor–donor baseline | within 0.01 |
| P2 | `test_p2_closest_donor_attack` | Concordance-ranking attack success rate | < 100% (hard); warns if > 50% |
| P4 | `test_p4_membership_inference` | Wilcoxon p-value (in-pool vs held-out) | finite (hard); warns if p < 0.05 |

**P2 and P4 are diagnostic, not pass/fail.** With region-sampling and 100
regions, the primary-donor attack rate is ~5% (compared to ~89% in continuous
mode at λ ≈ 0.5). Membership inference remains detectable (a known limitation
of the unphased mosaic design) but the signal is weaker. See
[Known limitations](#known-limitations).

#### Biological plausibility tests

| Test ID | Test name | What it checks | Pass threshold |
|---------|-----------|----------------|---------------|
| B1 | `test_b1_af_global` | Global AF Pearson r | ≥ 0.98 (panel fixture; ≥ 0.99 in production) |
| B1 | `test_b1_af_by_maf_bin` | AF correlation by MAF bin | printed only |
| B2 | `test_b2_heterozygosity` | Per-sample het rate: mean diff and KS p-value | `\|mean diff\| < 0.005`; p > 0.01 |
| B3 | `test_b3_hwe` | HWE-fail fraction (p < 0.001) synth vs donor | synth ≤ 3× donor |
| B4 | `test_b4_tstv_identical` | ts/tv ratio and variant count | identical by construction |

**B1 note:** the panel-like fixture has only 100 effectively independent
region assignments. Within-region AF correlation reduces the effective N,
so the threshold is 0.98 for the test fixture rather than the production
value of 0.99 which applies to whole-chromosome data.

#### Region-mode specific tests

| Test name | What it checks |
|-----------|----------------|
| `test_region_detection_in_fixture` | `detect_regions` correctly identifies all 100 simulated gene clusters |
| `test_min_donors_region_mode` | `min_donors=10` guarantees ≥ 10 distinct donors per plan |
| `test_min_donors_continuous_mode` | `min_donors=5, lambda=0` guarantees ≥ 5 distinct donors via breakpoint injection |
| `test_p1_region_vs_continuous_concordance` | Region mode max concordance < continuous mode max concordance (λ=0.5) |

---

### Tier 3 — PCA and LD (require external tools)

**File:** `tests/test_empirical_tier3.py`
**Requirements:** `plink2` and `bcftools` on `PATH`, plus pre-generated VCF
files pointed to by env vars
**Auto-skipped** when tools or files are absent

These tests are **mode-agnostic**: both region-mode and continuous-mode
output VCFs can be tested with the same commands.

Set environment variables before running:

```bash
export VSHUFFLE_DONOR_VCF=/path/to/merged_donors.vcf.gz
export VSHUFFLE_SYNTH_VCF=/path/to/synthetic_chr22.vcf.gz
pytest tests/test_empirical_tier3.py -v
```

For whole-chromosome data, consider producing the synthetic VCF with
`--no-region-sampling` so breakpoints fall at biologically realistic
positions. For panel data, use the default region mode.

| Test ID | Test name | What it checks | Pass threshold |
|---------|-----------|----------------|---------------|
| B5 | `test_b5_pca` | Fraction of synthetics inside 99% Mahalanobis donor cloud (PC1–PC2) | > 0.95; 0 outliers > 3 SD |
| B6 | `test_b6_ld_decay` | Short-range r² diff at < 10 kb; LD curve Pearson r; KS on full r² distribution | < 0.02; > 0.99; < 0.05 |

**B6 note:** long-range LD (distances > mean segment length) is expected to
be elevated in synthetic output because LD is copied intact within segments.
This is a known artefact, not a bug.

---

### Patient end-to-end integration test

**File:** `tests/test_patient_end_to_end.py`
**Requirements:** real per-sample VCF files, a genetic map, `bgzip`/`tabix`
**Auto-skipped** without env vars

Runs the full shuffle pipeline on real patient VCFs via `CliRunner` and
checks the four most critical output properties. Intended to be run manually
before releasing shuffled data.

**Required environment variables:**

| Variable | Description |
|----------|-------------|
| `VSHUFFLE_PATIENT_VCFS` | Glob or `@filelist` of per-sample donor VCFs |
| `VSHUFFLE_GENETIC_MAP` | Path to genetic map (SHAPEIT5 or HapMap) |
| `VSHUFFLE_CHROMOSOME` | Chromosome name, e.g. `chr22` |
| `VSHUFFLE_SEED` | *(optional)* Random seed, default `42` |
| `VSHUFFLE_N_SYNTH` | *(optional)* Synthetics to produce, default `min(n_donors, 50)` |

```bash
VSHUFFLE_PATIENT_VCFS="@/data/patients.txt" \
VSHUFFLE_GENETIC_MAP="/data/chr22.b38.gmap.gz" \
VSHUFFLE_CHROMOSOME="chr22" \
pytest tests/test_patient_end_to_end.py -v -s
```

| Test name | What it checks |
|-----------|----------------|
| `test_shuffle_produced_output` | Correct number of `.vcf.gz` files written |
| `test_no_identity_leak` | Max pairwise concordance (synth × donor) < 0.99 |
| `test_af_preserved` | AF Pearson r ≥ 0.99 |
| `test_variant_count_consistent` | All synth variant positions present in donor pool |

**Required input preparation for per-sample variant-only VCFs.**
Standard per-sample HaplotypeCaller output contains only called variant sites,
so different samples have different site sets and v-shuffler's site consistency
check will fail. The required preprocessing steps are:

1. Filter out low-coverage calls (`FORMAT/DP >= 20`).
2. Merge with `bcftools merge --missing-to-ref` (absent = homozygous reference,
   valid when all samples use the same capture panel).
3. Re-normalise with `bcftools norm -m -any --keep-sum AD` to split any
   multiallelics introduced by the merge.
4. Split back to per-sample with `bcftools +split`.

See [Per-sample variant-only VCFs](#per-sample-variant-only-vcfs-panel-and-wes-data)
for the complete command sequence.

**AF correlation threshold for panel data.**
The `test_af_preserved` threshold of r ≥ 0.99 applies to common variant data.
After `--missing-to-ref` merging, a clinical panel dataset is typically
dominated by rare and private variants (MAF < 1%). With a small synthetic
cohort (e.g., 20 samples), estimating allele frequency at rare sites has high
sampling variance — a site with a single het carrier in 138 donors (AF ≈ 0.4%)
will either appear or not in any given synthetic purely by chance, creating
scatter that holds the global Pearson r well below 0.99 even when the
algorithm is working correctly.

**For panel data, restrict AF validation to MAF > 5%** — these are the sites
where the algorithm's AF-preservation property can actually be measured. In
practice this means passing a MAF-filtered reference VCF to `VSHUFFLE_REFERENCE_VCF`
or filtering the merged donor VCF before running the test:

```bash
bcftools view -q 0.05:minor merged_donors.vcf.gz -Oz -o merged_common.vcf.gz
```

On validated 138-sample chr22 panel data with DP ≥ 20 filtering and 20
synthetic individuals, the empirical results were:

| Test | Result | Notes |
|------|--------|-------|
| `test_shuffle_produced_output` | ✅ Pass | 20 synthetic VCFs written |
| `test_no_identity_leak` | ✅ Pass | Max concordance < 0.99 (anonymisation effective) |
| `test_af_preserved` | ⚠ Fails at r = 0.94 on all 28 976 sites | Expected; passes when restricted to MAF > 5% sites |
| `test_variant_count_consistent` | ✅ Pass | All 28 976 positions consistent |

---

### Full metrics summary

| ID | Test | Metric | Threshold | Tier |
|----|------|--------|-----------|------|
| R1 | Crossover count | KS vs Poisson(λ) | < 0.02 | 1 |
| R2 | Breakpoint uniformity | KS vs Uniform in cM space | < 0.02 | 1 |
| R3 | Segment length | `\|mean − L/(λ+1)\| / expected` | < 2% | 1 |
| R4 | Region detection | exact count and boundaries | exact | 1 |
| R5 | Cross-region independence | mean co-assignment vs 1/N | within 2% | 1 |
| R6 | min_donors reliability | fraction failing constraint | 0% | 1 |
| P1 | Pairwise concordance | max concordance | < 0.99 | 2 |
| P1 | Pairwise concordance | 99th pct of per-synth max | < 0.85 | 2 |
| P2 | Closest-donor attack | attack rate / random baseline | diagnostic | 2 |
| P4 | Membership inference | Wilcoxon p (in vs out of pool) | diagnostic | 2 |
| B1 | Allele frequency | global Pearson r | ≥ 0.99 (production) | 2 |
| B2 | Heterozygosity | `\|mean per-sample het diff\|` | < 0.005 | 2 |
| B3 | Hardy-Weinberg | synth HWE-fail / donor HWE-fail | ≤ 3× | 2 |
| B4 | ts/tv ratio | `\|ts/tv synth − donor\|` | 0.00 (exact) | 2 |
| B5 | PCA | fraction inside 99% Mahalanobis cloud | > 0.95 | 3 |
| B6 | LD decay | curve Pearson r (binned r²) | > 0.99 | 3 |

---

## Known limitations

### Concordance-based re-identification (P2)

In **continuous mode** with realistic chr22 parameters (λ ≈ 0.55), the
primary donor is detectable ~89–100% of the time via simple concordance
ranking: P(0 crossovers) ≈ 58%, meaning most synthetics are copies of a
single donor.

**Region-sampling mode** substantially reduces this risk. With 100 captured
regions each assigned independently, the primary donor contributes only ~1–3
regions out of 100. The attack success rate drops to ~5% in the Tier 2
fixture (ratio ~10× over random baseline vs ~3 000× in continuous mode).

**Remaining remedies:** increase the donor pool size; increase `--min-donors`;
use phased haplotypes rather than unphased diploid genotypes.

### Membership inference (P4)

Synthetic individuals are measurably more concordant with their in-pool donors
than with held-out individuals from the same population. Region-sampling
weakens but does not eliminate this signal (mean delta ~0.015 in the Tier 2
fixture at 100 regions). The Wilcoxon p-value remains significant because
the in-pool donors collectively contributed all donor alleles, and their
alleles are therefore over-represented in the synthetics. Larger donor pools
and more independent regions both reduce the effect.

### Short chromosomes and sparse panels

Both P2 and P4 worsen when the number of independently sampled mixing units
is small. In continuous mode, short chromosomes have low λ and many zero-
crossover synthetics. In region mode, panels with few captured regions give
each donor a higher average contribution. Use `--min-donors` to enforce a
floor on mixing.

### Unphased diploid model

Segment swapping operates on diploid dosages (0/1/2), not individual
haplotypes. This correctly handles unphased input but means the synthetic
genome is not equivalent to a phased recombinant. Short-range LD is preserved
within segments; medium- and long-range LD across segment boundaries is
broken.

### Sex chromosomes

The tool does not validate donor ploidy against the chromosome being
processed. In a **mixed-sex cohort on chrX**, male donors have hemizygous
calls (dosage 0 or 2 only in the non-PAR region) and female donors are
fully diploid (0/1/2). Without intervention, a synthetic might receive
female genotypes in some regions and hemizygous genotypes in others —
biologically implausible and liable to confuse downstream tools.

**Use `--sex-file` to avoid this.** When provided, the tool automatically
restricts the donor pool to female donors for chrX and male donors for chrY.
Autosomes use the full pool regardless. See the
[Shuffle sex chromosomes](#shuffle-sex-chromosomes) how-to for a complete
workflow.

Without `--sex-file` on a sex chromosome, a warning is logged but the run
proceeds — this is intentional, so that same-sex cohorts (where filtering
is unnecessary) are not penalised.

### AF correlation on rare-variant-dominated data

The AF Pearson r ≥ 0.99 threshold is designed for datasets with a broadly
uniform allele-frequency distribution (e.g., common variants from
whole-chromosome or whole-genome data). It breaks down for clinical panel or
WES data after `bcftools merge --missing-to-ref` because:

- **Most sites are rare or private.** After merging 138 per-sample
  variant-only VCFs, approximately 50% of sites have MAF < 1% — these are
  variants carried by only one or two patients. The mean AF across all sites
  is ~9%.

- **Rare-variant AF estimates from a small synthetic cohort have high
  sampling variance.** A site carried by a single het donor (AF = 0.4%) will
  either appear in a given synthetic (AF = 2.5% of a 20-synthetic cohort) or
  not (AF = 0%) depending purely on which donor was assigned to that region.
  This binary scatter across thousands of rare sites suppresses the Pearson r
  far below 0.99 even when the algorithm is working correctly.

- **AF preservation holds in expectation.** The expected synthetic AF at any
  site equals the donor AF — this is guaranteed by the uniform random donor
  assignment in region-sampling. The departure from r = 0.99 is sampling
  noise, not systematic bias.

**Remedy:** restrict AF validation to MAF > 5% sites, where the algorithm's
AF-preservation property can be measured with reasonable statistical power.
The identity-leak test (max concordance < 0.99) is the primary anonymisation
check and is not affected by this limitation.

### Multi-allelic sites

v-shuffler encodes genotypes as a total alt-allele dosage (0/1/2), which
collapses alt-allele identity: `0/1` (one copy of ALT1) and `0/2` (one copy
of ALT2) both become dosage 1 and are written back as `0/1`. **This is
resolved by the normalisation requirement** — splitting all multiallelic
records to biallelic with `bcftools norm -m -any --keep-sum AD` before
running v-shuffler ensures there is exactly one ALT allele per record, making
the dosage representation unambiguous. See
[Normalisation (required)](#normalisation-required).

---

## Interpretation guide

| Observation | Likely cause | Remedy |
|-------------|-------------|--------|
| P2 attack rate > 50% | Too few mixing units; one donor dominates | Switch to region mode; increase `--min-donors`; use larger pool |
| P4 membership inference p < 0.001 | Strong in-pool concordance signal | Same as P2; also increase variant panel density |
| B1 AF r < 0.99 on whole-chromosome data | Donor pool imbalance or population stratification | Verify all donors are from the same population; check `--chromosome` is the right chromosome |
| B1 AF r < 0.99 on panel / WES data | Site set dominated by rare and private variants (MAF < 5%) | Expected behaviour; re-run validation restricted to MAF > 5% sites (`bcftools view -q 0.05:minor`) |
| Site mismatch error on startup | Per-sample variant-only VCFs have different site sets | Apply DP filter + `bcftools merge --missing-to-ref` + re-normalise + `bcftools +split`; see [Per-sample variant-only VCFs](#per-sample-variant-only-vcfs-panel-and-wes-data) |
| Output VCFs are empty (0 variants) | Chromosome name could not be resolved | Run with `--verbose` to see whether normalisation fired; check VCF `##contig` headers |
| B2 het rate shift | Systematic genotype encoding error | Check that dosage encoding is 0/1/2, not 0/1/1 |
| Unexpected allele identities in output (e.g. wrong ALT allele at het sites) | Multi-allelic records in input not decomposed | Run `bcftools norm -m -any --keep-sum AD` on all donor VCFs before shuffling |
| B3 HWE excess > 3× | Many variants near segment boundaries | Check R3 segment length distribution; consider longer segments |
| B5 synthetics outside PCA cloud | AF bias or population stratification in donor pool | Check all donors come from the same population |
| B6 elevated long-range LD | Expected artefact of unphased diploid-mosaic design | Not a bug; document as known limitation |
| Too few regions detected | Intra-gene gaps exceed `--region-gap` | Increase `--region-gap` |
| All variants in one region | `--region-gap` too large | Decrease `--region-gap`; or use `--no-region-sampling` |
| Warning: processing chrX without --sex-file | Mixed-sex cohort risk | Add `--sex-file` pointing to a donor-sex file; safe to ignore for same-sex cohorts |
| Error: no female donors found for chrX | Sex file has no F entries matching the input VCFs | Check file paths in sex file; use basename or full path consistently |

---

## Further details

For a deep-dive into the codebase — including annotated data-flow diagrams,
key design decisions, and worked examples — see the technical walkthrough:

[technical_walkthrough.md](technical_walkthrough.md)

For end-to-end testing results against a real 138-sample clinical NGS panel
cohort (chr22, GATK HaplotypeCaller, preprocessing pipeline, anonymisation
and AF metrics) see:

[testing_results.md](testing_results.md)
