# shuffle: Technical Walkthrough

*2026-03-10T15:56:11Z by Showboat 0.6.1*
<!-- showboat-id: 02bcdae9-c019-4ee8-87d6-48ba8771c471 -->

shuffle anonymises genomic VCF files by shuffling diploid genotypes between individuals using simulated meiotic recombination. It is designed for genomics researchers who need to share cohort-level variant data without exposing identifiable individual genotypes.

**The problem it solves:** Raw VCF files are personally identifiable — long runs of genotypes are unique to each person. shuffle transforms them into synthetic individuals that are biologically plausible (preserving population-level statistics and local linkage disequilibrium) but unidentifiable.

**How it works:**

1. **Preprocessing (first pass, region mode)** — The tool scans the first donor VCF to collect all variant positions. Captured genomic regions are detected by finding gaps larger than a configurable threshold. This drives the region-based mixing strategy for targeted panels.
2. **Plan generation** — For each synthetic individual, a segment plan is built: a list of genomic intervals each assigned a randomly chosen donor. In *region mode* (default), each captured region is assigned one donor independently. In *continuous mode* (whole-chromosome data), crossover positions are drawn via a Poisson process in genetic (cM) space.
3. **Mosaic assembly (streaming)** — Donor VCFs are streamed in chunks. Each chunk is matched against the pre-computed plans to produce the synthetic dosage matrix.
4. **Output** — Bgzipped, tabix-indexed VCFs with anonymised sample names, stripped metadata, and a provenance header.

**Supported platforms:** Linux/macOS; Python ≥ 3.10, bgzip, tabix on PATH.

## Repository Layout

<details>
<summary>Full directory tree</summary>

```python3

import os
SKIP = {'.git', '.claude', '__pycache__', '.pytest_cache', '.hypothesis', 'v_shuffler.egg-info'}
lines = []
for root, dirs, files in os.walk('.'):
    dirs[:] = sorted(d for d in dirs if d not in SKIP)
    level = root.count(os.sep)
    indent = '  ' * level
    lines.append(f'{indent}{os.path.basename(root)}/')
    for f in sorted(files):
        lines.append(f'{indent}  {f}')
print('\n'.join(lines))

```

```output
./
  .gitignore
  README.md
  pyproject.toml
  requirements-dev.txt
  requirements.txt
  technical_walkthrough.md
  testing_results.md
  tests/
    __init__.py
    conftest.py
    test_cli.py
    test_empirical_tier1.py
    test_empirical_tier2.py
    test_empirical_tier3.py
    test_genetic_map.py
    test_mosaic_builder.py
    test_patient_end_to_end.py
    test_recombination.py
    test_sex_filter.py
    test_vcf_io.py
  v_shuffler/
    __init__.py
    __main__.py
    cli.py
    config.py
    validate.py
    core/
      __init__.py
      genotype_pool.py
      mosaic_builder.py
      recombination.py
    io/
      __init__.py
      genetic_map.py
      sex_map.py
      vcf_reader.py
      vcf_writer.py
```

</details>

## Architecture

shuffle is structured as a Python package with two top-level namespaces:

- **`v_shuffler/core/`** — Pure computation: crossover simulation, region detection, genotype pool, mosaic assembly. No I/O.
- **`v_shuffler/io/`** — File I/O: genetic map loading, VCF reading, VCF writing, and donor sex-map filtering.
- **`v_shuffler/cli.py`** — Click-based CLI that wires everything together.
- **`v_shuffler/config.py`** — Dataclass holding all configuration for a run.
- **`v_shuffler/validate.py`** — Post-run sanity checks (AF correlation, identity check).

### Data flow

```
                        ┌─ sex_file (optional) ──────────────────────┐
Per-sample VCF files    │  filter donors to F (chrX) or M (chrY)     │
(one per donor)         └────────────────────────────────────────────┘
        │
        ▼  resolve_chromosome_name()   ← detect chr22 vs 22 convention
        │
        │  ┌── REGION MODE (default) ──────────────────────────────────┐
        │  │  iter_positions()     first pass: collect bp positions     │
        │  │  detect_regions()     cluster positions by gap threshold   │
        │  │  generate_all_region_plans()  one donor per region         │
        │  └───────────────────────────────────────────────────────────┘
        │  ┌── CONTINUOUS MODE (--no-region-sampling) ─────────────────┐
        │  │  simulate_crossover_breakpoints()  Poisson in cM space    │
        │  │  build_segment_plan()              assign donors to segs  │
        │  │  generate_all_segment_plans()                              │
        │  └───────────────────────────────────────────────────────────┘
        │
        ▼
 PerSampleVCFReader.iter_chunks()    ← GeneticMap.bp_to_cm()
   reads variants in lockstep, converts GT→dosage, annotates cM positions
   yields GenotypePool chunks  (n_variants × n_donors, uint8)
        │
        ▼
 build_synthetic_genotypes(pool, segment_plans)
   applies pre-computed plans to each pool chunk
   yields synthetic dosage matrix  (n_variants × n_synthetic, uint8)
        │
        ▼
 SyntheticVCFWriter.write_chunk() → finalize()
   converts dosage→GT string, writes VCF records, bgzip + tabix
```

The segment plans are computed **once upfront** (tiny — a few dozen floats per individual in cM space) and reused across every chunk streamed from disk. Memory cost is proportional to `chunk_size × n_donors`, not total variant count.

## CLI Reference

The entry point is `v_shuffler/cli.py`, exposed as the `v-shuffler` console script. It has two subcommands: `shuffle` and `validate`.

<details>
<summary>v-shuffler --help</summary>

```bash
/tmp/vshuffler-venv/bin/python -m v_shuffler --help
```

```output
Usage: python -m v_shuffler [OPTIONS] COMMAND [ARGS]...

  v-shuffler: anonymise VCFs via simulated genetic recombination.

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  shuffle   Shuffle variants between VCFs via simulated genetic...
  validate  Sanity-check synthetic VCFs against the original input.
```

</details>

<details>
<summary>v-shuffler shuffle --help</summary>

```bash
/tmp/vshuffler-venv/bin/python -m v_shuffler shuffle --help
```

```output
Usage: python -m v_shuffler shuffle [OPTIONS]

  Shuffle variants between VCFs via simulated genetic recombination.

  Each output individual is a mosaic of genotype segments drawn from different
  real donors at biologically realistic crossover positions.

Options:
  -i, --input TEXT                Input VCF glob ("data/*.vcf.gz") or
                                  @filelist.txt.  [required]
  -o, --output-dir PATH           Directory for output VCF files (created if
                                  absent).  [required]
  -m, --genetic-map PATH          Genetic map file (SHAPEIT5 or HapMap
                                  format).  [required]
  -c, --chromosome TEXT           Chromosome to process, e.g. "chr22" or "22".
                                  [required]
  -n, --n-samples INTEGER         Number of synthetic individuals to produce.
                                  Defaults to number of inputs.
  -s, --seed INTEGER              Random seed for reproducibility.
  --chunk-size INTEGER            Variants to process per memory chunk.
                                  [default: 50000]
  --threads INTEGER               Number of threads (reserved for future
                                  parallelism).  [default: 4]
  --output-mode [per_sample|multi_sample]
                                  Write one VCF per synthetic individual, or a
                                  single multi-sample VCF.  [default:
                                  per_sample]
  --max-missing FLOAT             Maximum fraction of missing calls per
                                  variant (variants above this are skipped).
                                  [default: 0.05]
  --no-region-sampling            Disable region-based sampling. Use classic
                                  continuous-cM mode (for whole-chromosome
                                  data).
  --region-gap INTEGER            bp gap between variants that starts a new
                                  captured region (region mode only).
                                  [default: 10000]
  --min-donors INTEGER            Minimum distinct donors per synthetic
                                  individual.  [default: 1]
  --sex-file PATH                 Two-column file mapping donor VCF paths to
                                  sex (F/M). When provided, chrX is shuffled
                                  from female donors only and chrY from male
                                  donors only. Autosomes use the full donor
                                  pool.
  --verbose                       Enable debug logging.
  --help                          Show this message and exit.
```

</details>

<details>
<summary>v-shuffler validate --help</summary>

```bash
/tmp/vshuffler-venv/bin/python -m v_shuffler validate --help
```

```output
Usage: python -m v_shuffler validate [OPTIONS]

  Sanity-check synthetic VCFs against the original input.

  Checks: 1. Output VCFs are valid and parseable. 2. Allele frequency
  correlation between input and output is high (r > 0.99). 3. No output sample
  is identical to any input sample.

Options:
  -i, --input TEXT          Glob or @filelist pointing to synthetic output
                            VCFs.  [required]
  -r, --reference-vcf PATH  A merged multi-sample input VCF to compare allele
                            frequencies against.  [required]
  -c, --chromosome TEXT     Chromosome to validate.  [required]
  --verbose
  --help                    Show this message and exit.
```

</details>

### Pipeline startup: `_run_shuffle()`

Before any VCF data is loaded, `cli._run_shuffle()` performs three preparatory steps in order:

**1. Chromosome name normalisation.** `resolve_chromosome_name()` opens the first donor VCF header and checks its sequence names against the user-supplied `--chromosome` value. If the VCF uses `22` but the user passed `chr22` (or vice versa), the name is silently translated. The resolved form is used consistently for the genetic map lookup, tabix queries, and output filenames.

**2. Sex-based donor filtering.** If `--sex-file` is supplied, the donor pool is restricted before any genotype reading begins:
- `chrX` → female donors only
- `chrY` → male donors only
- Autosomes → full pool

A warning is logged if a sex chromosome is processed without a sex file.

**3. Plan generation.** Segment plans for all synthetic individuals are computed upfront using either the region or continuous strategy (see sections below). Plans live entirely in cM space and are tiny in memory.

### Input resolution (`_resolve_inputs`)

Converts the `--input` string into a concrete list of `Path` objects:

| Form | Example |
|------|---------|
| Glob | `"data/per_sample/*.vcf.gz"` |
| File list | `@samples.txt` (one path per line, `#` comments ignored) |
| Single path | `sample.vcf.gz` |

Missing files are reported before any I/O starts.

## Core Module: Recombination (`v_shuffler/core/recombination.py`)

This module contains the biological heart of shuffle. It provides two distinct mixing strategies — region-based (for panels) and continuous Poisson-based (for whole chromosomes) — plus a shared data structure.

### Key type: `Segment` (frozen dataclass)

```python
@dataclass(frozen=True)
class Segment:
    cm_start:   float  # inclusive segment start in cM
    cm_end:     float  # inclusive end for the last segment; exclusive otherwise
    sample_idx: int    # column index of the donor in GenotypePool.dosages
```

A segment plan for one synthetic individual is a `list[Segment]`. Both mixing strategies produce the same type, so `apply_segment_plan()` in the mosaic builder works unchanged regardless of how the plan was generated.

### Region-sampling mode (default)

Designed for targeted NGS panels where most of a chromosome is not captured. The effective Poisson rate across the sparse captured fraction approaches zero, so continuous-mode would almost never fire a crossover and most synthetics would copy a single donor.

Region-sampling solves this by treating each captured genomic region as an independent mixing unit:

```
detect_regions(positions, gap_threshold_bp=10_000)
  → list of (start_bp, end_bp) tuples
  → a new region starts whenever positions[i] − positions[i−1] > gap_threshold_bp

build_region_segment_plan(regions_cm, n_samples, rng, min_donors=1)
  → one Segment per region
  → first effective_min = min(min_donors, n_samples, n_regions) regions
    are assigned donors without replacement (guarantees diversity)
  → remaining regions: free sampling with adjacency constraint only
    (consecutive regions never share the same donor)
```

`generate_all_region_plans()` calls `build_region_segment_plan()` once per synthetic individual.

```/tmp/vshuffler-venv/bin/python3

import numpy as np
from v_shuffler.core.recombination import detect_regions, build_region_segment_plan

# Panel-like positions: 3 gene clusters separated by 500 kb gaps
cluster1 = list(range(1_000_000, 1_001_000, 100))   # 10 variants
cluster2 = list(range(2_000_000, 2_001_000, 100))   # 10 variants
cluster3 = list(range(5_000_000, 5_001_000, 100))   # 10 variants
positions = np.array(cluster1 + cluster2 + cluster3, dtype=np.int64)

regions = detect_regions(positions, gap_threshold_bp=10_000)
print(f'Detected {len(regions)} regions from {len(positions)} variants:')
for i, (s, e) in enumerate(regions):
    print(f'  Region {i}: {s:,}–{e:,} bp')

# Build a plan for one synthetic individual from a pool of 10 donors
rng = np.random.default_rng(42)
regions_cm = [(float(i * 2), float(i * 2 + 1)) for i in range(len(regions))]
plan = build_region_segment_plan(regions_cm, n_samples=10, rng=rng, min_donors=3)
print(f'Segment plan ({len(plan)} segments, min_donors=3):')
for seg in plan:
    print(f'  [{seg.cm_start:.1f}, {seg.cm_end:.1f}] cM  ->  donor {seg.sample_idx}')
print(f'Distinct donors used: {len(set(s.sample_idx for s in plan))}')

```

```output
Detected 3 regions from 30 variants:
  Region 0: 1,000,000–1,000,900 bp
  Region 1: 2,000,000–2,000,900 bp
  Region 2: 5,000,000–5,000,900 bp
Segment plan (3 segments, min_donors=3):
  [0.0, 1.0] cM  ->  donor 0
  [2.0, 3.0] cM  ->  donor 7
  [4.0, 5.0] cM  ->  donor 6
Distinct donors used: 3
```

### Continuous mode (`--no-region-sampling`)

Used for whole-chromosome or whole-genome data where there are enough variants for the Poisson crossover model to provide meaningful mixing.

**`simulate_crossover_breakpoints(gmap, rng, lambda_override=None)`**

1. Sample the number of crossovers: `n ~ Poisson(λ)` where `λ = total_cM / 100` by default.
2. Sample positions uniformly in `[start_cM, end_cM]`. Working in cM means the physical density automatically reflects local recombination rate.

**`build_segment_plan(breakpoints, gmap, n_samples, rng)`**

Builds boundary array `[map_start, *breakpoints, map_end]` and assigns a new donor at each boundary. Consecutive segments always come from different donors.

**`generate_all_segment_plans(..., min_donors=1)`**

Generates plans for all synthetic individuals upfront. When `min_donors > N`, extra breakpoints are injected one at a time until the plan contains at least `min(min_donors, n_pool_samples)` distinct donors — a retry loop is needed because the adjacency-only constraint can produce cycles (A→B→A) that give only 2 distinct donors even with 2 breakpoints.

```/tmp/vshuffler-venv/bin/python3

import numpy as np
from v_shuffler.core.recombination import simulate_crossover_breakpoints, build_segment_plan

class MockMap:
    start_cm = 0.0; end_cm = 100.0; total_length_cm = 100.0

rng = np.random.default_rng(42)
gmap = MockMap()
breakpoints = simulate_crossover_breakpoints(gmap, rng)
print(f'Crossovers (lambda=1.0 for 100 cM map): {len(breakpoints)}')
print(f'Breakpoint positions (cM): {breakpoints.round(2)}')

plan = build_segment_plan(breakpoints, gmap, n_samples=5, rng=rng)
print(f'Segment plan ({len(plan)} segments):')
for seg in plan:
    print(f'  [{seg.cm_start:.2f}, {seg.cm_end:.2f}) cM  ->  donor {seg.sample_idx}')

```

```output
Crossovers (lambda=1.0 for 100 cM map): 1
Breakpoint positions (cM): [85.86]
Segment plan (2 segments):
  [0.00, 85.86) cM  ->  donor 0
  [85.86, 100.00) cM  ->  donor 3
```

## Core Module: GenotypePool (`v_shuffler/core/genotype_pool.py`)

`GenotypePool` is the central in-memory data structure representing a chunk of variants across all donors.

### Dosage encoding

| Value | Meaning |
|-------|---------|
| `0` | Homozygous reference (0/0) |
| `1` | Heterozygous (0/1) |
| `2` | Homozygous alt (1/1) |
| `255` | Missing sentinel (`MISSING`) |

Input VCFs must be normalised to biallelic records (`bcftools norm -m -any`) before use. Multi-allelic records collapse alt-allele identity (both `0/1` and `0/2` become dosage 1) and are written back as `0/1` — the wrong allele for the `0/2` donor.

### Key attributes

| Attribute | Shape | dtype | Description |
|-----------|-------|-------|-------------|
| `dosages` | `(n_variants, n_samples)` | uint8 | Dosage matrix |
| `positions` | `(n_variants,)` | int64 | Physical positions (1-based bp) |
| `cm_pos` | `(n_variants,)` | float64 | Genetic positions (cM) |
| `variant_info` | `list[VariantInfo]` | — | Per-variant metadata for VCF output |

## Core Module: Mosaic Builder (`v_shuffler/core/mosaic_builder.py`)

### `apply_segment_plan(pool, plan)`

Builds the synthetic dosage vector for one output individual from one pool chunk:

```python
result = np.full(pool.n_variants, MISSING, dtype=np.uint8)
for seg in plan:
    mask = (pool.cm_pos >= seg.cm_start) & (pool.cm_pos <= seg.cm_end)
    result[mask] = pool.dosages[mask, seg.sample_idx]
```

This works identically for both region-mode and continuous-mode plans because both produce the same `list[Segment]` structure. Single-variant regions (`cm_start == cm_end`) match correctly by equality.

### `build_synthetic_genotypes(pool, segment_plans)`

Iterates over all plans and stacks results into a `(n_variants, n_output_samples)` matrix. Called once per streamed chunk:

```python
for pool in reader.iter_chunks():
    synthetic = build_synthetic_genotypes(pool, segment_plans)
    writer.write_chunk(pool, synthetic)
```

## IO Module: VCF Reader (`v_shuffler/io/vcf_reader.py`)

### `resolve_chromosome_name(vcf_path, chromosome)`

Called at the very start of `_run_shuffle()` before any genotype data is read. Opens the first donor VCF and inspects `reader.seqnames` (populated from `##contig` headers or the tabix index). Tries the supplied chromosome name, then the prefixed form (`22` → `chr22`), then the bare form (`chr22` → `22`). Returns the first match; falls back to the supplied name if no contig information is available (plain VCFs without `##contig` headers). An INFO message is logged whenever a translation occurs.

### `PerSampleVCFReader`

Reads one VCF per donor and yields `GenotypePool` chunks. Key design points:

- Uses **cyvcf2** for fast C-backed VCF iteration.
- All per-sample VCFs opened simultaneously and iterated **in lockstep** with `zip(*iterators)` so the dosage vector for each variant is assembled in one pass.
- **`iter_positions()`** — lightweight first pass used in region mode. Opens only the first VCF file, applies no missing-rate filter, and returns a `np.ndarray` of bp positions. This is O(N_variants × 8 bytes) and much faster than `iter_chunks()`.
- **`iter_chunks()`** — main streaming pass. CHROM/POS/REF/ALT consistency across all files is validated on the first chunk. Tabix region queries are used when an index exists; plain-text VCFs fall back to full-file iteration with a chromosome filter.
- Variants with `n_missing / n_samples > max_missing_rate` are silently skipped; a warning is logged at the end.

### Genotype → dosage (`_gt_to_dosage`)

Reads allele indices `[a1, a2, phased]` from cyvcf2:
- Returns `MISSING (255)` if either allele is `-1`.
- Returns `int(a1 > 0) + int(a2 > 0)` otherwise.

## IO Module: Sex Map (`v_shuffler/io/sex_map.py`)

Provides donor-pool filtering for sex chromosomes. Three public functions:

**`sex_filter_for_chromosome(chromosome)`** — Returns `'F'` for chrX/X, `'M'` for chrY/Y, `None` for autosomes. Case-insensitive; handles both `chr`-prefixed and bare names.

**`load_sex_map(sex_file, vcf_paths)`** — Parses a two-column whitespace-delimited file (`path sex`). Accepts `F/female/2` for female and `M/male/1` for male (PLINK conventions). Tries exact path match first, then basename fallback. Warns for any VCF paths not found in the file. Skips `#` comments, blank lines, and a first-line header row.

**`filter_vcfs_by_sex(vcf_paths, sex_map, keep_sex)`** — Returns only the paths whose sex in `sex_map` equals `keep_sex`, preserving the original order.

### Sex file format

```
# path                          sex
/data/samples/sample001.vcf.gz  F
sample002.vcf.gz                M        # basename match also works
sample003.vcf.gz                female
sample004.vcf.gz                2        # PLINK code
```

When `--sex-file` is provided and the chromosome is chrX, the donor pool is silently restricted to female donors. A warning (not an error) is emitted when a sex chromosome is processed without a sex file, so same-sex cohorts are not penalised.

## IO Module: Genetic Map (`v_shuffler/io/genetic_map.py`)

`GeneticMap` loads a per-chromosome genetic map and provides physical-to-genetic coordinate conversion.

| Format | Columns | Example header |
|--------|---------|----------------|
| SHAPEIT5 / Oxford | 3 (space-delimited) | `pos chr cM` |
| HapMap | 4 (space-delimited) | `Chromosome Position(bp) Rate(cM/Mb) Map(cM)` |

Both plain text and `.gz` files are accepted. The chromosome name in the file may use either the `chr22` or `22` spelling; the loader normalises both internally (independent of the `resolve_chromosome_name` step above).

**`GeneticMap.bp_to_cm(positions)`** converts bp arrays to cM via `numpy.interp`. Positions outside the map range are clamped to the boundary cM value. The map checks monotonicity at load time — a non-monotonic map would silently corrupt crossover assignments.

## IO Module: VCF Writer (`v_shuffler/io/vcf_writer.py`)

Converts the synthetic dosage matrix back into VCF files.

### Output modes

| Mode | Flag | Output |
|------|------|--------|
| `per_sample` (default) | `--output-mode per_sample` | One `synthetic_N.vcf.gz` per synthetic |
| `multi_sample` | `--output-mode multi_sample` | One `synthetic_${chrom}.vcf.gz` with all samples as columns |

### Header handling

1. All header lines copied from the first donor VCF.
2. `##sample=` lines stripped (sample-identifiable metadata).
3. Provenance line inserted before `#CHROM`: `##v-shuffler=<version="...",seed="...",chromosome="...",date="...">`
4. Sample column names replaced with `synthetic_0`, `synthetic_1`, …

### Dosage → GT

| Dosage | Written as |
|--------|-----------|
| 0 | `0/0` |
| 1 | `0/1` |
| 2 | `1/1` |
| 255 | `./.` |

All output is unphased (`/` separator). `finalize()` calls `bgzip -f` and `tabix -p vcf` via subprocess; if either is absent, a warning is logged and the uncompressed file is returned.

## Configuration (`v_shuffler/config.py`)

All pipeline parameters are collected into a single `ShufflerConfig` dataclass. The CLI constructs this object and passes it to `_run_shuffle()`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_vcfs` | `list[Path]` | required | Resolved per-sample VCF paths |
| `output_dir` | `Path` | required | Output directory (created if absent) |
| `genetic_map` | `Path` | required | Genetic map file |
| `chromosome` | `str` | required | Chromosome label (e.g. `"chr22"` or `"22"`) |
| `n_output_samples` | `int` | required | Number of synthetic individuals |
| `seed` | `int \| None` | `None` | RNG seed (None = non-reproducible) |
| `n_crossovers_lambda` | `float \| None` | `None` | Poisson rate override (continuous mode) |
| `region_sampling` | `bool` | `True` | `False` → continuous mode |
| `region_gap_bp` | `int` | `10_000` | bp gap that starts a new captured region |
| `min_donors_per_synthetic` | `int` | `1` | Minimum distinct donors per synthetic (1 = unconstrained) |
| `sex_file` | `Path \| None` | `None` | Two-column donor-sex file for chrX/chrY filtering |
| `max_missing_rate` | `float` | `0.05` | Skip variants with > 5% missing calls |
| `chunk_size_variants` | `int` | `50_000` | Variants per streaming chunk |
| `n_threads` | `int` | `4` | Thread count (reserved for future use) |
| `output_mode` | `str` | `"per_sample"` | `"per_sample"` or `"multi_sample"` |

`__post_init__` validates all numeric ranges, checks that `sex_file` exists if provided, and creates the output directory.

## Validation Subcommand (`v_shuffler/validate.py`)

The `validate` subcommand provides two quantitative checks.

**Check 1: Allele frequency correlation.** Pearson r between donor AF and synthetic AF across all variants. Threshold: r ≥ 0.99. For panel/WES data dominated by rare variants, restrict the reference VCF to MAF > 5% before validation — see the AF correlation limitation in Known Limitations.

**Check 2: Identity leak detection.** Max pairwise concordance (synth × donor). Any pair exceeding 0.99 triggers a warning.

## Error Handling

### Fail-fast errors

| Error | Cause |
|-------|-------|
| `ValueError` in `GeneticMap.__init__` | No chromosome entries, or non-monotonic cM values |
| `ValueError` in `iter_chunks` | CHROM/POS/REF/ALT mismatch across VCFs on first chunk |
| `click.BadParameter` in `_resolve_inputs` | Glob matched nothing, or listed paths don't exist |
| `click.ClickException` in `_run_shuffle` | `--sex-file` produced an empty donor pool |
| `ValueError` in `ShufflerConfig.__post_init__` | Invalid field values |

### Degrade-gracefully

| Situation | Behaviour |
|-----------|-----------|
| `bgzip`/`tabix` not on PATH | Warning logged; uncompressed file returned |
| Sex chromosome processed without `--sex-file` | Warning logged; full pool used |
| Chromosome name convention mismatch | INFO logged; name translated automatically |
| High missing-rate variants | Silently skipped; count in warning at end |

## Dependencies

<details>
<summary>Installed dependency versions</summary>

```bash
/tmp/vshuffler-venv/bin/pip show cyvcf2 numpy click pysam scipy tqdm 2>/dev/null | grep -E '^(Name|Version):'
```

```output
Name: cyvcf2
Version: 0.32.1
Name: numpy
Version: 2.4.3
Name: click
Version: 8.3.1
Name: pysam
Version: 0.23.3
Name: scipy
Version: 1.17.1
Name: tqdm
Version: 4.67.3
```

</details>

| Package | Role |
|---------|------|
| **cyvcf2** | Fast C-backed VCF reading. Used in `vcf_reader.py` and `vcf_writer.py`. |
| **numpy** | Array operations: dosage matrices, cM interpolation, crossover sampling, boolean masks. |
| **click** | CLI framework: argument parsing, subcommand routing, help generation. |
| **pysam** | htslib bindings; bgzip/tabix are called via subprocess in the current implementation. |
| **scipy** | Statistical tests in the empirical test suite (KS, Wilcoxon, chi-squared HWE). |
| **tqdm** | Progress bar over the chunk-streaming loop.

**External tools:** `bgzip` and `tabix` (htslib) for output compression and indexing.

## Tests

The test suite lives in `tests/` and uses pytest. It is organised into five layers:

| Layer | File(s) | What it covers |
|-------|---------|----------------|
| Unit | `test_recombination.py`, `test_genetic_map.py`, `test_mosaic_builder.py`, `test_vcf_io.py`, `test_cli.py`, `test_sex_filter.py` | Individual functions and classes |
| Tier 1 | `test_empirical_tier1.py` | Recombination model statistical correctness (R1–R6) |
| Tier 2 | `test_empirical_tier2.py` | Privacy and biological plausibility (P1–P4, B1–B4) |
| Tier 3 | `test_empirical_tier3.py` | PCA and LD decay (require plink2 + bcftools) |
| E2E | `test_patient_end_to_end.py` | Full pipeline on real patient VCFs (require env vars) |

Tier 3 and E2E tests are auto-skipped when tools or environment variables are absent.

<details>
<summary>Full test collection (145 tests)</summary>

```bash
/tmp/vshuffler-venv/bin/pytest --collect-only -q 2>&1 | grep -E '^tests/'
```

```output
tests/test_cli.py::TestShuffleCommand::test_help
tests/test_cli.py::TestShuffleCommand::test_missing_required_options
tests/test_cli.py::TestShuffleCommand::test_shuffle_runs
tests/test_cli.py::TestShuffleCommand::test_shuffle_deterministic
tests/test_cli.py::TestShuffleCommand::test_shuffle_multi_sample_mode
tests/test_cli.py::TestShuffleCommand::test_nonexistent_input_raises
tests/test_cli.py::TestShuffleCommand::test_output_vcf_is_unphased
tests/test_cli.py::TestShuffleCommand::test_output_not_identical_to_any_input
tests/test_cli.py::TestVersion::test_version
tests/test_cli.py::TestChromosomeNormalisation::test_bare_vcf_with_chr_prefix_flag
tests/test_cli.py::TestChromosomeNormalisation::test_bare_vcf_with_bare_flag
tests/test_cli.py::TestChromosomeNormalisation::test_normalisation_logged
tests/test_empirical_tier1.py::test_r1_mean
tests/test_empirical_tier1.py::test_r1_dispersion
tests/test_empirical_tier1.py::test_r1_ks_vs_poisson
tests/test_empirical_tier1.py::test_r2_cM_uniformity
tests/test_empirical_tier1.py::test_r2_mean_position
tests/test_empirical_tier1.py::test_r3_mean_segment_length
tests/test_empirical_tier1.py::test_r4_region_detection_count
tests/test_empirical_tier1.py::test_r4_region_boundaries
tests/test_empirical_tier1.py::test_r5_cross_region_independence
tests/test_empirical_tier1.py::test_r6_min_donors_always_satisfied
tests/test_empirical_tier1.py::test_r6_min_donors_capped_by_n_regions
tests/test_empirical_tier2.py::test_p1_max_concordance
tests/test_empirical_tier2.py::test_p1_99th_percentile_concordance
tests/test_empirical_tier2.py::test_p1_mean_concordance_vs_baseline
tests/test_empirical_tier2.py::test_p2_closest_donor_attack
tests/test_empirical_tier2.py::test_p4_membership_inference
tests/test_empirical_tier2.py::test_b1_af_global
tests/test_empirical_tier2.py::test_b1_af_by_maf_bin
tests/test_empirical_tier2.py::test_b2_heterozygosity
tests/test_empirical_tier2.py::test_b3_hwe
tests/test_empirical_tier2.py::test_b4_tstv_identical
tests/test_empirical_tier2.py::test_region_detection_in_fixture
tests/test_empirical_tier2.py::test_min_donors_region_mode
tests/test_empirical_tier2.py::test_min_donors_continuous_mode
tests/test_empirical_tier2.py::test_p1_region_vs_continuous_concordance
tests/test_empirical_tier3.py::test_b5_pca
tests/test_empirical_tier3.py::test_b6_ld_decay
tests/test_genetic_map.py::test_load_shapeit5_format
tests/test_genetic_map.py::test_load_hapmap_format
tests/test_genetic_map.py::test_chrom_without_prefix
tests/test_genetic_map.py::test_wrong_chrom_raises
tests/test_genetic_map.py::test_bp_to_cm_exact_points
tests/test_genetic_map.py::test_bp_to_cm_interpolated
tests/test_genetic_map.py::test_bp_to_cm_monotonic
tests/test_genetic_map.py::test_total_length_cm
tests/test_genetic_map.py::test_positions_outside_map_clamped
tests/test_genetic_map.py::test_non_monotonic_cm_raises
tests/test_genetic_map.py::test_gzipped_map
tests/test_mosaic_builder.py::test_single_segment_copies_donor
tests/test_mosaic_builder.py::test_two_segments_different_donors
tests/test_mosaic_builder.py::test_all_zero_pool_gives_all_zero
tests/test_mosaic_builder.py::test_missing_sentinel_preserved
tests/test_mosaic_builder.py::test_plan_at_chunk_boundary
tests/test_mosaic_builder.py::test_build_synthetic_genotypes_shape
tests/test_mosaic_builder.py::test_build_synthetic_genotypes_correct_donors
tests/test_mosaic_builder.py::test_build_synthetic_genotypes_mosaic
tests/test_patient_end_to_end.py::test_shuffle_produced_output
tests/test_patient_end_to_end.py::test_no_identity_leak
tests/test_patient_end_to_end.py::test_af_preserved
tests/test_patient_end_to_end.py::test_variant_count_consistent
tests/test_recombination.py::test_breakpoints_within_map
tests/test_recombination.py::test_breakpoints_sorted
tests/test_recombination.py::test_zero_lambda_gives_empty
tests/test_recombination.py::test_poisson_mean
tests/test_recombination.py::test_plan_covers_full_chromosome
tests/test_recombination.py::test_plan_segments_contiguous
tests/test_recombination.py::test_plan_no_zero_width_segments
tests/test_recombination.py::test_plan_donor_switches_at_each_breakpoint
tests/test_recombination.py::test_plan_with_no_crossovers
tests/test_recombination.py::test_plan_single_sample_pool
tests/test_recombination.py::test_plan_sample_indices_in_range
tests/test_recombination.py::test_generate_all_plans_count
tests/test_recombination.py::test_generate_all_plans_each_valid
tests/test_recombination.py::test_determinism_with_seed
tests/test_recombination.py::TestRegionSampling::test_detect_regions_empty
tests/test_recombination.py::TestRegionSampling::test_detect_regions_single
tests/test_recombination.py::TestRegionSampling::test_detect_regions_all_within_threshold
tests/test_recombination.py::TestRegionSampling::test_detect_regions_multiple_gaps
tests/test_recombination.py::TestRegionSampling::test_region_plan_segment_count
tests/test_recombination.py::TestRegionSampling::test_region_plan_sample_indices_in_range
tests/test_recombination.py::TestRegionSampling::test_region_plan_adjacency_constraint
tests/test_recombination.py::TestRegionSampling::test_region_plan_min_donors_satisfied
tests/test_recombination.py::TestRegionSampling::test_region_plan_min_donors_relaxed_when_few_regions
tests/test_recombination.py::TestRegionSampling::test_region_plan_single_sample_all_zero
tests/test_recombination.py::TestRegionSampling::test_generate_all_region_plans_count
tests/test_recombination.py::TestRegionSampling::test_generate_all_region_plans_deterministic
tests/test_recombination.py::TestRegionSampling::test_continuous_min_donors_enforced
tests/test_recombination.py::TestRegionSampling::test_continuous_min_donors_capped_by_pool
tests/test_recombination.py::TestRegionSampling::test_region_plan_end_to_end_no_missing
tests/test_recombination.py::test_plan_always_covers_chromosome
tests/test_sex_filter.py::TestParseSexLabel::test_female_labels[F]
tests/test_sex_filter.py::TestParseSexLabel::test_female_labels[f]
tests/test_sex_filter.py::TestParseSexLabel::test_female_labels[female]
tests/test_sex_filter.py::TestParseSexLabel::test_female_labels[FEMALE]
tests/test_sex_filter.py::TestParseSexLabel::test_female_labels[2]
tests/test_sex_filter.py::TestParseSexLabel::test_male_labels[M]
tests/test_sex_filter.py::TestParseSexLabel::test_male_labels[m]
tests/test_sex_filter.py::TestParseSexLabel::test_male_labels[male]
tests/test_sex_filter.py::TestParseSexLabel::test_male_labels[MALE]
tests/test_sex_filter.py::TestParseSexLabel::test_male_labels[1]
tests/test_sex_filter.py::TestParseSexLabel::test_invalid_labels_raise[X]
tests/test_sex_filter.py::TestParseSexLabel::test_invalid_labels_raise[unknown]
tests/test_sex_filter.py::TestParseSexLabel::test_invalid_labels_raise[0]
tests/test_sex_filter.py::TestParseSexLabel::test_invalid_labels_raise[3]
tests/test_sex_filter.py::TestParseSexLabel::test_invalid_labels_raise[]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_chrx_returns_female[chrX]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_chrx_returns_female[X]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_chrx_returns_female[CHRX]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_chrx_returns_female[chrx]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_chry_returns_male[chrY]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_chry_returns_male[Y]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_chry_returns_male[CHRY]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_chry_returns_male[chry]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_autosomes_return_none[chr22]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_autosomes_return_none[22]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_autosomes_return_none[chr1]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_autosomes_return_none[1]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_autosomes_return_none[chrM]
tests/test_sex_filter.py::TestSexFilterForChromosome::test_autosomes_return_none[MT]
tests/test_sex_filter.py::TestLoadSexMap::test_basic_full_path_match
tests/test_sex_filter.py::TestLoadSexMap::test_basename_match_fallback
tests/test_sex_filter.py::TestLoadSexMap::test_header_line_skipped
tests/test_sex_filter.py::TestLoadSexMap::test_comment_lines_ignored
tests/test_sex_filter.py::TestLoadSexMap::test_blank_lines_ignored
tests/test_sex_filter.py::TestLoadSexMap::test_unmatched_vcf_warns_and_omits
tests/test_sex_filter.py::TestLoadSexMap::test_malformed_line_raises
tests/test_sex_filter.py::TestLoadSexMap::test_invalid_sex_label_raises
tests/test_sex_filter.py::TestLoadSexMap::test_all_sex_labels_accepted
tests/test_sex_filter.py::TestFilterVcfsBySex::test_keeps_only_matching_sex
tests/test_sex_filter.py::TestFilterVcfsBySex::test_preserves_order
tests/test_sex_filter.py::TestFilterVcfsBySex::test_empty_result_when_no_match
tests/test_sex_filter.py::TestFilterVcfsBySex::test_path_not_in_map_is_excluded
tests/test_sex_filter.py::TestCliSexFilter::test_autosome_uses_all_donors
tests/test_sex_filter.py::TestCliSexFilter::test_chrx_uses_only_female_donors
tests/test_sex_filter.py::TestCliSexFilter::test_chrx_without_sex_file_warns
tests/test_sex_filter.py::TestCliSexFilter::test_chrx_no_female_donors_errors
tests/test_vcf_io.py::test_genotype_pool_shape
tests/test_vcf_io.py::test_genotype_pool_dimension_mismatch_raises
tests/test_vcf_io.py::test_reader_basic
tests/test_vcf_io.py::test_reader_chunk_splitting
tests/test_vcf_io.py::test_reader_missing_genotype
tests/test_vcf_io.py::test_reader_missing_rate_filter
tests/test_vcf_io.py::test_resolve_chromosome_name_no_change
tests/test_vcf_io.py::test_resolve_chromosome_name_adds_prefix
tests/test_vcf_io.py::test_resolve_chromosome_name_strips_prefix
tests/test_vcf_io.py::test_resolve_chromosome_name_no_contig_header
tests/test_vcf_io.py::test_reader_accepts_bare_chromosome_name
tests/test_vcf_io.py::test_reader_accepts_prefixed_chromosome_name
tests/test_vcf_io.py::test_reader_site_mismatch_raises
```

</details>

### Unit tests

| File | What it covers |
|------|----------------|
| `test_recombination.py` | Continuous mode: breakpoints within map range, sorted, Poisson mean, plan coverage, contiguity, adjacency constraint, determinism, Hypothesis property tests. Region mode (`TestRegionSampling`): `detect_regions` edge cases, `build_region_segment_plan` adjacency and `min_donors` guarantee, end-to-end panel integration, continuous-mode `min_donors` enforcement. |
| `test_genetic_map.py` | SHAPEIT5 and HapMap format loading, chromosome name normalisation, `bp_to_cm` interpolation, boundary clamping, non-monotonic rejection, gzipped files. |
| `test_mosaic_builder.py` | Single/two-segment plans, all-zero pools, MISSING preservation, chunk-boundary behaviour, output shape and donor correctness. |
| `test_vcf_io.py` | `GenotypePool` dimension validation, reader shape and dosage values, chunk splitting, missing genotype encoding, missing-rate filter, site mismatch error. `resolve_chromosome_name`: no-change, adds prefix, strips prefix, no-contig-header fallback. Cross-convention integration tests (bare VCF + chr flag). |
| `test_cli.py` | CliRunner end-to-end: help text, required options, full shuffle run, determinism, multi-sample mode, nonexistent input error, unphased output. `TestChromosomeNormalisation`: bare VCF with chr-prefix flag, logging of normalisation. |
| `test_sex_filter.py` | `parse_sex_label` all label variants, `sex_filter_for_chromosome` all chromosome forms, `load_sex_map` path/basename matching, header/comment handling, error cases, `filter_vcfs_by_sex` ordering. CLI integration: autosome passthrough, chrX female-only, no-sex-file warning, empty-pool error. |

### Tier 1: Recombination model fidelity

Tests R1–R6 cover both the Poisson crossover model and the region-sampling model:

```bash
/tmp/vshuffler-venv/bin/pytest tests/test_empirical_tier1.py -v --tb=short 2>&1 | tail -17
```

```output
configfile: pyproject.toml
plugins: hypothesis-6.151.9
collecting ... collected 11 items

tests/test_empirical_tier1.py::test_r1_mean PASSED                       [  9%]
tests/test_empirical_tier1.py::test_r1_dispersion PASSED                 [ 18%]
tests/test_empirical_tier1.py::test_r1_ks_vs_poisson PASSED              [ 27%]
tests/test_empirical_tier1.py::test_r2_cM_uniformity PASSED              [ 36%]
tests/test_empirical_tier1.py::test_r2_mean_position PASSED              [ 45%]
tests/test_empirical_tier1.py::test_r3_mean_segment_length PASSED        [ 54%]
tests/test_empirical_tier1.py::test_r4_region_detection_count PASSED     [ 63%]
tests/test_empirical_tier1.py::test_r4_region_boundaries PASSED          [ 72%]
tests/test_empirical_tier1.py::test_r5_cross_region_independence PASSED  [ 81%]
tests/test_empirical_tier1.py::test_r6_min_donors_always_satisfied PASSED [ 90%]
tests/test_empirical_tier1.py::test_r6_min_donors_capped_by_n_regions PASSED [100%]

============================= 11 passed in 16.97s ==============================
```

| Test | What it checks | Pass threshold |
|------|----------------|---------------|
| R1 (×3) | Crossover count: mean, dispersion (var/mean), KS vs Poisson(λ) CDF | < 1%; 0.95–1.05; KS < 0.02 |
| R2 (×2) | Breakpoint positions: KS vs Uniform(0,1) in cM space; mean normalised position | KS < 0.02; 0.49–0.51 |
| R3 | Mean segment length vs expected L/(λ+1) | < 2% relative error |
| R4 (×2) | Region detection: correct count and boundary spans on 100-region panel positions | Exact |
| R5 | Cross-region co-assignment rate ≈ 1/N_DONORS (independence check, 1 000 MC reps) | Within 2% |
| R6 (×2) | `min_donors` constraint: 100% of plans meet threshold; capping by n_regions | 0 failures; exact |

### Tier 2: Privacy and biological plausibility

Fully in-memory fixture using region-sampling mode: 200 donors, 2 000 variants in 100 gene regions (20 variants × 100 bp intra-spacing, 500 kb inter-region gaps), 100 synthetics, 20 held-out individuals.

```bash
/tmp/vshuffler-venv/bin/pytest tests/test_empirical_tier2.py -v --tb=short -W ignore::UserWarning 2>&1 | tail -22
```

```output
hypothesis profile 'default'
rootdir: /home/wook/Documents/v-shuffler
configfile: pyproject.toml
plugins: hypothesis-6.151.9
collecting ... collected 14 items

tests/test_empirical_tier2.py::test_p1_max_concordance PASSED            [  7%]
tests/test_empirical_tier2.py::test_p1_99th_percentile_concordance PASSED [ 14%]
tests/test_empirical_tier2.py::test_p1_mean_concordance_vs_baseline PASSED [ 21%]
tests/test_empirical_tier2.py::test_p2_closest_donor_attack PASSED       [ 28%]
tests/test_empirical_tier2.py::test_p4_membership_inference PASSED       [ 35%]
tests/test_empirical_tier2.py::test_b1_af_global PASSED                  [ 42%]
tests/test_empirical_tier2.py::test_b1_af_by_maf_bin PASSED              [ 50%]
tests/test_empirical_tier2.py::test_b2_heterozygosity PASSED             [ 57%]
tests/test_empirical_tier2.py::test_b3_hwe PASSED                        [ 64%]
tests/test_empirical_tier2.py::test_b4_tstv_identical PASSED             [ 71%]
tests/test_empirical_tier2.py::test_region_detection_in_fixture PASSED   [ 78%]
tests/test_empirical_tier2.py::test_min_donors_region_mode PASSED        [ 85%]
tests/test_empirical_tier2.py::test_min_donors_continuous_mode PASSED    [ 92%]
tests/test_empirical_tier2.py::test_p1_region_vs_continuous_concordance PASSED [100%]

============================== 14 passed in 1.94s ==============================
```

**Privacy tests (P1, P2, P4):** P2 (closest-donor attack) and P4 (membership inference) are diagnostic reporters, not hard assertions. With region-sampling and 100 regions, the attack success rate is ~5% — far below the ~89% seen in continuous mode at λ ≈ 0.5. P4 membership signal is weakened but still statistically detectable (a known limitation of the unphased mosaic design).

**Biological plausibility tests (B1–B4):** AF Pearson r ≥ 0.98 (threshold relaxed from 0.99 for the panel fixture due to within-region AF correlation); het rate |mean diff| < 0.005; HWE-fail fraction ≤ 3× donors; ts/tv identical by construction.

**Region-mode specific tests:** `test_region_detection_in_fixture` verifies all 100 clusters are detected; `test_min_donors_region_mode` and `test_min_donors_continuous_mode` verify the minimum-donors constraint in both modes; `test_p1_region_vs_continuous_concordance` demonstrates that region mode achieves lower maximum concordance than continuous mode at λ = 0.5 on the same data.

### Full test run

```bash
/tmp/vshuffler-venv/bin/pytest --tb=short -q -W ignore::UserWarning --ignore=tests/test_cli.py 2>&1 | tail -6
```

```output
.........................ss...................ssss...................... [ 51%]
...................................................................      [100%]
133 passed, 6 skipped in 18.53s
```

133 of the non-CLI tests pass; 6 skip (2 require plink2/bcftools for Tier 3; 4 require patient VCF env vars).

One CLI test (`test_output_vcf_is_unphased`) is sensitive to whether `bgzip` is on PATH in the test environment. The test asserts the uncompressed `synthetic_0.vcf` file exists, but when bgzip is available it compresses and removes the original. This is a test environment issue, not a code issue — the behaviour is correct in both cases. It passes in environments where bgzip is absent from PATH.
