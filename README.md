# v-shuffler

Anonymise genomic VCF files by shuffling variants between individuals using
simulated genetic recombination.  Each output individual is a mosaic of
genotype segments from different real donors, stitched together at crossover
positions drawn from a real genetic map.  This preserves local linkage
disequilibrium (haplotype block structure) while making individuals
unidentifiable.

## How it works

1. **Segment plan generation** — For each synthetic output individual, the
   tool simulates meiotic crossovers along the chromosome using a Poisson
   process scaled to the genetic map length (≈1 crossover per Morgan per
   meiosis).  Crossover positions are drawn uniformly in cM space, so
   their distribution in physical coordinates automatically reflects the
   local recombination rate.

2. **Genotype mosaic assembly** — The chromosome is divided into segments
   by the crossover positions.  Each segment is assigned a randomly chosen
   donor individual.  The tool streams variants in chunks and copies
   the donor's diploid genotype calls (0/0, 0/1, 1/1) for each segment.

3. **Output** — New VCF files are written with anonymised sample names
   (`synthetic_0`, `synthetic_1`, …).  Sample-identifying header metadata
   is stripped.  A provenance line records the tool version and random seed.

**Why this is biologically plausible:** Each segment is a contiguous stretch
of real genotypes from one real person, so short-range LD is perfectly
preserved.  The recombination positions are drawn from a validated genetic
map (e.g. SHAPEIT5 GRCh38), matching the true recombination landscape.

**Why this is unidentifiable:** Each output individual is a mosaic of many
donors.  With a typical chromosome of 100 cM and ≈1–2 crossovers, an output
individual typically has 2–3 segments from 2–3 different donors, so no
long run of genotypes uniquely identifies any original sample.

---

## Installation

```bash
pip install -e .
# or for development:
pip install -e ".[dev]"
```

Dependencies: `cyvcf2`, `numpy`, `click`, `pysam`, `scipy`, `tqdm`.

---

## Genetic map

v-shuffler uses HapMap-format genetic maps.  SHAPEIT5 provides high-quality
GRCh38 maps:

```
https://github.com/odelaneau/shapeit5/tree/main/maps
```

Download the per-chromosome map files (e.g. `chr22.b38.gmap.gz`) and pass
the path to `--genetic-map`.

---

## Usage

### Shuffle

```bash
v-shuffler shuffle \
    --input "data/per_sample/*.vcf.gz" \
    --output-dir shuffled/ \
    --genetic-map chr22.b38.gmap.gz \
    --chromosome chr22 \
    --n-samples 1000 \
    --seed 42
```

Or use a file list:

```bash
ls data/per_sample/*.vcf.gz > samples.txt
v-shuffler shuffle \
    --input @samples.txt \
    --output-dir shuffled/ \
    --genetic-map chr22.b38.gmap.gz \
    --chromosome chr22 \
    --n-samples 1000 \
    --seed 42
```

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--input` | required | Glob or `@filelist.txt` of per-sample VCFs |
| `--output-dir` | required | Output directory |
| `--genetic-map` | required | SHAPEIT5 / HapMap format map file |
| `--chromosome` | required | Chromosome to process, e.g. `chr22` |
| `--n-samples` | = n inputs | Number of synthetic individuals to produce |
| `--seed` | None | Random seed (for reproducibility) |
| `--output-mode` | `per_sample` | `per_sample` or `multi_sample` |
| `--chunk-size` | 50000 | Variants per memory chunk |
| `--max-missing` | 0.05 | Skip variants with >N% missing calls |

### Validate

After shuffling, run the validate command to check biological plausibility:

```bash
v-shuffler validate \
    --input "shuffled/synthetic_*.vcf.gz" \
    --reference-vcf merged_input.vcf.gz \
    --chromosome chr22
```

Checks:
- Allele frequency correlation between input and output (should be r > 0.99)
- No output sample is near-identical to any single input sample

---

## Input requirements

- **One VCF per sample** — approximately 1000 files.
- VCFs need **not** be phased (unphased `0/1` calls are used directly).
- All VCFs must cover **the same set of variants** (same CHROM/POS/REF/ALT
  in the same order).  Use `bcftools isec` or `bcftools view --regions` to
  harmonise sites if needed.
- VCFs should be sorted, bgzipped and tabix-indexed for random access.

### Processing multiple chromosomes

Run v-shuffler once per chromosome:

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
```

Then concatenate the per-chromosome per-sample VCFs with `bcftools concat`.

---

## Limitations

- Sex chromosomes (chrX, chrY) are not supported; only autosomes.
- Multi-allelic sites are encoded as total alt allele count (dosage); the
  specific alt allele identity is not preserved if dosage > 0.
- The tool performs **genotype-level** (diploid) segment swapping, not true
  haplotype recombination.  This is correct for unphased input and
  preserves local LD, but the output is not equivalent to phased
  recombination.

---

## Running tests

```bash
pip install -r requirements-dev.txt
pytest -v
```

---

## Biological validation (manual)

After shuffling, check that the output is biologically plausible:

**PCA check** — compute principal components of input + output jointly.
Output individuals should fall within the cloud of input individuals.

```bash
# Using plink2
plink2 --vcf merged_input_and_output.vcf.gz --pca 10 --out pca
```

**LD check** — compare LD (r²) in sliding windows between input and output.
Values should be similar.

**IBD check** — scan for long identical-by-state (IBS) runs between each
output individual and all input individuals.  No run longer than a few cM
should be shared, confirming that the recombination is working as expected.
