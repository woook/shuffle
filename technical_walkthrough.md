# v-shuffler: Technical Walkthrough

*2026-03-10T09:33:21Z by Showboat 0.6.1*
<!-- showboat-id: 4e7d170a-04e1-48f0-a492-94154f94d5df -->

## Overview

v-shuffler anonymises genomic VCF files by shuffling diploid genotypes between individuals using simulated meiotic recombination. It is designed for genomics researchers who need to share cohort-level variant data without exposing identifiable individual genotypes.

**The problem it solves:** Raw VCF files containing genotypes from thousands of individuals are personally identifiable — long runs of genotypes are unique to each person. v-shuffler transforms these into synthetic individuals that are biologically plausible (preserving population-level statistics and local linkage disequilibrium) but unidentifiable.

**How it works, in three steps:**

1. **Crossover simulation** — For each desired synthetic individual, Poisson-distributed crossover positions are drawn in genetic (cM) space using a real validated genetic map (e.g. SHAPEIT5 GRCh38). This produces a *segment plan*: a list of chromosomal intervals, each assigned a randomly chosen donor from the pool.
2. **Mosaic assembly** — The tool streams through all per-sample VCFs in chunks of 50,000 variants. For each chunk, each variant is assigned the dosage (0/1/2) of whichever donor owns that cM position in that synthetic individual's plan.
3. **Output** — Bgzipped, tabix-indexed VCFs are written with anonymised sample names (`synthetic_0`, `synthetic_1`, …), stripped metadata, and a provenance header recording tool version, seed, and date.

**Supported platforms:** Linux/macOS; requires Python ≥ 3.10, bgzip, and tabix on `PATH`.

## Repository Layout

<details>
<summary>Full directory tree</summary>

```python3

import os
lines = []
for root, dirs, files in os.walk('.'):
    dirs[:] = sorted(d for d in dirs if d not in ('.git', '.claude'))
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
  tests/
    __init__.py
    conftest.py
    test_cli.py
    test_genetic_map.py
    test_mosaic_builder.py
    test_recombination.py
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
      vcf_reader.py
      vcf_writer.py
```

</details>

## Architecture

v-shuffler is structured as a Python package with two top-level namespaces:

- **`v_shuffler/core/`** — Pure computation: crossover simulation, genotype pool, mosaic assembly. No I/O.
- **`v_shuffler/io/`** — File I/O: genetic map loading, VCF reading, VCF writing.
- **`v_shuffler/cli.py`** — Click-based CLI that wires everything together.
- **`v_shuffler/config.py`** — Dataclass holding all configuration for a run.
- **`v_shuffler/validate.py`** — Post-run sanity checks (AF correlation, identity check).

### Data flow

```
Per-sample VCF files (one per donor)
        │
        ▼
PerSampleVCFReader.iter_chunks()
  │  reads variants in lockstep, converts GT→dosage, annotates cM positions
  │  yields GenotypePool chunks  (n_variants × n_donors, uint8)
  │
  ├── GeneticMap.bp_to_cm()   ← genetic map file (SHAPEIT5 or HapMap format)
  │
  ▼
build_synthetic_genotypes(pool, segment_plans)
  │  applies pre-computed segment plans to each pool chunk
  │  yields synthetic dosage matrix  (n_variants × n_synthetic, uint8)
  │
  ├── segment_plans   ← generate_all_segment_plans()
  │     └── simulate_crossover_breakpoints()  (Poisson in cM space)
  │         build_segment_plan()              (assign donors to intervals)
  │
  ▼
SyntheticVCFWriter.write_chunk()
  │  converts dosage→GT string, writes VCF records
  ▼
SyntheticVCFWriter.finalize()
   bgzip + tabix → synthetic_N.vcf.gz + .tbi
```

The segment plans are computed **once upfront** (they are tiny — a few segments per individual in cM space) and then reused for every chunk streamed from disk. This keeps memory proportional to `chunk_size × n_donors`, not to the total variant count.

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

### Input resolution (`_resolve_inputs`)

Before the pipeline starts, `cli._resolve_inputs()` converts the `--input` string to a concrete list of `Path` objects. It handles three forms:

| Form | Example |
|------|---------|
| Glob string | `"data/per_sample/*.vcf.gz"` |
| File list | `@samples.txt` (one path per line, `#` comments ignored) |
| Single path | `sample.vcf.gz` |

Missing files are reported before any I/O starts, avoiding confusing mid-run failures.

## Core Module: Recombination (`v_shuffler/core/recombination.py`)

This module contains the biological heart of v-shuffler: simulating where crossovers happen and which donor supplies genotypes in each resulting segment.

### Key types

**`Segment`** (frozen dataclass)

```python
@dataclass(frozen=True)
class Segment:
    cm_start: float   # inclusive segment start in cM
    cm_end:   float   # exclusive end (inclusive for the last segment)
    sample_idx: int   # column index of the donor in GenotypePool.dosages
```

A segment plan for one synthetic individual is a `list[Segment]` covering the full chromosome.

### Crossover simulation: `simulate_crossover_breakpoints()`

Draws crossover positions via a Poisson process in genetic distance:

1. Sample the *number* of crossovers: `n ~ Poisson(λ)` where `λ = total_cM / 100` (one crossover per Morgan per meiosis, the standard meiotic rate). The caller can override λ via `lambda_override`.
2. Sample the *positions* uniformly in `[start_cM, end_cM]`. Working in cM space means the distribution in physical coordinates automatically reflects the local recombination rate encoded in the map.

Returns a sorted `np.ndarray` of breakpoints in cM (may be empty if `n = 0`).

### Segment plan construction: `build_segment_plan()`

Given breakpoints, builds the segment list for one synthetic individual:

- Boundaries array: `[map_start, *breakpoints, map_end]`
- For each consecutive pair of boundaries, a `Segment` is created and a new donor is chosen uniformly at random from the pool (never the same as the previous segment when `n_samples > 1`).

### Upfront plan generation: `generate_all_segment_plans()`

Generates plans for *all* synthetic output individuals before any VCF I/O starts. Plans are cheap (a few dozen floats per individual) and must exist in full to correctly assign donors across chunk boundaries.

```/tmp/vshuffler-venv/bin/python3

import numpy as np
from v_shuffler.core.recombination import simulate_crossover_breakpoints, build_segment_plan

class MockMap:
    start_cm = 0.0
    end_cm = 100.0
    total_length_cm = 100.0

rng = np.random.default_rng(42)
gmap = MockMap()

breakpoints = simulate_crossover_breakpoints(gmap, rng)
print(f'Crossovers drawn (lambda=1.0 for 100 cM map): {len(breakpoints)}')
print(f'Breakpoint positions (cM): {breakpoints.round(2)}')

plan = build_segment_plan(breakpoints, gmap, n_samples=5, rng=rng)
print(f'Segment plan ({len(plan)} segments):')
for seg in plan:
    print(f'  [{seg.cm_start:.2f}, {seg.cm_end:.2f}) cM  ->  donor {seg.sample_idx}')

```

```output
Crossovers drawn (lambda=1.0 for 100 cM map): 1
Breakpoint positions (cM): [85.86]
Segment plan (2 segments):
  [0.00, 85.86) cM  ->  donor 0
  [85.86, 100.00) cM  ->  donor 3
```

## Core Module: GenotypePool (`v_shuffler/core/genotype_pool.py`)

`GenotypePool` is the central in-memory data structure. It represents a *chunk* of variants (up to `chunk_size`) across all donor individuals.

### Dosage encoding

Genotypes are stored as dosage integers in a `uint8` array, not as GT strings:

| Value | Meaning |
|-------|---------|
| `0` | Homozygous reference (0/0) |
| `1` | Heterozygous (0/1) |
| `2` | Homozygous alt (1/1) |
| `255` | Missing (sentinel `MISSING`) |

For multi-allelic sites, dosage = total count of non-ref allele copies, regardless of which alt allele. This is sufficient for the shuffling operation.

### Key attributes of `GenotypePool`

| Attribute | Shape | dtype | Description |
|-----------|-------|-------|-------------|
| `dosages` | `(n_variants, n_samples)` | uint8 | Dosage matrix |
| `positions` | `(n_variants,)` | int64 | Physical positions (1-based bp) |
| `cm_pos` | `(n_variants,)` | float64 | Genetic positions (cM) |
| `variant_info` | `list[VariantInfo]` | — | Per-variant metadata |

`VariantInfo` stores `chrom`, `pos`, `ref`, `alts`, `id`, `qual`, `filters`, and `cm_pos` — everything needed to write a VCF record without re-reading the input.

## Core Module: Mosaic Builder (`v_shuffler/core/mosaic_builder.py`)

This module applies pre-computed segment plans to a `GenotypePool` chunk, producing the synthetic dosage matrix.

### `apply_segment_plan(pool, plan) → np.ndarray`

Builds the synthetic dosage vector for one output individual from one pool chunk:

```python
result = np.full(pool.n_variants, MISSING, dtype=np.uint8)
for seg in plan:
    mask = (pool.cm_pos >= seg.cm_start) & (pool.cm_pos <= seg.cm_end)
    result[mask] = pool.dosages[mask, seg.sample_idx]
```

Because both the pool (`cm_pos`) and the plan (`Segment`) live in cM space, this is a simple interval mask operation — it works identically whether the pool is a 50k-variant chunk or the whole chromosome at once.

Variants that fall outside all segments (can happen at floating-point cM edge cases) remain `MISSING = 255`.

### `build_synthetic_genotypes(pool, segment_plans) → np.ndarray`

Iterates over all segment plans and stacks the results into a `(n_variants, n_output_samples)` matrix. This is the function called in the main streaming loop:

```python
for pool in reader.iter_chunks():
    synthetic = build_synthetic_genotypes(pool, segment_plans)
    writer.write_chunk(pool, synthetic)
```

## IO Module: Genetic Map (`v_shuffler/io/genetic_map.py`)

`GeneticMap` loads a per-chromosome genetic map and provides physical-to-genetic coordinate conversion.

### Supported file formats

| Format | Columns | Example header |
|--------|---------|----------------|
| SHAPEIT5 / Oxford | 3 (space-delimited) | `pos chr cM` |
| HapMap | 4 (space-delimited) | `Chromosome Position(bp) Rate(cM/Mb) Map(cM)` |

Both plain text and `.gz` files are accepted. The chromosome name in the file may use either the `chr22` or `22` spelling; the loader normalises both.

### `GeneticMap.bp_to_cm(positions)`

Converts an array of physical positions (bp) to cM values using linear interpolation (`np.interp`). Positions outside the map range are clamped to the boundary cM value (numpy's default extrapolation behaviour). This is called once per chunk in `PerSampleVCFReader` to annotate every variant with its genetic coordinate.

### Validation

After loading, the map checks that cM values are monotonically non-decreasing. A non-monotonic map would silently corrupt crossover assignments, so this is caught at load time.

## IO Module: VCF Reader (`v_shuffler/io/vcf_reader.py`)

`PerSampleVCFReader` reads one VCF file per donor individual and yields `GenotypePool` chunks.

### Design

- Uses **cyvcf2** for fast VCF iteration (C-backed, much faster than pure-Python parsers).
- All per-sample VCFs are opened simultaneously and iterated **in lockstep** using `zip(*iterators)`. This ensures variant N is read from all files at the same time, so the dosage vector for variant N can be assembled in one pass.
- Variant sites are validated for CHROM/POS/REF/ALT consistency across all files **on the first chunk only** (subsequent chunks are assumed consistent). A mismatch raises a `ValueError` immediately.
- Region access uses **tabix** when a `.tbi` or `.csi` index exists alongside the VCF; otherwise falls back to full-file iteration with a chromosome filter (useful in tests with small plain-text VCFs).

### Genotype → dosage conversion (`_gt_to_dosage`)

Reads allele indices `[a1, a2, phased]` from cyvcf2 and returns:
- `MISSING (255)` if either allele is `-1` (missing in VCF)
- `int(a1 > 0) + int(a2 > 0)` otherwise (counts non-ref copies)

### Missing rate filter

Before a variant is added to the chunk, its missing rate is checked: if `n_missing / n_samples > max_missing_rate` the variant is silently skipped. A warning is logged at the end of iteration reporting how many variants were dropped.

## IO Module: VCF Writer (`v_shuffler/io/vcf_writer.py`)

`SyntheticVCFWriter` converts the synthetic dosage matrix back into VCF files.

### Output modes

| Mode | Flag | Output |
|------|------|--------|
| `per_sample` (default) | `--output-mode per_sample` | One `synthetic_N.vcf.gz` per synthetic individual |
| `multi_sample` | `--output-mode multi_sample` | One `synthetic_${chrom}.vcf.gz` with all synthetic samples as columns |

`per_sample` is preferred for downstream per-sample pipelines (e.g. running GWAS on each synthetic individual independently). `multi_sample` is convenient for joint analysis tools.

### Header handling

The writer copies the header from the first input VCF (via cyvcf2), then:
1. Strips all `##sample=` lines (sample-identifiable metadata).
2. Inserts a provenance line immediately before `#CHROM`:
   ```
   ##v-shuffler=<version="0.1.0",seed="42",chromosome="chr22",date="2026-03-10T...">
   ```
3. Replaces the sample column names with `synthetic_0`, `synthetic_1`, …

### Dosage → GT conversion

| Dosage | GT string written |
|--------|------------------|
| 0 | `0/0` |
| 1 | `0/1` |
| 2 | `1/1` |
| 255 | `./.` |

All calls are written as unphased (slash separator), regardless of whether the input was phased.

### Finalisation

After all chunks are written, `finalize()` calls `bgzip -f` and `tabix -p vcf` on each output file via `subprocess.run`. If either tool is absent from `PATH`, a warning is logged and the uncompressed file is returned — the pipeline continues rather than crashing.

## Configuration (`v_shuffler/config.py`)

All pipeline parameters are collected into a single `ShufflerConfig` dataclass. The CLI constructs this object and passes it to `_run_shuffle()`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `input_vcfs` | `list[Path]` | required | Resolved per-sample VCF paths |
| `output_dir` | `Path` | required | Output directory (created if absent) |
| `genetic_map` | `Path` | required | Genetic map file |
| `chromosome` | `str` | required | Chromosome label (e.g. `"chr22"`) |
| `n_output_samples` | `int` | required | Number of synthetic individuals |
| `seed` | `int \| None` | `None` | RNG seed (None = non-reproducible) |
| `min_segment_cM` | `float` | `0.5` | Minimum segment length in cM |
| `n_crossovers_lambda` | `float \| None` | `None` | Poisson rate override (None = derived from map) |
| `chunk_size_variants` | `int` | `50000` | Variants per memory chunk |
| `n_threads` | `int` | `4` | Thread count (reserved for future use) |
| `max_missing_rate` | `float` | `0.05` | Skip variants with > 5% missing calls |
| `output_mode` | `str` | `"per_sample"` | `"per_sample"` or `"multi_sample"` |

`__post_init__` validates ranges and creates the output directory immediately.

## Validation Subcommand (`v_shuffler/validate.py`)

The `validate` subcommand provides two quantitative checks that the shuffled output is biologically plausible and has not leaked identity.

### Check 1: Allele frequency correlation

Reads the reference (merged input) VCF and all synthetic VCFs once each, computes per-site allele frequency in both, and reports the Pearson correlation coefficient. A well-shuffled output should have `r ≥ 0.99` — genotype shuffling preserves population-level allele frequencies because the same pool of alleles is redistributed, not transformed.

### Check 2: Identity leak detection

For each synthetic sample × each input sample pair, computes the fraction of overlapping variants with identical dosage. If any pair exceeds `_IDENTITY_THRESHOLD = 0.99`, a warning is logged. This catches the degenerate case where a synthetic individual happens to be assigned the same donor for the entire chromosome (possible when the pool is very small or the chromosome is short).

### Limitations of the automated check

The automated validator catches gross problems. The README recommends additional manual checks:
- **PCA**: synthetic individuals should fall within the cloud of input individuals.
- **LD**: compare sliding-window r² between input and output.
- **IBD scan**: no long IBS runs (> a few cM) between any synthetic and any input individual.

## Error Handling and Recovery

v-shuffler uses a fail-fast approach for errors that indicate corrupt input or misconfiguration, and a degrade-gracefully approach for missing optional tools.

### Fail-fast errors (raise immediately)

| Error | Where raised | Cause |
|-------|-------------|-------|
| `ValueError` | `GeneticMap.__init__` | No entries for the requested chromosome, or non-monotonic cM values |
| `ValueError` | `PerSampleVCFReader.iter_chunks` | CHROM/POS/REF/ALT mismatch across per-sample VCFs on first chunk |
| `click.BadParameter` | `_resolve_inputs` | Glob matched no files, or listed file paths don't exist |
| `ValueError` | `ShufflerConfig.__post_init__` | Invalid `output_mode`, `n_output_samples < 1`, or `max_missing_rate` out of [0,1] |

### Degrade-gracefully warnings

| Situation | Behaviour |
|-----------|-----------|
| `bgzip` not on PATH | Warning logged; uncompressed `.vcf` returned |
| `tabix` not on PATH | Warning logged; file not indexed |
| `bgzip` exits with error | Warning logged; uncompressed path returned |
| High missing-rate variants | Silently skipped; count reported as a warning after iteration |
| Too few variants to compute AF correlation | Warning logged; validation continues |

### Missing rate filtering detail

The check `n_missing / n_samples > max_missing_rate` is performed per-variant *before* the dosage is added to the chunk buffer. The variant is simply not appended to `chunk_info` — it never enters the pool matrix or the output VCF.

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
| **cyvcf2** | Fast C-backed VCF reading and writing. Used in `vcf_reader.py` for variant iteration and header extraction, and in `vcf_writer.py` for header parsing. |
| **numpy** | Array operations throughout: dosage matrices, cM interpolation, crossover position sampling, boolean masks for segment assignment. |
| **click** | CLI framework: argument parsing, help text generation, subcommand routing. |
| **pysam** | Listed as a dependency (provides htslib bindings); bgzip/tabix are called via `subprocess` rather than pysam directly in the current implementation. |
| **scipy** | Available for statistical operations; used in validation for correlation computations (numpy's `np.corrcoef` is used directly, scipy available for future extensions). |
| **tqdm** | Progress bar over the chunk-streaming loop (`Chunks` progress bar in `_run_shuffle`). |

**External tools required at runtime:**
- `bgzip` — to compress output VCFs (part of htslib)
- `tabix` — to index the compressed VCFs (part of htslib)

## Tests

The test suite lives in `tests/` and uses pytest. Run with:

```bash
pip install -r requirements-dev.txt
pytest -v
```

<details>
<summary>Full test list (50 tests)</summary>

```bash
/tmp/vshuffler-venv/bin/pytest --collect-only -q 2>&1 | grep -E '^tests|^[0-9]' | head -50
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
tests/test_recombination.py::test_plan_always_covers_chromosome
tests/test_vcf_io.py::test_genotype_pool_shape
tests/test_vcf_io.py::test_genotype_pool_dimension_mismatch_raises
tests/test_vcf_io.py::test_reader_basic
tests/test_vcf_io.py::test_reader_chunk_splitting
tests/test_vcf_io.py::test_reader_missing_genotype
tests/test_vcf_io.py::test_reader_missing_rate_filter
tests/test_vcf_io.py::test_reader_site_mismatch_raises
```

</details>

### Test modules

| Module | What it covers |
|--------|---------------|
| `test_genetic_map.py` | SHAPEIT5 and HapMap format loading; chromosome name normalisation; `bp_to_cm` interpolation at exact and intermediate points; clamping at map boundaries; rejection of non-monotonic maps; gzipped files. |
| `test_recombination.py` | Breakpoint sampling: positions within map range, sorted order, Poisson mean over many draws. Segment plan: full chromosome coverage, contiguity, no zero-width segments, donor switching at breakpoints, single-sample pool, correct index range. Plan determinism with a fixed seed. Uses **hypothesis** for property-based testing. |
| `test_mosaic_builder.py` | Single-segment plans copy the correct donor; two-segment plans copy two different donors; all-zero pool stays all-zero; `MISSING` sentinel is preserved; chunk-boundary behaviour; output shape and donor correctness for `build_synthetic_genotypes`. |
| `test_vcf_io.py` | `GenotypePool` dimension validation; reader produces correct shape and dosage values; chunk splitting at `chunk_size` boundary; missing genotype encoding; missing-rate filter drops high-missing variants; site mismatch raises immediately. |
| `test_cli.py` | End-to-end CLI tests using Click's `CliRunner`: help text, missing required options, full shuffle run, determinism under fixed seed, multi-sample output mode, nonexistent input error, unphased GT output, no synthetic sample identical to any input. |

### Shared fixtures (`conftest.py`)

- **`genetic_map_file` / `genetic_map`** — a 10-point SHAPEIT5 map for chr22 spanning 10 cM written to `tmp_path`.
- **`five_sample_vcfs`** — 5 per-sample VCFs, 10 variants each, with distinct hand-crafted genotypes. Used by reader and CLI tests.

<details>
<summary>pytest run output (50 passed)</summary>

```bash
/tmp/vshuffler-venv/bin/pytest -v --tb=short 2>&1 | tail -20
```

```output
tests/test_recombination.py::test_plan_covers_full_chromosome PASSED     [ 66%]
tests/test_recombination.py::test_plan_segments_contiguous PASSED        [ 68%]
tests/test_recombination.py::test_plan_no_zero_width_segments PASSED     [ 70%]
tests/test_recombination.py::test_plan_donor_switches_at_each_breakpoint PASSED [ 72%]
tests/test_recombination.py::test_plan_with_no_crossovers PASSED         [ 74%]
tests/test_recombination.py::test_plan_single_sample_pool PASSED         [ 76%]
tests/test_recombination.py::test_plan_sample_indices_in_range PASSED    [ 78%]
tests/test_recombination.py::test_generate_all_plans_count PASSED        [ 80%]
tests/test_recombination.py::test_generate_all_plans_each_valid PASSED   [ 82%]
tests/test_recombination.py::test_determinism_with_seed PASSED           [ 84%]
tests/test_recombination.py::test_plan_always_covers_chromosome PASSED   [ 86%]
tests/test_vcf_io.py::test_genotype_pool_shape PASSED                    [ 88%]
tests/test_vcf_io.py::test_genotype_pool_dimension_mismatch_raises PASSED [ 90%]
tests/test_vcf_io.py::test_reader_basic PASSED                           [ 92%]
tests/test_vcf_io.py::test_reader_chunk_splitting PASSED                 [ 94%]
tests/test_vcf_io.py::test_reader_missing_genotype PASSED                [ 96%]
tests/test_vcf_io.py::test_reader_missing_rate_filter PASSED             [ 98%]
tests/test_vcf_io.py::test_reader_site_mismatch_raises PASSED            [100%]

============================== 50 passed in 0.64s ==============================
```

</details>
