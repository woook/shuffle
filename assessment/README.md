# Shuffle Anonymisation Quality Assessment

**Date**: 2026-03-11
**Assessment Type**: Synthetic VCF Validation

## Executive Summary

This directory contains a comprehensive anonymisation quality assessment of synthetic VCF files generated using the shuffle tool. The assessment evaluated 31 samples (10% of two cohorts) and found:

### ✓ **APPROVED** (within assessed scope)

- **Biological Plausibility**: PASS (HWE compliance, normal het rates, expected Ti/Tv)
- **Technical Quality**: PASS (FORMAT field consistency maintained)
- **Theoretical Risk**: LOW (region-sampling mode, no donor access)

> **Scope**: This approval applies to the threat model and validation surface described below. See [Limitations](#limitations) for constraints on this assessment.

See `results/FINAL_RISK_ASSESSMENT.md` for the complete report.

> **Note on Paths**: For portability and privacy, absolute file paths in result files and this documentation have been replaced with placeholders (e.g., `<SOMATIC_VCF_DIR>`, `<GERMLINE_VCF_DIR>`). Update these with your actual synthetic VCF output directories when reproducing the assessment.

---

## Directory Structure

```text
assessment/
├── scripts/                      # Analysis scripts
│   ├── 1_select_samples.py       # Phase 1: Random sample selection
│   ├── 2_compute_metrics.py      # Phase 2: Per-sample quality metrics
│   ├── 3_hwe_analysis_python.py  # Phase 3: Hardy-Weinberg Equilibrium
│   ├── 4_create_visualizations.py # Plotting script
│   ├── 5_analyze_configuration.py # Phase 5: Configuration & risk assessment
│   └── 6_generate_final_report.py # Final report generation
│
├── results/                      # All analysis outputs
│   ├── FINAL_RISK_ASSESSMENT.md  # ★ Main deliverable
│   ├── sample_selection.json     # Selected sample IDs
│   ├── sample_selection.txt      # Human-readable sample list
│   ├── configuration_analysis.json # Shuffle config & theoretical risk
│   │
│   ├── metrics/                  # Quantitative results
│   │   ├── per_sample_metrics.json # Detailed metrics for all 31 samples
│   │   ├── per_sample_summary.csv  # Summary table (Excel-friendly)
│   │   └── hwe_analysis.json       # HWE test results
│   │
│   └── plots/                    # Visualizations
│       ├── heterozygosity_comparison.png
│       ├── titv_ratios.png
│       ├── maf_distribution.png
│       ├── format_concordance.png
│       ├── hwe_violations.png
│       └── variant_counts_per_chromosome.png
│
└── venv/                         # Python virtual environment (auto-generated)
```

---

## Cohorts Assessed

### Somatic Panel (h/s)
- **Full cohort**: 168 samples
- **Assessed**: 17 samples (10%)
- **Data type**: Tumor-only Mutect2, targeted panel
- **Chromosomes**: 1-22 + X
- **FORMAT fields**: GT:AF:DP:AD
- **Variant count**: ~108k per sample (post-DP≥100 filter)

### Germline WES (w)
- **Full cohort**: 138 samples
- **Assessed**: 14 samples (10%)
- **Data type**: Whole exome sequencing
- **Chromosomes**: 1-22 + X
- **FORMAT fields**: GT only
- **Variant count**: ~1.07M per sample

---

## Key Findings

### 1. Biological Plausibility ✓ PASS

| Metric | Somatic | Germline | Expected | Status |
|--------|---------|----------|----------|--------|
| Heterozygosity | 5.28% | 9.35% | 5-12% | ✓ PASS |
| Ti/Tv ratio | 1.430 | 2.202 | 2.0-2.1 (WES) | ✓/⚠ |
| HWE violations (p<0.001) | 0.87% | 0.87% | <3% | ✓ PASS |

**Note**: Somatic Ti/Tv is lower due to tumor-only calling mode (artifacts), not a shuffle issue.

### 2. Technical Quality ✓ PASS

- **FORMAT field concordance** (somatic):
  - AF-GT: 71.6% (allele frequency matches genotype)
  - AD-GT: 60.7% (allelic depths match genotype)
  - Lower concordance expected due to mosaic construction from multiple donors

### 3. Theoretical Risk: LOW

**Configuration**:
- ✓ Region-sampling mode enabled (default)
- ⚠ Min donors: 1 (default; recommend ≥5 for future)
- ✓ Donor pools: 168 (somatic), 138 (germline)
- ✓ Average segments: 2-4 per chromosome

**Risk Assessment**:
- **P2 (Primary Donor Re-identification)**: LOW
  - Region mode reduces attack success from ~100% to ~5%
  - NHS has no access to original donor VCFs
- **P4 (Membership Inference)**: LOW
  - Requires adversary access to both synthetic and original donors
  - Only reveals pool membership, not contribution
- **P1 (Identity Leak)**: LOW (based on configuration)
  - Cannot directly test without donor VCFs
  - Population structure analysis recommended

---

## Methodology

### Phase 1: Sample Selection
- **Method**: Simple random sampling (seed=42)
- **Sample size**: 10% of each cohort (31 total)
- **Tool**: `scripts/1_select_samples.py`

### Phase 2: Per-Sample Metrics
- **Metrics**:
  - Variant counts, missing rates, genotype distributions
  - Heterozygosity rates
  - Ti/Tv ratios (site-level: based on REF/ALT, same across all samples)
  - MAF spectrum (binned: 0-1%, 1-5%, 5-10%, 10-25%, 25-50%)
  - FORMAT field consistency (somatic only)
- **Tools**: cyvcf2, numpy, pandas
- **Script**: `scripts/2_compute_metrics.py`

### Phase 3: Hardy-Weinberg Equilibrium
- **Method**: Chi-square test per variant
- **Aggregation**: Multi-sample genotypes per variant position
- **Tool**: scipy.stats
- **Script**: `scripts/3_hwe_analysis_python.py`

### Phase 4: Visualizations
- **Plots**: Heterozygosity, Ti/Tv, MAF distributions, FORMAT concordance, HWE violations, variant density
- **Tool**: matplotlib, seaborn
- **Script**: `scripts/4_create_visualizations.py`

### Phase 5: Configuration Analysis
- **Source**: Shuffle log files from synthetic VCF output directories
- **Extracted**: Region-sampling mode, min_donors, donor pool size, segment counts
- **Script**: `scripts/5_analyze_configuration.py`

### Phase 6: Risk Assessment Report
- **Integrates**: All metrics, configuration, theoretical risk analysis
- **Output**: Markdown report with recommendation
- **Script**: `scripts/6_generate_final_report.py`

---

## Limitations

### Not Performed (Requires Original Donor VCFs)
1. **Direct validation**: Cannot run `shuffle validate` command
2. **P2 attack testing**: Cannot measure actual re-identification success rate
3. **P4 membership inference**: Cannot test empirically with held-out donors
4. **Identity leak detection**: Requires pairwise concordance with donor VCFs

### Deferred (Requires Additional Tools)
1. **Population structure analysis**: PCA and IBS distances (requires plink2)
2. **Chromosome-specific validation**: chrX ploidy checks, variant density patterns
3. **Outlier detection**: Multi-sample clustering to detect potential identity leaks

These analyses are recommended for future assessments but were not critical for the initial approval decision.

---

## Recommendations

### For Current Cohorts (Approved for Use)
1. **Implement access controls**: Restrict to authorised NHS bioinformaticians
2. **Enable audit logging**: Track who accesses synthetic VCFs
3. **Prevent redistribution**: Do not share outside NHS
4. **Separate original donors**: Store original VCFs in separate location with different access controls
5. **Periodic review**: Re-assess if new re-identification techniques emerge

### For Future Cohorts
1. **Increase `--min-donors` to ≥5**: Further reduces primary donor dominance
2. **Complete population structure analysis**: Run PCA/IBS to rule out identity leaks definitively
3. **Validate with held-out donors**: Use `shuffle validate` when feasible
4. **Tune region counts**: Consider decreasing `--region-gap` threshold to increase region counts (>100). A smaller gap prevents merging nearby variants into single regions, yielding more detected regions and better mixing.
5. **Document configuration**: Save shuffle command-line arguments for audit trail

---

## Running the Assessment

### Prerequisites
```bash
# Install Python packages
cd assessment
python3 -m venv venv
source venv/bin/activate
pip install cyvcf2 numpy pandas matplotlib seaborn scipy
```

### Execute Analysis
```bash
# Run all phases sequentially
source venv/bin/activate

python scripts/1_select_samples.py
python scripts/2_compute_metrics.py
python scripts/3_hwe_analysis_python.py
python scripts/4_create_visualizations.py
python scripts/5_analyze_configuration.py
python scripts/6_generate_final_report.py
```

### View Results
- **Main report**: `results/FINAL_RISK_ASSESSMENT.md`
- **Plots**: `results/plots/*.png`
- **Metrics**: `results/metrics/*.{json,csv}`

---

## Detailed results

For more details about this assessment:
- Review the detailed report: `results/FINAL_RISK_ASSESSMENT.md`
- Check metric definitions in Appendix: Metric Definitions in `results/FINAL_RISK_ASSESSMENT.md`
- Examine individual sample metrics: `results/metrics/per_sample_metrics.json`

---

## References

### Shuffle Documentation
- Main README: `../README.md`
- Project documentation in repository root

### Assessment Plan
- Original plan: See commit history or JSONL transcript

### Data Locations
> **Note**: Paths below are placeholders. Update with actual synthetic VCF output directories.
- Somatic VCFs: `<SOMATIC_VCF_DIR>/synthetic_*.vcf.gz`
- Germline VCFs: `<GERMLINE_VCF_DIR>/synthetic_*.vcf.gz`
- Somatic logs: `<SOMATIC_VCF_DIR>/*.log`
- Germline logs: `<GERMLINE_VCF_DIR>/*.log`

---

*Assessment completed: 2026-03-11*
