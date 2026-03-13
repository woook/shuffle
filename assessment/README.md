# Shuffle Anonymisation Quality Assessment

## Summary

This directory contains a comprehensive anonymisation quality assessment of synthetic VCF files generated using the shuffle tool.

### **Assessment** 

- **Biological Plausibility**: HWE compliance, normal het rates, expected Ti/Tv
- **Technical Quality**: FORMAT field consistency maintained
- **Theoretical Identification Risk**: Effectiveness of region-sampling mode with no access to original genotypes

> **Scope**: This assessment applies to the threat model and validation surface described below. See [Limitations](#limitations) for constraints on this assessment.

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

## Methodology

### Phase 1: Sample Selection
- **Method**: Use all samples or use simple random sampling (seed=42)
- **Sample size**: As defined
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
  - Check metric definitions in Appendix: Metric Definitions in `results/FINAL_RISK_ASSESSMENT.md`
  - Examine individual sample metrics: `results/metrics/per_sample_metrics.json`

---

## References

### Shuffle Documentation
- Main README: `../README.md`
- Project documentation in repository root

### Data Locations
> **Note**: Paths below are placeholders. Update with actual synthetic VCF output directories.
- Somatic VCFs: `<SOMATIC_VCF_DIR>/synthetic_*.vcf.gz`
- Germline VCFs: `<GERMLINE_VCF_DIR>/synthetic_*.vcf.gz`
- Somatic logs: `<SOMATIC_VCF_DIR>/*.log`
- Germline logs: `<GERMLINE_VCF_DIR>/*.log`
