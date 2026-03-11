#!/usr/bin/env python3
"""
Phase 2/3: Create Visualizations for Assessment

Generates plots for:
- Heterozygosity rates
- Ti/Tv ratios
- MAF distributions
- FORMAT field concordance
- HWE violation rates
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Configuration
RESULTS_DIR = Path("/home/wook/Documents/shuffle/assessment/results")
METRICS_DIR = RESULTS_DIR / "metrics"
PLOTS_DIR = RESULTS_DIR / "plots"
PLOTS_DIR.mkdir(exist_ok=True)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 150


def load_data():
    """Load all metrics data."""
    with open(METRICS_DIR / "per_sample_metrics.json") as f:
        metrics = json.load(f)

    with open(METRICS_DIR / "hwe_analysis.json") as f:
        hwe = json.load(f)

    return metrics, hwe


def plot_heterozygosity_comparison(metrics: list):
    """Plot heterozygosity rates comparison between cohorts."""
    somatic = [m for m in metrics if m["cohort"] == "somatic"]
    germline = [m for m in metrics if m["cohort"] == "germline"]

    somatic_het = [m["heterozygosity_rate"] for m in somatic]
    germline_het = [m["heterozygosity_rate"] for m in germline]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Histogram
    ax1.hist(somatic_het, bins=15, alpha=0.6, label='Somatic (n=17)', color='#e74c3c')
    ax1.hist(germline_het, bins=15, alpha=0.6, label='Germline (n=14)', color='#3498db')
    ax1.set_xlabel('Heterozygosity Rate')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of Heterozygosity Rates')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Box plot
    data = pd.DataFrame({
        'Heterozygosity': somatic_het + germline_het,
        'Cohort': ['Somatic'] * len(somatic_het) + ['Germline'] * len(germline_het)
    })
    sns.boxplot(data=data, x='Cohort', y='Heterozygosity', ax=ax2, palette=['#e74c3c', '#3498db'])
    ax2.set_title('Heterozygosity by Cohort')
    ax2.set_ylabel('Heterozygosity Rate')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "heterozygosity_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Created: heterozygosity_comparison.png")


def plot_titv_ratios(metrics: list):
    """Plot Ti/Tv ratios."""
    df = pd.DataFrame([
        {
            'sample': m['sample_name'],
            'cohort': m['cohort'],
            'ti_tv': m['ti_tv_ratio']
        }
        for m in metrics
    ])

    fig, ax = plt.subplots(figsize=(10, 6))

    # Bar plot
    somatic_df = df[df['cohort'] == 'somatic']
    germline_df = df[df['cohort'] == 'germline']

    x_pos_somatic = np.arange(len(somatic_df))
    x_pos_germline = np.arange(len(germline_df)) + len(somatic_df) + 1

    ax.bar(x_pos_somatic, somatic_df['ti_tv'], color='#e74c3c', alpha=0.7, label='Somatic')
    ax.bar(x_pos_germline, germline_df['ti_tv'], color='#3498db', alpha=0.7, label='Germline')

    # Expected ranges
    ax.axhline(y=2.0, color='green', linestyle='--', linewidth=2, alpha=0.5, label='Expected WES (2.0-2.1)')
    ax.axhline(y=2.1, color='green', linestyle='--', linewidth=2, alpha=0.5)

    ax.set_xlabel('Sample')
    ax.set_ylabel('Ti/Tv Ratio')
    ax.set_title('Ti/Tv Ratios Across Samples')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "titv_ratios.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Created: titv_ratios.png")


def plot_maf_distributions(metrics: list):
    """Plot MAF distributions."""
    maf_bins = ["0-1%", "1-5%", "5-10%", "10-25%", "25-50%"]

    # Aggregate by cohort
    somatic_maf = np.zeros(5)
    germline_maf = np.zeros(5)

    for m in metrics:
        maf_dist = m["allele_frequency"]["maf_distribution"]
        for i, bin_name in enumerate(maf_bins):
            if m["cohort"] == "somatic":
                somatic_maf[i] += maf_dist[bin_name]
            else:
                germline_maf[i] += maf_dist[bin_name]

    # Normalize to fractions
    somatic_maf = somatic_maf / somatic_maf.sum()
    germline_maf = germline_maf / germline_maf.sum()

    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(maf_bins))
    width = 0.35

    ax.bar(x - width/2, somatic_maf, width, label='Somatic', color='#e74c3c', alpha=0.7)
    ax.bar(x + width/2, germline_maf, width, label='Germline', color='#3498db', alpha=0.7)

    ax.set_xlabel('MAF Bin')
    ax.set_ylabel('Fraction of Variants')
    ax.set_title('Minor Allele Frequency Distribution')
    ax.set_xticks(x)
    ax.set_xticklabels(maf_bins)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "maf_distribution.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Created: maf_distribution.png")


def plot_format_concordance(metrics: list):
    """Plot FORMAT field concordance for somatic cohort."""
    somatic = [m for m in metrics if m["cohort"] == "somatic" and "format_consistency" in m]

    af_concordance = [m["format_consistency"]["concordance"]["af_gt_concordance_rate"] for m in somatic]
    ad_concordance = [m["format_consistency"]["concordance"]["ad_gt_concordance_rate"] for m in somatic]

    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(somatic))
    width = 0.35

    ax.bar(x - width/2, af_concordance, width, label='AF-GT Concordance', color='#9b59b6', alpha=0.7)
    ax.bar(x + width/2, ad_concordance, width, label='AD-GT Concordance', color='#e67e22', alpha=0.7)

    ax.axhline(y=0.7, color='green', linestyle='--', alpha=0.5, label='70% threshold')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Concordance Rate')
    ax.set_title('FORMAT Field Concordance (Somatic Cohort)')
    ax.set_ylim(0, 1)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "format_concordance.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Created: format_concordance.png")


def plot_hwe_summary(hwe: dict):
    """Plot HWE violation rates."""
    cohorts = ['Somatic', 'Germline']
    rates_001 = [hwe['somatic']['violation_rate_p001'], hwe['germline']['violation_rate_p001']]
    rates_0001 = [hwe['somatic']['violation_rate_p0001'], hwe['germline']['violation_rate_p0001']]

    fig, ax = plt.subplots(figsize=(8, 6))

    x = np.arange(len(cohorts))
    width = 0.35

    ax.bar(x - width/2, rates_001, width, label='p < 0.001', color='#e74c3c', alpha=0.7)
    ax.bar(x + width/2, rates_0001, width, label='p < 0.0001', color='#c0392b', alpha=0.7)

    ax.axhline(y=0.03, color='orange', linestyle='--', linewidth=2, alpha=0.7, label='Expected (<3%)')
    ax.axhline(y=0.10, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Elevated (>10%)')

    ax.set_xlabel('Cohort')
    ax.set_ylabel('HWE Violation Rate')
    ax.set_title('Hardy-Weinberg Equilibrium Violation Rates')
    ax.set_xticks(x)
    ax.set_xticklabels(cohorts)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "hwe_violations.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Created: hwe_violations.png")


def plot_variant_counts(metrics: list):
    """Plot variant counts per chromosome."""
    # Get chromosome variant counts from first sample of each cohort
    somatic_sample = [m for m in metrics if m["cohort"] == "somatic"][0]
    germline_sample = [m for m in metrics if m["cohort"] == "germline"][0]

    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
              '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']

    somatic_counts = [somatic_sample["variants_per_chromosome"].get(c, 0) for c in chroms]
    germline_counts = [germline_sample["variants_per_chromosome"].get(c, 0) for c in chroms]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Somatic
    ax1.bar(range(len(chroms)), somatic_counts, color='#e74c3c', alpha=0.7)
    ax1.set_xlabel('Chromosome')
    ax1.set_ylabel('Variant Count')
    ax1.set_title('Variant Counts per Chromosome (Somatic)')
    ax1.set_xticks(range(len(chroms)))
    ax1.set_xticklabels(chroms)
    ax1.grid(True, alpha=0.3, axis='y')

    # Germline
    ax2.bar(range(len(chroms)), germline_counts, color='#3498db', alpha=0.7)
    ax2.set_xlabel('Chromosome')
    ax2.set_ylabel('Variant Count')
    ax2.set_title('Variant Counts per Chromosome (Germline)')
    ax2.set_xticks(range(len(chroms)))
    ax2.set_xticklabels(chroms)
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "variant_counts_per_chromosome.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Created: variant_counts_per_chromosome.png")


def main():
    """Main visualization creation."""
    print("=" * 80)
    print("CREATING VISUALIZATIONS")
    print("=" * 80)

    # Load data
    print("\nLoading data...")
    metrics, hwe = load_data()

    # Create plots
    print("\nGenerating plots...")
    plot_heterozygosity_comparison(metrics)
    plot_titv_ratios(metrics)
    plot_maf_distributions(metrics)
    plot_format_concordance(metrics)
    plot_hwe_summary(hwe)
    plot_variant_counts(metrics)

    print("\n" + "=" * 80)
    print("VISUALIZATIONS COMPLETE")
    print("=" * 80)
    print(f"\nAll plots saved to: {PLOTS_DIR}")


if __name__ == "__main__":
    main()
