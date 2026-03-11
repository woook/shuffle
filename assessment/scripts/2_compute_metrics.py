#!/usr/bin/env python3
"""
Phase 2: Per-Sample Quality Metrics

Computes comprehensive quality metrics for each selected sample:
- Basic statistics (variant count, missing rate, genotype distribution, het rate)
- Ti/Tv ratio
- Allele frequency distribution
- FORMAT field consistency (somatic only)
"""

import json
import subprocess
import sys
from contextlib import closing
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict, Counter
import re

try:
    from cyvcf2 import VCF
    import numpy as np
    import pandas as pd
except ImportError:
    print("ERROR: Required packages not installed")
    print("Please install: pip install cyvcf2 numpy pandas")
    sys.exit(1)

# Configuration
RESULTS_DIR = Path("/home/wook/Documents/shuffle/assessment/results")
SELECTION_FILE = RESULTS_DIR / "sample_selection.json"
METRICS_DIR = RESULTS_DIR / "metrics"
METRICS_DIR.mkdir(exist_ok=True)


# Transition and transversion definitions
TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
TRANSVERSIONS = {
    ('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'),
    ('G', 'C'), ('C', 'G'), ('G', 'T'), ('T', 'G')
}


def load_sample_selection() -> Dict:
    """Load sample selection from Phase 1."""
    with open(SELECTION_FILE) as f:
        return json.load(f)


def classify_mutation(ref: str, alt: str) -> str:
    """Classify mutation as transition or transversion."""
    if len(ref) != 1 or len(alt) != 1:
        return "other"  # Indel or complex
    pair = (ref.upper(), alt.upper())
    if pair in TRANSITIONS:
        return "transition"
    elif pair in TRANSVERSIONS:
        return "transversion"
    return "other"


def compute_basic_stats(vcf_path: Path, sample_name: str) -> Dict:
    """
    Compute basic statistics for a single sample.

    Returns:
        - Total variant count
        - Missing genotype rate
        - Genotype distribution (0/0, 0/1, 1/1)
        - Observed heterozygosity
        - Ti/Tv ratio
    """
    total_variants = 0
    missing_count = 0
    gt_counts = Counter()
    ti_count = 0
    tv_count = 0
    chr_variants = defaultdict(int)

    with closing(VCF(str(vcf_path))) as vcf:
        for variant in vcf:
            total_variants += 1
            chr_variants[variant.CHROM] += 1

            # Get genotype for first (only) sample
            gt = variant.genotypes[0]  # [allele1, allele2, phased]

            # Count genotype
            if gt[0] == -1 or gt[1] == -1:
                missing_count += 1
                gt_counts["missing"] += 1
            elif gt[0] == 0 and gt[1] == 0:
                gt_counts["0/0"] += 1
            elif gt[0] != gt[1]:
                gt_counts["0/1"] += 1
            elif gt[0] == gt[1] and gt[0] > 0:
                gt_counts["1/1"] += 1

            # Ti/Tv classification (only for SNPs)
            if variant.is_snp and len(variant.ALT) == 1:
                ref = variant.REF
                alt = variant.ALT[0]
                mut_type = classify_mutation(ref, alt)
                if mut_type == "transition":
                    ti_count += 1
                elif mut_type == "transversion":
                    tv_count += 1

    # Calculate rates
    missing_rate = missing_count / total_variants if total_variants > 0 else 0
    het_count = gt_counts["0/1"]
    hom_ref = gt_counts["0/0"]
    hom_alt = gt_counts["1/1"]
    called_sites = hom_ref + het_count + hom_alt
    het_rate = het_count / called_sites if called_sites > 0 else 0
    ti_tv_ratio = ti_count / tv_count if tv_count > 0 else 0

    return {
        "sample_name": sample_name,
        "total_variants": total_variants,
        "missing_count": missing_count,
        "missing_rate": missing_rate,
        "genotype_counts": dict(gt_counts),
        "heterozygosity_rate": het_rate,
        "ti_count": ti_count,
        "tv_count": tv_count,
        "ti_tv_ratio": ti_tv_ratio,
        "variants_per_chromosome": dict(chr_variants)
    }


def compute_allele_frequency_dist(vcf_path: Path) -> Dict:
    """
    Compute allele frequency distribution (MAF spectrum).

    Bins: 0-1%, 1-5%, 5-10%, 10-25%, 25-50%
    Note: For single-sample VCFs, AF is computed from GT
    """
    maf_bins = {
        "0-1%": 0,
        "1-5%": 0,
        "5-10%": 0,
        "10-25%": 0,
        "25-50%": 0
    }

    singleton_count = 0
    total_variants = 0

    with closing(VCF(str(vcf_path))) as vcf:
        for variant in vcf:
            total_variants += 1
            gt = variant.genotypes[0]

            # Calculate AF from genotype (for single sample)
            if gt[0] == -1 or gt[1] == -1:
                continue  # Skip missing

            allele_count = gt[0] + gt[1]  # 0, 1, or 2
            af = allele_count / 2.0  # AF for this sample

            # Convert to MAF (minor allele frequency)
            maf = min(af, 1 - af)

            # Bin the MAF
            if maf < 0.01:
                maf_bins["0-1%"] += 1
            elif maf < 0.05:
                maf_bins["1-5%"] += 1
            elif maf < 0.10:
                maf_bins["5-10%"] += 1
            elif maf < 0.25:
                maf_bins["10-25%"] += 1
            elif maf <= 0.50:
                maf_bins["25-50%"] += 1

            # Count singletons (heterozygous sites in single sample context)
            if allele_count == 1:
                singleton_count += 1

    return {
        "maf_distribution": maf_bins,
        "singleton_count": singleton_count,
        "total_variants": total_variants
    }


def compute_format_consistency(vcf_path: Path, sample_name: str) -> Dict:
    """
    Check FORMAT field consistency for somatic samples (GT:AF:DP:AD).

    Validates:
    - AF vs GT concordance
    - AD vs GT concordance
    - DP distribution
    - Missing value rates
    """
    # Storage for validation
    af_values = []
    dp_values = []
    missing_af = 0
    missing_dp = 0
    missing_ad = 0

    # Concordance checks
    af_gt_concordant = 0
    af_gt_discordant = 0
    ad_gt_concordant = 0
    ad_gt_discordant = 0

    total_variants = 0

    with closing(VCF(str(vcf_path))) as vcf:
        for variant in vcf:
            total_variants += 1

            # Get genotype
            gt = variant.genotypes[0]
            if gt[0] == -1 or gt[1] == -1:
                continue  # Skip missing genotypes

            gt_dosage = gt[0] + gt[1]  # 0, 1, or 2

            # Check AF
            try:
                af_field = variant.format('AF')
                if af_field is not None and len(af_field) > 0:
                    af_val = af_field[0][0]  # First sample, first alt allele
                    if af_val is not None and af_val >= 0:
                        af_values.append(af_val)

                        # Check AF vs GT concordance
                        if gt_dosage == 0 and af_val < 0.1:
                            af_gt_concordant += 1
                        elif gt_dosage == 1 and 0.3 <= af_val <= 0.7:
                            af_gt_concordant += 1
                        elif gt_dosage == 2 and af_val > 0.9:
                            af_gt_concordant += 1
                        else:
                            af_gt_discordant += 1
                    else:
                        missing_af += 1
            except (AttributeError, IndexError, KeyError, TypeError):
                missing_af += 1

            # Check DP
            try:
                dp_field = variant.format('DP')
                if dp_field is not None and len(dp_field) > 0:
                    dp_val = dp_field[0][0]  # First sample
                    if dp_val is not None and dp_val >= 0:
                        dp_values.append(dp_val)
                    else:
                        missing_dp += 1
            except (AttributeError, IndexError, KeyError, TypeError):
                missing_dp += 1

            # Check AD
            try:
                ad_field = variant.format('AD')
                if ad_field is not None and len(ad_field) > 0:
                    ad_vals = ad_field[0]  # First sample: [ref_depth, alt_depth]
                    if ad_vals is not None and len(ad_vals) >= 2:
                        ref_depth = ad_vals[0] if ad_vals[0] >= 0 else 0
                        alt_depth = ad_vals[1] if ad_vals[1] >= 0 else 0
                        total_depth = ref_depth + alt_depth

                        if total_depth > 0:
                            alt_ratio = alt_depth / total_depth

                            # Check AD vs GT concordance
                            if gt_dosage == 0 and alt_ratio < 0.2:
                                ad_gt_concordant += 1
                            elif gt_dosage == 1 and 0.3 <= alt_ratio <= 0.7:
                                ad_gt_concordant += 1
                            elif gt_dosage == 2 and alt_ratio > 0.8:
                                ad_gt_concordant += 1
                            else:
                                ad_gt_discordant += 1
                    else:
                        missing_ad += 1
            except (AttributeError, IndexError, KeyError, TypeError):
                missing_ad += 1

    # Calculate statistics
    af_missing_rate = missing_af / total_variants if total_variants > 0 else 0
    dp_missing_rate = missing_dp / total_variants if total_variants > 0 else 0
    ad_missing_rate = missing_ad / total_variants if total_variants > 0 else 0

    total_af_checked = af_gt_concordant + af_gt_discordant
    af_concordance_rate = af_gt_concordant / total_af_checked if total_af_checked > 0 else 0

    total_ad_checked = ad_gt_concordant + ad_gt_discordant
    ad_concordance_rate = ad_gt_concordant / total_ad_checked if total_ad_checked > 0 else 0

    return {
        "sample_name": sample_name,
        "format_fields_checked": ["AF", "DP", "AD"],
        "af_stats": {
            "values_count": len(af_values),
            "missing_count": missing_af,
            "missing_rate": af_missing_rate,
            "mean": float(np.mean(af_values)) if af_values else None,
            "median": float(np.median(af_values)) if af_values else None,
            "std": float(np.std(af_values)) if af_values else None
        },
        "dp_stats": {
            "values_count": len(dp_values),
            "missing_count": missing_dp,
            "missing_rate": dp_missing_rate,
            "mean": float(np.mean(dp_values)) if dp_values else None,
            "median": float(np.median(dp_values)) if dp_values else None,
            "q25": float(np.percentile(dp_values, 25)) if dp_values else None,
            "q75": float(np.percentile(dp_values, 75)) if dp_values else None
        },
        "ad_stats": {
            "missing_count": missing_ad,
            "missing_rate": ad_missing_rate
        },
        "concordance": {
            "af_gt_concordance_rate": af_concordance_rate,
            "af_gt_concordant_count": af_gt_concordant,
            "af_gt_discordant_count": af_gt_discordant,
            "ad_gt_concordance_rate": ad_concordance_rate,
            "ad_gt_concordant_count": ad_gt_concordant,
            "ad_gt_discordant_count": ad_gt_discordant
        }
    }


def process_sample(vcf_path: Path, sample_id: int, cohort_type: str) -> Dict:
    """Process a single sample and compute all metrics."""
    sample_name = f"{cohort_type}_synthetic_{sample_id}"

    print(f"  Processing: {sample_name}")

    # Basic statistics
    basic_stats = compute_basic_stats(vcf_path, sample_name)

    # Allele frequency distribution
    af_dist = compute_allele_frequency_dist(vcf_path)

    # Combine results
    results = {
        "sample_id": sample_id,
        "sample_name": sample_name,
        "cohort": cohort_type,
        "vcf_path": str(vcf_path),
        **basic_stats,
        "allele_frequency": af_dist
    }

    # FORMAT consistency (somatic only)
    if cohort_type == "somatic":
        format_consistency = compute_format_consistency(vcf_path, sample_name)
        results["format_consistency"] = format_consistency

    return results


def main():
    """Main processing logic."""
    print("=" * 80)
    print("PHASE 2: PER-SAMPLE QUALITY METRICS")
    print("=" * 80)

    # Load sample selection
    selection = load_sample_selection()

    all_metrics = []

    # Process somatic samples
    print("\nProcessing SOMATIC samples...")
    somatic_dir = Path(selection["somatic"]["vcf_directory"])
    for sample_id in selection["somatic"]["selected_ids"]:
        vcf_path = somatic_dir / f"synthetic_{sample_id}.vcf.gz"
        metrics = process_sample(vcf_path, sample_id, "somatic")
        all_metrics.append(metrics)

    # Process germline samples
    print("\nProcessing GERMLINE samples...")
    germline_dir = Path(selection["germline"]["vcf_directory"])
    for sample_id in selection["germline"]["selected_ids"]:
        vcf_path = germline_dir / f"synthetic_{sample_id}.vcf.gz"
        metrics = process_sample(vcf_path, sample_id, "germline")
        all_metrics.append(metrics)

    # Save detailed JSON results
    output_json = METRICS_DIR / "per_sample_metrics.json"
    with open(output_json, 'w') as f:
        json.dump(all_metrics, f, indent=2)

    print(f"\n✓ Detailed metrics saved to: {output_json}")

    # Create summary CSV
    create_summary_csv(all_metrics)

    print("\n" + "=" * 80)
    print("PHASE 2 COMPLETE")
    print("=" * 80)


def create_summary_csv(metrics_list: List[Dict]):
    """Create a summary CSV file for easy analysis."""
    rows = []

    for m in metrics_list:
        row = {
            "sample_name": m["sample_name"],
            "cohort": m["cohort"],
            "total_variants": m["total_variants"],
            "missing_rate": f"{m['missing_rate']:.4f}",
            "het_rate": f"{m['heterozygosity_rate']:.4f}",
            "ti_tv_ratio": f"{m['ti_tv_ratio']:.3f}",
            "hom_ref": m["genotype_counts"].get("0/0", 0),
            "het": m["genotype_counts"].get("0/1", 0),
            "hom_alt": m["genotype_counts"].get("1/1", 0),
            "singletons": m["allele_frequency"]["singleton_count"],
            "maf_0-1%": m["allele_frequency"]["maf_distribution"]["0-1%"],
            "maf_1-5%": m["allele_frequency"]["maf_distribution"]["1-5%"],
            "maf_5-10%": m["allele_frequency"]["maf_distribution"]["5-10%"],
            "maf_10-25%": m["allele_frequency"]["maf_distribution"]["10-25%"],
            "maf_25-50%": m["allele_frequency"]["maf_distribution"]["25-50%"]
        }

        # Add FORMAT consistency metrics for somatic
        if m["cohort"] == "somatic" and "format_consistency" in m:
            fc = m["format_consistency"]
            row["af_mean"] = f"{fc['af_stats']['mean']:.3f}" if fc['af_stats']['mean'] is not None else "N/A"
            row["dp_median"] = f"{fc['dp_stats']['median']:.0f}" if fc['dp_stats']['median'] is not None else "N/A"
            row["af_gt_concordance"] = f"{fc['concordance']['af_gt_concordance_rate']:.3f}"
            row["ad_gt_concordance"] = f"{fc['concordance']['ad_gt_concordance_rate']:.3f}"

        rows.append(row)

    df = pd.DataFrame(rows)

    # Save CSV
    output_csv = METRICS_DIR / "per_sample_summary.csv"
    df.to_csv(output_csv, index=False)

    print(f"✓ Summary CSV saved to: {output_csv}")

    # Print summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)

    for cohort in ["somatic", "germline"]:
        cohort_df = df[df["cohort"] == cohort]
        print(f"\n{cohort.upper()} Cohort (n={len(cohort_df)}):")
        print(f"  Variant count: {cohort_df['total_variants'].mean():.0f} ± {cohort_df['total_variants'].std():.0f}")
        print(f"  Heterozygosity: {cohort_df['het_rate'].astype(float).mean():.4f} ± {cohort_df['het_rate'].astype(float).std():.4f}")
        print(f"  Ti/Tv ratio: {cohort_df['ti_tv_ratio'].astype(float).mean():.3f} ± {cohort_df['ti_tv_ratio'].astype(float).std():.3f}")


if __name__ == "__main__":
    main()
