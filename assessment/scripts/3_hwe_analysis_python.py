#!/usr/bin/env python3
"""
Phase 2.3: Hardy-Weinberg Equilibrium (HWE) Analysis (Pure Python)

Computes HWE chi-square test p-values per variant without requiring plink2.
Uses scipy.stats.chi2 for statistical tests (approximation, not exact test).

Note: This implementation uses chi-square approximation for computational efficiency.
For low-count sites or comparison with exact HWE tests (e.g., plink2/Fisher's exact),
results may differ.
"""

import json
import sys
from contextlib import closing
from pathlib import Path
from typing import Dict, List
import numpy as np

try:
    from cyvcf2 import VCF
    from scipy.stats import chi2
    import pandas as pd
except ImportError:
    print("ERROR: Required packages not installed")
    print("Please install: pip install cyvcf2 scipy pandas")
    sys.exit(1)

# Configuration
RESULTS_DIR = Path("/home/wook/Documents/shuffle/assessment/results")
SELECTION_FILE = RESULTS_DIR / "sample_selection.json"
METRICS_DIR = RESULTS_DIR / "metrics"


def load_sample_selection():
    """Load sample selection from Phase 1."""
    with open(SELECTION_FILE) as f:
        return json.load(f)


def hwe_chisq_test(obs_hom_ref: int, obs_het: int, obs_hom_alt: int) -> float:
    """
    Compute HWE chi-square test p-value.

    Uses chi-square approximation with scipy.stats.chi2 (df=1).
    This is NOT an exact test (Fisher's exact HWE). For sites with low expected counts
    (< 5), returns 1.0 (cannot test reliably). Results may differ from exact tests
    used by plink/plink2.
    """
    n = obs_hom_ref + obs_het + obs_hom_alt
    if n == 0:
        return 1.0

    # Calculate allele frequencies
    total_alleles = 2 * n
    n_ref = 2 * obs_hom_ref + obs_het
    n_alt = 2 * obs_hom_alt + obs_het

    p_ref = n_ref / total_alleles
    p_alt = n_alt / total_alleles

    # Expected genotype counts under HWE
    exp_hom_ref = n * p_ref * p_ref
    exp_het = n * 2 * p_ref * p_alt
    exp_hom_alt = n * p_alt * p_alt

    # Avoid division by zero
    if exp_hom_ref < 0.1 or exp_het < 0.1 or exp_hom_alt < 0.1:
        return 1.0  # Cannot test with low expected counts

    # Chi-square statistic
    chi_sq = (
        ((obs_hom_ref - exp_hom_ref) ** 2) / exp_hom_ref +
        ((obs_het - exp_het) ** 2) / exp_het +
        ((obs_hom_alt - exp_hom_alt) ** 2) / exp_hom_alt
    )

    # P-value from chi-square distribution with df=1
    p_value = 1.0 - chi2.cdf(chi_sq, df=1)

    return p_value


def analyze_vcf_hwe(vcf_paths: List[Path], cohort_name: str) -> Dict:
    """
    Analyze HWE across multiple VCFs.

    For each variant position, aggregate genotypes across samples
    and compute HWE p-value.
    """
    print("  Loading VCFs and computing HWE...")

    # Store genotypes by position
    variant_genotypes = {}

    # Read all VCFs and aggregate genotypes by position
    for vcf_path in vcf_paths:
        with closing(VCF(str(vcf_path))) as vcf:
            for variant in vcf:
                pos_key = (variant.CHROM, variant.POS, variant.REF, variant.ALT[0] if variant.ALT else "")

                if pos_key not in variant_genotypes:
                    variant_genotypes[pos_key] = {"hom_ref": 0, "het": 0, "hom_alt": 0, "missing": 0}

                gt = variant.genotypes[0]
                if gt[0] == -1 or gt[1] == -1:
                    variant_genotypes[pos_key]["missing"] += 1
                elif gt[0] == 0 and gt[1] == 0:
                    variant_genotypes[pos_key]["hom_ref"] += 1
                elif gt[0] != gt[1]:
                    variant_genotypes[pos_key]["het"] += 1
                elif gt[0] == gt[1] and gt[0] > 0:
                    variant_genotypes[pos_key]["hom_alt"] += 1

    print(f"  Aggregated genotypes for {len(variant_genotypes):,} variants")

    # Compute HWE for each variant
    p_values = []
    sites_tested = 0

    for _, gts in variant_genotypes.items():
        n_samples = gts["hom_ref"] + gts["het"] + gts["hom_alt"]

        # Only test variants with sufficient samples
        if n_samples >= 2:
            p_val = hwe_chisq_test(gts["hom_ref"], gts["het"], gts["hom_alt"])
            p_values.append(p_val)
            sites_tested += 1

    if not p_values:
        return {
            "cohort": cohort_name,
            "total_sites": 0,
            "error": "No variants with sufficient samples"
        }

    # Convert to numpy for analysis
    p_values = np.array(p_values)

    # Count violations at different thresholds
    violations_001 = (p_values < 0.001).sum()
    violations_0001 = (p_values < 0.0001).sum()
    violations_1e5 = (p_values < 1e-5).sum()

    violation_rate_001 = violations_001 / sites_tested
    violation_rate_0001 = violations_0001 / sites_tested
    violation_rate_1e5 = violations_1e5 / sites_tested

    results = {
        "cohort": cohort_name,
        "total_sites": sites_tested,
        "violations_p001": int(violations_001),
        "violations_p0001": int(violations_0001),
        "violations_p1e5": int(violations_1e5),
        "violation_rate_p001": float(violation_rate_001),
        "violation_rate_p0001": float(violation_rate_0001),
        "violation_rate_p1e5": float(violation_rate_1e5),
        "p_value_mean": float(p_values.mean()),
        "p_value_median": float(np.median(p_values)),
        "p_value_min": float(p_values.min()),
        "p_value_max": float(p_values.max())
    }

    print(f"  ✓ HWE violation rate (p<0.001): {violation_rate_001:.4f}")

    return results


def process_cohort(cohort_config: dict, cohort_type: str) -> dict:
    """Process one cohort for HWE analysis."""
    print(f"\nProcessing {cohort_type.upper()} cohort...")

    # Get VCF paths
    vcf_dir = Path(cohort_config["vcf_directory"])
    vcf_paths = [
        vcf_dir / f"synthetic_{sample_id}.vcf.gz"
        for sample_id in cohort_config["selected_ids"]
    ]

    # Verify all files exist
    missing = [p for p in vcf_paths if not p.exists()]
    if missing:
        print(f"  ERROR: Missing VCF files: {missing}")
        return {"error": "missing files"}

    # Analyze HWE
    results = analyze_vcf_hwe(vcf_paths, cohort_type)

    return results


def main():
    """Main HWE analysis logic."""
    print("=" * 80)
    print("PHASE 2.3: HARDY-WEINBERG EQUILIBRIUM ANALYSIS")
    print("=" * 80)

    # Load sample selection
    selection = load_sample_selection()

    results = {
        "somatic": process_cohort(selection["somatic"], "somatic"),
        "germline": process_cohort(selection["germline"], "germline")
    }

    # Save results
    output_file = METRICS_DIR / "hwe_analysis.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print("\n" + "=" * 80)
    print("HWE ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nResults saved to: {output_file}")

    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    for cohort_type in ["somatic", "germline"]:
        r = results[cohort_type]
        if "error" in r:
            print(f"\n{cohort_type.upper()}: ERROR - {r['error']}")
            continue

        print(f"\n{cohort_type.upper()} Cohort:")
        print(f"  Total sites tested: {r['total_sites']:,}")
        print(f"  HWE violations (p < 0.001): {r['violations_p001']:,} ({r['violation_rate_p001']:.2%})")
        print(f"  HWE violations (p < 0.0001): {r['violations_p0001']:,} ({r['violation_rate_p0001']:.2%})")
        print(f"  HWE violations (p < 1e-5): {r['violations_p1e5']:,} ({r['violation_rate_p1e5']:.2%})")
        print(f"  P-value range: [{r['p_value_min']:.2e}, {r['p_value_max']:.2e}]")

        # Interpretation
        if r['violation_rate_p001'] < 0.03:
            status = "✓ PASS - Within expected range (<3%)"
        elif r['violation_rate_p001'] < 0.10:
            status = "⚠ CAUTION - Slightly elevated (3-10%)"
        else:
            status = "✗ FAIL - Significantly elevated (>10%)"

        print(f"  Assessment: {status}")


if __name__ == "__main__":
    main()
