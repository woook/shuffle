#!/usr/bin/env python3
"""
Phase 2.3: Hardy-Weinberg Equilibrium (HWE) Analysis

Computes HWE exact test p-values per variant and assesses violation rates.
This is done using plink2 for computational efficiency.
"""

import json
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np

# Configuration
RESULTS_DIR = Path("/home/wook/Documents/shuffle/assessment/results")
SELECTION_FILE = RESULTS_DIR / "sample_selection.json"
METRICS_DIR = RESULTS_DIR / "metrics"
TEMP_DIR = RESULTS_DIR / "temp_hwe"
TEMP_DIR.mkdir(exist_ok=True)


def load_sample_selection():
    """Load sample selection from Phase 1."""
    with open(SELECTION_FILE) as f:
        return json.load(f)


def merge_cohort_vcfs(vcf_paths: list, output_vcf: Path, cohort_name: str):
    """Merge multiple single-sample VCFs for HWE analysis."""
    print(f"  Merging {len(vcf_paths)} VCFs...")

    # Create file list
    file_list = TEMP_DIR / f"{cohort_name}_vcf_list.txt"
    with open(file_list, 'w') as f:
        for vcf_path in vcf_paths:
            f.write(f"{vcf_path}\n")

    # Merge VCFs using bcftools
    cmd = [
        "bcftools", "merge",
        "-l", str(file_list),
        "-m", "none",  # Don't merge variant records
        "-Oz",
        "-o", str(output_vcf)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR merging VCFs: {result.stderr}")
        return False

    # Index the merged VCF
    subprocess.run(["bcftools", "index", "-t", str(output_vcf)],
                  capture_output=True)

    print(f"  ✓ Merged VCF: {output_vcf}")
    return True


def run_plink_hwe(merged_vcf: Path, cohort_name: str) -> Path:
    """Run plink2 HWE test on merged VCF."""
    print(f"  Running plink2 --hardy...")

    output_prefix = TEMP_DIR / f"{cohort_name}_plink"

    # Convert VCF to plink format and run HWE
    cmd = [
        "plink2",
        "--vcf", str(merged_vcf),
        "--hardy",
        "--allow-extra-chr",
        "--out", str(output_prefix)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running plink2: {result.stderr}")
        return None

    hardy_file = Path(str(output_prefix) + ".hardy")
    if hardy_file.exists():
        print(f"  ✓ HWE results: {hardy_file}")
        return hardy_file
    else:
        print(f"  ERROR: Hardy file not created")
        return None


def analyze_hwe_results(hardy_file: Path, cohort_name: str) -> dict:
    """Analyze HWE test results."""
    print(f"  Analyzing HWE results...")

    # Read plink2 hardy output
    # Format: #CHROM  POS     ID      REF     ALT     A1      AX      HOM_A1_CT       HET_A1_CT       TWO_AX_CT       O(HET_A1)       E(HET_A1)       P
    df = pd.read_csv(hardy_file, sep='\t', comment='#')

    # Filter for valid p-values
    valid_df = df[df['P'].notna()]

    total_sites = len(valid_df)
    if total_sites == 0:
        return {
            "cohort": cohort_name,
            "total_sites": 0,
            "error": "No valid HWE results"
        }

    # Count violations at different thresholds
    violations_001 = (valid_df['P'] < 0.001).sum()
    violations_0001 = (valid_df['P'] < 0.0001).sum()
    violations_1e5 = (valid_df['P'] < 1e-5).sum()

    violation_rate_001 = violations_001 / total_sites
    violation_rate_0001 = violations_0001 / total_sites
    violation_rate_1e5 = violations_1e5 / total_sites

    # Statistics on p-values
    p_mean = valid_df['P'].mean()
    p_median = valid_df['P'].median()
    p_min = valid_df['P'].min()

    results = {
        "cohort": cohort_name,
        "total_sites": total_sites,
        "violations_p001": int(violations_001),
        "violations_p0001": int(violations_0001),
        "violations_p1e5": int(violations_1e5),
        "violation_rate_p001": float(violation_rate_001),
        "violation_rate_p0001": float(violation_rate_0001),
        "violation_rate_p1e5": float(violation_rate_1e5),
        "p_value_mean": float(p_mean),
        "p_value_median": float(p_median),
        "p_value_min": float(p_min)
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

    # Merge VCFs
    merged_vcf = TEMP_DIR / f"{cohort_type}_merged.vcf.gz"
    if not merge_cohort_vcfs(vcf_paths, merged_vcf, cohort_type):
        return {"error": "merge failed"}

    # Run HWE test
    hardy_file = run_plink_hwe(merged_vcf, cohort_type)
    if hardy_file is None:
        return {"error": "plink failed"}

    # Analyze results
    results = analyze_hwe_results(hardy_file, cohort_type)

    return results


def main():
    """Main HWE analysis logic."""
    print("=" * 80)
    print("PHASE 2.3: HARDY-WEINBERG EQUILIBRIUM ANALYSIS")
    print("=" * 80)

    # Check if plink2 is available
    result = subprocess.run(["which", "plink2"], capture_output=True)
    if result.returncode != 0:
        print("\nERROR: plink2 not found in PATH")
        print("Please install plink2: https://www.cog-genomics.org/plink/2.0/")
        return

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

        # Interpretation
        if r['violation_rate_p001'] < 0.03:
            status = "PASS - Within expected range"
        elif r['violation_rate_p001'] < 0.10:
            status = "CAUTION - Slightly elevated"
        else:
            status = "FAIL - Significantly elevated"

        print(f"  Assessment: {status}")


if __name__ == "__main__":
    main()
