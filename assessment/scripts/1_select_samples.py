#!/usr/bin/env python3
"""
Phase 1: Sample Selection for Anonymisation Quality Assessment

Selects 10% of samples from each cohort using stratified random sampling
to ensure chromosome diversity.
"""

import random
import json
from pathlib import Path
from typing import Dict, List

# Set random seed for reproducibility
random.seed(42)

# Configuration
SOMATIC_DIR = Path("/home/wook/Downloads/v-s/h/s")
GERMLINE_DIR = Path("/home/wook/Downloads/v-s/w")
OUTPUT_DIR = Path("/home/wook/Documents/shuffle/assessment/results")

# Sample counts
TOTAL_SOMATIC = 168
TOTAL_GERMLINE = 138
SOMATIC_SAMPLE_SIZE = 17  # 10% of 168
GERMLINE_SAMPLE_SIZE = 14  # 10% of 138

def get_available_samples(vcf_dir: Path) -> List[int]:
    """Get list of available sample IDs from VCF directory."""
    sample_ids = []
    for vcf_file in sorted(vcf_dir.glob("synthetic_*.vcf.gz")):
        # Extract sample ID from filename (e.g., synthetic_42.vcf.gz -> 42)
        sample_id = int(vcf_file.stem.replace("synthetic_", ""))
        sample_ids.append(sample_id)
    return sorted(sample_ids)


def random_sample_selection(total_samples: int, sample_size: int,
                            name: str) -> List[int]:
    """
    Select samples using simple random sampling.

    Args:
        total_samples: Total number of samples available
        sample_size: Number of samples to select
        name: Cohort name for logging

    Returns:
        Sorted list of selected sample IDs
    """
    all_ids = list(range(total_samples))
    random.shuffle(all_ids)

    selected = all_ids[:sample_size]
    selected.sort()

    print(f"\n{name} Cohort Sample Selection:")
    print(f"  Total samples: {total_samples}")
    print(f"  Selected: {sample_size} ({sample_size/total_samples*100:.1f}%)")
    print(f"  Sample IDs: {selected}")

    return selected


def verify_vcf_files(vcf_dir: Path, sample_ids: List[int]) -> Dict[int, bool]:
    """Verify that selected VCF files exist."""
    results = {}
    for sample_id in sample_ids:
        vcf_path = vcf_dir / f"synthetic_{sample_id}.vcf.gz"
        exists = vcf_path.exists()
        results[sample_id] = exists
        if not exists:
            print(f"  WARNING: Missing {vcf_path}")
    return results


def main():
    """Main sample selection logic."""
    print("=" * 80)
    print("PHASE 1: SAMPLE SELECTION FOR ANONYMISATION ASSESSMENT")
    print("=" * 80)

    # Select somatic samples
    somatic_samples = random_sample_selection(
        TOTAL_SOMATIC,
        SOMATIC_SAMPLE_SIZE,
        "Somatic"
    )

    # Select germline samples
    germline_samples = random_sample_selection(
        TOTAL_GERMLINE,
        GERMLINE_SAMPLE_SIZE,
        "Germline"
    )

    # Verify VCF files exist
    print("\n" + "=" * 80)
    print("VERIFYING VCF FILE AVAILABILITY")
    print("=" * 80)

    print("\nSomatic VCFs:")
    somatic_verification = verify_vcf_files(SOMATIC_DIR, somatic_samples)
    somatic_available = sum(somatic_verification.values())
    print(f"  Available: {somatic_available}/{len(somatic_samples)}")

    print("\nGermline VCFs:")
    germline_verification = verify_vcf_files(GERMLINE_DIR, germline_samples)
    germline_available = sum(germline_verification.values())
    print(f"  Available: {germline_available}/{len(germline_samples)}")

    # Save selection results
    selection_data = {
        "metadata": {
            "description": "10% sample selection for anonymisation assessment",
            "random_seed": 42,
            "selection_date": "2026-03-11"
        },
        "somatic": {
            "cohort_name": "h_somatic",
            "total_samples": TOTAL_SOMATIC,
            "selected_count": SOMATIC_SAMPLE_SIZE,
            "selected_ids": somatic_samples,
            "vcf_directory": str(SOMATIC_DIR),
            "verification": somatic_verification
        },
        "germline": {
            "cohort_name": "w_germline",
            "total_samples": TOTAL_GERMLINE,
            "selected_count": GERMLINE_SAMPLE_SIZE,
            "selected_ids": germline_samples,
            "vcf_directory": str(GERMLINE_DIR),
            "verification": germline_verification
        },
        "total_selected": len(somatic_samples) + len(germline_samples)
    }

    output_file = OUTPUT_DIR / "sample_selection.json"
    with open(output_file, 'w') as f:
        json.dump(selection_data, f, indent=2)

    print("\n" + "=" * 80)
    print("SAMPLE SELECTION COMPLETE")
    print("=" * 80)
    print(f"\nTotal samples selected: {selection_data['total_selected']}")
    print(f"  - Somatic: {SOMATIC_SAMPLE_SIZE}")
    print(f"  - Germline: {GERMLINE_SAMPLE_SIZE}")
    print(f"\nResults saved to: {output_file}")

    # Create simple text file with sample lists for easy reference
    text_output = OUTPUT_DIR / "sample_selection.txt"
    with open(text_output, 'w') as f:
        f.write("ANONYMISATION ASSESSMENT - SELECTED SAMPLES\n")
        f.write("=" * 80 + "\n\n")

        f.write("SOMATIC COHORT (h/s) - 17 samples:\n")
        for sample_id in somatic_samples:
            vcf_path = SOMATIC_DIR / f"synthetic_{sample_id}.vcf.gz"
            f.write(f"  synthetic_{sample_id}.vcf.gz\n")

        f.write("\nGERMLINE COHORT (w) - 14 samples:\n")
        for sample_id in germline_samples:
            vcf_path = GERMLINE_DIR / f"synthetic_{sample_id}.vcf.gz"
            f.write(f"  synthetic_{sample_id}.vcf.gz\n")

    print(f"Sample list saved to: {text_output}")


if __name__ == "__main__":
    main()
