#!/usr/bin/env python3
"""
Phase 5: Analyze Shuffle Configuration

Extracts configuration parameters from shuffle log files to assess
anonymisation strength and theoretical re-identification risk.
"""

import re
import json
from pathlib import Path
from typing import Dict, List

# Configuration
RESULTS_DIR = Path("/home/wook/Documents/shuffle/assessment/results")
SOMATIC_LOG_DIR = Path("/home/wook/Downloads/v-s/h/s")
GERMLINE_LOG_DIR = Path("/home/wook/Downloads/v-s/w")


def parse_log_file(log_path: Path) -> Dict:
    """Parse a shuffle log file to extract configuration parameters."""
    if not log_path.exists():
        return {"error": "file not found"}

    with open(log_path, 'r') as f:
        content = f.read()

    config = {}

    # Extract region sampling mode
    if "--no-region-sampling" in content or "region_sampling: false" in content.lower():
        config["region_sampling"] = False
    else:
        config["region_sampling"] = True  # Default

    # Extract region count
    region_match = re.search(r'(\d+)\s+regions\s+detected', content, re.IGNORECASE)
    if region_match:
        config["region_count"] = int(region_match.group(1))

    # Extract min_donors setting
    min_donors_match = re.search(r'min[_-]donors[:\s=]+(\d+)', content, re.IGNORECASE)
    if min_donors_match:
        config["min_donors"] = int(min_donors_match.group(1))
    else:
        config["min_donors"] = 1  # Default

    # Extract region gap threshold
    gap_match = re.search(r'region[_-]gap[:\s=]+(\d+)', content, re.IGNORECASE)
    if gap_match:
        config["region_gap_bp"] = int(gap_match.group(1))

    # Extract donor pool size
    donor_match = re.search(r'(\d+)\s+donors?\s+in\s+pool', content, re.IGNORECASE)
    if donor_match:
        config["donor_pool_size"] = int(donor_match.group(1))

    # Check for sex filtering
    if "sex-file" in content.lower() or "female donors" in content.lower():
        config["sex_filtering"] = True
    else:
        config["sex_filtering"] = False

    # Check for FORMAT fields
    format_match = re.search(r'carry[_-]format[_-]fields[:\s=]+([A-Z,]+)', content, re.IGNORECASE)
    if format_match:
        config["format_fields"] = format_match.group(1).split(',')
    else:
        config["format_fields"] = []

    # Extract chromosome info
    chr_match = re.search(r'chromosome[:\s]+([0-9XY]+)', content, re.IGNORECASE)
    if chr_match:
        config["chromosome"] = chr_match.group(1)

    return config


def analyze_cohort_logs(log_dir: Path, cohort_name: str) -> Dict:
    """Analyze all log files in a directory."""
    print(f"\nAnalyzing {cohort_name} logs...")

    log_files = list(log_dir.glob("*.log"))

    if not log_files:
        print(f"  WARNING: No log files found in {log_dir}")
        return {"error": "no log files"}

    # Parse all logs
    configs = []
    for log_file in log_files:
        config = parse_log_file(log_file)
        if "error" not in config:
            configs.append(config)

    if not configs:
        print(f"  WARNING: Could not parse any log files")
        return {"error": "parsing failed"}

    # Aggregate configuration
    aggregated = {
        "cohort": cohort_name,
        "log_files_analyzed": len(configs),
        "region_sampling_enabled": any(c.get("region_sampling", False) for c in configs),
        "region_counts": [c.get("region_count") for c in configs if "region_count" in c],
        "min_donors_values": [c.get("min_donors") for c in configs if "min_donors" in c],
        "donor_pool_sizes": [c.get("donor_pool_size") for c in configs if "donor_pool_size" in c],
        "sex_filtering": any(c.get("sex_filtering", False) for c in configs),
        "format_fields": list(set(field for c in configs if "format_fields" in c for field in c["format_fields"])),
        "chromosomes": list(set(c.get("chromosome") for c in configs if "chromosome" in c))
    }

    # Summary statistics
    if aggregated["region_counts"]:
        aggregated["region_count_mean"] = sum(aggregated["region_counts"]) / len(aggregated["region_counts"])
        aggregated["region_count_min"] = min(aggregated["region_counts"])
        aggregated["region_count_max"] = max(aggregated["region_counts"])

    if aggregated["min_donors_values"]:
        aggregated["min_donors_mode"] = max(set(aggregated["min_donors_values"]),
                                           key=aggregated["min_donors_values"].count)

    if aggregated["donor_pool_sizes"]:
        aggregated["donor_pool_size_mode"] = max(set(aggregated["donor_pool_sizes"]),
                                                key=aggregated["donor_pool_sizes"].count)

    print(f"  ✓ Analyzed {len(configs)} log files")
    print(f"  Region sampling: {'ENABLED' if aggregated['region_sampling_enabled'] else 'DISABLED'}")
    if aggregated.get("region_count_mean"):
        print(f"  Mean regions per chromosome: {aggregated['region_count_mean']:.0f}")
    if aggregated.get("min_donors_mode"):
        print(f"  Min donors per synthetic: {aggregated['min_donors_mode']}")
    if aggregated.get("donor_pool_size_mode"):
        print(f"  Donor pool size: {aggregated['donor_pool_size_mode']}")

    return aggregated


def assess_theoretical_risk(config: Dict) -> Dict:
    """Assess theoretical re-identification risk based on configuration."""
    risk_assessment = {
        "overall_risk": "UNKNOWN",
        "risk_factors": [],
        "protective_factors": [],
        "recommendations": []
    }

    # Region sampling protection
    if config.get("region_sampling_enabled"):
        risk_assessment["protective_factors"].append(
            "Region-sampling mode enabled (reduces P2 attack rate from ~100% to ~5%)"
        )
    else:
        risk_assessment["risk_factors"].append(
            "CRITICAL: Continuous mode used (high P2 re-identification risk)"
        )

    # Min donors setting
    min_donors = config.get("min_donors_mode", 1)
    if min_donors >= 5:
        risk_assessment["protective_factors"].append(
            f"High min_donors setting ({min_donors}) increases donor diversity"
        )
    elif min_donors >= 3:
        risk_assessment["protective_factors"].append(
            f"Moderate min_donors setting ({min_donors})"
        )
    else:
        risk_assessment["risk_factors"].append(
            f"Low min_donors setting ({min_donors}) - allows primary donor dominance"
        )
        risk_assessment["recommendations"].append(
            "Consider increasing min_donors to ≥5 for future cohorts"
        )

    # Donor pool size
    donor_pool = config.get("donor_pool_size_mode")
    if donor_pool:
        if donor_pool >= 100:
            risk_assessment["protective_factors"].append(
                f"Large donor pool ({donor_pool}) reduces membership inference"
            )
        elif donor_pool >= 50:
            risk_assessment["protective_factors"].append(
                f"Moderate donor pool ({donor_pool})"
            )
        else:
            risk_assessment["risk_factors"].append(
                f"Small donor pool ({donor_pool}) increases membership inference risk"
            )

    # Region count
    region_count = config.get("region_count_mean")
    if region_count:
        if region_count >= 100:
            risk_assessment["protective_factors"].append(
                f"High region count ({region_count:.0f}) dilutes primary donor contribution"
            )
        elif region_count >= 50:
            risk_assessment["protective_factors"].append(
                f"Moderate region count ({region_count:.0f})"
            )
        else:
            risk_assessment["risk_factors"].append(
                f"Low region count ({region_count:.0f}) may allow donor tracking"
            )

    # Overall risk assessment
    if not config.get("region_sampling_enabled"):
        risk_assessment["overall_risk"] = "HIGH"
    elif len(risk_assessment["risk_factors"]) >= 2:
        risk_assessment["overall_risk"] = "MODERATE"
    else:
        risk_assessment["overall_risk"] = "LOW"

    return risk_assessment


def main():
    """Main configuration analysis."""
    print("=" * 80)
    print("PHASE 5: CONFIGURATION ANALYSIS")
    print("=" * 80)

    # Analyze somatic cohort
    somatic_config = analyze_cohort_logs(SOMATIC_LOG_DIR, "somatic")

    # Analyze germline cohort
    germline_config = analyze_cohort_logs(GERMLINE_LOG_DIR, "germline")

    # Assess risks
    print("\n" + "=" * 80)
    print("THEORETICAL RISK ASSESSMENT")
    print("=" * 80)

    somatic_risk = assess_theoretical_risk(somatic_config)
    germline_risk = assess_theoretical_risk(germline_config)

    results = {
        "somatic": {
            "configuration": somatic_config,
            "risk_assessment": somatic_risk
        },
        "germline": {
            "configuration": germline_config,
            "risk_assessment": germline_risk
        }
    }

    # Save results
    output_file = RESULTS_DIR / "configuration_analysis.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\n✓ Analysis saved to: {output_file}")

    # Print risk summaries
    for cohort_name in ["somatic", "germline"]:
        cohort_results = results[cohort_name]
        risk = cohort_results["risk_assessment"]

        print(f"\n{cohort_name.upper()} Cohort Risk: {risk['overall_risk']}")
        print(f"  Protective factors: {len(risk['protective_factors'])}")
        for factor in risk['protective_factors']:
            print(f"    ✓ {factor}")
        if risk['risk_factors']:
            print(f"  Risk factors: {len(risk['risk_factors'])}")
            for factor in risk['risk_factors']:
                print(f"    ✗ {factor}")

    print("\n" + "=" * 80)
    print("CONFIGURATION ANALYSIS COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
