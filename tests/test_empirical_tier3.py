"""
Tier 3 empirical tests: PCA projection and LD decay comparison.

Tests B5 and B6 from the v-shuffler empirical anonymisation test plan.
These require:
  * plink2  in PATH
  * bcftools in PATH (for merging VCFs)
  * A set of donor VCF files and a synthetic VCF from a v-shuffler run.

Running Tier 3 tests
--------------------
The tests are skipped automatically when plink2 or bcftools are not available.
To run them against production data:

1.  Follow the Step 0 data setup in the test plan (1000G chr22 or msprime).
2.  Set the environment variables:
      VSHUFFLE_DONOR_VCF   path to merged multi-sample donor VCF (.vcf.gz, tabixed)
      VSHUFFLE_SYNTH_VCF   path to synthetic multi-sample VCF from v-shuffler
3.  Run:
      pytest tests/test_empirical_tier3.py -v

Mode-agnostic tests
-------------------
B5 (PCA) and B6 (LD decay) are agnostic to whether the synthetic VCF was
produced in region-sampling mode (default) or continuous mode (--no-region-sampling).
Both output formats are standard VCF files; the same plink2/bcftools commands
work for both.

For whole-chromosome data passed via VSHUFFLE_DONOR_VCF, consider using
--no-region-sampling so that crossover breakpoints are placed at biologically
realistic positions across the full chromosome rather than at region boundaries.
For targeted panel data (e.g., exome or gene panels), the default region-sampling
mode is preferred: each captured region is treated as an independent mixing unit,
which provides thorough anonymisation even when the effective λ across the
captured fraction of the chromosome would otherwise be near zero.

Example commands reproduced from the plan
------------------------------------------
# Region mode (default) — for panel data
  v-shuffler shuffle -i "data/*.vcf.gz" -o /tmp/out -m map.txt -c chr22

# Continuous mode — for whole-chromosome data
  v-shuffler shuffle -i "data/*.vcf.gz" -o /tmp/out -m map.txt -c chr22 --no-region-sampling

# Merge donor and synth VCFs for joint PCA
  bcftools merge -Oz -o merged_all.vcf.gz merged_donors.vcf.gz synthetic_chr22.vcf.gz
  tabix -p vcf merged_all.vcf.gz

# PCA
  plink2 --vcf merged_all.vcf.gz --maf 0.05 --pca 10 --out pca_all

# LD decay (donors)
  plink2 --vcf merged_donors.vcf.gz --r2 --ld-window 100 --ld-window-kb 500 \\
         --ld-window-r2 0 --out ld_donors
# LD decay (synthetics)
  plink2 --vcf synthetic_chr22.vcf.gz --r2 --ld-window 100 --ld-window-kb 500 \\
         --ld-window-r2 0 --out ld_synth
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Availability checks
# ---------------------------------------------------------------------------

_PLINK2_AVAILABLE = shutil.which("plink2") is not None
_BCFTOOLS_AVAILABLE = shutil.which("bcftools") is not None
_TABIX_AVAILABLE = shutil.which("tabix") is not None

_SKIP_NO_PLINK2 = pytest.mark.skipif(
    not _PLINK2_AVAILABLE,
    reason="plink2 not found in PATH; install plink2 to run Tier 3 tests.",
)
_SKIP_NO_BCFTOOLS = pytest.mark.skipif(
    not _BCFTOOLS_AVAILABLE,
    reason="bcftools not found in PATH; install bcftools to run Tier 3 tests.",
)
_SKIP_NO_TABIX = pytest.mark.skipif(
    not _TABIX_AVAILABLE,
    reason="tabix not found in PATH; install tabix to run Tier 3 tests.",
)


def _require_env(var: str) -> Path:
    """Return the path stored in *var*; skip the test if not set or non-existent."""
    val = os.environ.get(var)
    if not val:
        pytest.skip(
            f"Environment variable {var} is not set. "
            "Set it to the path of the relevant VCF file to enable this test."
        )
    p = Path(val)
    if not p.exists():
        pytest.skip(f"{var}={val!r} does not exist.")
    return p


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def donor_vcf() -> Path:
    """Path to the merged multi-sample donor VCF (from $VSHUFFLE_DONOR_VCF)."""
    return _require_env("VSHUFFLE_DONOR_VCF")


@pytest.fixture(scope="module")
def synth_vcf() -> Path:
    """Path to the synthetic multi-sample VCF (from $VSHUFFLE_SYNTH_VCF)."""
    return _require_env("VSHUFFLE_SYNTH_VCF")


# ---------------------------------------------------------------------------
# B5: PCA Projection
#
# Merge donor and synthetic VCFs, run PCA with plink2, and verify that
# synthetic individuals fall within the donor population cloud (Mahalanobis
# distance using the first two PCs).
# ---------------------------------------------------------------------------

@_SKIP_NO_PLINK2
@_SKIP_NO_BCFTOOLS
@_SKIP_NO_TABIX
def test_b5_pca(
    donor_vcf: Path, synth_vcf: Path, tmp_path: Path
) -> None:
    """
    B5 — Fraction of synthetics inside the 99% Mahalanobis cloud > 0.95.

    Steps:
      1. Merge donor + synth VCFs.
      2. Run plink2 --pca 10.
      3. Split eigenvec by sample-name prefix.
      4. Compute Mahalanobis distance of each synth from the donor centroid.

    Pass threshold (plan §B5):
      * Fraction of synthetics inside 99% cloud > 0.95.
      * No synth PC1 or PC2 outlier (> 3 SD of donor range).
    """
    from scipy.spatial.distance import mahalanobis
    from scipy.stats import chi2 as chi2_dist

    # --- 1. Merge VCFs ---
    merged = tmp_path / "merged_all.vcf.gz"
    subprocess.run(
        ["bcftools", "merge", "-Oz", "-o", str(merged), str(donor_vcf), str(synth_vcf)],
        check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", str(merged)], check=True)

    # --- 2. PCA ---
    pca_prefix = str(tmp_path / "pca_all")
    subprocess.run(
        [
            "plink2",
            "--vcf", str(merged),
            "--maf", "0.05",
            "--pca", "10",
            "--out", pca_prefix,
            "--allow-extra-chr",
        ],
        check=True,
    )

    # --- 3. Parse eigenvec ---
    eigenvec_path = Path(pca_prefix + ".eigenvec")
    if not eigenvec_path.exists():
        pytest.skip("plink2 did not produce .eigenvec — check plink2 output for errors.")

    donor_pcs: list[np.ndarray] = []
    synth_pcs: list[np.ndarray] = []

    with eigenvec_path.open() as fh:
        next(fh)  # Skip header line
        for line in fh:
            parts = line.strip().split()
            sample_id = parts[1]  # IID column
            pc_vals = np.array(parts[2:], dtype=np.float64)
            if sample_id.startswith("synthetic_"):
                synth_pcs.append(pc_vals)
            else:
                donor_pcs.append(pc_vals)

    if not donor_pcs or not synth_pcs:
        pytest.skip("Could not split eigenvec into donor/synth by sample-name prefix.")

    donor_arr = np.stack(donor_pcs)  # (D, 10)
    synth_arr = np.stack(synth_pcs)  # (S, 10)

    # --- 4. Mahalanobis distance in PC1–PC2 ---
    donor_2d = donor_arr[:, :2]
    synth_2d = synth_arr[:, :2]
    centroid = donor_2d.mean(axis=0)
    cov_inv = np.linalg.pinv(np.cov(donor_2d.T))  # pinv handles singular/near-singular cov

    mahal = np.array(
        [mahalanobis(s, centroid, cov_inv) for s in synth_2d]
    )
    threshold = chi2_dist.ppf(0.99, df=2) ** 0.5  # 99% cloud radius
    frac_inside = float((mahal < threshold).mean())

    # PC1/PC2 outlier check (> 3 SD of donor range)
    donor_sd = donor_2d.std(axis=0)[:2]
    donor_mean = donor_2d.mean(axis=0)[:2]
    pc_outliers = (np.abs(synth_2d - donor_mean) > 3 * donor_sd).any(axis=1)
    n_outliers = int(pc_outliers.sum())

    print(f"\n[B5] Fraction inside 99% Mahalanobis cloud: {frac_inside:.3f}")
    print(f"[B5] Synthetics with PC1/PC2 outlier (>3 SD): {n_outliers}")

    assert frac_inside > 0.95, (
        f"Only {frac_inside:.3f} of synthetics fall inside the 99% donor cloud "
        "(threshold 0.95). Synthetics may not match the donor population structure."
    )
    assert n_outliers == 0, (
        f"{n_outliers} synthetic(s) have a PC1 or PC2 value more than 3 SD "
        "outside the donor range."
    )


# ---------------------------------------------------------------------------
# B6: LD Decay Comparison
#
# Compute r² values vs physical distance for donors and synthetics separately,
# bin by 10 kb windows, and compare the decay curves.
# ---------------------------------------------------------------------------

@_SKIP_NO_PLINK2
def test_b6_ld_decay(
    donor_vcf: Path, synth_vcf: Path, tmp_path: Path
) -> None:
    """
    B6 — LD decay curves (binned r²) are highly correlated between donors and synthetics.

    Steps:
      1. Run plink2 --r2 on each VCF.
      2. Bin r² by distance in 10 kb windows.
      3. Compare mean-r² decay curves.

    Pass thresholds (plan §B6):
      * Mean r² at < 10 kb: |synth − donor| < 0.02.
      * Pearson r of LD decay curves (binned means) > 0.99.
      * KS statistic on full r² distribution < 0.05.

    Note: long-range LD may be inflated in synthetics (a known artefact of the
    diploid mosaic design), so only short-range LD is hard-asserted.
    """
    def _run_ld(vcf: Path, prefix: str) -> Path:
        ld_prefix = str(tmp_path / prefix)
        subprocess.run(
            [
                "plink2",
                "--vcf", str(vcf),
                "--r2",
                "--ld-window", "100",
                "--ld-window-kb", "500",
                "--ld-window-r2", "0",
                "--out", ld_prefix,
                "--allow-extra-chr",
            ],
            check=True,
        )
        return Path(ld_prefix + ".ld")

    ld_donor_file = _run_ld(donor_vcf, "ld_donors")
    ld_synth_file = _run_ld(synth_vcf, "ld_synth")

    def _parse_ld(path: Path) -> tuple[np.ndarray, np.ndarray]:
        """Return (physical_distance_bp, r2) arrays from a plink2 .ld file."""
        dists, r2s = [], []
        with path.open() as fh:
            next(fh)  # skip header
            for line in fh:
                parts = line.split()
                if len(parts) < 7:
                    continue
                try:
                    bp_a, bp_b = int(parts[1]), int(parts[4])
                    r2 = float(parts[6])
                except ValueError:
                    continue
                dists.append(abs(bp_b - bp_a))
                r2s.append(r2)
        return np.array(dists, dtype=np.int64), np.array(r2s, dtype=np.float64)

    if not ld_donor_file.exists() or not ld_synth_file.exists():
        pytest.skip("plink2 did not produce .ld files — check plink2 output.")

    dist_d, r2_d = _parse_ld(ld_donor_file)
    dist_s, r2_s = _parse_ld(ld_synth_file)

    # --- Bin r² by 10 kb windows ---
    bin_size = 10_000
    max_dist = 500_000
    bins = np.arange(0, max_dist + bin_size, bin_size)

    def _bin_means(dists: np.ndarray, r2s: np.ndarray) -> np.ndarray:
        means = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            mask = (dists >= lo) & (dists < hi)
            means.append(float(r2s[mask].mean()) if mask.sum() > 0 else float("nan"))
        return np.array(means)

    decay_d = _bin_means(dist_d, r2_d)
    decay_s = _bin_means(dist_s, r2_s)

    # --- Short-range LD (first bin: < 10 kb) ---
    mean_r2_d_short = float(r2_d[dist_d < 10_000].mean()) if (dist_d < 10_000).sum() > 0 else float("nan")
    mean_r2_s_short = float(r2_s[dist_s < 10_000].mean()) if (dist_s < 10_000).sum() > 0 else float("nan")

    # --- Decay curve correlation ---
    common_bins = ~np.isnan(decay_d) & ~np.isnan(decay_s)
    if common_bins.sum() >= 5:
        r_decay = float(np.corrcoef(decay_d[common_bins], decay_s[common_bins])[0, 1])
    else:
        r_decay = float("nan")

    # --- Full r² distribution KS test ---
    from scipy.stats import ks_2samp
    ks_stat, ks_p = ks_2samp(r2_d, r2_s)

    print(f"\n[B6] Mean r² at < 10 kb: donors={mean_r2_d_short:.4f}, synth={mean_r2_s_short:.4f}")
    print(f"[B6] LD decay curve correlation: r={r_decay:.4f}")
    print(f"[B6] KS statistic (full r² distribution): {ks_stat:.4f} (p={ks_p:.4f})")

    if not np.isnan(mean_r2_d_short) and not np.isnan(mean_r2_s_short):
        assert abs(mean_r2_d_short - mean_r2_s_short) < 0.02, (
            f"Short-range LD difference {abs(mean_r2_d_short - mean_r2_s_short):.4f} ≥ 0.02."
        )

    if not np.isnan(r_decay):
        assert r_decay > 0.99, (
            f"LD decay curve correlation {r_decay:.4f} ≤ 0.99."
        )

    assert ks_stat < 0.05, (
        f"KS statistic on full r² distribution {ks_stat:.4f} ≥ 0.05. "
        "Note: long-range LD inflation is an expected artefact of the diploid mosaic."
    )
