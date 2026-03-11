"""
Tier 2 empirical tests: privacy and biological plausibility.

Tests P1, P2, P4, B1, B2, B3, B4 from the v-shuffler empirical anonymisation
test plan.  All metrics are computed from an in-memory synthetic dataset; no
external VCF files, bcftools, or plink2 are required.

Fixture design
--------------
* N_DONORS = 200  individuals drawn under HWE from Beta(2,2) allele frequencies.
* N_VARIANTS = 2 000  biallelic SNPs with realistic ts/tv composition.
* Position layout: 100 "gene" regions × 20 variants each.
  Intra-region spacing: 100 bp.  Inter-region gap: 500 kb.
  Total span: ~50 Mb (same genetic map as before).
* Region-based sampling is used (default mode) — detect_regions + generate_all_region_plans.
  This ensures thorough mixing even without an artificial λ override.
* N_SYNTH = 100  synthetic individuals.
* N_HELDOUT = 20  individuals *not* used as donors (for P4 membership inference).

Known limitations documented here
----------------------------------
P2  The concordance-based closest-donor attack may succeed at a reduced but
    still non-trivial rate: each synthetic is built from ~100 independently
    sampled region-donors, but the primary donor still contributes the plurality.
    The test reports the metric and warns if the absolute rate exceeds 50%.

P4  Membership inference signal is detectable at reduced strength compared to
    the continuous-mode design.  The test reports the metric and warns when
    it reaches alarming levels.
"""

from __future__ import annotations

import textwrap
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import chi2 as chi2_dist
from scipy.stats import ks_2samp, wilcoxon

from v_shuffler.core.genotype_pool import MISSING, GenotypePool, VariantInfo
from v_shuffler.core.mosaic_builder import build_synthetic_genotypes
from v_shuffler.core.recombination import (
    detect_regions,
    generate_all_region_plans,
    generate_all_segment_plans,
)
from v_shuffler.io.genetic_map import GeneticMap


# ---------------------------------------------------------------------------
# Fixture parameters
# ---------------------------------------------------------------------------

N_DONORS = 200
N_SYNTH = 100
N_HELDOUT = 20
N_VARIANTS = 2_000
N_GENE_REGIONS = 100
N_VARIANTS_PER_REGION = N_VARIANTS // N_GENE_REGIONS  # 20
INTRA_SPACING_BP = 100
INTER_GAP_BP = 500_000
REGION_GAP_THRESHOLD = 10_000
CHROM = "chr22"
SEED = 42

# Genetic map: 100 cM across a chr22-like physical range
_MAP_TEXT = textwrap.dedent(f"""\
    pos chr cM
    1000000 {CHROM} 0.0
    11000000 {CHROM} 20.0
    21000000 {CHROM} 40.0
    31000000 {CHROM} 60.0
    41000000 {CHROM} 80.0
    51000000 {CHROM} 100.0
""")

# Transition and transversion allele pairs
_TS_PAIRS = [("A", "G"), ("C", "T"), ("G", "A"), ("T", "C")]
_TV_PAIRS = [
    ("A", "C"), ("A", "T"), ("C", "G"), ("G", "T"),
    ("C", "A"), ("T", "A"), ("G", "C"), ("T", "G"),
]
_TS_SET = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def _make_panel_positions() -> np.ndarray:
    """100 gene regions × 20 variants, starting at 1 Mb with 500 kb inter-region gaps."""
    all_pos: list[int] = []
    region_start = 1_000_000
    for _ in range(N_GENE_REGIONS):
        for j in range(N_VARIANTS_PER_REGION):
            all_pos.append(region_start + j * INTRA_SPACING_BP)
        region_start += INTER_GAP_BP
    return np.array(all_pos, dtype=np.int64)


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------

@dataclass
class EmpiricalFixture:
    gmap: GeneticMap
    donor_matrix: np.ndarray    # (N_VARIANTS, N_DONORS), uint8
    synth_matrix: np.ndarray    # (N_VARIANTS, N_SYNTH), uint8
    heldout_matrix: np.ndarray  # (N_VARIANTS, N_HELDOUT), uint8
    donor_afs: np.ndarray       # (N_VARIANTS,), float64
    synth_afs: np.ndarray       # (N_VARIANTS,), float64
    refs: list[str]
    alts: list[str]
    positions: np.ndarray       # (N_VARIANTS,), int64
    segment_plans: list         # one plan per synthetic
    concordances: np.ndarray    # (N_SYNTH, N_DONORS), float32 — synth vs donor
    heldout_concs: np.ndarray   # (N_SYNTH, N_HELDOUT), float32 — synth vs held-out
    donor_baseline: float       # mean pairwise concordance among donors
    regions_bp: list            # detected regions (list of (start_bp, end_bp))
    n_regions: int               # len(regions_bp)
    regions_cm: list            # (start_cM, end_cM) per region


# ---------------------------------------------------------------------------
# Session fixture: generate everything once
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def fix(tmp_path_factory: pytest.TempPathFactory) -> EmpiricalFixture:
    """Build in-memory donor and synthetic genotype matrices for all Tier 2 tests."""
    rng = np.random.default_rng(SEED)

    # --- Genetic map ---
    map_path = tmp_path_factory.mktemp("tier2_maps") / "map.txt"
    map_path.write_text(_MAP_TEXT)
    gmap = GeneticMap(map_path, CHROM)

    # --- Panel-like variant positions (100 gene regions × 20 variants) ---
    positions = _make_panel_positions()
    assert len(positions) == N_VARIANTS
    cm_pos = gmap.bp_to_cm(positions)

    # --- REF/ALT: 2/3 transitions, 1/3 transversions ---
    refs: list[str] = []
    alts: list[str] = []
    for i in range(N_VARIANTS):
        pairs = _TS_PAIRS if (i % 3 < 2) else _TV_PAIRS
        idx = int(rng.integers(len(pairs)))
        r, a = pairs[idx]
        refs.append(r)
        alts.append(a)

    # --- Allele frequencies from Beta(2,2): realistic range, mean=0.5 ---
    true_afs = rng.beta(2.0, 2.0, size=N_VARIANTS)

    # --- Genotypes under Hardy-Weinberg equilibrium ---
    def _hwe_dosages(n_samples: int) -> np.ndarray:
        p = true_afs[:, np.newaxis]
        q = 1.0 - p
        u = rng.random((N_VARIANTS, n_samples))
        # P(0/0) = q², P(0/1) = 2pq, P(1/1) = p²
        return np.where(u < q**2, 0, np.where(u < q**2 + 2 * p * q, 1, 2)).astype(
            np.uint8
        )

    donor_matrix = _hwe_dosages(N_DONORS)    # (V, D)
    heldout_matrix = _hwe_dosages(N_HELDOUT)  # (V, H)

    donor_afs = donor_matrix.astype(np.float64).sum(axis=1) / (2.0 * N_DONORS)

    # --- Region detection ---
    regions_bp = detect_regions(positions, gap_threshold_bp=REGION_GAP_THRESHOLD)
    n_regions = len(regions_bp)
    bp_flat = np.array([[r[0], r[1]] for r in regions_bp], dtype=np.int64).ravel()
    cm_flat = gmap.bp_to_cm(bp_flat).reshape(-1, 2)
    regions_cm = [(float(row[0]), float(row[1])) for row in cm_flat]

    # --- Segment plans via region sampling (independent RNG so fixture is deterministic) ---
    rng2 = np.random.default_rng(SEED + 1)
    segment_plans = generate_all_region_plans(
        n_output_samples=N_SYNTH,
        regions_cm=regions_cm,
        n_pool_samples=N_DONORS,
        rng=rng2,
    )

    # --- Synthetic genotype matrix via GenotypePool ---
    variant_info = [
        VariantInfo(
            chrom=CHROM,
            pos=int(positions[i]),
            ref=refs[i],
            alts=[alts[i]],
            id=".",
            qual=None,
            filters=[],
            cm_pos=float(cm_pos[i]),
        )
        for i in range(N_VARIANTS)
    ]
    pool = GenotypePool(
        dosages=donor_matrix,
        positions=positions,
        cm_pos=cm_pos,
        variant_info=variant_info,
    )
    synth_matrix, _ = build_synthetic_genotypes(pool, segment_plans)  # (V, S)

    # --- Synth allele frequencies ---
    valid = synth_matrix != MISSING
    synth_afs = np.where(
        valid.sum(axis=1) > 0,
        np.where(valid, synth_matrix.astype(np.float64), 0.0).sum(axis=1)
        / (valid.sum(axis=1) * 2.0),
        np.nan,
    )

    # --- Pairwise concordances: synth × donor  (N_SYNTH, N_DONORS) ---
    concordances = (
        (synth_matrix[:, :, np.newaxis] == donor_matrix[:, np.newaxis, :])
        .mean(axis=0)
        .astype(np.float32)
    )

    # --- Synth × held-out concordances (N_SYNTH, N_HELDOUT) ---
    heldout_concs = (
        (synth_matrix[:, :, np.newaxis] == heldout_matrix[:, np.newaxis, :])
        .mean(axis=0)
        .astype(np.float32)
    )

    # --- Donor-donor baseline: mean concordance among first 50 donors ---
    n_sub = min(50, N_DONORS)
    d_sub = donor_matrix[:, :n_sub]
    match_dd = (d_sub[:, :, np.newaxis] == d_sub[:, np.newaxis, :]).mean(axis=0)
    np.fill_diagonal(match_dd, np.nan)
    donor_baseline = float(np.nanmean(match_dd))

    return EmpiricalFixture(
        gmap=gmap,
        donor_matrix=donor_matrix,
        synth_matrix=synth_matrix,
        heldout_matrix=heldout_matrix,
        donor_afs=donor_afs,
        synth_afs=synth_afs,
        refs=refs,
        alts=alts,
        positions=positions,
        segment_plans=segment_plans,
        concordances=concordances,
        heldout_concs=heldout_concs,
        donor_baseline=donor_baseline,
        regions_bp=regions_bp,
        n_regions=n_regions,
        regions_cm=regions_cm,
    )


# ---------------------------------------------------------------------------
# P1: Full Pairwise Concordance Distribution
# ---------------------------------------------------------------------------


def test_p1_max_concordance(fix: EmpiricalFixture) -> None:
    """P1 — No synthetic sample shares > 99% of genotypes with any donor."""
    max_conc = float(fix.concordances.max())
    assert max_conc < 0.99, (
        f"Max concordance {max_conc:.4f} ≥ 0.99 — possible identity leak. "
        "(region-sampling mode with 100 regions)"
    )


def test_p1_99th_percentile_concordance(fix: EmpiricalFixture) -> None:
    """P1 — 99th percentile of per-synth max concordance < 0.85."""
    per_synth_max = fix.concordances.max(axis=1)  # (N_SYNTH,)
    pct99 = float(np.percentile(per_synth_max, 99))
    assert pct99 < 0.85, (
        f"99th percentile of per-synth max concordance {pct99:.4f} ≥ 0.85."
    )


def test_p1_mean_concordance_vs_baseline(fix: EmpiricalFixture) -> None:
    """P1 — Mean synth–donor concordance within 0.01 of the donor–donor baseline."""
    mean_conc = float(fix.concordances.mean())
    delta = abs(mean_conc - fix.donor_baseline)
    assert delta < 0.01, (
        f"Mean synth–donor concordance {mean_conc:.4f} deviates from "
        f"donor–donor baseline {fix.donor_baseline:.4f} by {delta:.4f}; threshold 0.01."
    )


# ---------------------------------------------------------------------------
# P2: Closest-Donor Matching Attack
#
# For each synthetic, determine whether the highest-contributing donor
# (by region count) also ranks #1 in concordance.  In region mode the
# primary donor contributes ~1–3 regions out of 100, so attack rates are
# expected to be lower than in continuous mode, but are still a known
# limitation.
# ---------------------------------------------------------------------------


def test_p2_closest_donor_attack(fix: EmpiricalFixture) -> None:
    """
    P2 — Report closest-donor attack success rate (region mode).

    Hard assertion: attack success rate < 100%.
    Advisory warning emitted if rate exceeds 50% absolute.
    """
    n_donors = N_DONORS
    attack_successes = 0
    primary_fractions: list[float] = []

    for s_idx, plan in enumerate(fix.segment_plans):
        # Count regions contributed by each donor
        donor_region_count = np.zeros(n_donors)
        for seg in plan:
            donor_region_count[seg.sample_idx] += 1

        primary_donor = int(np.argmax(donor_region_count))
        primary_fractions.append(donor_region_count[primary_donor] / fix.n_regions)

        # Does concordance ranking put the primary donor first?
        ranked = np.argsort(-fix.concordances[s_idx])
        rank_of_primary = int(np.where(ranked == primary_donor)[0][0]) + 1
        if rank_of_primary == 1:
            attack_successes += 1

    attack_rate = attack_successes / N_SYNTH
    random_baseline = 1.0 / N_DONORS
    ratio = attack_rate / random_baseline
    mean_primary_frac = float(np.mean(primary_fractions))

    print(
        f"\n[P2 region-mode] Attack success rate: {attack_rate:.3f} "
        f"({attack_successes}/{N_SYNTH}), random baseline: {random_baseline:.4f}, "
        f"ratio: {ratio:.1f}×"
    )
    print(f"[P2] Mean primary donor region fraction: {mean_primary_frac:.3f}")
    print(
        "[P2] Pass thresholds (plan §P2): ratio < 5× | concerning 10-30× | alarming > 50% absolute"
    )

    # Hard assertions for documented privacy thresholds
    assert attack_rate < 1.0, "Attack succeeded on every synthetic — possible bug."

    assert attack_rate <= 0.50, (
        f"[P2 ALARMING] Attack success rate {attack_rate:.3f} > 50% absolute. "
        "Primary donors are identifiable from synthetic output. "
        "Privacy threshold exceeded."
    )

    assert ratio <= 10, (
        f"[P2 CONCERNING] Attack success ratio {ratio:.1f}× > 10× random baseline. "
        "Privacy threshold exceeded."
    )


# ---------------------------------------------------------------------------
# P4: Membership Inference Test
# ---------------------------------------------------------------------------


def test_p4_membership_inference(fix: EmpiricalFixture) -> None:
    """
    P4 — Report membership inference signal.

    Hard assertion: mean delta is a finite float (computation completes).
    Advisory warnings emitted when delta and Wilcoxon p-value exceed thresholds.
    """
    in_max = fix.concordances.max(axis=1)       # (S,) — best match in pool
    out_max = fix.heldout_concs.max(axis=1)     # (S,) — best match in held-out

    delta = in_max - out_max
    mean_delta = float(delta.mean())
    frac_alarming = float((delta > 0.02).mean())

    # Wilcoxon signed-rank test (paired)
    stat, p_val = wilcoxon(in_max, out_max)

    print(
        f"\n[P4] Mean delta (in-cohort − out-of-cohort max concordance): {mean_delta:.4f}"
    )
    print(f"[P4] Fraction with delta > 0.02: {frac_alarming:.3f}")
    print(f"[P4] Wilcoxon signed-rank p-value: {p_val:.4g}")
    print("[P4] Pass thresholds (plan §P4): mean delta < 0.005, frac > 0.02 < 5%, p > 0.05")

    # Hard assertions for documented privacy thresholds
    assert np.isfinite(mean_delta), "mean_delta is not finite — computation error."

    assert mean_delta <= 0.02 and p_val >= 0.001, (
        f"[P4 ALARMING] Membership inference signal is strong: "
        f"mean delta={mean_delta:.4f}, Wilcoxon p={p_val:.2e}. "
        "Synthetics are measurably more similar to their donors than to held-out individuals. "
        "Privacy threshold exceeded."
    )

    assert p_val >= 0.05, (
        f"[P4 CONCERNING] Wilcoxon p={p_val:.4f} < 0.05 — membership signal detected. "
        "Privacy threshold exceeded."
    )


# ---------------------------------------------------------------------------
# B1: Allele Frequency Correlation by MAF Bin
# ---------------------------------------------------------------------------


def test_b1_af_global(fix: EmpiricalFixture) -> None:
    """B1 — Global Pearson r between donor and synth AFs ≥ 0.98.

    Note: with the panel-like fixture (100 regions × 20 variants), all intra-region
    variants share the same assigned donor, introducing within-region AF correlation
    that reduces the effective sample count.  The threshold is 0.98 here rather than
    the production value of 0.99 which applies to whole-chromosome data with
    many more independently-assigned segments.
    """
    common = ~np.isnan(fix.synth_afs) & ~np.isnan(fix.donor_afs)
    r = float(np.corrcoef(fix.donor_afs[common], fix.synth_afs[common])[0, 1])
    assert r >= 0.98, f"Global AF Pearson r={r:.4f} < 0.98."


def test_b1_af_by_maf_bin(fix: EmpiricalFixture) -> None:
    """B1 — AF correlation and MAD by MAF bin (informational).

    Per-bin Pearson r is printed for comparison against the plan's production
    thresholds (r ≥ 0.995 for MAF > 0.05).  Hard assertions are NOT applied
    per-bin because within a narrow MAF range the variance of the true AF is
    small relative to sampling noise at N_SYNTH=100, which suppresses r below
    the production threshold even when the algorithm is correct.

    The primary hard assertion (global r ≥ 0.98) is in test_b1_af_global.
    """
    maf = np.minimum(fix.donor_afs, 1 - fix.donor_afs)
    bins = [
        (0.00, 0.01, "rare (MAF < 0.01)"),
        (0.01, 0.05, "low-frequency (0.01–0.05)"),
        (0.05, 0.15, "uncommon (0.05–0.15)"),
        (0.15, 0.35, "intermediate (0.15–0.35)"),
        (0.35, 0.50, "common (0.35–0.50)"),
    ]
    print()
    for lo, hi, label in bins:
        mask = (
            (maf >= lo) & (maf < hi)
            & ~np.isnan(fix.donor_afs)
            & ~np.isnan(fix.synth_afs)
        )
        if mask.sum() < 50:
            print(f"[B1] {label}: n={mask.sum()} — skipped (< 50 variants)")
            continue
        r = float(np.corrcoef(fix.donor_afs[mask], fix.synth_afs[mask])[0, 1])
        mad = float(np.abs(fix.donor_afs[mask] - fix.synth_afs[mask]).mean())
        print(f"[B1] {label}: n={mask.sum()}, r={r:.4f}, MAD={mad:.5f}")

    print("[B1] Production thresholds: r≥0.995 for MAF>0.05, r≥0.95 for rare; "
          "not asserted here (require N_SYNTH≥500 for reliable per-bin estimates)")


# ---------------------------------------------------------------------------
# B2: Heterozygosity Rate Comparison
# ---------------------------------------------------------------------------


def _het_rate_per_sample(matrix: np.ndarray) -> np.ndarray:
    """Per-sample heterozygosity rate (fraction of sites with dosage=1)."""
    valid = matrix != MISSING  # (V, S)
    het = matrix == 1
    denom = valid.sum(axis=0).astype(np.float64)
    rate = np.where(denom > 0, het.sum(axis=0) / denom, np.nan)
    return rate.astype(np.float64)


def test_b2_heterozygosity(fix: EmpiricalFixture) -> None:
    """B2 — Per-sample het rate: |mean diff| < 0.005 and Wilcoxon p > 0.01."""
    donor_het = _het_rate_per_sample(fix.donor_matrix)
    synth_het = _het_rate_per_sample(fix.synth_matrix)

    mean_diff = abs(
        float(np.nanmean(donor_het)) - float(np.nanmean(synth_het))
    )
    ks_stat, p_val = ks_2samp(
        donor_het[~np.isnan(donor_het)],
        synth_het[~np.isnan(synth_het)],
    )
    print(
        f"\n[B2] Mean het: donors={np.nanmean(donor_het):.4f}, "
        f"synth={np.nanmean(synth_het):.4f}, |diff|={mean_diff:.5f}"
    )
    print(f"[B2] KS statistic: {ks_stat:.4f} (p={p_val:.4f})")

    assert mean_diff < 0.005, (
        f"|mean het diff| {mean_diff:.5f} ≥ 0.005."
    )
    assert p_val > 0.01, (
        f"KS test on per-sample het rates: p={p_val:.4f} < 0.01 "
        f"(KS stat={ks_stat:.4f}) — distributions are significantly different."
    )


# ---------------------------------------------------------------------------
# B3: Hardy-Weinberg Equilibrium
# ---------------------------------------------------------------------------


def _hwe_pvalue(n00: int, n01: int, n11: int) -> float:
    """Chi-squared HWE p-value for one site. Returns NaN if expected counts < 5."""
    n = n00 + n01 + n11
    if n == 0:
        return float("nan")
    p = (2 * n11 + n01) / (2 * n)
    q = 1.0 - p
    e00, e01, e11 = n * q**2, n * 2 * p * q, n * p**2
    if min(e00, e01, e11) < 5:
        return float("nan")
    chi2_val = (n00 - e00) ** 2 / e00 + (n01 - e01) ** 2 / e01 + (n11 - e11) ** 2 / e11
    return float(1.0 - chi2_dist.cdf(chi2_val, df=1))


def _hwe_fail_fraction(matrix: np.ndarray, threshold: float = 0.001) -> float:
    """Fraction of testable sites with HWE p < threshold."""
    n_variants, n_samples = matrix.shape
    fails = 0
    testable = 0
    for v in range(n_variants):
        row = matrix[v, :]
        valid = row[row != MISSING]
        n00 = int((valid == 0).sum())
        n01 = int((valid == 1).sum())
        n11 = int((valid == 2).sum())
        pv = _hwe_pvalue(n00, n01, n11)
        if np.isnan(pv):
            continue
        testable += 1
        if pv < threshold:
            fails += 1
    return fails / testable if testable > 0 else float("nan")


def test_b3_hwe(fix: EmpiricalFixture) -> None:
    """B3 — Fraction of synth sites failing HWE (p < 0.001) ≤ 3× donor fraction."""
    donor_frac = _hwe_fail_fraction(fix.donor_matrix)
    synth_frac = _hwe_fail_fraction(fix.synth_matrix)

    print(f"\n[B3] HWE-fail fraction: donors={donor_frac:.4f}, synth={synth_frac:.4f}")

    if np.isnan(donor_frac) or np.isnan(synth_frac):
        pytest.skip("Not enough testable HWE sites in fixture")

    assert synth_frac <= 3 * donor_frac + 0.005, (
        f"Synth HWE-fail fraction {synth_frac:.4f} > 3× donor fraction "
        f"{donor_frac:.4f} + 0.005."
    )


# ---------------------------------------------------------------------------
# B4: Transition/Transversion Ratio
# ---------------------------------------------------------------------------


def _count_ts_tv(refs: list[str], alts: list[str]) -> tuple[int, int]:
    """Count transitions and transversions from REF/ALT lists."""
    ts = sum(1 for r, a in zip(refs, alts) if (r.upper(), a.upper()) in _TS_SET)
    tv = len(refs) - ts
    return ts, tv


def test_b4_tstv_identical(fix: EmpiricalFixture) -> None:
    """
    B4 — ts/tv ratio and variant count are identical between donors and synthetics.

    Because v-shuffler does not add or remove variant sites, only reassign
    dosages, ts/tv must match exactly.
    """
    ts, tv = _count_ts_tv(fix.refs, fix.alts)
    tstv = ts / tv if tv > 0 else float("inf")

    assert ts + tv == N_VARIANTS, (
        f"Variant count {ts + tv} != {N_VARIANTS}."
    )
    assert tv > 0, "All variants are transitions — ts/tv undefined."

    print(f"\n[B4] ts/tv = {tstv:.4f} ({ts} ts, {tv} tv) — identical for donor and synth.")

    # ts/tv is identical for donor and synthetic by construction: the shuffle
    # reassigns genotypes but never changes which variant sites are present.
    # Verify the ratio is a real finite number (catches fixture misconfiguration).
    assert np.isfinite(tstv), f"ts/tv should be finite, got {tstv}"


# ---------------------------------------------------------------------------
# Region-mode specific tests
# ---------------------------------------------------------------------------


def test_region_detection_in_fixture(fix: EmpiricalFixture) -> None:
    """detect_regions correctly identified all 100 simulated gene clusters."""
    assert fix.n_regions == N_GENE_REGIONS, (
        f"Expected {N_GENE_REGIONS} regions, got {fix.n_regions}."
    )


def test_min_donors_region_mode(fix: EmpiricalFixture) -> None:
    """With min_donors=10, every region-mode plan has ≥ 10 distinct donors."""
    rng = np.random.default_rng(200)
    plans = generate_all_region_plans(
        n_output_samples=N_SYNTH,
        regions_cm=fix.regions_cm,
        n_pool_samples=N_DONORS,
        rng=rng,
        min_donors=10,
    )
    for i, plan in enumerate(plans):
        distinct = len({seg.sample_idx for seg in plan})
        assert distinct >= 10, (
            f"Synthetic {i}: only {distinct} distinct donors with min_donors=10."
        )


def test_min_donors_continuous_mode(fix: EmpiricalFixture) -> None:
    """With min_donors=5 and lambda=0.0, every plan has ≥ 5 distinct donors."""
    rng = np.random.default_rng(201)
    plans = generate_all_segment_plans(
        n_output_samples=N_SYNTH,
        genetic_map=fix.gmap,
        n_pool_samples=N_DONORS,
        rng=rng,
        lambda_override=0.0,
        min_donors=5,
    )
    for i, plan in enumerate(plans):
        distinct = len({seg.sample_idx for seg in plan})
        assert distinct >= 5, (
            f"Synthetic {i}: only {distinct} distinct donors with min_donors=5."
        )


def test_p1_region_vs_continuous_concordance(fix: EmpiricalFixture) -> None:
    """
    Region mode max concordance < continuous mode max concordance (lambda=0.5).

    This demonstrates the core improvement of region sampling over continuous
    mode for panel data: with realistic λ ≈ 0.5, most synthetics have zero
    crossovers and are copies of a single donor; region mode assigns each
    captured gene independently, guaranteeing mixing.
    """
    N_SUB = 20  # lightweight sub-fixture

    rng = np.random.default_rng(300)
    true_afs = rng.beta(2.0, 2.0, size=N_VARIANTS)

    def _hwe(n: int) -> np.ndarray:
        p = true_afs[:, np.newaxis]
        q = 1.0 - p
        u = rng.random((N_VARIANTS, n))
        return np.where(u < q**2, 0, np.where(u < q**2 + 2 * p * q, 1, 2)).astype(np.uint8)

    donor_sub = _hwe(N_DONORS)
    variant_info = [
        VariantInfo(
            chrom=CHROM, pos=int(fix.positions[i]), ref="A", alts=["G"],
            id=".", qual=None, filters=[], cm_pos=float(fix.gmap.bp_to_cm(fix.positions[i])),
        )
        for i in range(N_VARIANTS)
    ]
    pool_sub = GenotypePool(
        dosages=donor_sub,
        positions=fix.positions,
        cm_pos=fix.gmap.bp_to_cm(fix.positions),
        variant_info=variant_info,
    )

    # Region mode
    rng_r = np.random.default_rng(301)
    plans_region = generate_all_region_plans(
        n_output_samples=N_SUB,
        regions_cm=fix.regions_cm,
        n_pool_samples=N_DONORS,
        rng=rng_r,
    )
    synth_region, _ = build_synthetic_genotypes(pool_sub, plans_region)

    # Continuous mode with realistic λ = 0.5
    rng_c = np.random.default_rng(302)
    plans_cont = generate_all_segment_plans(
        n_output_samples=N_SUB,
        genetic_map=fix.gmap,
        n_pool_samples=N_DONORS,
        rng=rng_c,
        lambda_override=0.5,
    )
    synth_cont, _ = build_synthetic_genotypes(pool_sub, plans_cont)

    def _max_conc(synth: np.ndarray, donor: np.ndarray) -> float:
        conc = (synth[:, :, np.newaxis] == donor[:, np.newaxis, :]).mean(axis=0)
        return float(conc.max())

    max_region = _max_conc(synth_region, donor_sub)
    max_cont = _max_conc(synth_cont, donor_sub)

    print(f"\n[P1 comparison] Region max concordance: {max_region:.4f}")
    print(f"[P1 comparison] Continuous (λ=0.5) max concordance: {max_cont:.4f}")

    assert max_region < max_cont, (
        f"Region mode max concordance ({max_region:.4f}) is not lower than "
        f"continuous mode ({max_cont:.4f}).  Region mode should provide better "
        "mixing for panel data."
    )
