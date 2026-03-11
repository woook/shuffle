"""
Tier 1 empirical tests: recombination model statistical correctness.

Tests R1, R2, R3 from the v-shuffler empirical anonymisation test plan.
No external tools or VCF data required — these operate purely on v-shuffler
recombination internals + scipy/numpy.

Run time: ~30 s (10 000 simulations each).
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import kstest, poisson

from v_shuffler.core.recombination import (
    build_region_segment_plan,
    build_segment_plan,
    detect_regions,
    generate_all_region_plans,
    simulate_crossover_breakpoints,
)
from v_shuffler.io.genetic_map import GeneticMap


# ---------------------------------------------------------------------------
# Shared map: 0–150 cM across a 15 Mb region
# ---------------------------------------------------------------------------

_MAP_TEXT = textwrap.dedent("""\
    pos chr cM
    1000 chr1 0.0
    5000000 chr1 50.0
    10000000 chr1 100.0
    15000000 chr1 150.0
""")

N_REPS = 10_000  # Monte Carlo repetitions for R1–R3

# ---------------------------------------------------------------------------
# Panel-like positions fixture for R4–R6
# 100 regions × 20 variants, 100 bp intra-region spacing, 500 kb inter-region gaps
# ---------------------------------------------------------------------------

N_REGIONS = 100
N_VARIANTS_PER_REGION = 20
INTRA_SPACING_BP = 100
INTER_GAP_BP = 500_000
REGION_GAP_THRESHOLD = 10_000


def _make_panel_positions() -> np.ndarray:
    """Build the 100-region panel position array (deterministic)."""
    all_pos: list[int] = []
    region_start = 1_000_000
    for _ in range(N_REGIONS):
        for j in range(N_VARIANTS_PER_REGION):
            all_pos.append(region_start + j * INTRA_SPACING_BP)
        region_start += INTER_GAP_BP
    return np.array(all_pos, dtype=np.int64)


_PANEL_POSITIONS = _make_panel_positions()


@pytest.fixture(scope="module")
def gmap_150(tmp_path_factory: pytest.TempPathFactory) -> GeneticMap:
    p = tmp_path_factory.mktemp("tier1_maps") / "map_150cM.txt"
    p.write_text(_MAP_TEXT)
    return GeneticMap(p, "chr1")


# ---------------------------------------------------------------------------
# Pre-computed crossover counts (shared by R1 tests)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def crossover_sample(
    gmap_150: GeneticMap,
) -> tuple[np.ndarray, float]:
    """Draw N_REPS crossover counts; return (counts, lambda)."""
    lam = gmap_150.total_length_cm / 100.0
    rng = np.random.default_rng(42)
    counts = np.array(
        [len(simulate_crossover_breakpoints(gmap_150, rng)) for _ in range(N_REPS)],
        dtype=np.int64,
    )
    return counts, lam


# ---------------------------------------------------------------------------
# R1: Crossover Count Distribution
#
# Number of crossovers per meiosis should follow Poisson(λ),
# λ = total_cM / 100.  We verify the mean, dispersion, and the
# cumulative KS statistic against the Poisson CDF.
# ---------------------------------------------------------------------------


def test_r1_mean(crossover_sample: tuple[np.ndarray, float]) -> None:
    """R1 — |empirical mean − λ| / λ < 1%."""
    counts, lam = crossover_sample
    rel_error = abs(counts.mean() - lam) / lam
    assert rel_error < 0.01, (
        f"Crossover mean {counts.mean():.4f} deviates from λ={lam:.4f} "
        f"by {rel_error:.2%}; threshold is 1%."
    )


def test_r1_dispersion(crossover_sample: tuple[np.ndarray, float]) -> None:
    """R1 — Poisson dispersion (var/mean) in [0.95, 1.05]."""
    counts, _ = crossover_sample
    dispersion = counts.var() / counts.mean()
    assert 0.95 <= dispersion <= 1.05, (
        f"Dispersion {dispersion:.4f} outside [0.95, 1.05]; "
        "expected ~1.0 for a Poisson process."
    )


def test_r1_ks_vs_poisson(crossover_sample: tuple[np.ndarray, float]) -> None:
    """R1 — Discrete KS statistic vs Poisson(λ) CDF < 0.02."""
    counts, lam = crossover_sample
    max_k = int(counts.max()) + 1
    obs_cdf = np.cumsum(np.bincount(counts, minlength=max_k) / len(counts))
    theo_cdf = np.cumsum([poisson.pmf(k, lam) for k in range(max_k)])
    ks_stat = float(np.abs(obs_cdf - theo_cdf).max())
    assert ks_stat < 0.02, (
        f"KS statistic {ks_stat:.4f} vs Poisson({lam:.2f}) exceeds 0.02."
    )


# ---------------------------------------------------------------------------
# R2: Breakpoint Position Uniformity in cM Space
#
# Crossover positions are drawn as uniform[start_cM, end_cM], so after
# normalising to [0,1] they should be indistinguishable from Uniform(0,1).
# Physical space is NOT expected to be uniform (recombination hotspots).
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def breakpoints_normalised(gmap_150: GeneticMap) -> np.ndarray:
    """
    Collect all crossover positions from N_REPS meioses,
    normalised to [0, 1] in cM space.
    """
    rng = np.random.default_rng(42)
    cm_positions: list[float] = []
    for _ in range(N_REPS):
        bps = simulate_crossover_breakpoints(gmap_150, rng)
        cm_positions.extend(bps.tolist())
    arr = np.array(cm_positions)
    return (arr - gmap_150.start_cm) / gmap_150.total_length_cm


def test_r2_cM_uniformity(breakpoints_normalised: np.ndarray) -> None:
    """R2 — KS statistic vs Uniform(0,1) in cM space < 0.02."""
    stat, p = kstest(breakpoints_normalised, "uniform")
    assert stat < 0.02, (
        f"KS stat {stat:.4f} (p={p:.4f}) vs Uniform(0,1) in cM space exceeds 0.02. "
        "Breakpoints should be drawn uniformly in genetic distance."
    )


def test_r2_mean_position(breakpoints_normalised: np.ndarray) -> None:
    """R2 — Mean normalised position in [0.49, 0.51]."""
    mean_pos = float(breakpoints_normalised.mean())
    assert 0.49 <= mean_pos <= 0.51, (
        f"Mean normalised cM position {mean_pos:.4f} outside [0.49, 0.51]."
    )


# ---------------------------------------------------------------------------
# R3: Segment Length Distribution
#
# Expected mean segment length = L_cM / (λ + 1).
# We test this to within 2% relative error.
# ---------------------------------------------------------------------------


def test_r3_mean_segment_length(gmap_150: GeneticMap) -> None:
    """R3 — |empirical mean segment length − expected| / expected < 2%."""
    lam = gmap_150.total_length_cm / 100.0
    expected_mean = gmap_150.total_length_cm / (lam + 1)

    rng = np.random.default_rng(0)
    all_lengths: list[float] = []
    for _ in range(N_REPS):
        bps = simulate_crossover_breakpoints(gmap_150, rng)
        plan = build_segment_plan(bps, gmap_150, n_samples=200, rng=rng)
        all_lengths.extend(seg.cm_end - seg.cm_start for seg in plan)

    empirical_mean = float(np.mean(all_lengths))
    rel_error = abs(empirical_mean - expected_mean) / expected_mean
    assert rel_error < 0.02, (
        f"Mean segment length {empirical_mean:.3f} cM deviates from "
        f"expected {expected_mean:.3f} cM by {rel_error:.2%}; threshold is 2%."
    )


# ---------------------------------------------------------------------------
# R4: Region Detection Accuracy
# ---------------------------------------------------------------------------


def test_r4_region_detection_count() -> None:
    """R4 — detect_regions returns exactly N_REGIONS regions from the panel positions."""
    regions = detect_regions(_PANEL_POSITIONS, gap_threshold_bp=REGION_GAP_THRESHOLD)
    assert len(regions) == N_REGIONS, (
        f"Expected {N_REGIONS} regions, got {len(regions)}."
    )


def test_r4_region_boundaries() -> None:
    """R4 — Each region spans exactly N_VARIANTS_PER_REGION variants."""
    regions = detect_regions(_PANEL_POSITIONS, gap_threshold_bp=REGION_GAP_THRESHOLD)
    expected_span = (N_VARIANTS_PER_REGION - 1) * INTRA_SPACING_BP
    for i, (start, end) in enumerate(regions):
        span = end - start
        assert span == expected_span, (
            f"Region {i}: span={span} bp, expected {expected_span} bp."
        )


# ---------------------------------------------------------------------------
# R5: Cross-Region Donor Independence
# ---------------------------------------------------------------------------

_N_REPS_R5 = 1_000
_N_DONORS_R5 = 200


def test_r5_cross_region_independence(gmap_150: GeneticMap) -> None:
    """
    R5 — Cross-region co-assignment frequency is within 2% of 1/N_DONORS.

    Under independence, the probability that the same donor is assigned to
    region i and region j equals 1/N_DONORS.  We verify the mean co-assignment
    rate over all region pairs is close to this expectation.
    """
    regions = detect_regions(_PANEL_POSITIONS, gap_threshold_bp=REGION_GAP_THRESHOLD)
    bp_flat = np.array([[r[0], r[1]] for r in regions], dtype=np.int64).ravel()
    # Synthesise a trivial uniform cM map for region coordinates
    cm_flat = (bp_flat / 1_000_000.0)  # 1 Mb → 1 cM (arbitrary scale)
    regions_cm = [(float(cm_flat[2 * i]), float(cm_flat[2 * i + 1])) for i in range(len(regions))]

    rng = np.random.default_rng(42)
    n_regions = len(regions_cm)

    # Accumulate co-assignment counts across Monte Carlo reps
    # co_count[i, j] = number of reps where region i and j share the same donor
    co_count = np.zeros((n_regions, n_regions), dtype=np.int64)

    for _ in range(_N_REPS_R5):
        plan = build_region_segment_plan(regions_cm, n_samples=_N_DONORS_R5, rng=rng)
        donors = np.array([seg.sample_idx for seg in plan], dtype=np.int64)
        # Broadcast comparison: co_assign[i,j] = (donors[i] == donors[j])
        co_count += (donors[:, np.newaxis] == donors[np.newaxis, :]).astype(np.int64)

    # Exclude diagonal (self-comparison)
    np.fill_diagonal(co_count, 0)
    n_off_diag = n_regions * (n_regions - 1)
    mean_co_rate = float(co_count.sum()) / (_N_REPS_R5 * n_off_diag)
    expected_rate = 1.0 / _N_DONORS_R5

    assert abs(mean_co_rate - expected_rate) < 0.02 * expected_rate + 0.0005, (
        f"Mean cross-region co-assignment rate {mean_co_rate:.6f} deviates from "
        f"expected {expected_rate:.6f} by more than 2% + 0.0005.  "
        "This suggests regions are NOT assigned independently."
    )


# ---------------------------------------------------------------------------
# R6: min_donors Constraint Reliability
# ---------------------------------------------------------------------------

_N_REPS_R6 = 1_000


def test_r6_min_donors_always_satisfied() -> None:
    """R6 — With min_donors=10 and 200 donors, 100% of plans have ≥ 10 distinct donors."""
    regions = detect_regions(_PANEL_POSITIONS, gap_threshold_bp=REGION_GAP_THRESHOLD)
    bp_flat = np.array([[r[0], r[1]] for r in regions], dtype=np.int64).ravel()
    cm_flat = bp_flat / 1_000_000.0
    regions_cm = [(float(cm_flat[2 * i]), float(cm_flat[2 * i + 1])) for i in range(len(regions))]

    rng = np.random.default_rng(7)
    min_donors = 10
    n_samples = 200

    fails = 0
    for _ in range(_N_REPS_R6):
        plan = build_region_segment_plan(regions_cm, n_samples=n_samples, rng=rng, min_donors=min_donors)
        if len({seg.sample_idx for seg in plan}) < min_donors:
            fails += 1

    assert fails == 0, (
        f"{fails}/{_N_REPS_R6} plans had fewer than {min_donors} distinct donors."
    )


def test_r6_min_donors_capped_by_n_regions() -> None:
    """R6 — With min_donors=200 and 100 regions, all plans have exactly 100 distinct donors."""
    regions = detect_regions(_PANEL_POSITIONS, gap_threshold_bp=REGION_GAP_THRESHOLD)
    bp_flat = np.array([[r[0], r[1]] for r in regions], dtype=np.int64).ravel()
    cm_flat = bp_flat / 1_000_000.0
    regions_cm = [(float(cm_flat[2 * i]), float(cm_flat[2 * i + 1])) for i in range(len(regions))]

    rng = np.random.default_rng(8)
    n_samples = 200

    for _ in range(_N_REPS_R6):
        plan = build_region_segment_plan(regions_cm, n_samples=n_samples, rng=rng, min_donors=200)
        distinct = len({seg.sample_idx for seg in plan})
        # effective_min is capped by min(min_donors, n_samples, n_regions) = 100
        assert distinct == N_REGIONS, (
            f"Expected exactly {N_REGIONS} distinct donors (capped by n_regions), got {distinct}."
        )
