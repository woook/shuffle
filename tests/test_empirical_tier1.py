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
    build_segment_plan,
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

N_REPS = 10_000  # Monte Carlo repetitions for each statistical test


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
