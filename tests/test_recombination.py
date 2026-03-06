"""Tests for v_shuffler.core.recombination."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from v_shuffler.io.genetic_map import GeneticMap
from v_shuffler.core.recombination import (
    Segment,
    build_segment_plan,
    simulate_crossover_breakpoints,
    generate_all_segment_plans,
)


# ---------------------------------------------------------------------------
# Fixture: a small genetic map covering 0–150 cM
# ---------------------------------------------------------------------------

MAP_TEXT = textwrap.dedent("""\
    pos chr cM
    1000 chr1 0.0
    5000000 chr1 50.0
    10000000 chr1 100.0
    15000000 chr1 150.0
""")


@pytest.fixture
def gmap(tmp_path: Path) -> GeneticMap:
    p = tmp_path / "map.txt"
    p.write_text(MAP_TEXT)
    return GeneticMap(p, "chr1")


def make_rng(seed: int = 42) -> np.random.Generator:
    return np.random.default_rng(seed)


# ---------------------------------------------------------------------------
# simulate_crossover_breakpoints
# ---------------------------------------------------------------------------

def test_breakpoints_within_map(gmap: GeneticMap) -> None:
    rng = make_rng(1)
    for _ in range(50):
        bps = simulate_crossover_breakpoints(gmap, rng)
        assert np.all(bps >= gmap.start_cm)
        assert np.all(bps <= gmap.end_cm)


def test_breakpoints_sorted(gmap: GeneticMap) -> None:
    rng = make_rng(2)
    bps = simulate_crossover_breakpoints(gmap, rng, lambda_override=10.0)
    assert np.all(np.diff(bps) >= 0)


def test_zero_lambda_gives_empty(gmap: GeneticMap) -> None:
    """With lambda=0 the Poisson draw is always 0."""
    rng = make_rng(3)
    for _ in range(20):
        bps = simulate_crossover_breakpoints(gmap, rng, lambda_override=0.0)
        assert len(bps) == 0


def test_poisson_mean(gmap: GeneticMap) -> None:
    """Average number of crossovers should be close to lambda."""
    rng = make_rng(99)
    counts = [
        len(simulate_crossover_breakpoints(gmap, rng))
        for _ in range(5000)
    ]
    # Default lambda = 150/100 = 1.5
    assert abs(np.mean(counts) - 1.5) < 0.1


# ---------------------------------------------------------------------------
# build_segment_plan
# ---------------------------------------------------------------------------

def test_plan_covers_full_chromosome(gmap: GeneticMap) -> None:
    rng = make_rng(10)
    bps = simulate_crossover_breakpoints(gmap, rng, lambda_override=5.0)
    plan = build_segment_plan(bps, gmap, n_samples=10, rng=rng)

    assert plan[0].cm_start == pytest.approx(gmap.start_cm)
    assert plan[-1].cm_end == pytest.approx(gmap.end_cm)


def test_plan_segments_contiguous(gmap: GeneticMap) -> None:
    rng = make_rng(11)
    bps = simulate_crossover_breakpoints(gmap, rng, lambda_override=8.0)
    plan = build_segment_plan(bps, gmap, n_samples=20, rng=rng)

    for i in range(len(plan) - 1):
        assert plan[i].cm_end == pytest.approx(plan[i + 1].cm_start)


def test_plan_no_zero_width_segments(gmap: GeneticMap) -> None:
    rng = make_rng(12)
    bps = np.array([10.0, 10.0, 50.0])  # duplicate breakpoint
    plan = build_segment_plan(bps, gmap, n_samples=5, rng=rng)
    for seg in plan:
        assert seg.cm_end > seg.cm_start


def test_plan_donor_switches_at_each_breakpoint(gmap: GeneticMap) -> None:
    """Consecutive segments should never share the same donor (when n_samples > 1)."""
    rng = make_rng(13)
    bps = simulate_crossover_breakpoints(gmap, rng, lambda_override=10.0)
    plan = build_segment_plan(bps, gmap, n_samples=50, rng=rng)

    for i in range(len(plan) - 1):
        assert plan[i].sample_idx != plan[i + 1].sample_idx, (
            f"Segment {i} and {i+1} share the same donor"
        )


def test_plan_with_no_crossovers(gmap: GeneticMap) -> None:
    """Zero crossovers → single segment covering the whole chromosome."""
    rng = make_rng(14)
    plan = build_segment_plan(np.array([]), gmap, n_samples=5, rng=rng)
    assert len(plan) == 1
    assert plan[0].cm_start == pytest.approx(gmap.start_cm)
    assert plan[0].cm_end == pytest.approx(gmap.end_cm)


def test_plan_single_sample_pool(gmap: GeneticMap) -> None:
    """With only one donor, donor index is always 0."""
    rng = make_rng(15)
    bps = np.array([25.0, 75.0])
    plan = build_segment_plan(bps, gmap, n_samples=1, rng=rng)
    for seg in plan:
        assert seg.sample_idx == 0


def test_plan_sample_indices_in_range(gmap: GeneticMap) -> None:
    rng = make_rng(16)
    n_samples = 100
    bps = simulate_crossover_breakpoints(gmap, rng, lambda_override=10.0)
    plan = build_segment_plan(bps, gmap, n_samples=n_samples, rng=rng)
    for seg in plan:
        assert 0 <= seg.sample_idx < n_samples


# ---------------------------------------------------------------------------
# generate_all_segment_plans
# ---------------------------------------------------------------------------

def test_generate_all_plans_count(gmap: GeneticMap) -> None:
    rng = make_rng(20)
    plans = generate_all_segment_plans(
        n_output_samples=10, genetic_map=gmap, n_pool_samples=50, rng=rng
    )
    assert len(plans) == 10


def test_generate_all_plans_each_valid(gmap: GeneticMap) -> None:
    rng = make_rng(21)
    plans = generate_all_segment_plans(
        n_output_samples=5, genetic_map=gmap, n_pool_samples=20, rng=rng
    )
    for plan in plans:
        assert plan[0].cm_start == pytest.approx(gmap.start_cm)
        assert plan[-1].cm_end == pytest.approx(gmap.end_cm)


def test_determinism_with_seed(gmap: GeneticMap) -> None:
    plans_a = generate_all_segment_plans(
        n_output_samples=5, genetic_map=gmap, n_pool_samples=20,
        rng=np.random.default_rng(99)
    )
    plans_b = generate_all_segment_plans(
        n_output_samples=5, genetic_map=gmap, n_pool_samples=20,
        rng=np.random.default_rng(99)
    )
    for pa, pb in zip(plans_a, plans_b):
        for sa, sb in zip(pa, pb):
            assert sa == sb


# ---------------------------------------------------------------------------
# Property-based tests (Hypothesis)
# ---------------------------------------------------------------------------

@given(
    n_bps=st.integers(min_value=0, max_value=50),
    n_samples=st.integers(min_value=1, max_value=200),
    seed=st.integers(min_value=0, max_value=2**32 - 1),
)
@settings(max_examples=200, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_plan_always_covers_chromosome(
    n_bps: int, n_samples: int, seed: int, gmap: GeneticMap
) -> None:
    rng = np.random.default_rng(seed)
    bps = np.sort(rng.uniform(gmap.start_cm, gmap.end_cm, size=n_bps))
    plan = build_segment_plan(bps, gmap, n_samples=n_samples, rng=rng)

    assert len(plan) >= 1
    assert plan[0].cm_start == pytest.approx(gmap.start_cm)
    assert plan[-1].cm_end == pytest.approx(gmap.end_cm)

    # Contiguous
    for i in range(len(plan) - 1):
        assert plan[i].cm_end == pytest.approx(plan[i + 1].cm_start)

    # No zero-width
    for seg in plan:
        assert seg.cm_end > seg.cm_start
