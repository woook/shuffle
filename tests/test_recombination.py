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
    build_region_segment_plan,
    detect_regions,
    generate_all_region_plans,
    generate_all_segment_plans,
    simulate_crossover_breakpoints,
)
from v_shuffler.core.genotype_pool import MISSING, GenotypePool, VariantInfo
from v_shuffler.core.mosaic_builder import build_synthetic_genotypes


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
# TestRegionSampling
# ---------------------------------------------------------------------------


class TestRegionSampling:
    """Tests for detect_regions, build_region_segment_plan, generate_all_region_plans."""

    # --- detect_regions ---

    def test_detect_regions_empty(self) -> None:
        """Empty input → empty list."""
        assert detect_regions(np.array([], dtype=np.int64)) == []

    def test_detect_regions_single(self) -> None:
        """Single position → one region spanning that position."""
        result = detect_regions(np.array([5000], dtype=np.int64))
        assert result == [(5000, 5000)]

    def test_detect_regions_all_within_threshold(self) -> None:
        """All positions within the gap threshold → single region."""
        positions = np.array([1000, 2000, 3000, 4000], dtype=np.int64)
        result = detect_regions(positions, gap_threshold_bp=10_000)
        assert result == [(1000, 4000)]

    def test_detect_regions_multiple_gaps(self) -> None:
        """Three clusters separated by large gaps → three regions."""
        positions = np.array(
            [1000, 1500, 2000,
             100_000, 100_500,
             500_000, 500_100, 500_200],
            dtype=np.int64,
        )
        result = detect_regions(positions, gap_threshold_bp=10_000)
        assert result == [(1000, 2000), (100_000, 100_500), (500_000, 500_200)]

    # --- build_region_segment_plan ---

    def test_region_plan_segment_count(self, gmap: GeneticMap) -> None:
        """Number of segments equals number of regions."""
        regions_cm = [(0.0, 10.0), (20.0, 30.0), (50.0, 60.0)]
        plan = build_region_segment_plan(regions_cm, n_samples=10, rng=make_rng(1))
        assert len(plan) == 3

    def test_region_plan_sample_indices_in_range(self, gmap: GeneticMap) -> None:
        """All sample_idx in [0, n_samples)."""
        n_samples = 20
        regions_cm = [(float(i * 10), float(i * 10 + 5)) for i in range(10)]
        plan = build_region_segment_plan(regions_cm, n_samples=n_samples, rng=make_rng(2))
        for seg in plan:
            assert 0 <= seg.sample_idx < n_samples

    def test_region_plan_adjacency_constraint(self, gmap: GeneticMap) -> None:
        """No two consecutive segments share the same donor (n_samples > 1)."""
        regions_cm = [(float(i * 5), float(i * 5 + 4)) for i in range(20)]
        plan = build_region_segment_plan(regions_cm, n_samples=10, rng=make_rng(3))
        for i in range(len(plan) - 1):
            assert plan[i].sample_idx != plan[i + 1].sample_idx

    def test_region_plan_min_donors_satisfied(self, gmap: GeneticMap) -> None:
        """Plan has at least min_donors distinct donors."""
        regions_cm = [(float(i * 5), float(i * 5 + 4)) for i in range(20)]
        min_donors = 5
        plan = build_region_segment_plan(
            regions_cm, n_samples=50, rng=make_rng(4), min_donors=min_donors
        )
        assert len({seg.sample_idx for seg in plan}) >= min_donors

    def test_region_plan_min_donors_relaxed_when_few_regions(self) -> None:
        """min_donors gracefully relaxed when n_regions < min_donors."""
        regions_cm = [(0.0, 5.0), (10.0, 15.0)]  # only 2 regions
        plan = build_region_segment_plan(
            regions_cm, n_samples=50, rng=make_rng(5), min_donors=10
        )
        # Should succeed without error; distinct donors ≤ 2
        assert len(plan) == 2

    def test_region_plan_single_sample_all_zero(self) -> None:
        """n_samples == 1 → all segments get sample_idx=0."""
        regions_cm = [(0.0, 5.0), (10.0, 15.0), (20.0, 25.0)]
        plan = build_region_segment_plan(regions_cm, n_samples=1, rng=make_rng(6))
        assert all(seg.sample_idx == 0 for seg in plan)

    # --- generate_all_region_plans ---

    def test_generate_all_region_plans_count(self, gmap: GeneticMap) -> None:
        """Output list length equals n_output_samples."""
        regions_cm = [(float(i * 10), float(i * 10 + 8)) for i in range(5)]
        plans = generate_all_region_plans(
            n_output_samples=7, regions_cm=regions_cm, n_pool_samples=30, rng=make_rng(10)
        )
        assert len(plans) == 7

    def test_generate_all_region_plans_deterministic(self, gmap: GeneticMap) -> None:
        """Same seed → identical plans."""
        regions_cm = [(float(i * 10), float(i * 10 + 8)) for i in range(5)]
        plans_a = generate_all_region_plans(
            n_output_samples=4, regions_cm=regions_cm, n_pool_samples=20,
            rng=np.random.default_rng(77),
        )
        plans_b = generate_all_region_plans(
            n_output_samples=4, regions_cm=regions_cm, n_pool_samples=20,
            rng=np.random.default_rng(77),
        )
        for pa, pb in zip(plans_a, plans_b):
            for sa, sb in zip(pa, pb):
                assert sa == sb

    # --- Continuous-mode min_donors ---

    def test_continuous_min_donors_enforced(self, gmap: GeneticMap) -> None:
        """lambda=0 but min_donors=3 → plan has ≥ 3 distinct donors."""
        rng = make_rng(20)
        plans = generate_all_segment_plans(
            n_output_samples=20, genetic_map=gmap, n_pool_samples=10,
            rng=rng, lambda_override=0.0, min_donors=3,
        )
        for plan in plans:
            assert len({seg.sample_idx for seg in plan}) >= 3

    def test_continuous_min_donors_capped_by_pool(self, gmap: GeneticMap) -> None:
        """min_donors=100 with n_pool_samples=3 → no error, ≤ 3 distinct donors."""
        rng = make_rng(21)
        plans = generate_all_segment_plans(
            n_output_samples=10, genetic_map=gmap, n_pool_samples=3,
            rng=rng, lambda_override=0.0, min_donors=100,
        )
        for plan in plans:
            distinct = len({seg.sample_idx for seg in plan})
            assert distinct <= 3

    # --- Integration test ---

    def test_region_plan_end_to_end_no_missing(self, gmap: GeneticMap) -> None:
        """Panel-like fixture: region plans + build_synthetic_genotypes → no MISSING."""
        # Three variant clusters separated by large gaps
        cluster1 = np.arange(1_000_000, 1_000_300, 100, dtype=np.int64)   # 3 variants
        cluster2 = np.arange(5_000_000, 5_000_300, 100, dtype=np.int64)   # 3 variants
        cluster3 = np.arange(10_000_000, 10_000_300, 100, dtype=np.int64) # 3 variants
        positions = np.concatenate([cluster1, cluster2, cluster3])
        cm_pos = gmap.bp_to_cm(positions)

        n_variants = len(positions)
        n_donors = 5
        rng = make_rng(30)
        dosages = rng.integers(0, 3, size=(n_variants, n_donors), dtype=np.uint8)

        variant_info = [
            VariantInfo(
                chrom="chr1", pos=int(positions[i]), ref="A", alts=["G"],
                id=".", qual=None, filters=[], cm_pos=float(cm_pos[i]),
            )
            for i in range(n_variants)
        ]
        pool = GenotypePool(
            dosages=dosages,
            positions=positions,
            cm_pos=cm_pos,
            variant_info=variant_info,
        )

        regions_bp = detect_regions(positions, gap_threshold_bp=10_000)
        assert len(regions_bp) == 3

        bp_flat = np.array([[r[0], r[1]] for r in regions_bp], dtype=np.int64).ravel()
        cm_flat = gmap.bp_to_cm(bp_flat).reshape(-1, 2)
        regions_cm = [(float(row[0]), float(row[1])) for row in cm_flat]

        plans = generate_all_region_plans(
            n_output_samples=3, regions_cm=regions_cm, n_pool_samples=n_donors, rng=make_rng(31)
        )
        synth, _ = build_synthetic_genotypes(pool, plans)
        assert synth.shape == (n_variants, 3)
        assert not np.any(synth == MISSING)


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
