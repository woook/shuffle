"""Tests for v_shuffler.core.mosaic_builder."""

from __future__ import annotations

import numpy as np
import pytest

from v_shuffler.core.genotype_pool import MISSING, GenotypePool, VariantInfo
from v_shuffler.core.mosaic_builder import apply_segment_plan, build_synthetic_genotypes
from v_shuffler.core.recombination import Segment


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_pool(dosages: np.ndarray, cm_pos: np.ndarray) -> GenotypePool:
    n_var = dosages.shape[0]
    return GenotypePool(
        dosages=dosages,
        positions=np.arange(1, n_var + 1, dtype=np.int64) * 100,
        cm_pos=cm_pos,
        variant_info=[
            VariantInfo("chr1", (i + 1) * 100, "A", ["T"], ".", None, [], float(cm_pos[i]))
            for i in range(n_var)
        ],
    )


# ---------------------------------------------------------------------------
# apply_segment_plan
# ---------------------------------------------------------------------------

def test_single_segment_copies_donor() -> None:
    """Single segment → all variants come from that donor."""
    # 5 variants, 3 samples; each sample has a unique constant dosage
    dosages = np.array([
        [0, 1, 2],
        [0, 1, 2],
        [0, 1, 2],
        [0, 1, 2],
        [0, 1, 2],
    ], dtype=np.uint8)
    cm_pos = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    pool = make_pool(dosages, cm_pos)

    plan = [Segment(cm_start=0.0, cm_end=4.0, sample_idx=1)]
    result = apply_segment_plan(pool, plan)

    # All values should be dosage of sample 1 = 1
    np.testing.assert_array_equal(result, np.ones(5, dtype=np.uint8))


def test_two_segments_different_donors() -> None:
    """Two segments from two different donors."""
    # variants 0-2 in segment A (sample 0, dosage 0)
    # variants 3-4 in segment B (sample 1, dosage 2)
    dosages = np.array([
        [0, 2],
        [0, 2],
        [0, 2],
        [0, 2],
        [0, 2],
    ], dtype=np.uint8)
    cm_pos = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    pool = make_pool(dosages, cm_pos)

    plan = [
        Segment(cm_start=0.0, cm_end=2.5, sample_idx=0),
        Segment(cm_start=2.5, cm_end=4.0, sample_idx=1),
    ]
    result = apply_segment_plan(pool, plan)

    expected = np.array([0, 0, 0, 2, 2], dtype=np.uint8)
    np.testing.assert_array_equal(result, expected)


def test_all_zero_pool_gives_all_zero() -> None:
    """Constant zero pool → all-zero output."""
    dosages = np.zeros((10, 5), dtype=np.uint8)
    cm_pos = np.linspace(0.0, 9.0, 10)
    pool = make_pool(dosages, cm_pos)

    plan = [Segment(cm_start=0.0, cm_end=9.0, sample_idx=3)]
    result = apply_segment_plan(pool, plan)
    assert np.all(result == 0)


def test_missing_sentinel_preserved() -> None:
    """MISSING values in the pool are copied through to the output."""
    dosages = np.full((5, 2), MISSING, dtype=np.uint8)
    dosages[2, 0] = 1  # one non-missing call
    cm_pos = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    pool = make_pool(dosages, cm_pos)

    plan = [Segment(cm_start=0.0, cm_end=4.0, sample_idx=0)]
    result = apply_segment_plan(pool, plan)

    assert result[2] == 1
    assert result[0] == MISSING


def test_plan_at_chunk_boundary() -> None:
    """
    Variants at exact segment boundaries are covered by the segment that
    claims them (>=cm_start and <=cm_end inclusive logic).
    """
    dosages = np.array([[0, 1], [0, 1]], dtype=np.uint8)
    cm_pos = np.array([2.5, 2.5])  # both variants sit exactly on a boundary
    pool = make_pool(dosages, cm_pos)

    plan = [
        Segment(cm_start=0.0, cm_end=2.5, sample_idx=0),
        Segment(cm_start=2.5, cm_end=5.0, sample_idx=1),
    ]
    result = apply_segment_plan(pool, plan)
    # Both match both segments due to overlapping inclusive bounds;
    # the last matching segment wins (loop order).
    # The important thing is no MISSING entries.
    assert np.all(result != MISSING)


# ---------------------------------------------------------------------------
# build_synthetic_genotypes
# ---------------------------------------------------------------------------

def test_build_synthetic_genotypes_shape() -> None:
    dosages = np.zeros((8, 4), dtype=np.uint8)
    cm_pos = np.linspace(0.0, 7.0, 8)
    pool = make_pool(dosages, cm_pos)

    plans = [
        [Segment(0.0, 7.0, 0)],
        [Segment(0.0, 7.0, 1)],
        [Segment(0.0, 7.0, 2)],
    ]
    result = build_synthetic_genotypes(pool, plans)
    assert result.shape == (8, 3)


def test_build_synthetic_genotypes_correct_donors() -> None:
    """Each output column should match the assigned donor column."""
    n_var = 6
    # Sample i has constant dosage i
    dosages = np.zeros((n_var, 3), dtype=np.uint8)
    for i in range(3):
        dosages[:, i] = i
    cm_pos = np.linspace(0.0, 5.0, n_var)
    pool = make_pool(dosages, cm_pos)

    plans = [
        [Segment(0.0, 5.0, 2)],  # output 0 → donor 2 → all 2s
        [Segment(0.0, 5.0, 0)],  # output 1 → donor 0 → all 0s
        [Segment(0.0, 5.0, 1)],  # output 2 → donor 1 → all 1s
    ]
    result = build_synthetic_genotypes(pool, plans)

    np.testing.assert_array_equal(result[:, 0], np.full(n_var, 2, dtype=np.uint8))
    np.testing.assert_array_equal(result[:, 1], np.zeros(n_var, dtype=np.uint8))
    np.testing.assert_array_equal(result[:, 2], np.ones(n_var, dtype=np.uint8))


def test_build_synthetic_genotypes_mosaic() -> None:
    """
    Multi-segment plan: first half from donor 0, second half from donor 1.
    """
    n_var = 10
    dosages = np.zeros((n_var, 2), dtype=np.uint8)
    dosages[:, 0] = 0  # donor 0: all 0
    dosages[:, 1] = 2  # donor 1: all 2
    cm_pos = np.linspace(0.0, 9.0, n_var)  # 0,1,...,9
    pool = make_pool(dosages, cm_pos)

    plan = [
        Segment(cm_start=0.0, cm_end=4.5, sample_idx=0),
        Segment(cm_start=4.5, cm_end=9.0, sample_idx=1),
    ]
    result = build_synthetic_genotypes(pool, [plan])[:, 0]

    # First 5 variants (cm 0-4) → donor 0 → dosage 0
    # Last 5 variants (cm 5-9) → donor 1 → dosage 2
    assert np.all(result[:5] == 0)
    assert np.all(result[5:] == 2)
