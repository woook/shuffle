"""Tests for v_shuffler.io.genetic_map."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest

from v_shuffler.io.genetic_map import GeneticMap


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SHAPEIT5_MAP = textwrap.dedent("""\
    pos chr cM
    100 chr22 0.0
    500 chr22 0.4
    1000 chr22 1.0
    5000 chr22 2.5
    10000 chr22 4.0
""")

HAPMAP_MAP = textwrap.dedent("""\
    Chromosome Position(bp) Rate(cM/Mb) Map(cM)
    chr22 100 0.0 0.0
    chr22 500 1000.0 0.4
    chr22 1000 1200.0 1.0
    chr22 5000 375.0 2.5
    chr22 10000 300.0 4.0
""")


@pytest.fixture
def shapeit5_map_file(tmp_path: Path) -> Path:
    p = tmp_path / "map_shapeit5.txt"
    p.write_text(SHAPEIT5_MAP)
    return p


@pytest.fixture
def hapmap_map_file(tmp_path: Path) -> Path:
    p = tmp_path / "map_hapmap.txt"
    p.write_text(HAPMAP_MAP)
    return p


# ---------------------------------------------------------------------------
# Format detection and loading
# ---------------------------------------------------------------------------

def test_load_shapeit5_format(shapeit5_map_file: Path) -> None:
    gmap = GeneticMap(shapeit5_map_file, "chr22")
    assert len(gmap.positions) == 5
    assert gmap.positions[0] == 100
    assert gmap.positions[-1] == 10000
    np.testing.assert_allclose(gmap.cm_values[0], 0.0)
    np.testing.assert_allclose(gmap.cm_values[-1], 4.0)


def test_load_hapmap_format(hapmap_map_file: Path) -> None:
    gmap = GeneticMap(hapmap_map_file, "chr22")
    assert len(gmap.positions) == 5
    np.testing.assert_allclose(gmap.cm_values[-1], 4.0)


def test_chrom_without_prefix(shapeit5_map_file: Path) -> None:
    """Should match chr22 even when caller passes '22'."""
    gmap = GeneticMap(shapeit5_map_file, "22")
    assert len(gmap.positions) == 5


def test_wrong_chrom_raises(shapeit5_map_file: Path) -> None:
    with pytest.raises(ValueError, match="No entries found"):
        GeneticMap(shapeit5_map_file, "chr1")


# ---------------------------------------------------------------------------
# bp_to_cm
# ---------------------------------------------------------------------------

def test_bp_to_cm_exact_points(shapeit5_map_file: Path) -> None:
    gmap = GeneticMap(shapeit5_map_file, "chr22")
    positions = np.array([100, 500, 1000, 5000, 10000])
    expected = np.array([0.0, 0.4, 1.0, 2.5, 4.0])
    np.testing.assert_allclose(gmap.bp_to_cm(positions), expected)


def test_bp_to_cm_interpolated(shapeit5_map_file: Path) -> None:
    gmap = GeneticMap(shapeit5_map_file, "chr22")
    # Midpoint between 100 (0.0 cM) and 500 (0.4 cM) => pos 300 => ~0.2 cM
    result = gmap.bp_to_cm(np.array([300]))
    assert 0.19 < result[0] < 0.21


def test_bp_to_cm_monotonic(shapeit5_map_file: Path) -> None:
    gmap = GeneticMap(shapeit5_map_file, "chr22")
    positions = np.linspace(100, 10000, 200).astype(np.int64)
    cm = gmap.bp_to_cm(positions)
    assert np.all(np.diff(cm) >= 0), "bp_to_cm output is not monotonically non-decreasing"


# ---------------------------------------------------------------------------
# total_length_cm
# ---------------------------------------------------------------------------

def test_total_length_cm(shapeit5_map_file: Path) -> None:
    gmap = GeneticMap(shapeit5_map_file, "chr22")
    np.testing.assert_allclose(gmap.total_length_cm, 4.0)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_positions_outside_map_clamped(shapeit5_map_file: Path) -> None:
    """Positions before/after map should return boundary cM values."""
    gmap = GeneticMap(shapeit5_map_file, "chr22")
    before = gmap.bp_to_cm(np.array([1]))    # before first entry
    after = gmap.bp_to_cm(np.array([99999])) # after last entry
    np.testing.assert_allclose(before, [0.0])
    np.testing.assert_allclose(after, [4.0])


def test_non_monotonic_cm_raises(tmp_path: Path) -> None:
    bad = textwrap.dedent("""\
        pos chr cM
        100 chr22 0.0
        500 chr22 2.0
        1000 chr22 1.0
    """)
    p = tmp_path / "bad_map.txt"
    p.write_text(bad)
    with pytest.raises(ValueError, match="non-monotonic"):
        GeneticMap(p, "chr22")


def test_gzipped_map(tmp_path: Path) -> None:
    import gzip
    p = tmp_path / "map.txt.gz"
    with gzip.open(p, "wt") as f:
        f.write(SHAPEIT5_MAP)
    gmap = GeneticMap(p, "chr22")
    assert len(gmap.positions) == 5
