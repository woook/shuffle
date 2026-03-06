"""Tests for v_shuffler.io.vcf_reader and v_shuffler.core.genotype_pool."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest

from v_shuffler.core.genotype_pool import MISSING, GenotypePool, VariantInfo
from v_shuffler.io.genetic_map import GeneticMap


# ---------------------------------------------------------------------------
# Helpers: build minimal VCF content
# ---------------------------------------------------------------------------

VCF_HEADER = textwrap.dedent("""\
    ##fileformat=VCFv4.1
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##contig=<ID=chr22,length=51304566>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
""")

# 5 variants, one sample
def make_vcf_lines(sample: str, genotypes: list[str]) -> str:
    """genotypes: list of GT strings like '0/0', '0/1', '1/1'"""
    positions = [100, 200, 300, 400, 500]
    alts = ["A", "G", "C", "T", "A"]
    refs = ["T", "C", "G", "A", "G"]
    header = VCF_HEADER.format(sample=sample)
    lines = []
    for pos, ref, alt, gt in zip(positions, refs, alts, genotypes):
        lines.append(
            f"chr22\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}"
        )
    return header + "\n".join(lines) + "\n"


MAP_TEXT = textwrap.dedent("""\
    pos chr cM
    1 chr22 0.0
    1000 chr22 1.0
""")


@pytest.fixture
def gmap(tmp_path: Path) -> GeneticMap:
    p = tmp_path / "map.txt"
    p.write_text(MAP_TEXT)
    return GeneticMap(p, "chr22")


def write_vcfs(tmp_path: Path, all_gts: list[list[str]]) -> list[Path]:
    """Write one VCF per sample and return paths."""
    paths = []
    for i, gts in enumerate(all_gts):
        content = make_vcf_lines(f"sample{i}", gts)
        p = tmp_path / f"sample{i}.vcf"
        p.write_text(content)
        paths.append(p)
    return paths


# ---------------------------------------------------------------------------
# GenotypePool unit tests (no VCF needed)
# ---------------------------------------------------------------------------

def test_genotype_pool_shape() -> None:
    n_var, n_samp = 10, 5
    pool = GenotypePool(
        dosages=np.zeros((n_var, n_samp), dtype=np.uint8),
        positions=np.arange(n_var, dtype=np.int64),
        cm_pos=np.linspace(0, 1, n_var),
        variant_info=[
            VariantInfo("chr22", i, "A", ["T"], ".", None, [], 0.0)
            for i in range(n_var)
        ],
    )
    assert pool.n_variants == n_var
    assert pool.n_samples == n_samp


def test_genotype_pool_dimension_mismatch_raises() -> None:
    with pytest.raises(AssertionError):
        GenotypePool(
            dosages=np.zeros((5, 3), dtype=np.uint8),
            positions=np.zeros(4, dtype=np.int64),   # wrong length
            cm_pos=np.zeros(4),
            variant_info=[
                VariantInfo("chr22", i, "A", ["T"], ".", None, [], 0.0)
                for i in range(4)
            ],
        )


# ---------------------------------------------------------------------------
# PerSampleVCFReader integration tests
# ---------------------------------------------------------------------------

def test_reader_basic(tmp_path: Path, gmap: GeneticMap) -> None:
    """Reader yields correct dosages for a simple 3-sample case."""
    from v_shuffler.io.vcf_reader import PerSampleVCFReader

    all_gts = [
        ["0/0", "0/1", "1/1", "0/0", "0/1"],  # sample 0
        ["0/1", "1/1", "0/0", "0/1", "0/0"],  # sample 1
        ["1/1", "0/0", "0/1", "1/1", "0/1"],  # sample 2
    ]
    paths = write_vcfs(tmp_path, all_gts)

    reader = PerSampleVCFReader(
        vcf_paths=paths,
        chromosome="chr22",
        genetic_map=gmap,
        chunk_size=10,
        max_missing_rate=0.05,
    )
    chunks = list(reader.iter_chunks())
    assert len(chunks) == 1
    pool = chunks[0]

    assert pool.n_variants == 5
    assert pool.n_samples == 3

    # Check dosage for variant 0: 0/0=0, 0/1=1, 1/1=2
    expected_v0 = np.array([0, 1, 2], dtype=np.uint8)
    np.testing.assert_array_equal(pool.dosages[0], expected_v0)


def test_reader_chunk_splitting(tmp_path: Path, gmap: GeneticMap) -> None:
    """With chunk_size=2, 5 variants → 3 chunks (2+2+1)."""
    from v_shuffler.io.vcf_reader import PerSampleVCFReader

    all_gts = [["0/0", "0/1", "1/1", "0/0", "0/1"]]
    paths = write_vcfs(tmp_path, all_gts)

    reader = PerSampleVCFReader(
        vcf_paths=paths,
        chromosome="chr22",
        genetic_map=gmap,
        chunk_size=2,
    )
    chunks = list(reader.iter_chunks())
    assert len(chunks) == 3
    assert chunks[0].n_variants == 2
    assert chunks[1].n_variants == 2
    assert chunks[2].n_variants == 1


def test_reader_missing_genotype(tmp_path: Path, gmap: GeneticMap) -> None:
    """Missing genotype '.' is encoded as MISSING=255."""
    from v_shuffler.io.vcf_reader import PerSampleVCFReader

    # Sample 0 has a missing call at variant index 2
    all_gts = [
        ["0/0", "0/1", "./.", "0/0", "0/1"],
        ["0/1", "1/1", "0/0", "0/1", "0/0"],
    ]
    paths = write_vcfs(tmp_path, all_gts)

    reader = PerSampleVCFReader(
        vcf_paths=paths,
        chromosome="chr22",
        genetic_map=gmap,
        chunk_size=10,
        max_missing_rate=1.0,  # don't filter any missing
    )
    chunks = list(reader.iter_chunks())
    pool = chunks[0]
    # variant index 2, sample 0 should be MISSING
    assert pool.dosages[2, 0] == MISSING


def test_reader_missing_rate_filter(tmp_path: Path, gmap: GeneticMap) -> None:
    """Variants with >50% missing are dropped when max_missing_rate=0.5."""
    from v_shuffler.io.vcf_reader import PerSampleVCFReader

    # 2 samples; variant index 1 has both missing → 100% missing → dropped
    all_gts = [
        ["0/0", "./.", "1/1", "0/0", "0/1"],
        ["0/1", "./.", "0/0", "0/1", "0/0"],
    ]
    paths = write_vcfs(tmp_path, all_gts)

    reader = PerSampleVCFReader(
        vcf_paths=paths,
        chromosome="chr22",
        genetic_map=gmap,
        chunk_size=10,
        max_missing_rate=0.5,
    )
    chunks = list(reader.iter_chunks())
    pool = chunks[0]
    # 5 variants - 1 dropped = 4
    assert pool.n_variants == 4
    # POS 200 should be absent
    assert 200 not in pool.positions


def test_reader_site_mismatch_raises(tmp_path: Path, gmap: GeneticMap) -> None:
    """Site mismatch between VCFs raises ValueError."""
    from v_shuffler.io.vcf_reader import PerSampleVCFReader

    # sample0 has normal sites; sample1 has a different POS at variant 0
    vcf0 = make_vcf_lines("sample0", ["0/0", "0/1", "1/1", "0/0", "0/1"])
    # Manually swap the first variant position to something different
    vcf1_header = VCF_HEADER.format(sample="sample1")
    vcf1_body = "\n".join([
        "chr22\t999\t.\tT\tA\t.\tPASS\t.\tGT\t0/1",  # different POS
        "chr22\t200\t.\tC\tG\t.\tPASS\t.\tGT\t1/1",
        "chr22\t300\t.\tG\tC\t.\tPASS\t.\tGT\t0/0",
        "chr22\t400\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
        "chr22\t500\t.\tG\tA\t.\tPASS\t.\tGT\t0/0",
    ]) + "\n"

    p0 = tmp_path / "s0.vcf"
    p1 = tmp_path / "s1.vcf"
    p0.write_text(vcf0)
    p1.write_text(vcf1_header + vcf1_body)

    reader = PerSampleVCFReader(
        vcf_paths=[p0, p1],
        chromosome="chr22",
        genetic_map=gmap,
        chunk_size=10,
    )
    with pytest.raises(ValueError, match="Site mismatch"):
        list(reader.iter_chunks())
