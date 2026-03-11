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
    with pytest.raises((AssertionError, ValueError)):
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


# ---------------------------------------------------------------------------
# carry_format_fields
# ---------------------------------------------------------------------------


def test_carry_format_fields_reads_af(tmp_path: Path) -> None:
    """carry_format_fields=['AF'] populates pool.format_fields with a float matrix."""
    import textwrap
    from v_shuffler.io.vcf_reader import PerSampleVCFReader

    content = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##contig=<ID=chr22,length=51304566>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
        chr22\t100\t.\tT\tA\t.\tPASS\t.\tGT:AF\t0/1:0.42
        chr22\t200\t.\tC\tG\t.\tPASS\t.\tGT:AF\t1/1:0.91
        chr22\t300\t.\tG\tC\t.\tPASS\t.\tGT:AF\t0/1:.
    """)
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(content)

    map_path = tmp_path / "map.txt"
    map_path.write_text("pos chr cM\n1 chr22 0.0\n1000 chr22 1.0\n")
    gmap = GeneticMap(map_path, "chr22")

    reader = PerSampleVCFReader([vcf], "chr22", gmap, carry_format_fields=("AF",))
    chunks = list(reader.iter_chunks())
    assert len(chunks) == 1
    pool = chunks[0]

    assert "AF" in pool.format_fields
    af = pool.format_fields["AF"]
    assert af.shape == (3, 1)
    assert af.dtype == object
    assert af[0, 0] == "0.42"
    assert af[1, 0] == "0.91"
    assert af[2, 0] == "."   # missing → "."


def test_carry_format_fields_propagates_to_output(tmp_path: Path) -> None:
    """AF carried through build_synthetic_genotypes appears in the output matrix."""
    import textwrap
    import numpy as np
    from v_shuffler.core.genotype_pool import GenotypePool, VariantInfo, MISSING
    from v_shuffler.core.mosaic_builder import build_synthetic_genotypes
    from v_shuffler.core.recombination import Segment

    n_variants = 4
    n_donors = 2
    dosages = np.array([[1, 2], [0, 1], [2, 0], [1, 1]], dtype=np.uint8)
    af_vals = np.array([["0.45", "0.9"], ["0.0", "0.5"], ["0.88", "0.0"], ["0.3", "0.7"]],
                       dtype=object)
    cm_pos = np.array([0.0, 1.0, 2.0, 3.0])
    variant_info = [
        VariantInfo("chr22", i * 100, "A", ["G"], ".", None, [], float(i))
        for i in range(n_variants)
    ]
    pool = GenotypePool(
        dosages=dosages,
        positions=np.array([100, 200, 300, 400], dtype=np.int64),
        cm_pos=cm_pos,
        variant_info=variant_info,
        format_fields={"AF": af_vals},
    )

    # One synthetic: first half from donor 0, second half from donor 1
    plan = [
        Segment(cm_start=0.0, cm_end=1.5, sample_idx=0),
        Segment(cm_start=1.5, cm_end=3.0, sample_idx=1),
    ]
    synth_dosages, synth_fields = build_synthetic_genotypes(pool, [plan])

    assert "AF" in synth_fields
    af_out = synth_fields["AF"][:, 0]

    # Variants 0-1 (cm 0.0, 1.0) → donor 0
    assert af_out[0] == "0.45"
    assert af_out[1] == "0.0"
    # Variants 2-3 (cm 2.0, 3.0) → donor 1
    assert af_out[2] == "0.0"
    assert af_out[3] == "0.7"


def test_no_format_fields_returns_empty_dict(tmp_path: Path) -> None:
    """Without carry_format_fields, build_synthetic_genotypes returns an empty dict."""
    import numpy as np
    from v_shuffler.core.genotype_pool import GenotypePool, VariantInfo
    from v_shuffler.core.mosaic_builder import build_synthetic_genotypes
    from v_shuffler.core.recombination import Segment

    dosages = np.array([[1], [0]], dtype=np.uint8)
    pool = GenotypePool(
        dosages=dosages,
        positions=np.array([100, 200], dtype=np.int64),
        cm_pos=np.array([0.0, 1.0]),
        variant_info=[
            VariantInfo("chr22", 100, "A", ["G"], ".", None, [], 0.0),
            VariantInfo("chr22", 200, "C", ["T"], ".", None, [], 1.0),
        ],
    )
    plan = [Segment(cm_start=0.0, cm_end=1.0, sample_idx=0)]
    dosage_mat, fields = build_synthetic_genotypes(pool, [plan])
    assert fields == {}
    assert dosage_mat.shape == (2, 1)


# ---------------------------------------------------------------------------
# resolve_chromosome_name
# ---------------------------------------------------------------------------

def _make_vcf_with_contig(tmp_path: Path, contig_id: str, chrom_in_data: str) -> Path:
    """Write a minimal single-sample VCF whose ##contig line and data use *contig_id*."""
    content = textwrap.dedent(f"""\
        ##fileformat=VCFv4.1
        ##contig=<ID={contig_id},length=51304566>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
        {chrom_in_data}\t100\t.\tT\tA\t.\tPASS\t.\tGT\t0/1
    """)
    p = tmp_path / f"vcf_{contig_id}.vcf"
    p.write_text(content)
    return p


def test_resolve_chromosome_name_no_change(tmp_path: Path) -> None:
    """User passes chr22, VCF uses chr22 → unchanged."""
    from v_shuffler.io.vcf_reader import resolve_chromosome_name
    vcf = _make_vcf_with_contig(tmp_path, "chr22", "chr22")
    assert resolve_chromosome_name(vcf, "chr22") == "chr22"


def test_resolve_chromosome_name_adds_prefix(tmp_path: Path) -> None:
    """User passes 22, VCF uses chr22 → normalised to chr22."""
    from v_shuffler.io.vcf_reader import resolve_chromosome_name
    vcf = _make_vcf_with_contig(tmp_path, "chr22", "chr22")
    assert resolve_chromosome_name(vcf, "22") == "chr22"


def test_resolve_chromosome_name_strips_prefix(tmp_path: Path) -> None:
    """User passes chr22, VCF uses 22 → normalised to 22."""
    from v_shuffler.io.vcf_reader import resolve_chromosome_name
    vcf = _make_vcf_with_contig(tmp_path, "22", "22")
    assert resolve_chromosome_name(vcf, "chr22") == "22"


def test_resolve_chromosome_name_no_contig_header(tmp_path: Path) -> None:
    """VCF has no ##contig lines → original name returned unchanged."""
    from v_shuffler.io.vcf_reader import resolve_chromosome_name
    content = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
        chr22\t100\t.\tT\tA\t.\tPASS\t.\tGT\t0/1
    """)
    vcf = tmp_path / "no_contig.vcf"
    vcf.write_text(content)
    assert resolve_chromosome_name(vcf, "chr22") == "chr22"
    assert resolve_chromosome_name(vcf, "22") == "22"


def test_reader_accepts_bare_chromosome_name(tmp_path: Path) -> None:
    """Reader with chromosome='22' reads a VCF whose CHROM column says 'chr22'."""
    from v_shuffler.io.vcf_reader import PerSampleVCFReader

    # VCF uses chr22 in both the ##contig header and data rows.
    content = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##contig=<ID=chr22,length=51304566>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
        chr22\t100\t.\tT\tA\t.\tPASS\t.\tGT\t0/1
        chr22\t200\t.\tC\tG\t.\tPASS\t.\tGT\t1/1
        chr22\t300\t.\tG\tC\t.\tPASS\t.\tGT\t0/0
    """)
    vcf = tmp_path / "chr22.vcf"
    vcf.write_text(content)

    map_text = "pos chr cM\n1 chr22 0.0\n1000 chr22 1.0\n"
    gmap_path = tmp_path / "map.txt"
    gmap_path.write_text(map_text)
    gmap = GeneticMap(gmap_path, "chr22")

    # Pass the bare form; the plain-VCF fallback resolves both forms.
    reader = PerSampleVCFReader([vcf], chromosome="22", genetic_map=gmap)
    chunks = list(reader.iter_chunks())
    assert sum(p.n_variants for p in chunks) == 3


def test_reader_accepts_prefixed_chromosome_name(tmp_path: Path) -> None:
    """Reader with chromosome='chr22' reads a VCF whose CHROM column says '22'."""
    from v_shuffler.io.vcf_reader import PerSampleVCFReader

    content = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##contig=<ID=22,length=51304566>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
        22\t100\t.\tT\tA\t.\tPASS\t.\tGT\t0/1
        22\t200\t.\tC\tG\t.\tPASS\t.\tGT\t1/1
        22\t300\t.\tG\tC\t.\tPASS\t.\tGT\t0/0
    """)
    vcf = tmp_path / "22.vcf"
    vcf.write_text(content)

    map_text = "pos chr cM\n1 22 0.0\n1000 22 1.0\n"
    gmap_path = tmp_path / "map.txt"
    gmap_path.write_text(map_text)
    gmap = GeneticMap(gmap_path, "22")

    # Pass the prefixed form; the plain-VCF fallback resolves both forms.
    reader = PerSampleVCFReader([vcf], chromosome="chr22", genetic_map=gmap)
    chunks = list(reader.iter_chunks())
    assert sum(p.n_variants for p in chunks) == 3


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


def test_iter_positions_respects_missing_rate_filter(tmp_path: Path, gmap: GeneticMap) -> None:
    """iter_positions() filters high-missing variants like iter_chunks()."""
    from v_shuffler.io.vcf_reader import PerSampleVCFReader

    # Create 3 samples with 5 variants
    # Variant at position 300: 2/3 samples missing (66.7% missing rate)
    # Should be filtered with max_missing_rate=0.5
    vcf0_content = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##contig=<ID=chr22,length=51304566>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE0
        chr22\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
        chr22\t200\t.\tC\tT\t.\tPASS\t.\tGT\t1/1
        chr22\t300\t.\tG\tA\t.\tPASS\t.\tGT\t./.
        chr22\t400\t.\tT\tC\t.\tPASS\t.\tGT\t0/1
        chr22\t500\t.\tA\tT\t.\tPASS\t.\tGT\t1/1
    """)

    vcf1_content = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##contig=<ID=chr22,length=51304566>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
        chr22\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0
        chr22\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/1
        chr22\t300\t.\tG\tA\t.\tPASS\t.\tGT\t./.
        chr22\t400\t.\tT\tC\t.\tPASS\t.\tGT\t1/1
        chr22\t500\t.\tA\tT\t.\tPASS\t.\tGT\t0/1
    """)

    vcf2_content = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##contig=<ID=chr22,length=51304566>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE2
        chr22\t100\t.\tA\tG\t.\tPASS\t.\tGT\t1/1
        chr22\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/0
        chr22\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/1
        chr22\t400\t.\tT\tC\t.\tPASS\t.\tGT\t0/0
        chr22\t500\t.\tA\tT\t.\tPASS\t.\tGT\t0/0
    """)

    vcf0 = tmp_path / "s0.vcf"
    vcf1 = tmp_path / "s1.vcf"
    vcf2 = tmp_path / "s2.vcf"
    vcf0.write_text(vcf0_content)
    vcf1.write_text(vcf1_content)
    vcf2.write_text(vcf2_content)

    reader = PerSampleVCFReader(
        vcf_paths=[vcf0, vcf1, vcf2],
        chromosome="chr22",
        genetic_map=gmap,
        max_missing_rate=0.5,  # 50% threshold
    )

    # iter_positions() should filter position 300 (66.7% missing)
    positions = reader.iter_positions()
    assert list(positions) == [100, 200, 400, 500]


def test_iter_positions_matches_iter_chunks_filtering(tmp_path: Path, gmap: GeneticMap) -> None:
    """Positions from iter_positions() exactly match those in iter_chunks()."""
    from v_shuffler.io.vcf_reader import PerSampleVCFReader

    # Create 2 samples with varying missing rates
    vcf0_content = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##contig=<ID=chr22,length=51304566>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE0
        chr22\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
        chr22\t200\t.\tC\tT\t.\tPASS\t.\tGT\t./.
        chr22\t300\t.\tG\tA\t.\tPASS\t.\tGT\t1/1
        chr22\t400\t.\tT\tC\t.\tPASS\t.\tGT\t0/1
    """)

    vcf1_content = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##contig=<ID=chr22,length=51304566>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
        chr22\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0
        chr22\t200\t.\tC\tT\t.\tPASS\t.\tGT\t./.
        chr22\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/1
        chr22\t400\t.\tT\tC\t.\tPASS\t.\tGT\t1/1
    """)

    vcf0 = tmp_path / "s0.vcf"
    vcf1 = tmp_path / "s1.vcf"
    vcf0.write_text(vcf0_content)
    vcf1.write_text(vcf1_content)

    reader = PerSampleVCFReader(
        vcf_paths=[vcf0, vcf1],
        chromosome="chr22",
        genetic_map=gmap,
        max_missing_rate=0.3,
    )

    # Get positions from both methods
    positions_from_iter = reader.iter_positions()
    chunks = list(reader.iter_chunks())
    positions_from_chunks = np.concatenate([chunk.positions for chunk in chunks])

    # Must match exactly (variant 200 filtered: 100% missing)
    assert np.array_equal(positions_from_iter, positions_from_chunks)
    assert list(positions_from_iter) == [100, 300, 400]
