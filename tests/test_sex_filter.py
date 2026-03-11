"""Tests for v_shuffler.io.sex_map and sex-chromosome donor filtering."""

from __future__ import annotations

import textwrap
import warnings
from pathlib import Path

import pytest
from click.testing import CliRunner

from v_shuffler.cli import main
from v_shuffler.io.sex_map import (
    filter_vcfs_by_sex,
    load_sex_map,
    parse_sex_label,
    sex_filter_for_chromosome,
)


# ---------------------------------------------------------------------------
# parse_sex_label
# ---------------------------------------------------------------------------

class TestParseSexLabel:
    @pytest.mark.parametrize("label", ["F", "f", "female", "FEMALE", "2"])
    def test_female_labels(self, label: str) -> None:
        assert parse_sex_label(label) == "F"

    @pytest.mark.parametrize("label", ["M", "m", "male", "MALE", "1"])
    def test_male_labels(self, label: str) -> None:
        assert parse_sex_label(label) == "M"

    @pytest.mark.parametrize("label", ["X", "unknown", "0", "3", ""])
    def test_invalid_labels_raise(self, label: str) -> None:
        with pytest.raises(ValueError, match="Unknown sex label"):
            parse_sex_label(label)


# ---------------------------------------------------------------------------
# sex_filter_for_chromosome
# ---------------------------------------------------------------------------

class TestSexFilterForChromosome:
    @pytest.mark.parametrize("chrom", ["chrX", "X", "CHRX", "chrx"])
    def test_chrx_returns_female(self, chrom: str) -> None:
        assert sex_filter_for_chromosome(chrom) == "F"

    @pytest.mark.parametrize("chrom", ["chrY", "Y", "CHRY", "chry"])
    def test_chry_returns_male(self, chrom: str) -> None:
        assert sex_filter_for_chromosome(chrom) == "M"

    @pytest.mark.parametrize("chrom", ["chr22", "22", "chr1", "1", "chrM", "MT"])
    def test_autosomes_return_none(self, chrom: str) -> None:
        assert sex_filter_for_chromosome(chrom) is None


# ---------------------------------------------------------------------------
# load_sex_map
# ---------------------------------------------------------------------------

class TestLoadSexMap:
    def _make_vcfs(self, tmp_path: Path, names: list[str]) -> list[Path]:
        return [tmp_path / name for name in names]

    def test_basic_full_path_match(self, tmp_path: Path) -> None:
        vcfs = self._make_vcfs(tmp_path, ["a.vcf.gz", "b.vcf.gz"])
        sex_file = tmp_path / "sex.txt"
        sex_file.write_text(f"{vcfs[0]}  F\n{vcfs[1]}  M\n")
        result = load_sex_map(sex_file, vcfs)
        assert result == {vcfs[0]: "F", vcfs[1]: "M"}

    def test_basename_match_fallback(self, tmp_path: Path) -> None:
        vcfs = self._make_vcfs(tmp_path, ["sample1.vcf.gz", "sample2.vcf.gz"])
        sex_file = tmp_path / "sex.txt"
        # Use only basenames in the file
        sex_file.write_text("sample1.vcf.gz  F\nsample2.vcf.gz  M\n")
        result = load_sex_map(sex_file, vcfs)
        assert result == {vcfs[0]: "F", vcfs[1]: "M"}

    def test_header_line_skipped(self, tmp_path: Path) -> None:
        vcfs = self._make_vcfs(tmp_path, ["s1.vcf.gz"])
        sex_file = tmp_path / "sex.txt"
        sex_file.write_text(f"path sex\n{vcfs[0]}  F\n")
        result = load_sex_map(sex_file, vcfs)
        assert result == {vcfs[0]: "F"}

    def test_comment_lines_ignored(self, tmp_path: Path) -> None:
        vcfs = self._make_vcfs(tmp_path, ["s1.vcf.gz"])
        sex_file = tmp_path / "sex.txt"
        sex_file.write_text(f"# this is a comment\n{vcfs[0]}  F\n")
        result = load_sex_map(sex_file, vcfs)
        assert result == {vcfs[0]: "F"}

    def test_blank_lines_ignored(self, tmp_path: Path) -> None:
        vcfs = self._make_vcfs(tmp_path, ["s1.vcf.gz"])
        sex_file = tmp_path / "sex.txt"
        sex_file.write_text(f"\n{vcfs[0]}  F\n\n")
        result = load_sex_map(sex_file, vcfs)
        assert result == {vcfs[0]: "F"}

    def test_unmatched_vcf_warns_and_omits(self, tmp_path: Path) -> None:
        vcfs = self._make_vcfs(tmp_path, ["s1.vcf.gz", "s2.vcf.gz"])
        sex_file = tmp_path / "sex.txt"
        sex_file.write_text(f"{vcfs[0]}  F\n")  # s2 not listed
        result = load_sex_map(sex_file, vcfs)  # should not raise
        assert vcfs[0] in result
        assert vcfs[1] not in result

    def test_malformed_line_raises(self, tmp_path: Path) -> None:
        vcfs = self._make_vcfs(tmp_path, ["s1.vcf.gz"])
        sex_file = tmp_path / "sex.txt"
        sex_file.write_text("only_one_column\n")
        with pytest.raises(ValueError, match="expected at least 2 columns"):
            load_sex_map(sex_file, vcfs)

    def test_invalid_sex_label_raises(self, tmp_path: Path) -> None:
        # Put the invalid label on line 2 so line 1 is not silently treated as a header.
        vcfs = self._make_vcfs(tmp_path, ["s1.vcf.gz", "s2.vcf.gz"])
        sex_file = tmp_path / "sex.txt"
        sex_file.write_text(f"{vcfs[0]}  F\n{vcfs[1]}  X\n")
        with pytest.raises(ValueError):
            load_sex_map(sex_file, vcfs)

    def test_all_sex_labels_accepted(self, tmp_path: Path) -> None:
        labels = [("F", "F"), ("female", "F"), ("2", "F"),
                  ("M", "M"), ("male", "M"), ("1", "M")]
        for raw, expected in labels:
            vcfs = [tmp_path / "s.vcf.gz"]
            sex_file = tmp_path / f"sex_{raw}.txt"
            sex_file.write_text(f"{vcfs[0]}  {raw}\n")
            result = load_sex_map(sex_file, vcfs)
            assert result[vcfs[0]] == expected, f"label {raw!r} should map to {expected}"


# ---------------------------------------------------------------------------
# filter_vcfs_by_sex
# ---------------------------------------------------------------------------

class TestFilterVcfsBySex:
    def test_keeps_only_matching_sex(self, tmp_path: Path) -> None:
        vcfs = [tmp_path / f"s{i}.vcf.gz" for i in range(4)]
        sex_map = {vcfs[0]: "F", vcfs[1]: "M", vcfs[2]: "F", vcfs[3]: "M"}
        assert filter_vcfs_by_sex(vcfs, sex_map, "F") == [vcfs[0], vcfs[2]]
        assert filter_vcfs_by_sex(vcfs, sex_map, "M") == [vcfs[1], vcfs[3]]

    def test_preserves_order(self, tmp_path: Path) -> None:
        vcfs = [tmp_path / f"s{i}.vcf.gz" for i in range(3)]
        sex_map = {vcfs[0]: "F", vcfs[1]: "F", vcfs[2]: "F"}
        assert filter_vcfs_by_sex(vcfs, sex_map, "F") == vcfs

    def test_empty_result_when_no_match(self, tmp_path: Path) -> None:
        vcfs = [tmp_path / "s0.vcf.gz"]
        sex_map = {vcfs[0]: "M"}
        assert filter_vcfs_by_sex(vcfs, sex_map, "F") == []

    def test_path_not_in_map_is_excluded(self, tmp_path: Path) -> None:
        vcfs = [tmp_path / "s0.vcf.gz", tmp_path / "s1.vcf.gz"]
        sex_map = {vcfs[0]: "F"}  # s1 not in map
        assert filter_vcfs_by_sex(vcfs, sex_map, "F") == [vcfs[0]]


# ---------------------------------------------------------------------------
# CLI integration tests
# ---------------------------------------------------------------------------

# Reuse conftest fixtures via imports: five_sample_vcfs, genetic_map_file.
# We need a chrX-flavoured set of fixtures too.

_CHRX_MAP_TEXT = textwrap.dedent("""\
    pos chr cM
    100 chrX 0.0
    100000 chrX 1.0
    500000 chrX 3.0
    1000000 chrX 5.0
    5000000 chrX 8.0
    10000000 chrX 10.0
""")

_CHRX_VCF_HEADER = textwrap.dedent("""\
    ##fileformat=VCFv4.1
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##contig=<ID=chrX,length=155270560>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
""")

_VARIANT_POSITIONS = [100, 5000, 50000, 200000, 400000, 600000, 800000, 2000000, 4000000, 9000000]
_REFS = ["T", "C", "G", "A", "T", "C", "G", "A", "T", "C"]
_ALTS = ["A", "G", "C", "T", "A", "G", "C", "T", "A", "G"]


def _make_chrx_vcf(sample: str, genotypes: list[str]) -> str:
    header = _CHRX_VCF_HEADER.format(sample=sample)
    lines = [
        f"chrX\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}"
        for pos, ref, alt, gt in zip(_VARIANT_POSITIONS, _REFS, _ALTS, genotypes)
    ]
    return header + "\n".join(lines) + "\n"


@pytest.fixture
def chrx_genetic_map(tmp_path: Path) -> Path:
    p = tmp_path / "map_chrX.txt"
    p.write_text(_CHRX_MAP_TEXT)
    return p


@pytest.fixture
def chrx_vcfs(tmp_path: Path) -> list[Path]:
    """6 per-sample VCFs on chrX: samples 0-3 female (diploid), 4-5 male (hom-only)."""
    female_gts = [
        ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0"],
        ["0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"],
        ["1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1"],
        ["0/0", "0/1", "0/1", "1/1", "0/0", "0/1", "0/1", "1/1", "0/0", "0/1"],
    ]
    male_gts = [
        ["0/0", "0/0", "1/1", "0/0", "1/1", "0/0", "1/1", "0/0", "0/0", "1/1"],
        ["1/1", "1/1", "0/0", "1/1", "0/0", "1/1", "0/0", "1/1", "1/1", "0/0"],
    ]
    paths = []
    for i, gts in enumerate(female_gts):
        p = tmp_path / f"female_{i}.vcf"
        p.write_text(_make_chrx_vcf(f"female_{i}", gts))
        paths.append(p)
    for i, gts in enumerate(male_gts):
        p = tmp_path / f"male_{i}.vcf"
        p.write_text(_make_chrx_vcf(f"male_{i}", gts))
        paths.append(p)
    return paths  # [female_0..3, male_0..1]


@pytest.fixture
def sex_file_for_chrx(tmp_path: Path, chrx_vcfs: list[Path]) -> Path:
    """Sex file mapping the 6 chrX VCFs: first 4 F, last 2 M."""
    lines = []
    for i, vcf in enumerate(chrx_vcfs):
        sex = "F" if i < 4 else "M"
        lines.append(f"{vcf}  {sex}")
    p = tmp_path / "sex.txt"
    p.write_text("\n".join(lines) + "\n")
    return p


class TestCliSexFilter:
    def test_autosome_uses_all_donors(
        self,
        tmp_path: Path,
        five_sample_vcfs: list[Path],
        genetic_map_file: Path,
    ) -> None:
        """With --sex-file on an autosome (chr22), all donors are used."""
        # Sex file marks 3 female, 2 male
        sex_content = "\n".join(
            f"{p}  {'F' if i < 3 else 'M'}"
            for i, p in enumerate(five_sample_vcfs)
        )
        sex_file = tmp_path / "sex.txt"
        sex_file.write_text(sex_content + "\n")

        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in five_sample_vcfs) + "\n")
        out_dir = tmp_path / "out"

        result = CliRunner().invoke(main, [
            "shuffle",
            "--input", f"@{filelist}",
            "--output-dir", str(out_dir),
            "--genetic-map", str(genetic_map_file),
            "--chromosome", "chr22",
            "--n-samples", "3",
            "--seed", "42",
            "--sex-file", str(sex_file),
        ])
        assert result.exit_code == 0, result.output
        # Output should be produced normally
        assert len(list(out_dir.glob("synthetic_*.vcf*"))) >= 3

    def test_chrx_uses_only_female_donors(
        self,
        tmp_path: Path,
        chrx_vcfs: list[Path],
        chrx_genetic_map: Path,
        sex_file_for_chrx: Path,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """With --sex-file on chrX, only the 4 female donors are used."""
        import logging
        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in chrx_vcfs) + "\n")
        out_dir = tmp_path / "out"

        with caplog.at_level(logging.INFO, logger="v_shuffler"):
            result = CliRunner().invoke(main, [
                "shuffle",
                "--input", f"@{filelist}",
                "--output-dir", str(out_dir),
                "--genetic-map", str(chrx_genetic_map),
                "--chromosome", "chrX",
                "--n-samples", "3",
                "--seed", "42",
                "--sex-file", str(sex_file_for_chrx),
            ])
        assert result.exit_code == 0, result.output
        assert len(list(out_dir.glob("synthetic_*.vcf*"))) >= 3
        # Log should mention 4 female donors of 6 total
        assert any("4 female donors" in r.message for r in caplog.records)

    def test_chrx_without_sex_file_warns(
        self,
        tmp_path: Path,
        chrx_vcfs: list[Path],
        chrx_genetic_map: Path,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Running chrX without --sex-file should log a warning but still complete."""
        import logging
        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in chrx_vcfs) + "\n")
        out_dir = tmp_path / "out"

        with caplog.at_level(logging.WARNING, logger="v_shuffler"):
            result = CliRunner().invoke(main, [
                "shuffle",
                "--input", f"@{filelist}",
                "--output-dir", str(out_dir),
                "--genetic-map", str(chrx_genetic_map),
                "--chromosome", "chrX",
                "--n-samples", "2",
                "--seed", "42",
            ])
        assert result.exit_code == 0, result.output
        assert any("without --sex-file" in r.message for r in caplog.records)

    def test_chrx_no_female_donors_errors(
        self,
        tmp_path: Path,
        chrx_vcfs: list[Path],
        chrx_genetic_map: Path,
    ) -> None:
        """If the sex file has no female donors, the run should exit with an error."""
        # Mark all donors as male
        all_male = "\n".join(f"{p}  M" for p in chrx_vcfs) + "\n"
        sex_file = tmp_path / "sex_all_male.txt"
        sex_file.write_text(all_male)

        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in chrx_vcfs) + "\n")
        out_dir = tmp_path / "out"

        result = CliRunner().invoke(main, [
            "shuffle",
            "--input", f"@{filelist}",
            "--output-dir", str(out_dir),
            "--genetic-map", str(chrx_genetic_map),
            "--chromosome", "chrX",
            "--n-samples", "2",
            "--seed", "42",
            "--sex-file", str(sex_file),
        ])
        assert result.exit_code != 0
        assert "No female donors" in result.output
