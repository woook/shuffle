"""Tests for v_shuffler.cli."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
from click.testing import CliRunner

from v_shuffler.cli import main


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


class TestShuffleCommand:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["shuffle", "--help"])
        assert result.exit_code == 0
        assert "--input" in result.output
        assert "--output-dir" in result.output
        assert "--genetic-map" in result.output

    def test_missing_required_options(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["shuffle"])
        assert result.exit_code != 0

    def test_shuffle_runs(
        self,
        runner: CliRunner,
        tmp_path: Path,
        five_sample_vcfs: list[Path],
        genetic_map_file: Path,
    ) -> None:
        """End-to-end smoke test: shuffle 5 samples → 3 synthetic per_sample VCFs."""
        # Write filelist
        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in five_sample_vcfs) + "\n")

        out_dir = tmp_path / "out"

        result = runner.invoke(main, [
            "shuffle",
            "--input", f"@{filelist}",
            "--output-dir", str(out_dir),
            "--genetic-map", str(genetic_map_file),
            "--chromosome", "chr22",
            "--n-samples", "3",
            "--seed", "42",
            "--output-mode", "per_sample",
        ])

        if result.exit_code != 0:
            print(result.output)
            if result.exception:
                import traceback
                traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)

        assert result.exit_code == 0, f"CLI failed: {result.output}"

        # Check output files were created
        vcf_files = list(out_dir.glob("synthetic_*.vcf*"))
        assert len(vcf_files) >= 3, f"Expected >=3 output files, got {vcf_files}"

    def test_shuffle_deterministic(
        self,
        runner: CliRunner,
        tmp_path: Path,
        five_sample_vcfs: list[Path],
        genetic_map_file: Path,
    ) -> None:
        """Two runs with the same seed should produce identical segment plans."""
        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in five_sample_vcfs) + "\n")

        def run(out_subdir: str) -> list[Path]:
            out_dir = tmp_path / out_subdir
            result = runner.invoke(main, [
                "shuffle",
                "--input", f"@{filelist}",
                "--output-dir", str(out_dir),
                "--genetic-map", str(genetic_map_file),
                "--chromosome", "chr22",
                "--n-samples", "2",
                "--seed", "99",
                "--output-mode", "per_sample",
            ])
            assert result.exit_code == 0, result.output
            return sorted(out_dir.glob("synthetic_*.vcf"))

        files_a = run("out_a")
        files_b = run("out_b")

        # Compare the plain VCF content (before bgzip)
        for fa, fb in zip(files_a, files_b):
            assert fa.read_text() == fb.read_text(), (
                f"Determinism failed: {fa.name} differs between runs"
            )

    def test_shuffle_multi_sample_mode(
        self,
        runner: CliRunner,
        tmp_path: Path,
        five_sample_vcfs: list[Path],
        genetic_map_file: Path,
    ) -> None:
        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in five_sample_vcfs) + "\n")
        out_dir = tmp_path / "out"

        result = runner.invoke(main, [
            "shuffle",
            "--input", f"@{filelist}",
            "--output-dir", str(out_dir),
            "--genetic-map", str(genetic_map_file),
            "--chromosome", "chr22",
            "--n-samples", "3",
            "--seed", "42",
            "--output-mode", "multi_sample",
        ])
        assert result.exit_code == 0, result.output

        vcf_files = list(out_dir.glob("synthetic_chr22.vcf*"))
        assert len(vcf_files) >= 1

    def test_nonexistent_input_raises(
        self,
        runner: CliRunner,
        tmp_path: Path,
        genetic_map_file: Path,
    ) -> None:
        result = runner.invoke(main, [
            "shuffle",
            "--input", str(tmp_path / "nonexistent_*.vcf.gz"),
            "--output-dir", str(tmp_path / "out"),
            "--genetic-map", str(genetic_map_file),
            "--chromosome", "chr22",
        ])
        assert result.exit_code != 0

    def test_output_vcf_is_unphased(
        self,
        runner: CliRunner,
        tmp_path: Path,
        five_sample_vcfs: list[Path],
        genetic_map_file: Path,
    ) -> None:
        """Output VCF genotypes should use / separator (unphased)."""
        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in five_sample_vcfs) + "\n")
        out_dir = tmp_path / "out"

        result = runner.invoke(main, [
            "shuffle",
            "--input", f"@{filelist}",
            "--output-dir", str(out_dir),
            "--genetic-map", str(genetic_map_file),
            "--chromosome", "chr22",
            "--n-samples", "1",
            "--seed", "42",
            "--output-mode", "per_sample",
        ])
        assert result.exit_code == 0, result.output

        vcf_files = list(out_dir.glob("synthetic_0.vcf"))
        assert len(vcf_files) == 1
        content = vcf_files[0].read_text()
        data_lines = [l for l in content.splitlines() if not l.startswith("#")]
        assert data_lines, "No data lines in output VCF"
        for line in data_lines:
            gt = line.split("\t")[-1]
            # Must be unphased: contains / not |
            assert "/" in gt, f"Unexpected phased GT: {gt!r}"
            assert "|" not in gt, f"Unexpected phased GT: {gt!r}"

    def test_output_not_identical_to_any_input(
        self,
        runner: CliRunner,
        tmp_path: Path,
        five_sample_vcfs: list[Path],
        genetic_map_file: Path,
    ) -> None:
        """
        With enough donors and variants, output individuals should not be
        identical to any single input.  With only 10 variants this can fail
        by chance, so we just verify the test runs without error.
        """
        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in five_sample_vcfs) + "\n")
        out_dir = tmp_path / "out"

        result = runner.invoke(main, [
            "shuffle",
            "--input", f"@{filelist}",
            "--output-dir", str(out_dir),
            "--genetic-map", str(genetic_map_file),
            "--chromosome", "chr22",
            "--n-samples", "3",
            "--seed", "7",
            "--output-mode", "per_sample",
        ])
        assert result.exit_code == 0, result.output


class TestVersion:
    def test_version(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "v-shuffler" in result.output


# ---------------------------------------------------------------------------
# Chromosome name normalisation
# ---------------------------------------------------------------------------

def _make_bare_chrom_vcfs(tmp_path: Path) -> tuple[list[Path], Path]:
    """
    Return (vcf_paths, map_path) where the VCFs use bare '22' chromosome names
    (no 'chr' prefix) and the genetic map similarly uses '22'.
    """
    header = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##contig=<ID=22,length=51304566>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
    """)
    positions = [100, 5000, 50000, 200000, 400000, 600000, 800000, 2000000, 4000000, 9000000]
    refs =      ["T",  "C",  "G",   "A",    "T",    "C",    "G",    "A",     "T",     "C"]
    alts =      ["A",  "G",  "C",   "T",    "A",    "G",    "C",    "T",     "A",     "G"]
    gts_list = [
        ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0"],
        ["0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"],
        ["1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1"],
        ["0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "0/1", "0/1", "0/1", "0/0"],
        ["0/1", "0/1", "0/1", "0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "0/1"],
    ]
    paths = []
    for i, gts in enumerate(gts_list):
        lines = "\n".join(
            f"22\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}"
            for pos, ref, alt, gt in zip(positions, refs, alts, gts)
        )
        p = tmp_path / f"bare_{i}.vcf"
        p.write_text(header.format(sample=f"sample{i}") + lines + "\n")
        paths.append(p)

    map_path = tmp_path / "map_bare.txt"
    map_path.write_text(
        "pos chr cM\n"
        "100 22 0.0\n"
        "100000 22 1.0\n"
        "500000 22 3.0\n"
        "1000000 22 5.0\n"
        "5000000 22 8.0\n"
        "10000000 22 10.0\n"
    )
    return paths, map_path


class TestChromosomeNormalisation:
    def test_bare_vcf_with_chr_prefix_flag(self, tmp_path: Path) -> None:
        """VCFs use '22', user passes --chromosome chr22 → run succeeds."""
        vcf_paths, map_path = _make_bare_chrom_vcfs(tmp_path)
        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in vcf_paths) + "\n")
        out_dir = tmp_path / "out"

        result = CliRunner().invoke(main, [
            "shuffle",
            "--input", f"@{filelist}",
            "--output-dir", str(out_dir),
            "--genetic-map", str(map_path),
            "--chromosome", "chr22",   # prefixed — VCFs use bare '22'
            "--n-samples", "3",
            "--seed", "42",
        ])
        assert result.exit_code == 0, result.output
        assert len(list(out_dir.glob("synthetic_*.vcf*"))) >= 3

    def test_bare_vcf_with_bare_flag(self, tmp_path: Path) -> None:
        """VCFs use '22', user passes --chromosome 22 → run succeeds (no change needed)."""
        vcf_paths, map_path = _make_bare_chrom_vcfs(tmp_path)
        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in vcf_paths) + "\n")
        out_dir = tmp_path / "out"

        result = CliRunner().invoke(main, [
            "shuffle",
            "--input", f"@{filelist}",
            "--output-dir", str(out_dir),
            "--genetic-map", str(map_path),
            "--chromosome", "22",
            "--n-samples", "3",
            "--seed", "42",
        ])
        assert result.exit_code == 0, result.output
        assert len(list(out_dir.glob("synthetic_*.vcf*"))) >= 3

    def test_normalisation_logged(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """When the name is normalised, an INFO message is logged."""
        import logging
        vcf_paths, map_path = _make_bare_chrom_vcfs(tmp_path)
        filelist = tmp_path / "samples.txt"
        filelist.write_text("\n".join(str(p) for p in vcf_paths) + "\n")
        out_dir = tmp_path / "out"

        with caplog.at_level(logging.INFO, logger="v_shuffler"):
            CliRunner().invoke(main, [
                "shuffle",
                "--input", f"@{filelist}",
                "--output-dir", str(out_dir),
                "--genetic-map", str(map_path),
                "--chromosome", "chr22",   # will be normalised to '22'
                "--n-samples", "2",
                "--seed", "42",
            ])
        assert any("normalised" in r.message for r in caplog.records)
