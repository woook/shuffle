"""Shared fixtures for v-shuffler tests."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from v_shuffler.io.genetic_map import GeneticMap


VCF_HEADER_TEMPLATE = textwrap.dedent("""\
    ##fileformat=VCFv4.1
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##contig=<ID=chr22,length=51304566>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
""")

# A genetic map for chr22 spanning 10 cM
MAP_TEXT = textwrap.dedent("""\
    pos chr cM
    100 chr22 0.0
    100000 chr22 1.0
    500000 chr22 3.0
    1000000 chr22 5.0
    5000000 chr22 8.0
    10000000 chr22 10.0
""")

# 10 variants spread across the map
VARIANT_POSITIONS = [100, 5000, 50000, 200000, 400000, 600000, 800000, 2000000, 4000000, 9000000]


def make_vcf_content(sample: str, genotypes: list[str]) -> str:
    assert len(genotypes) == len(VARIANT_POSITIONS)
    header = VCF_HEADER_TEMPLATE.format(sample=sample)
    lines = []
    refs = ["T", "C", "G", "A", "T", "C", "G", "A", "T", "C"]
    alts = ["A", "G", "C", "T", "A", "G", "C", "T", "A", "G"]
    for pos, ref, alt, gt in zip(VARIANT_POSITIONS, refs, alts, genotypes):
        lines.append(f"chr22\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}")
    return header + "\n".join(lines) + "\n"


@pytest.fixture
def genetic_map_file(tmp_path: Path) -> Path:
    p = tmp_path / "map_chr22.txt"
    p.write_text(MAP_TEXT)
    return p


@pytest.fixture
def genetic_map(genetic_map_file: Path) -> GeneticMap:
    return GeneticMap(genetic_map_file, "chr22")


@pytest.fixture
def five_sample_vcfs(tmp_path: Path) -> list[Path]:
    """5 per-sample VCFs, 10 variants each on chr22."""
    sample_gts = [
        ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0"],
        ["0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"],
        ["1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1"],
        ["0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "0/1", "0/1", "0/1", "0/0"],
        ["0/1", "0/1", "0/1", "0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "0/1"],
    ]
    paths = []
    for i, gts in enumerate(sample_gts):
        content = make_vcf_content(f"sample{i}", gts)
        p = tmp_path / f"sample{i}.vcf"
        p.write_text(content)
        paths.append(p)
    return paths
