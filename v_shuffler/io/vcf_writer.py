"""
VCF writer for v-shuffler.

Writes synthetic unphased genotype calls for all output individuals.

Two modes:
  - per_sample: one bgzipped VCF per synthetic individual.
  - multi_sample: one bgzipped multi-sample VCF.

The writer:
  - Copies all header lines from the first input VCF (via cyvcf2).
  - Strips sample-identifiable metadata (##sample=).
  - Adds a provenance ##v-shuffler= header line.
  - Writes GT as dosage → unphased format: 0→0/0, 1→0/1, 2→1/1, 255→./.
  - Calls bgzip + tabix (via subprocess) on each output file after writing.
"""

from __future__ import annotations

import logging
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

try:
    from cyvcf2 import VCF
except ImportError as exc:
    raise ImportError(
        "cyvcf2 is required. Install it with: pip install cyvcf2"
    ) from exc

from v_shuffler.core.genotype_pool import MISSING, GenotypePool

logger = logging.getLogger(__name__)

# Dosage → unphased GT string
_DOSAGE_TO_GT: dict[int, str] = {0: "0/0", 1: "0/1", 2: "1/1", MISSING: "./."}


def _dosage_to_gt_str(dosage: int) -> str:
    return _DOSAGE_TO_GT.get(int(dosage), "./.")


def _build_sample_str(dosage: int, field_vals: list[str]) -> str:
    """
    Build the per-sample column value, e.g. ``'0/1:0.4531:127'`` or
    ``'0/1:0.4531:3932:1904,3028'``.

    *field_vals* are VCF-ready strings already formatted by
    ``_get_format_str`` in the reader (single values like ``"0.4531"`` or
    multi-value like ``"1904,3028"``).  They are passed straight through.
    """
    parts = [_dosage_to_gt_str(dosage)] + list(field_vals)
    return ":".join(parts)


def _make_provenance_line(version: str, seed: int | None, chromosome: str) -> str:
    ts = datetime.now(tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    seed_str = str(seed) if seed is not None else "None"
    return (
        f'##v-shuffler=<version="{version}",seed="{seed_str}",'
        f'chromosome="{chromosome}",date="{ts}">'
    )


def _build_header_string(template_vcf: VCF, sample_names: list[str], provenance_line: str) -> str:
    """
    Build a VCF header string from a template VCF, replacing sample names
    and adding a provenance line.  Strips ##sample= lines.
    """
    header_lines = []
    for line in str(template_vcf.raw_header).splitlines():
        if line.startswith("##sample="):
            continue  # strip sample-identifiable metadata
        if line.startswith("#CHROM"):
            # Insert provenance before the column header line
            header_lines.append(provenance_line)
            # Rebuild the column header with new sample names
            cols = line.split("\t")[:9]  # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
            header_lines.append("\t".join(cols + sample_names))
            continue
        header_lines.append(line)
    return "\n".join(header_lines) + "\n"


class SyntheticVCFWriter:
    """
    Writes synthetic genotype data to VCF file(s).

    Parameters
    ----------
    output_dir : Path
    sample_names : list[str]
        Names for the synthetic output samples, e.g. ["synthetic_0", ...].
    template_vcf_path : Path
        A real input VCF; its header is used as a template.
    output_mode : str
        "per_sample" or "multi_sample".
    seed : int | None
    chromosome : str
    version : str
    """

    def __init__(
        self,
        output_dir: Path,
        sample_names: list[str],
        template_vcf_path: Path,
        output_mode: str,
        seed: int | None = None,
        chromosome: str = "",
        version: str = "0.1.0",
    ) -> None:
        self.output_dir = output_dir
        self.sample_names = sample_names
        self.output_mode = output_mode

        provenance = _make_provenance_line(version, seed, chromosome)

        template = VCF(str(template_vcf_path))
        self._template_vcf_path = template_vcf_path

        if output_mode == "multi_sample":
            out_path = output_dir / f"synthetic_{chromosome}.vcf"
            header_str = _build_header_string(template, sample_names, provenance)
            self._writers = [self._open_writer(out_path, header_str)]
            self._out_paths = [out_path]
        else:
            self._writers = []
            self._out_paths = []
            for name in sample_names:
                out_path = output_dir / f"{name}.vcf"
                header_str = _build_header_string(template, [name], provenance)
                self._writers.append(self._open_writer(out_path, header_str))
                self._out_paths.append(out_path)

        template.close()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def write_chunk(
        self,
        pool: GenotypePool,
        synthetic_dosages: np.ndarray,
        synthetic_fields: dict[str, np.ndarray] | None = None,
    ) -> None:
        """
        Write one chunk of synthetic genotypes.

        Parameters
        ----------
        pool : GenotypePool
            Source pool — used for variant metadata (CHROM, POS, REF, ALT, etc.).
        synthetic_dosages : np.ndarray, shape (n_variants, n_synthetic_samples), uint8
            Synthetic dosage matrix from mosaic_builder.build_synthetic_genotypes.
        synthetic_fields : dict[str, np.ndarray] or None
            Optional additional FORMAT fields to write alongside GT.  Each value
            must have shape ``(n_variants, n_synthetic_samples)``, dtype object
            (VCF-ready strings as produced by ``build_synthetic_genotypes``).
            Keys are written in iteration order, e.g. ``{"AF": ..., "DP": ...}``
            produces a FORMAT column of ``GT:AF:DP``.
        """
        n_out = len(self.sample_names)
        if synthetic_dosages.shape != (pool.n_variants, n_out):
            raise ValueError(
                f"synthetic_dosages shape {synthetic_dosages.shape} does not match "
                f"({pool.n_variants}, {n_out})"
            )

        fields = synthetic_fields or {}
        for fname, farr in fields.items():
            if farr.shape != (pool.n_variants, n_out):
                raise ValueError(
                    f"synthetic_fields[{fname!r}] has shape {farr.shape}, "
                    f"expected ({pool.n_variants}, {n_out})"
                )
        field_names = list(fields.keys())
        fmt_str = "GT" + (":" + ":".join(field_names) if field_names else "")

        if self.output_mode == "multi_sample":
            writer = self._writers[0]
            for v_idx, vi in enumerate(pool.variant_info):
                sample_cols = "\t".join(
                    _build_sample_str(
                        synthetic_dosages[v_idx, s],
                        [str(fields[f][v_idx, s]) for f in field_names],
                    )
                    for s in range(n_out)
                )
                writer.write(self._format_record(vi, sample_cols, fmt_str))
        else:
            for v_idx, vi in enumerate(pool.variant_info):
                for s_idx, writer in enumerate(self._writers):
                    sample_col = _build_sample_str(
                        synthetic_dosages[v_idx, s_idx],
                        [str(fields[f][v_idx, s_idx]) for f in field_names],
                    )
                    writer.write(self._format_record(vi, sample_col, fmt_str))

    def finalize(self) -> list[Path]:
        """
        Close all writers and compress + index the output files.

        Returns
        -------
        list[Path]
            Paths to the final .vcf.gz files.
        """
        for writer in self._writers:
            writer.close()

        final_paths = []
        for vcf_path in self._out_paths:
            gz_path = self._bgzip(vcf_path)
            self._tabix(gz_path)
            final_paths.append(gz_path)

        return final_paths

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _open_writer(path: Path, header_str: str):
        """Write a VCF header and return an open file handle (plain text)."""
        # We write plain VCF text; bgzip is called in finalize()
        fh = open(path, "w")
        fh.write(header_str)
        return fh

    @staticmethod
    def _format_record(vi, sample_cols: str, fmt_str: str = "GT") -> str:
        """Format a VCF data line for one variant."""
        qual = "." if vi.qual is None else str(vi.qual)
        filt = ";".join(vi.filters) if vi.filters else "PASS"
        alts = ",".join(vi.alts) if vi.alts else "."
        return (
            f"{vi.chrom}\t{vi.pos}\t{vi.id}\t{vi.ref}\t{alts}\t"
            f"{qual}\t{filt}\t.\t{fmt_str}\t{sample_cols}\n"
        )

    @staticmethod
    def _bgzip(path: Path) -> Path:
        """Compress a VCF with bgzip, return path to .gz file."""
        gz_path = path.with_suffix(".vcf.gz")
        try:
            result = subprocess.run(
                ["bgzip", "-f", str(path)],
                capture_output=True,
            )
        except FileNotFoundError:
            logger.warning("bgzip not found; output will not be compressed: %s", path)
            return path
        if result.returncode != 0:
            logger.warning(
                "bgzip failed for %s: %s", path, result.stderr.decode()
            )
            return path  # return uncompressed path if bgzip unavailable
        return gz_path

    @staticmethod
    def _tabix(path: Path) -> None:
        """Index a bgzipped VCF with tabix."""
        try:
            result = subprocess.run(
                ["tabix", "-p", "vcf", str(path)],
                capture_output=True,
            )
        except FileNotFoundError:
            logger.warning("tabix not found; output will not be indexed: %s", path)
            return
        if result.returncode != 0:
            logger.warning(
                "tabix failed for %s: %s", path, result.stderr.decode()
            )
