"""
VCF reader for v-shuffler.

Reads one VCF per sample (all covering the same chromosome and sites)
and yields GenotypePool chunks of shape (chunk_size, n_samples).

Design choices:
- Uses cyvcf2 for fast VCF iteration.
- Validates that all per-sample VCFs share the same sites (CHROM/POS/REF/ALT)
  on the first chunk; subsequent chunks are assumed consistent.
- Genotype encoding: dosage (0/1/2), with MISSING=255 sentinel.
- Multi-allelic sites: dosage = total number of non-ref allele copies
  (sum of all alt allele counts per individual).
- Variants with missing rate > max_missing_rate are dropped.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterator

import numpy as np

try:
    from cyvcf2 import VCF
except ImportError as exc:
    raise ImportError(
        "cyvcf2 is required for VCF reading. Install it with: pip install cyvcf2"
    ) from exc

from v_shuffler.core.genotype_pool import MISSING, GenotypePool, VariantInfo
from v_shuffler.io.genetic_map import GeneticMap

logger = logging.getLogger(__name__)


def resolve_chromosome_name(vcf_path: Path, chromosome: str) -> str:
    """
    Return the form of *chromosome* that the VCF file actually uses.

    Some pipelines write ``chr22``; others write ``22``.  This function
    opens the first VCF in the cohort (header + index only — no records
    are read), inspects the sequence names reported by cyvcf2, and returns
    the matching form.

    Resolution order (first match wins):
      1. *chromosome* as supplied by the user.
      2. With ``chr`` prefix added (e.g. ``22`` → ``chr22``).
      3. With ``chr`` prefix stripped (e.g. ``chr22`` → ``22``).

    Falls back to *chromosome* unchanged when the VCF has neither
    ``##contig`` header lines nor a tabix/CSI index (e.g. small plain-text
    test files without contig headers).  In that case ``_region_iter``'s
    plain-VCF fallback, which already checks all three forms, handles the
    iteration correctly.

    Parameters
    ----------
    vcf_path :
        Path to any one of the donor VCF files.  Only the header / index
        is inspected; no genotype records are read.
    chromosome :
        The chromosome name supplied by the user (e.g. from ``--chromosome``).

    Returns
    -------
    str
        The chromosome name as it appears in *vcf_path*, or *chromosome*
        unchanged if the VCF provides no sequence-name information.
    """
    reader = VCF(str(vcf_path))
    seqnames = set(reader.seqnames)
    reader.close()

    if not seqnames:
        return chromosome  # no contig info available; keep as-is

    bare = chromosome.lstrip("chr")
    prefixed = "chr" + bare

    for candidate in (chromosome, prefixed, bare):
        if candidate in seqnames:
            return candidate

    return chromosome  # fallback: let downstream code handle it


def _gt_to_dosage(gt: list) -> int:
    """
    Convert a cyvcf2 genotype list [allele1, allele2, phased] to a dosage.

    Returns
    -------
    int
        0, 1, or 2 for biallelic; sum of non-ref copies for multi-allelic;
        MISSING (255) if any allele is -1 (missing).
    """
    a1, a2 = gt[0], gt[1]
    if a1 < 0 or a2 < 0:
        return MISSING
    # Count non-ref copies (allele index > 0)
    return int(a1 > 0) + int(a2 > 0)



class PerSampleVCFReader:
    """
    Reads one VCF file per sample and merges them into GenotypePool chunks.

    All VCF files must cover the same variants in the same order.
    The reader validates site consistency across files for the first chunk.

    Parameters
    ----------
    vcf_paths : list[Path]
        Paths to per-sample VCF files (one sample column each).
    chromosome : str
        Chromosome to process (e.g. "chr22" or "22").
    genetic_map : GeneticMap
        Used to annotate variants with cM positions.
    chunk_size : int
        Number of variants per yielded chunk.
    max_missing_rate : float
        Variants with a higher missing rate are skipped.
    """

    def __init__(
        self,
        vcf_paths: list[Path],
        chromosome: str,
        genetic_map: GeneticMap,
        chunk_size: int = 50_000,
        max_missing_rate: float = 0.05,
    ) -> None:
        if not vcf_paths:
            raise ValueError("vcf_paths must not be empty")
        self.vcf_paths = vcf_paths
        self.chromosome = chromosome
        self.genetic_map = genetic_map
        self.chunk_size = chunk_size
        self.max_missing_rate = max_missing_rate
        self.n_samples = len(vcf_paths)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def iter_chunks(self) -> Iterator[GenotypePool]:
        """
        Yield GenotypePool objects, each covering up to chunk_size variants.

        Variants exceeding max_missing_rate are silently skipped.
        """
        readers = [VCF(str(p)) for p in self.vcf_paths]
        iterators = [
            self._region_iter(r, p, self.chromosome)
            for r, p in zip(readers, self.vcf_paths)
        ]

        chunk_dosages: list[np.ndarray] = []
        chunk_positions: list[int] = []
        chunk_cm: list[float] = []
        chunk_info: list[VariantInfo] = []

        skipped_missing = 0
        total_seen = 0
        first_chunk = True

        # Iterate in lockstep across all per-sample VCFs
        for variants in zip(*iterators):
            self._check_site_consistency(variants, first_chunk)
            first_chunk = False
            total_seen += 1

            # Build dosage vector for this variant
            dosages = np.empty(self.n_samples, dtype=np.uint8)
            for s_idx, v in enumerate(variants):
                dosages[s_idx] = _gt_to_dosage(v.genotypes[0])

            # Missing rate filter
            n_missing = int(np.sum(dosages == MISSING))
            if n_missing / self.n_samples > self.max_missing_rate:
                skipped_missing += 1
                continue

            v0 = variants[0]
            cm = float(self.genetic_map.bp_to_cm(v0.POS))

            chunk_dosages.append(dosages)
            chunk_positions.append(v0.POS)
            chunk_cm.append(cm)
            chunk_info.append(VariantInfo(
                chrom=v0.CHROM,
                pos=v0.POS,
                ref=v0.REF,
                alts=list(v0.ALT),
                id=v0.ID or ".",
                qual=v0.QUAL,
                filters=list(v0.FILTER) if v0.FILTER else [],
                cm_pos=cm,
            ))

            if len(chunk_dosages) == self.chunk_size:
                yield self._make_pool(chunk_dosages, chunk_positions, chunk_cm, chunk_info)
                chunk_dosages, chunk_positions, chunk_cm, chunk_info = [], [], [], []

        # Yield the final partial chunk
        if chunk_dosages:
            yield self._make_pool(chunk_dosages, chunk_positions, chunk_cm, chunk_info)

        for r in readers:
            r.close()

        if skipped_missing:
            logger.warning(
                "Skipped %d / %d variants due to high missing rate (>%.0f%%)",
                skipped_missing, total_seen, self.max_missing_rate * 100,
            )

    def iter_positions(self) -> np.ndarray:
        """
        Return all variant positions on self.chromosome from the first per-sample VCF.

        No missing-rate filter is applied; this is a lightweight first pass that
        collects only POS values for region detection.  The full genotype pass
        (with filtering) happens later in iter_chunks().

        Returns
        -------
        np.ndarray, int64, shape (n_variants,)
            Sorted array of base-pair positions.
        """
        reader = VCF(str(self.vcf_paths[0]))
        positions: list[int] = []
        for v in self._region_iter(reader, self.vcf_paths[0], self.chromosome):
            positions.append(v.POS)
        reader.close()
        return np.array(positions, dtype=np.int64)

    def get_header_vcf(self):
        """Return an open cyvcf2.VCF handle to the first file (for header extraction)."""
        return VCF(str(self.vcf_paths[0]))

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _make_pool(
        dosages_list: list[np.ndarray],
        positions: list[int],
        cm_list: list[float],
        variant_info: list[VariantInfo],
    ) -> GenotypePool:
        dosage_matrix = np.stack(dosages_list, axis=0)  # (n_variants, n_samples)
        return GenotypePool(
            dosages=dosage_matrix,
            positions=np.array(positions, dtype=np.int64),
            cm_pos=np.array(cm_list, dtype=np.float64),
            variant_info=variant_info,
        )

    @staticmethod
    def _region_iter(reader: "VCF", path: Path, chrom: str):
        """
        Return an iterator over variants on *chrom*.

        Uses tabix region query when the file is bgzipped + indexed (.tbi/.csi).
        Falls back to full-file iteration with chromosome filter for plain
        VCFs (useful in tests and small-scale usage).
        """
        import os
        fname = str(path)
        has_index = os.path.exists(fname + ".tbi") or os.path.exists(fname + ".csi")
        if has_index:
            return reader(chrom)
        # Fall back: iterate all records and filter by chromosome
        target_chroms = {chrom, chrom.lstrip("chr"), "chr" + chrom.lstrip("chr")}
        return (v for v in reader if v.CHROM in target_chroms)

    @staticmethod
    def _check_site_consistency(variants: tuple, first_chunk: bool) -> None:
        """Verify all per-sample VCFs agree on CHROM/POS/REF/ALT."""
        if not first_chunk:
            return
        v0 = variants[0]
        for v in variants[1:]:
            if v.CHROM != v0.CHROM or v.POS != v0.POS:
                raise ValueError(
                    f"Site mismatch between VCF files: "
                    f"expected {v0.CHROM}:{v0.POS}, got {v.CHROM}:{v.POS}. "
                    "All per-sample VCFs must have the same variants in the same order."
                )
            if v.REF != v0.REF or list(v.ALT) != list(v0.ALT):
                raise ValueError(
                    f"Allele mismatch at {v0.CHROM}:{v0.POS}: "
                    f"REF/ALT differ between VCF files."
                )
