"""
Validation utilities for v-shuffler output.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


def run_validate(
    synth_paths: list[Path],
    reference_vcf: Path,
    chromosome: str,
) -> None:
    """
    Validate synthetic VCFs against a reference multi-sample VCF.

    Parameters
    ----------
    synth_paths : list[Path]
        Per-sample synthetic VCF files.
    reference_vcf : Path
        Multi-sample input VCF (all original donors merged).
    chromosome : str
        Chromosome to validate.
    """
    try:
        from cyvcf2 import VCF
    except ImportError as exc:
        raise ImportError("cyvcf2 is required for validation") from exc

    logger.info("Validating %d synthetic VCFs against %s", len(synth_paths), reference_vcf)

    # --- Step 1: compute AF in reference ---
    logger.info("Computing allele frequencies in reference ...")
    ref_afs, ref_positions = _compute_af(reference_vcf, chromosome)

    if len(ref_positions) == 0:
        logger.error("No variants found in reference VCF for chromosome %s", chromosome)
        return

    # --- Step 2: compute AF in synthetic output ---
    logger.info("Computing allele frequencies in synthetic output ...")
    synth_dosages: list[np.ndarray] = []
    for p in synth_paths:
        _, _, dosages = _read_dosages(p, chromosome, ref_positions)
        synth_dosages.append(dosages)

    if not synth_dosages:
        logger.error("No synthetic dosage data could be read")
        return

    synth_matrix = np.stack(synth_dosages, axis=1)  # (n_variants, n_synth)
    valid = synth_matrix != 255
    synth_af = np.where(
        valid.sum(axis=1) > 0,
        np.where(valid, synth_matrix.astype(float), 0).sum(axis=1)
        / (valid.sum(axis=1) * 2.0 + 1e-9),
        np.nan,
    )

    # --- Step 3: AF correlation ---
    common = ~np.isnan(synth_af) & ~np.isnan(ref_afs)
    if common.sum() < 10:
        logger.warning("Too few overlapping variants to compute AF correlation")
    else:
        r = float(np.corrcoef(ref_afs[common], synth_af[common])[0, 1])
        logger.info("Allele frequency correlation (r): %.4f", r)
        if r < 0.99:
            logger.warning("AF correlation %.4f is below recommended threshold 0.99", r)
        else:
            logger.info("AF correlation OK (%.4f >= 0.99)", r)

    # --- Step 4: identity check ---
    logger.info("Checking that no synthetic sample is identical to any input ...")
    ref_sample_dosages = _read_all_samples(reference_vcf, chromosome, ref_positions)
    n_identical = 0
    for s_idx, synth_dos in enumerate(synth_dosages):
        for r_idx in range(ref_sample_dosages.shape[1]):
            ref_dos = ref_sample_dosages[:, r_idx]
            overlap = (synth_dos != 255) & (ref_dos != 255)
            if overlap.sum() == 0:
                continue
            frac_identical = (synth_dos[overlap] == ref_dos[overlap]).mean()
            if frac_identical > 0.99:
                n_identical += 1
                logger.warning(
                    "Synthetic sample %d shares %.1f%% of genotypes with input sample %d — "
                    "possible identity leak",
                    s_idx, frac_identical * 100, r_idx,
                )

    if n_identical == 0:
        logger.info("Identity check passed: no synthetic sample is near-identical to any input.")
    else:
        logger.warning(
            "%d potential identity leak(s) detected. This may occur with very small pools.", n_identical
        )

    logger.info("Validation complete.")


def _compute_af(vcf_path: Path, chromosome: str) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-site allele frequency from a multi-sample VCF."""
    from cyvcf2 import VCF
    vcf = VCF(str(vcf_path))
    afs = []
    positions = []
    for v in vcf(chromosome):
        gts = v.genotypes
        alleles = [a for gt in gts for a in gt[:2] if a >= 0]
        if not alleles:
            continue
        af = sum(a > 0 for a in alleles) / len(alleles)
        afs.append(af)
        positions.append(v.POS)
    vcf.close()
    return np.array(afs), np.array(positions, dtype=np.int64)


def _read_dosages(
    vcf_path: Path, chromosome: str, target_positions: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read dosages for a single-sample VCF at target_positions."""
    from cyvcf2 import VCF
    pos_set = set(target_positions.tolist())
    vcf = VCF(str(vcf_path))
    dosages = []
    positions = []
    for v in vcf(chromosome):
        if v.POS not in pos_set:
            continue
        gt = v.genotypes[0]
        a1, a2 = gt[0], gt[1]
        if a1 < 0 or a2 < 0:
            d = 255
        else:
            d = int(a1 > 0) + int(a2 > 0)
        dosages.append(d)
        positions.append(v.POS)
    vcf.close()
    return np.array(positions, dtype=np.int64), np.array([]), np.array(dosages, dtype=np.uint8)


def _read_all_samples(
    vcf_path: Path, chromosome: str, target_positions: np.ndarray
) -> np.ndarray:
    """Read dosages for all samples at target_positions. Returns (n_var, n_samples)."""
    from cyvcf2 import VCF
    pos_set = set(target_positions.tolist())
    vcf = VCF(str(vcf_path))
    rows = []
    for v in vcf(chromosome):
        if v.POS not in pos_set:
            continue
        row = []
        for gt in v.genotypes:
            a1, a2 = gt[0], gt[1]
            if a1 < 0 or a2 < 0:
                row.append(255)
            else:
                row.append(int(a1 > 0) + int(a2 > 0))
        rows.append(row)
    vcf.close()
    if not rows:
        return np.empty((0, 0), dtype=np.uint8)
    return np.array(rows, dtype=np.uint8)
