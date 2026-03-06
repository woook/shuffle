"""
Validation utilities for v-shuffler output.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from v_shuffler.core.genotype_pool import MISSING
from v_shuffler.io.vcf_reader import _gt_to_dosage

logger = logging.getLogger(__name__)

_IDENTITY_THRESHOLD = 0.99
_AF_CORRELATION_THRESHOLD = 0.99


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
        from cyvcf2 import VCF  # noqa: F401 — confirms cyvcf2 is importable
    except ImportError as exc:
        raise ImportError("cyvcf2 is required for validation") from exc

    logger.info("Validating %d synthetic VCFs against %s", len(synth_paths), reference_vcf)

    # --- Step 1: read reference once (AF + all-sample dosages) ---
    logger.info("Reading reference VCF ...")
    ref_afs, ref_positions, ref_sample_dosages = _read_reference(reference_vcf, chromosome)

    if len(ref_positions) == 0:
        logger.error("No variants found in reference VCF for chromosome %s", chromosome)
        return

    # --- Step 2: compute AF in synthetic output ---
    logger.info("Computing allele frequencies in synthetic output ...")
    synth_dosages: list[np.ndarray] = []
    for p in synth_paths:
        _, dosages = _read_dosages(p, chromosome, ref_positions)
        synth_dosages.append(dosages)

    if not synth_dosages:
        logger.error("No synthetic dosage data could be read")
        return

    synth_matrix = np.stack(synth_dosages, axis=1)  # (n_variants, n_synth)
    valid = synth_matrix != MISSING
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
        if r < _AF_CORRELATION_THRESHOLD:
            logger.warning(
                "AF correlation %.4f is below recommended threshold %.2f",
                r, _AF_CORRELATION_THRESHOLD,
            )
        else:
            logger.info("AF correlation OK (%.4f >= %.2f)", r, _AF_CORRELATION_THRESHOLD)

    # --- Step 4: identity check ---
    logger.info("Checking that no synthetic sample is identical to any input ...")
    n_identical = 0
    for s_idx, synth_dos in enumerate(synth_dosages):
        for r_idx in range(ref_sample_dosages.shape[1]):
            ref_dos = ref_sample_dosages[:, r_idx]
            overlap = (synth_dos != MISSING) & (ref_dos != MISSING)
            if overlap.sum() == 0:
                continue
            frac_identical = (synth_dos[overlap] == ref_dos[overlap]).mean()
            if frac_identical > _IDENTITY_THRESHOLD:
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


def _read_reference(
    vcf_path: Path, chromosome: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read reference VCF in a single pass.

    Returns
    -------
    ref_afs : np.ndarray, shape (n_variants,), float64
        Per-site allele frequency.
    ref_positions : np.ndarray, shape (n_variants,), int64
        Physical positions (1-based).
    ref_sample_dosages : np.ndarray, shape (n_variants, n_samples), uint8
        Dosage matrix for all donor samples.
    """
    from cyvcf2 import VCF
    vcf = VCF(str(vcf_path))
    afs = []
    positions = []
    rows = []
    for v in vcf(chromosome):
        gts = np.array(v.genotypes, dtype=np.int8)  # (n_samples, 3)
        a1, a2 = gts[:, 0], gts[:, 1]
        valid_mask = (a1 >= 0) & (a2 >= 0)
        n_valid_alleles = valid_mask.sum() * 2
        if n_valid_alleles == 0:
            continue
        af = float(((a1 > 0) & valid_mask).sum() + ((a2 > 0) & valid_mask).sum()) / n_valid_alleles
        dosages = np.where(valid_mask, (a1 > 0).astype(np.uint8) + (a2 > 0).astype(np.uint8), MISSING)
        afs.append(af)
        positions.append(v.POS)
        rows.append(dosages.astype(np.uint8))
    vcf.close()

    if not positions:
        return np.array([]), np.array([], dtype=np.int64), np.empty((0, 0), dtype=np.uint8)

    return (
        np.array(afs, dtype=np.float64),
        np.array(positions, dtype=np.int64),
        np.stack(rows, axis=0),  # (n_variants, n_samples)
    )


def _read_dosages(
    vcf_path: Path, chromosome: str, target_positions: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Read dosages for a single-sample VCF at target_positions."""
    from cyvcf2 import VCF
    pos_set = set(target_positions)
    vcf = VCF(str(vcf_path))
    dosages = []
    positions = []
    for v in vcf(chromosome):
        if v.POS not in pos_set:
            continue
        dosages.append(_gt_to_dosage(v.genotypes[0]))
        positions.append(v.POS)
    vcf.close()
    return np.array(positions, dtype=np.int64), np.array(dosages, dtype=np.uint8)
