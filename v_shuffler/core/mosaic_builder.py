"""
Mosaic builder for v-shuffler.

Assembles synthetic diploid genotypes by applying pre-computed segment plans
to a GenotypePool chunk.

Each synthetic individual is built from a single segment plan: a list of
Segment objects that together cover the full chromosome.  For each segment
the builder copies the dosage values from the assigned donor individual,
producing a synthetic genotype that is a mosaic of multiple real donors.

Because segment plans are defined in cM space and the pool carries cM
coordinates for every variant, the builder simply masks by cM interval —
this works identically whether the pool is a chunk of 50k variants or the
full chromosome in one shot.
"""

from __future__ import annotations

import numpy as np

from v_shuffler.core.genotype_pool import MISSING, GenotypePool
from v_shuffler.core.recombination import Segment


def apply_segment_plan(
    pool: GenotypePool,
    plan: list[Segment],
) -> np.ndarray:
    """
    Build a synthetic dosage vector for one output individual from a pool chunk.

    Variants not covered by any segment (which should not happen with a valid
    plan, but can arise at chunk edges due to floating-point cM interpolation)
    are filled with MISSING.

    Parameters
    ----------
    pool : GenotypePool
        A chunk of real genotype data.
    plan : list[Segment]
        Segment plan for one synthetic individual.

    Returns
    -------
    np.ndarray, shape (pool.n_variants,), dtype uint8
        Synthetic dosage values (0/1/2 or MISSING=255).
    """
    result = np.full(pool.n_variants, MISSING, dtype=np.uint8)
    for seg in plan:
        # Inclusive start, inclusive end for the last segment
        mask = (pool.cm_pos >= seg.cm_start) & (pool.cm_pos <= seg.cm_end)
        result[mask] = pool.dosages[mask, seg.sample_idx]
    return result


def build_synthetic_genotypes(
    pool: GenotypePool,
    segment_plans: list[list[Segment]],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """
    Apply all segment plans to one pool chunk, producing the full synthetic
    genotype matrix for that chunk.

    Parameters
    ----------
    pool : GenotypePool
        A chunk of real genotype data for n_pool_samples donor individuals.
    segment_plans : list[list[Segment]]
        One segment plan per synthetic output individual.

    Returns
    -------
    dosages : np.ndarray, shape (pool.n_variants, n_output_samples), dtype uint8
        Synthetic dosage matrix.  Column i corresponds to synthetic individual i.
    format_fields : dict[str, np.ndarray]
        For each field name in ``pool.format_fields``, a float32 array of shape
        ``(pool.n_variants, n_output_samples)`` containing the value copied from
        whichever donor was assigned to each variant's segment.  NaN where the
        source donor had a missing value.  Empty dict when no FORMAT fields were
        requested.
    """
    n_output = len(segment_plans)
    dosage_result = np.empty((pool.n_variants, n_output), dtype=np.uint8)

    # Initialise format-field output arrays ("." = not yet assigned).
    # Arrays use object dtype to carry VCF-ready strings for both single-value
    # fields (AF → "0.4531") and multi-value fields (AD → "1904,3028").
    field_results: dict[str, np.ndarray] = {
        name: np.full((pool.n_variants, n_output), ".", dtype=object)
        for name in pool.format_fields
    }

    for s_idx, plan in enumerate(segment_plans):
        dosage_result[:, s_idx] = apply_segment_plan(pool, plan)
        # Copy format-field values for each segment in this plan
        for seg in plan:
            mask = (pool.cm_pos >= seg.cm_start) & (pool.cm_pos <= seg.cm_end)
            for name, src in pool.format_fields.items():
                field_results[name][mask, s_idx] = src[mask, seg.sample_idx]

    return dosage_result, field_results
