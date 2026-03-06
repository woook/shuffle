"""
Recombination simulation for v-shuffler.

The core idea: simulate meiotic crossovers along a chromosome using a
Poisson process scaled to the genetic map, then assign a donor individual
to each inter-breakpoint segment.  Downstream code copies diploid genotypes
from the assigned donor for each segment, producing a synthetic individual
that is a mosaic of multiple donors at genetically realistic break-points.

No phasing is assumed or required: the unit being swapped is a full diploid
genotype call (dosage 0/1/2), not a single allele.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from v_shuffler.io.genetic_map import GeneticMap


@dataclass(frozen=True)
class Segment:
    """
    A contiguous chromosomal segment assigned to a specific donor individual.

    Parameters
    ----------
    cm_start : float
        Start of the segment in cM (inclusive).
    cm_end : float
        End of the segment in cM (exclusive, except for the last segment
        where it equals the map end and is treated as inclusive).
    sample_idx : int
        Column index of the donor individual in the GenotypePool matrix.
    """

    cm_start: float
    cm_end: float
    sample_idx: int


def simulate_crossover_breakpoints(
    genetic_map: GeneticMap,
    rng: np.random.Generator,
    lambda_override: float | None = None,
) -> np.ndarray:
    """
    Sample crossover positions along a chromosome via a Poisson process.

    The number of crossovers is drawn from Poisson(λ) where by default
    λ = total_length_cM / 100 (one crossover per Morgan per meiosis).
    Crossover positions are then drawn uniformly in [start_cM, end_cM],
    which is equivalent to a uniform process in genetic distance (the
    genetic map already encodes the spatially varying rate in physical
    coordinates).

    Parameters
    ----------
    genetic_map : GeneticMap
    rng : np.random.Generator
    lambda_override : float, optional
        Override the Poisson rate.  Useful for testing or for chromosomes
        where the default rate is inappropriate.

    Returns
    -------
    np.ndarray, float64, shape (k,)
        Sorted crossover positions in cM.  May be empty (no crossovers).
    """
    lam = lambda_override if lambda_override is not None else (
        genetic_map.total_length_cm / 100.0
    )
    n_crossovers: int = int(rng.poisson(lam))
    if n_crossovers == 0:
        return np.array([], dtype=np.float64)

    breakpoints = rng.uniform(
        genetic_map.start_cm,
        genetic_map.end_cm,
        size=n_crossovers,
    )
    return np.sort(breakpoints)


def build_segment_plan(
    breakpoints_cm: np.ndarray,
    genetic_map: GeneticMap,
    n_samples: int,
    rng: np.random.Generator,
    initial_sample: int | None = None,
) -> list[Segment]:
    """
    Build a segment plan for one synthetic (haploid-concept) track.

    The plan covers the full chromosome from map start to map end.
    At each crossover breakpoint the donor switches to a different
    individual chosen uniformly at random from the pool (never the
    same individual twice in a row).

    Parameters
    ----------
    breakpoints_cm : np.ndarray
        Sorted crossover positions in cM (may be empty).
    genetic_map : GeneticMap
    n_samples : int
        Total number of donor individuals available.
    rng : np.random.Generator
    initial_sample : int, optional
        Index of the first donor.  Drawn randomly if not provided.

    Returns
    -------
    list[Segment]
        Non-overlapping segments that together cover
        [genetic_map.start_cm, genetic_map.end_cm].
    """
    if n_samples < 1:
        raise ValueError("n_samples must be >= 1")

    start = genetic_map.start_cm
    end = genetic_map.end_cm

    # Choose initial donor
    current = (
        int(rng.integers(n_samples)) if initial_sample is None else initial_sample
    )

    boundaries = np.concatenate([[start], breakpoints_cm, [end]])
    segments: list[Segment] = []

    for i in range(len(boundaries) - 1):
        seg_start = float(boundaries[i])
        seg_end = float(boundaries[i + 1])

        # Skip zero-width segments that can arise from duplicate breakpoints
        if seg_end <= seg_start:
            continue

        segments.append(Segment(cm_start=seg_start, cm_end=seg_end, sample_idx=current))

        # Switch donor for the next segment (ensure it's different when possible)
        if n_samples > 1:
            candidates = np.arange(n_samples)
            candidates = candidates[candidates != current]
            current = int(rng.choice(candidates))
        # If only one sample in the pool, we can't switch — keep the same

    return segments


def generate_all_segment_plans(
    n_output_samples: int,
    genetic_map: GeneticMap,
    n_pool_samples: int,
    rng: np.random.Generator,
    lambda_override: float | None = None,
) -> list[list[Segment]]:
    """
    Pre-generate segment plans for all synthetic output samples.

    This is done upfront (before streaming variants) so the full
    chromosome plan is available when processing each variant chunk.
    The plans are tiny (a few dozen segments per individual at most)
    and live entirely in cM space.

    Parameters
    ----------
    n_output_samples : int
        Number of synthetic individuals to create.
    genetic_map : GeneticMap
    n_pool_samples : int
        Number of real donor individuals in the pool.
    rng : np.random.Generator
    lambda_override : float, optional

    Returns
    -------
    list[list[Segment]]
        One segment plan per output sample.
        segment_plans[i] is the list of Segment objects for synthetic individual i.
    """
    plans: list[list[Segment]] = []
    for _ in range(n_output_samples):
        breakpoints = simulate_crossover_breakpoints(genetic_map, rng, lambda_override)
        plan = build_segment_plan(breakpoints, genetic_map, n_pool_samples, rng)
        plans.append(plan)
    return plans
