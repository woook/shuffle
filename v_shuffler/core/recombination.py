"""
Recombination simulation for v-shuffler.

The core idea: simulate meiotic crossovers along a chromosome using a
Poisson process scaled to the genetic map, then assign a donor individual
to each inter-breakpoint segment.  Downstream code copies diploid genotypes
from the assigned donor for each segment, producing a synthetic individual
that is a mosaic of multiple donors at genetically realistic break-points.

No phasing is assumed or required: the unit being swapped is a full diploid
genotype call (dosage 0/1/2), not a single allele.

For targeted NGS panels, where the captured fraction of any chromosome is too
small for meaningful Poisson-based crossover placement, use the region-sampling
functions instead: detect_regions / build_region_segment_plan /
generate_all_region_plans.  Each captured region is treated as an independent
mixing unit and assigned a randomly-chosen donor, guaranteeing thorough mixing
regardless of how sparse the panel is.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from v_shuffler.io.genetic_map import GeneticMap

# Safety limit for the min_donors enforcement loop in generate_all_segment_plans().
# The adjacency constraint in build_segment_plan() can cause donor cycling
# (e.g. A→B→A→B...) that prevents reaching the requested min_donors count.
# If this limit is reached, behavior depends on strict mode:
#   - strict=True (default): Raises RuntimeError with diagnostic message
#   - strict=False: Logs warning and accepts the achieved donor count
# Value chosen to allow ~1000 extra crossover attempts before giving up.
MAX_MIN_DONOR_ITERATIONS = 1000


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
    min_donors: int = 1,
    strict: bool = True,
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
    min_donors : int
        Minimum number of distinct donors per synthetic.  If the Poisson draw
        produces fewer breakpoints than required, extra breakpoints are added
        uniformly at random within the map range.
    strict : bool
        If True (default), raises RuntimeError when min_donors cannot be achieved.
        If False, emits warning and continues with partial mixing.

    Returns
    -------
    list[list[Segment]]
        One segment plan per output sample.
        segment_plans[i] is the list of Segment objects for synthetic individual i.

    Raises
    ------
    RuntimeError
        If strict=True and min_donors cannot be achieved after MAX_MIN_DONOR_ITERATIONS.
    """
    effective_min = min(min_donors, n_pool_samples)
    plans: list[list[Segment]] = []
    failures: list[tuple[int, int]] = []  # (synthetic_idx, achieved_count)

    for synth_idx in range(n_output_samples):
        breakpoints = simulate_crossover_breakpoints(genetic_map, rng, lambda_override)
        # Ensure enough breakpoints for effective_min segments.
        required_bps = max(0, effective_min - 1)
        if len(breakpoints) < required_bps:
            n_extra = required_bps - len(breakpoints)
            extra = rng.uniform(genetic_map.start_cm, genetic_map.end_cm, size=n_extra)
            breakpoints = np.sort(np.concatenate([breakpoints, extra]))
        plan = build_segment_plan(breakpoints, genetic_map, n_pool_samples, rng)
        # The adjacency-only constraint in build_segment_plan can cause donor cycling
        # (e.g. A→B→A with 3 segments).  Add extra breakpoints one at a time until
        # the distinct-donor count reaches effective_min.
        iteration_count = 0
        while len({seg.sample_idx for seg in plan}) < effective_min:
            if iteration_count >= MAX_MIN_DONOR_ITERATIONS:
                achieved = len({seg.sample_idx for seg in plan})
                failures.append((synth_idx, achieved))
                break

            extra = rng.uniform(genetic_map.start_cm, genetic_map.end_cm, size=1)
            breakpoints = np.sort(np.concatenate([breakpoints, extra]))
            plan = build_segment_plan(breakpoints, genetic_map, n_pool_samples, rng)
            iteration_count += 1
        plans.append(plan)

    # Handle failures: raise exception or emit aggregate warning
    if failures:
        min_achieved = min(count for _, count in failures)
        max_achieved = max(count for _, count in failures)
        failure_msg = (
            f"{len(failures)}/{n_output_samples} synthetic individuals failed to achieve min_donors={effective_min}. "
            f"Achieved donors ranged from {min_achieved} to {max_achieved} "
            f"(pool size: {n_pool_samples}, iterations: {MAX_MIN_DONOR_ITERATIONS})."
        )

        if strict:
            raise RuntimeError(
                f"{failure_msg}\n"
                f"This indicates the requested min_donors is too high for the available pool size "
                f"or chromosome length. Please:\n"
                f"  1. Reduce --min-donors (try {max_achieved})\n"
                f"  2. Increase donor pool size (need at least {effective_min} donors)\n"
                f"  3. Use --allow-partial-mixing to emit warning instead (privacy trade-off)\n"
                f"  4. Use --no-region-sampling if currently in region mode"
            )
        else:
            logger = logging.getLogger(__name__)
            logger.error(
                "%s Output quality is COMPROMISED. "
                "Consider reducing --min-donors or increasing pool size.",
                failure_msg,
            )
    return plans


# ---------------------------------------------------------------------------
# Region-sampling functions (for targeted NGS panels)
# ---------------------------------------------------------------------------


def detect_regions(
    positions: np.ndarray,
    gap_threshold_bp: int = 10_000,
) -> list[tuple[int, int]]:
    """
    Identify contiguous captured regions from a sorted array of bp positions.

    A new region is started whenever the gap between consecutive positions
    exceeds *gap_threshold_bp*.

    Parameters
    ----------
    positions : np.ndarray
        Sorted array of int64 base-pair positions.
    gap_threshold_bp : int
        Minimum inter-variant gap (bp) that separates two regions.

    Returns
    -------
    list[tuple[int, int]]
        List of (start_bp, end_bp) tuples, one per region.
        Empty list if *positions* is empty.
    """
    if len(positions) == 0:
        return []

    regions: list[tuple[int, int]] = []
    region_start = int(positions[0])
    prev = int(positions[0])

    for i in range(1, len(positions)):
        cur = int(positions[i])
        if cur - prev > gap_threshold_bp:
            regions.append((region_start, prev))
            region_start = cur
        prev = cur

    regions.append((region_start, prev))
    return regions


def build_region_segment_plan(
    regions_cm: list[tuple[float, float]],
    n_samples: int,
    rng: np.random.Generator,
    min_donors: int = 1,
) -> list[Segment]:
    """
    Build a segment plan for one synthetic individual using region-based sampling.

    Each captured region is assigned one donor independently.  The first
    ``min_donors`` regions (up to the number of available donors and the
    number of regions) are filled with distinct donors; the rest are drawn
    freely subject only to the adjacency constraint (no two consecutive
    regions share the same donor).

    Parameters
    ----------
    regions_cm : list[tuple[float, float]]
        (start_cM, end_cM) for each captured region.
    n_samples : int
        Number of donor individuals available.
    rng : np.random.Generator
    min_donors : int
        Minimum number of distinct donors across all regions.

    Returns
    -------
    list[Segment]
        One Segment per region.  Empty list if *regions_cm* is empty.
    """
    if not regions_cm:
        return []
    if n_samples < 1:
        raise ValueError("n_samples must be >= 1")

    n_regions = len(regions_cm)

    if n_samples == 1:
        return [
            Segment(cm_start=float(s), cm_end=float(e), sample_idx=0)
            for s, e in regions_cm
        ]

    effective_min = min(min_donors, n_samples, n_regions)

    segments: list[Segment] = []
    used_donors: list[int] = []
    prev: int | None = None

    for i, (cm_start, cm_end) in enumerate(regions_cm):
        if i < effective_min:
            # Without-replacement phase: pick a donor not yet used and not prev
            unused = [d for d in range(n_samples) if d not in used_donors]
            candidates = [d for d in unused if d != prev]
            if not candidates:
                # Relax adjacency constraint (still without replacement)
                candidates = unused if unused else list(range(n_samples))
            chosen = int(rng.choice(candidates))
            used_donors.append(chosen)
        else:
            # Free-sampling phase: adjacency constraint only
            candidates_arr = np.arange(n_samples)
            if prev is not None:
                candidates_arr = candidates_arr[candidates_arr != prev]
            chosen = int(rng.choice(candidates_arr))

        segments.append(Segment(cm_start=float(cm_start), cm_end=float(cm_end), sample_idx=chosen))
        prev = chosen

    return segments


def generate_all_region_plans(
    n_output_samples: int,
    regions_cm: list[tuple[float, float]],
    n_pool_samples: int,
    rng: np.random.Generator,
    min_donors: int = 1,
) -> list[list[Segment]]:
    """
    Pre-generate region-based segment plans for all synthetic output samples.

    Direct analogue of ``generate_all_segment_plans`` for region mode.

    Parameters
    ----------
    n_output_samples : int
        Number of synthetic individuals to create.
    regions_cm : list[tuple[float, float]]
        (start_cM, end_cM) for each captured region (from detect_regions +
        bp_to_cm conversion).
    n_pool_samples : int
        Number of real donor individuals in the pool.
    rng : np.random.Generator
    min_donors : int
        Minimum number of distinct donors per synthetic individual.

    Returns
    -------
    list[list[Segment]]
        One segment plan (list of Segment) per output sample.
    """
    return [
        build_region_segment_plan(regions_cm, n_pool_samples, rng, min_donors=min_donors)
        for _ in range(n_output_samples)
    ]
