"""
GenotypePool: the in-memory representation of a chunk of variants
across all donor individuals.

The pool stores dosages (0/1/2 for biallelic sites) as a 2D uint8 array
of shape (n_variants, n_samples).  Missing genotypes are represented by
the sentinel value MISSING = 255.

Multi-allelic sites are handled by storing the sum of all non-ref allele
copies (i.e. number of alt alleles regardless of which alt), which is
the same as the biallelic dosage convention extended to multi-allelic
sites.  This is sufficient for the shuffling operation since we only need
to reproduce plausible genotype calls, not infer full allele composition
from dosage.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


MISSING: int = 255


@dataclass
class VariantInfo:
    """Metadata for a single variant site."""

    chrom: str
    pos: int         # 1-based
    ref: str
    alts: list[str]
    id: str
    qual: float | None
    filters: list[str]
    cm_pos: float    # computed from genetic map; 0.0 if not set


@dataclass
class GenotypePool:
    """
    A chunk of genotype data for all donor individuals at a set of variants.

    Attributes
    ----------
    dosages : np.ndarray, shape (n_variants, n_samples), dtype uint8
        Dosage values: 0 = hom-ref, 1 = het, 2 = hom-alt, 255 = missing.
    positions : np.ndarray, shape (n_variants,), dtype int64
        Physical positions in bp (1-based).
    cm_pos : np.ndarray, shape (n_variants,), dtype float64
        Genetic positions in cM.
    variant_info : list[VariantInfo]
        Per-variant metadata (same length as n_variants).
    """

    dosages: np.ndarray
    positions: np.ndarray
    cm_pos: np.ndarray
    variant_info: list[VariantInfo]

    def __post_init__(self) -> None:
        n = len(self.variant_info)
        assert self.dosages.shape[0] == n, (
            f"dosages first dimension {self.dosages.shape[0]} != "
            f"len(variant_info) {n}"
        )
        assert self.positions.shape == (n,)
        assert self.cm_pos.shape == (n,)

    @property
    def n_variants(self) -> int:
        return len(self.variant_info)

    @property
    def n_samples(self) -> int:
        return self.dosages.shape[1]
