"""
Genetic map parser for v-shuffler.

Supports two formats:
  1. SHAPEIT5 / Oxford format (3 columns, space-delimited, with header):
         pos  chr  cM
     Example:
         pos chr cM
         10177 chr22 0.0
         60000 chr22 0.012

  2. HapMap format (4 columns, with header):
         Chromosome  Position(bp)  Rate(cM/Mb)  Map(cM)
     Example:
         Chromosome Position(bp) Rate(cM/Mb) Map(cM)
         chr22 10177 0.0 0.0
         chr22 60000 1.2 0.012

Physical positions are 1-based (as in VCF).
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np


class GeneticMap:
    """
    Holds a genetic map for one chromosome.

    Attributes
    ----------
    positions : np.ndarray, shape (N,), int64
        Physical positions in bp (1-based).
    cm_values : np.ndarray, shape (N,), float64
        Cumulative genetic position in cM.
    """

    def __init__(self, path: Path, chrom: str) -> None:
        """
        Load the genetic map for *chrom* from *path*.

        Parameters
        ----------
        path : Path
            Path to the genetic map file (plain text or .gz).
        chrom : str
            Chromosome name to extract (e.g. "chr22" or "22").
            The loader accepts both "chr22" and "22" spellings and will
            match whichever the file uses.
        """
        self.path = path
        self.chrom = chrom
        self.positions, self.cm_values = self._load(path, chrom)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def bp_to_cm(self, positions: np.ndarray) -> np.ndarray:
        """
        Linearly interpolate bp positions to cM values.

        Positions outside the map range are extrapolated by numpy.interp
        (clamped to the boundary cM value).

        Parameters
        ----------
        positions : np.ndarray, int-like
            Physical positions in bp.

        Returns
        -------
        np.ndarray, float64, same shape as *positions*.
        """
        return np.interp(
            np.asarray(positions, dtype=np.float64),
            self.positions.astype(np.float64),
            self.cm_values,
        )

    @property
    def total_length_cm(self) -> float:
        """Total genetic length of the loaded region in cM."""
        return float(self.cm_values[-1] - self.cm_values[0])

    @property
    def start_cm(self) -> float:
        return float(self.cm_values[0])

    @property
    def end_cm(self) -> float:
        return float(self.cm_values[-1])

    # ------------------------------------------------------------------
    # Internal loading
    # ------------------------------------------------------------------

    @staticmethod
    def _normalise_chrom(name: str) -> str:
        """Strip 'chr' prefix for comparison."""
        return name.lower().lstrip("chr")

    @classmethod
    def _load(
        cls, path: Path, chrom: str
    ) -> tuple[np.ndarray, np.ndarray]:
        target = cls._normalise_chrom(chrom)

        opener = gzip.open if str(path).endswith(".gz") else open
        positions: list[int] = []
        cm_values: list[float] = []
        fmt: str | None = None  # "shapeit5" or "hapmap"

        with opener(path, "rt") as fh:
            for lineno, line in enumerate(fh):
                line = line.strip()
                if not line:
                    continue

                parts = line.split()

                # Detect format from header
                if lineno == 0:
                    fmt = cls._detect_format(parts)
                    continue  # skip header

                if fmt == "shapeit5":
                    # pos  chr  cM
                    pos_val, chr_val, cm_val = int(parts[0]), parts[1], float(parts[2])
                    if cls._normalise_chrom(chr_val) != target:
                        continue
                elif fmt == "hapmap":
                    # Chromosome  Position(bp)  Rate(cM/Mb)  Map(cM)
                    chr_val, pos_val, _, cm_val = (
                        parts[0], int(parts[1]), parts[2], float(parts[3])
                    )
                    if cls._normalise_chrom(chr_val) != target:
                        continue
                else:
                    raise ValueError(f"Unrecognised genetic map format in {path}")

                positions.append(pos_val)
                cm_values.append(cm_val)

        if not positions:
            raise ValueError(
                f"No entries found for chromosome {chrom!r} in {path}. "
                "Check that the chromosome name matches the file."
            )

        pos_arr = np.array(positions, dtype=np.int64)
        cm_arr = np.array(cm_values, dtype=np.float64)

        # Ensure sorted by position
        order = np.argsort(pos_arr)
        pos_arr = pos_arr[order]
        cm_arr = cm_arr[order]

        # Validate monotonicity of cM values
        if not np.all(np.diff(cm_arr) >= 0):
            raise ValueError(
                f"Genetic map for {chrom} has non-monotonic cM values. "
                "Check the file for errors."
            )

        return pos_arr, cm_arr

    @staticmethod
    def _detect_format(header_parts: list[str]) -> str:
        """Infer file format from header tokens."""
        h = [p.lower() for p in header_parts]
        if len(h) >= 3 and h[0] == "pos" and h[2] == "cm":
            return "shapeit5"
        if len(h) >= 4 and "position" in h[1] and "map" in h[3]:
            return "hapmap"
        # Fallback: try to detect by column count and content
        if len(h) == 3:
            return "shapeit5"
        if len(h) == 4:
            return "hapmap"
        raise ValueError(
            f"Cannot detect genetic map format from header: {header_parts!r}. "
            "Expected SHAPEIT5 (pos chr cM) or HapMap (Chromosome Position(bp) Rate(cM/Mb) Map(cM)) format."
        )
