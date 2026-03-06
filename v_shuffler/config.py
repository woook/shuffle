"""Configuration dataclass for v-shuffler."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ShufflerConfig:
    """All configuration for a single shuffle run."""

    # Input/output
    input_vcfs: list[Path]
    output_dir: Path
    genetic_map: Path
    chromosome: str

    # Synthetic sample parameters
    n_output_samples: int
    seed: int | None = None

    # Recombination parameters
    min_segment_cM: float = 0.5
    n_crossovers_lambda: float | None = None  # None = derived from map length

    # Processing
    chunk_size_variants: int = 50_000
    n_threads: int = 4
    max_missing_rate: float = 0.05  # skip variants with more missing calls

    # Output format
    output_mode: str = "per_sample"  # "per_sample" or "multi_sample"

    def __post_init__(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        if self.output_mode not in ("per_sample", "multi_sample"):
            raise ValueError(
                f"output_mode must be 'per_sample' or 'multi_sample', "
                f"got {self.output_mode!r}"
            )
        if self.n_output_samples < 1:
            raise ValueError("n_output_samples must be >= 1")
        if not 0.0 <= self.max_missing_rate <= 1.0:
            raise ValueError("max_missing_rate must be in [0, 1]")
