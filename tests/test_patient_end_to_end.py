"""
Patient end-to-end integration test for v-shuffler.

Exercises the full pipeline on real patient VCFs and verifies the three most
critical anonymity and plausibility properties:

  1. Output produced  — correct number of synthetic VCF files written.
  2. No identity leak — max pairwise concordance (synth × donor) < 0.99.
  3. AF preserved     — Pearson r ≥ 0.99 between donor and synthetic AFs.
  4. Variant count    — all synth positions are present in the donor pool.

Running these tests
-------------------
Set the following environment variables before invoking pytest:

  VSHUFFLE_PATIENT_VCFS   glob or @filelist of per-sample donor VCFs (required)
  VSHUFFLE_GENETIC_MAP    path to genetic map, SHAPEIT5 or HapMap format (required)
  VSHUFFLE_CHROMOSOME     chromosome name, e.g. "chr22" (required)
  VSHUFFLE_SEED           random seed for the shuffle run (optional, default 42)
  VSHUFFLE_N_SYNTH        number of synthetics to produce
                          (optional, default min(n_donors, 50))

Example:

  VSHUFFLE_PATIENT_VCFS="@/data/patients.txt" \\
  VSHUFFLE_GENETIC_MAP="/data/chr22.b38.gmap.gz" \\
  VSHUFFLE_CHROMOSOME="chr22" \\
  pytest tests/test_patient_end_to_end.py -v -s

Without the required env vars, all four tests are skipped cleanly.
"""

from __future__ import annotations

import dataclasses
import os
from pathlib import Path

import numpy as np
import pytest
from click import BadParameter
from click.testing import CliRunner

from v_shuffler.cli import _resolve_inputs, main
from v_shuffler.core.genotype_pool import MISSING
from v_shuffler.io.genetic_map import GeneticMap
from v_shuffler.io.vcf_reader import PerSampleVCFReader

# Maximum number of variants used for concordance computation.
# Subsampling keeps the (V, S, D) broadcast array to ≈ 100 MB for typical
# cohort sizes (S ≤ 50, D ≤ 200).
MAX_CHECK_VARIANTS = 10_000


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _require_env(var: str) -> str:
    """Return the value of *var*; skip the calling test/fixture if unset."""
    val = os.environ.get(var)
    if not val:
        pytest.skip(
            f"Environment variable {var} is not set. "
            "Set it to enable the patient end-to-end test."
        )
    return val


def _require_path_env(var: str) -> Path:
    """Return a Path from *var*; skip if unset or non-existent."""
    val = _require_env(var)
    p = Path(val)
    if not p.exists():
        pytest.skip(f"{var}={val!r} does not exist.")
    return p


def _load_full_matrix(
    vcf_paths: list[Path],
    chromosome: str,
    gmap: GeneticMap,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Stream all chunks from per-sample VCFs into a single (V, N) dosage matrix.

    Uses ``max_missing_rate=1.0`` so no variants are dropped — the caller is
    responsible for any position-level filtering.

    Returns
    -------
    matrix : np.ndarray, shape (V, N), uint8
    positions : np.ndarray, shape (V,), int64
    """
    reader = PerSampleVCFReader(
        vcf_paths=vcf_paths,
        chromosome=chromosome,
        genetic_map=gmap,
        chunk_size=50_000,
        max_missing_rate=1.0,
    )
    dosage_chunks: list[np.ndarray] = []
    position_chunks: list[np.ndarray] = []
    for pool in reader.iter_chunks():
        dosage_chunks.append(pool.dosages)      # (V_chunk, N)
        position_chunks.append(pool.positions)  # (V_chunk,)

    if not dosage_chunks:
        return np.empty((0, len(vcf_paths)), dtype=np.uint8), np.array([], dtype=np.int64)

    return (
        np.concatenate(dosage_chunks, axis=0),
        np.concatenate(position_chunks, axis=0),
    )


def _compute_afs(matrix: np.ndarray) -> np.ndarray:
    """
    Per-variant allele frequency, ignoring MISSING sentinels.

    Returns NaN for variants where all calls are missing.
    """
    valid = matrix != MISSING
    denom = valid.sum(axis=1) * 2.0
    numer = np.where(valid, matrix.astype(np.float64), 0.0).sum(axis=1)
    with np.errstate(invalid="ignore"):
        return np.where(denom > 0, numer / denom, np.nan)


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class PatientRunResult:
    synth_paths: list[Path]
    n_expected_synth: int
    donor_afs: np.ndarray       # (V_common,)
    synth_afs: np.ndarray       # (V_common,)
    concordances: np.ndarray    # (S, D) float32 — synth vs donor, subsampled
    n_synth_variants: int       # total variants in synth VCFs
    n_common_variants: int      # positions present in both donor and synth VCFs


# ---------------------------------------------------------------------------
# Module-scoped fixture: run the shuffler once, load matrices
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def patient_run(tmp_path_factory: pytest.TempPathFactory) -> PatientRunResult:
    """
    Full pipeline fixture.

    1. Resolve donor paths from VSHUFFLE_PATIENT_VCFS.
    2. Invoke the shuffle CLI via CliRunner.
    3. Load donor and synthetic dosage matrices with PerSampleVCFReader.
    4. Intersect positions; subsample to MAX_CHECK_VARIANTS for concordance.
    5. Return a PatientRunResult for all four sanity-check tests.
    """
    # --- Env vars ---
    vcf_spec = _require_env("VSHUFFLE_PATIENT_VCFS")
    gmap_path = _require_path_env("VSHUFFLE_GENETIC_MAP")
    chromosome = _require_env("VSHUFFLE_CHROMOSOME")
    seed = int(os.environ.get("VSHUFFLE_SEED", "42"))

    try:
        donor_paths = _resolve_inputs(vcf_spec)
    except (BadParameter, SystemExit) as exc:
        pytest.skip(f"Could not resolve VSHUFFLE_PATIENT_VCFS={vcf_spec!r}: {exc}")

    n_donors = len(donor_paths)
    n_synth = int(os.environ.get("VSHUFFLE_N_SYNTH", str(min(n_donors, 50))))

    # --- Run shuffler ---
    work_dir = tmp_path_factory.mktemp("patient_work")
    filelist = work_dir / "donors.txt"
    filelist.write_text("\n".join(str(p) for p in donor_paths) + "\n")
    synth_dir = work_dir / "synth"

    result = CliRunner().invoke(
        main,
        [
            "shuffle",
            "--input", f"@{filelist}",
            "--output-dir", str(synth_dir),
            "--genetic-map", str(gmap_path),
            "--chromosome", chromosome,
            "--n-samples", str(n_synth),
            "--seed", str(seed),
            "--output-mode", "per_sample",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, (
        f"v-shuffler shuffle failed (exit {result.exit_code}):\n{result.output}"
    )

    # --- Collect output paths (exclude .tbi / .csi index files) ---
    synth_paths = sorted(
        p for p in synth_dir.glob("synthetic_*.vcf*")
        if not p.name.endswith((".tbi", ".csi"))
    )

    # --- Load genetic map ---
    gmap = GeneticMap(gmap_path, chromosome)

    # --- Load donor matrix (all variants, no missing filter) ---
    donor_matrix, donor_positions = _load_full_matrix(donor_paths, chromosome, gmap)

    # --- Load synth matrix ---
    synth_matrix, synth_positions = _load_full_matrix(synth_paths, chromosome, gmap)

    n_synth_variants = len(synth_positions)

    # --- Align positions: synth is a subset of donor after missing-rate filter ---
    if len(donor_positions) == 0 or len(synth_positions) == 0:
        pytest.skip("No variants loaded — check VSHUFFLE_CHROMOSOME matches the VCF content.")

    lookup_idx = np.searchsorted(donor_positions, synth_positions)
    lookup_idx = lookup_idx.clip(0, len(donor_positions) - 1)
    valid_mask = donor_positions[lookup_idx] == synth_positions

    donor_aligned = donor_matrix[lookup_idx[valid_mask], :]   # (V_common, D)
    synth_aligned = synth_matrix[valid_mask, :]               # (V_common, S)
    n_common_variants = int(valid_mask.sum())

    # --- Allele frequencies on common variants ---
    donor_afs = _compute_afs(donor_aligned)
    synth_afs = _compute_afs(synth_aligned)

    # --- Subsample variants for concordance (memory-bounded) ---
    rng = np.random.default_rng(seed)
    n_check = min(MAX_CHECK_VARIANTS, n_common_variants)
    check_idx = np.sort(rng.choice(n_common_variants, size=n_check, replace=False))

    donor_sub = donor_aligned[check_idx, :]   # (n_check, D)
    synth_sub = synth_aligned[check_idx, :]   # (n_check, S)

    # --- Vectorised concordance (S, D), MISSING-aware ---
    # valid[v, s, d] = True when neither synth nor donor has MISSING at variant v
    valid_sd = (synth_sub[:, :, np.newaxis] != MISSING) & (donor_sub[:, np.newaxis, :] != MISSING)
    match_sd = (synth_sub[:, :, np.newaxis] == donor_sub[:, np.newaxis, :]) & valid_sd
    n_valid = valid_sd.sum(axis=0).clip(min=1)  # (S, D)
    concordances = (match_sd.sum(axis=0) / n_valid).astype(np.float32)  # (S, D)

    return PatientRunResult(
        synth_paths=synth_paths,
        n_expected_synth=n_synth,
        donor_afs=donor_afs,
        synth_afs=synth_afs,
        concordances=concordances,
        n_synth_variants=n_synth_variants,
        n_common_variants=n_common_variants,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_shuffle_produced_output(patient_run: PatientRunResult) -> None:
    """The shuffler must write exactly n_synth synthetic VCF files."""
    n_found = len(patient_run.synth_paths)
    assert n_found == patient_run.n_expected_synth, (
        f"Expected {patient_run.n_expected_synth} synthetic VCF(s), "
        f"found {n_found} in the output directory."
    )


def test_no_identity_leak(patient_run: PatientRunResult) -> None:
    """No synthetic individual may share ≥ 99% of genotypes with any donor (P1)."""
    max_conc = float(patient_run.concordances.max())
    assert max_conc < 0.99, (
        f"Max pairwise concordance {max_conc:.4f} ≥ 0.99 — possible identity leak. "
        "Consider increasing the donor pool size or using lambda_override."
    )


def test_af_preserved(patient_run: PatientRunResult) -> None:
    """Allele-frequency Pearson r between donor and synthetic output must be ≥ 0.99 (B1)."""
    common = ~np.isnan(patient_run.donor_afs) & ~np.isnan(patient_run.synth_afs)
    if common.sum() < 10:
        pytest.skip("Too few overlapping variants with non-NaN AFs to compute correlation.")

    r = float(np.corrcoef(patient_run.donor_afs[common], patient_run.synth_afs[common])[0, 1])
    print(
        f"\n[B1] AF Pearson r={r:.4f} over {common.sum()} variants "
        f"(donor mean={patient_run.donor_afs[common].mean():.4f}, "
        f"synth mean={patient_run.synth_afs[common].mean():.4f})"
    )
    assert r >= 0.99, (
        f"AF correlation r={r:.4f} < 0.99. "
        "The synthetic output does not preserve the donor allele-frequency spectrum."
    )


def test_variant_count_consistent(patient_run: PatientRunResult) -> None:
    """All variant positions in the synthetic output must be present in the donor pool."""
    assert patient_run.n_common_variants == patient_run.n_synth_variants, (
        f"{patient_run.n_synth_variants - patient_run.n_common_variants} synthetic "
        f"variant position(s) are absent from the donor VCFs "
        f"(synth={patient_run.n_synth_variants}, "
        f"common with donors={patient_run.n_common_variants}). "
        "This suggests the synthetic and donor VCFs were generated from different inputs."
    )
    print(
        f"\n[B4] Variant count consistent: {patient_run.n_common_variants} positions "
        "present in both donor and synthetic VCFs."
    )
