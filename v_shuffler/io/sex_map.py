"""
Sex-map loading and donor-pool filtering for sex chromosome handling.

Provides three things:
  1. ``load_sex_map``  — parse a two-column donor-sex file.
  2. ``filter_vcfs_by_sex``  — return only the VCF paths whose sex matches.
  3. ``sex_filter_for_chromosome``  — infer which sex to keep from a chromosome
     name (chrX → female only, chrY → male only, autosomes → no filter).
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Accepted chromosome names for X and Y (lower-cased for comparison).
_CHRX_NAMES = {"chrx", "x"}
_CHRY_NAMES = {"chry", "y"}

# Accepted sex-label strings (lower-cased).
# PLINK convention: 1 = male, 2 = female.
_FEMALE_LABELS = {"f", "female", "2"}
_MALE_LABELS = {"m", "male", "1"}


def parse_sex_label(raw: str) -> str:
    """
    Normalise a sex label to ``'F'`` (female) or ``'M'`` (male).

    Accepted values (case-insensitive):
      Female — ``F``, ``female``, ``2``
      Male   — ``M``, ``male``,   ``1``

    Raises ``ValueError`` for anything else.
    """
    normalised = raw.strip().lower()
    if normalised in _FEMALE_LABELS:
        return "F"
    if normalised in _MALE_LABELS:
        return "M"
    raise ValueError(
        f"Unknown sex label {raw!r}. "
        "Accepted values: F/female/2 (female), M/male/1 (male)."
    )


def load_sex_map(sex_file: Path, vcf_paths: list[Path]) -> dict[Path, str]:
    """
    Parse a donor-sex file and return a ``{vcf_path: 'F'|'M'}`` mapping.

    File format
    -----------
    Two whitespace-separated columns, one donor per line::

        /path/to/sample1.vcf.gz  F
        sample2.vcf.gz           M

    Lines that start with ``#`` or are blank are ignored.  An optional header
    line (e.g. ``path sex``) is automatically skipped if the second token is
    not a recognised sex label.

    Path matching
    -------------
    Exact string match is attempted first; if that fails, the file *basename*
    is compared.  VCF paths not found in the sex file are omitted from the
    result with a warning.

    Parameters
    ----------
    sex_file :
        Path to the two-column sex file.
    vcf_paths :
        The full list of donor VCF paths supplied to the run.  Used for
        matching and for emitting per-sample warnings.

    Returns
    -------
    dict[Path, str]
        Maps each recognised VCF path to ``'F'`` or ``'M'``.
    """
    # First pass: build raw string→sex lookup from the file.
    raw_map: dict[str, str] = {}
    with sex_file.open() as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(
                    f"{sex_file}:{lineno}: expected at least 2 columns, "
                    f"got {len(parts)} in line {line!r}"
                )
            path_str, sex_raw = parts[0], parts[1]
            # Skip an explicit header row (e.g. "path sex") on line 1 only, and
            # only when the path column doesn't look like a real VCF path.  A
            # typo in the sex label of a real donor row on line 1 would otherwise
            # be silently swallowed, producing a confusing "no donors found" error
            # later rather than pointing at the bad line.
            try:
                sex = parse_sex_label(sex_raw)
            except ValueError:
                looks_like_path = path_str.endswith(".vcf.gz") or "/" in path_str
                if lineno == 1 and not looks_like_path:
                    continue  # genuine header row, skip silently
                raise ValueError(
                    f"{sex_file}:{lineno}: unrecognised sex label {sex_raw!r}"
                ) from None
            raw_map[path_str] = sex

    # Second pass: match each supplied VCF path against the raw map.
    # Detect duplicate basenames up front so we can warn on ambiguous fallback.
    basenames = [p.name for p in vcf_paths]
    duplicate_basenames = {b for b in basenames if basenames.count(b) > 1}

    result: dict[Path, str] = {}
    for vcf_path in vcf_paths:
        full_key = str(vcf_path)
        base_key = vcf_path.name
        if full_key in raw_map:
            result[vcf_path] = raw_map[full_key]
        elif base_key in raw_map:
            if base_key in duplicate_basenames:
                logger.warning(
                    "Ambiguous basename match for %s — multiple donor VCFs share "
                    "this filename; use full paths in the sex file to avoid "
                    "incorrect sex assignment.",
                    base_key,
                )
            result[vcf_path] = raw_map[base_key]
        else:
            logger.warning(
                "Donor VCF %s not found in sex file %s — "
                "it will be excluded from sex-filtered chromosome runs.",
                vcf_path.name,
                sex_file.name,
            )

    n_matched = len(result)
    if n_matched < len(vcf_paths):
        logger.warning(
            "%d / %d donor VCFs matched in sex file %s.",
            n_matched, len(vcf_paths), sex_file.name,
        )

    return result


def filter_vcfs_by_sex(
    vcf_paths: list[Path],
    sex_map: dict[Path, str],
    keep_sex: str,
) -> list[Path]:
    """
    Return only those VCF paths whose sex in *sex_map* equals *keep_sex*.

    Preserves the original ordering of *vcf_paths*.

    Parameters
    ----------
    vcf_paths :
        Full donor VCF list.
    sex_map :
        Mapping returned by :func:`load_sex_map`.
    keep_sex :
        ``'F'`` to keep female donors, ``'M'`` to keep male donors.
    """
    return [p for p in vcf_paths if sex_map.get(p) == keep_sex]


def sex_filter_for_chromosome(chromosome: str) -> str | None:
    """
    Return the sex to keep for a given chromosome, or ``None`` for autosomes.

    Returns
    -------
    ``'F'``
        For chrX / X — use female donors only.
    ``'M'``
        For chrY / Y — use male donors only.
    ``None``
        For all other chromosomes — use the full donor pool.
    """
    normalised = chromosome.strip().lower()
    if normalised in _CHRX_NAMES:
        return "F"
    if normalised in _CHRY_NAMES:
        return "M"
    return None
