"""
Command-line interface for v-shuffler.

Subcommands:
  shuffle   -- main recombination shuffling pipeline
  validate  -- sanity-check output VCFs against input
"""

from __future__ import annotations

import glob
import logging
import sys
from pathlib import Path

import click
import numpy as np
from tqdm import tqdm

from v_shuffler import __version__
from v_shuffler.config import ShufflerConfig
from v_shuffler.core.mosaic_builder import build_synthetic_genotypes
from v_shuffler.core.recombination import generate_all_segment_plans
from v_shuffler.io.genetic_map import GeneticMap
from v_shuffler.io.vcf_reader import PerSampleVCFReader
from v_shuffler.io.vcf_writer import SyntheticVCFWriter

logger = logging.getLogger("v_shuffler")


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        level=level,
        stream=sys.stderr,
    )


def _resolve_inputs(input_spec: str) -> list[Path]:
    """
    Resolve --input to a list of VCF paths.

    Accepts:
      - A glob string:     "data/*.vcf.gz"
      - A @filelist path:  "@my_samples.txt"
      - A single path:     "sample.vcf.gz"
    """
    if input_spec.startswith("@"):
        list_file = Path(input_spec[1:])
        if not list_file.exists():
            raise click.BadParameter(f"File list {list_file} does not exist")
        paths = [
            Path(line.strip())
            for line in list_file.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]
    else:
        paths = [Path(p) for p in sorted(glob.glob(input_spec))]

    if not paths:
        raise click.BadParameter(
            f"No files matched: {input_spec!r}. "
            "Use quotes around glob patterns to prevent shell expansion."
        )

    missing = [p for p in paths if not p.exists()]
    if missing:
        raise click.BadParameter(
            f"{len(missing)} input file(s) do not exist, e.g.: {missing[0]}"
        )

    return paths


@click.group()
@click.version_option(version=__version__, prog_name="v-shuffler")
def main() -> None:
    """v-shuffler: anonymise VCFs via simulated genetic recombination."""


# ---------------------------------------------------------------------------
# shuffle subcommand
# ---------------------------------------------------------------------------

@main.command()
@click.option(
    "--input", "-i", "input_spec", required=True,
    help='Input VCF glob ("data/*.vcf.gz") or @filelist.txt.',
)
@click.option(
    "--output-dir", "-o", required=True, type=click.Path(),
    help="Directory for output VCF files (created if absent).",
)
@click.option(
    "--genetic-map", "-m", required=True, type=click.Path(exists=True),
    help="Genetic map file (SHAPEIT5 or HapMap format).",
)
@click.option(
    "--chromosome", "-c", required=True,
    help='Chromosome to process, e.g. "chr22" or "22".',
)
@click.option(
    "--n-samples", "-n", default=None, type=int,
    help="Number of synthetic individuals to produce. Defaults to number of inputs.",
)
@click.option("--seed", "-s", default=None, type=int, help="Random seed for reproducibility.")
@click.option("--chunk-size", default=50_000, show_default=True, type=int,
              help="Variants to process per memory chunk.")
@click.option("--threads", default=4, show_default=True, type=int,
              help="Number of threads (reserved for future parallelism).")
@click.option(
    "--output-mode", default="per_sample", show_default=True,
    type=click.Choice(["per_sample", "multi_sample"]),
    help="Write one VCF per synthetic individual, or a single multi-sample VCF.",
)
@click.option("--max-missing", default=0.05, show_default=True, type=float,
              help="Maximum fraction of missing calls per variant (variants above this are skipped).")
@click.option("--verbose", is_flag=True, help="Enable debug logging.")
def shuffle(
    input_spec: str,
    output_dir: str,
    genetic_map: str,
    chromosome: str,
    n_samples: int | None,
    seed: int | None,
    chunk_size: int,
    threads: int,
    output_mode: str,
    max_missing: float,
    verbose: bool,
) -> None:
    """
    Shuffle variants between VCFs via simulated genetic recombination.

    Each output individual is a mosaic of genotype segments drawn from
    different real donors at biologically realistic crossover positions.
    """
    _setup_logging(verbose)

    input_paths = _resolve_inputs(input_spec)
    logger.info("Found %d input VCF files", len(input_paths))

    out_dir = Path(output_dir)
    n_output = n_samples if n_samples is not None else len(input_paths)

    config = ShufflerConfig(
        input_vcfs=input_paths,
        output_dir=out_dir,
        genetic_map=Path(genetic_map),
        chromosome=chromosome,
        n_output_samples=n_output,
        seed=seed,
        chunk_size_variants=chunk_size,
        n_threads=threads,
        output_mode=output_mode,
        max_missing_rate=max_missing,
    )

    _run_shuffle(config)


def _run_shuffle(config: ShufflerConfig) -> None:
    """Main pipeline: generate plans → stream chunks → write output."""

    rng = np.random.default_rng(config.seed)

    logger.info("Loading genetic map for %s ...", config.chromosome)
    gmap = GeneticMap(config.genetic_map, config.chromosome)
    logger.info(
        "Map loaded: %.2f – %.2f cM (total %.2f cM)",
        gmap.start_cm, gmap.end_cm, gmap.total_length_cm,
    )

    n_pool_samples = len(config.input_vcfs)
    logger.info(
        "Generating segment plans for %d synthetic individuals from %d donors ...",
        config.n_output_samples, n_pool_samples,
    )

    segment_plans = generate_all_segment_plans(
        n_output_samples=config.n_output_samples,
        genetic_map=gmap,
        n_pool_samples=n_pool_samples,
        rng=rng,
        lambda_override=config.n_crossovers_lambda,
    )

    avg_segs = sum(len(p) for p in segment_plans) / len(segment_plans)
    logger.info("Average segments per synthetic individual: %.1f", avg_segs)

    sample_names = [f"synthetic_{i}" for i in range(config.n_output_samples)]

    writer = SyntheticVCFWriter(
        output_dir=config.output_dir,
        sample_names=sample_names,
        template_vcf_path=config.input_vcfs[0],
        output_mode=config.output_mode,
        seed=config.seed,
        chromosome=config.chromosome,
        version=__version__,
    )

    reader = PerSampleVCFReader(
        vcf_paths=config.input_vcfs,
        chromosome=config.chromosome,
        genetic_map=gmap,
        chunk_size=config.chunk_size_variants,
        max_missing_rate=config.max_missing_rate,
    )

    logger.info("Streaming variants and building synthetic genotypes ...")
    total_variants = 0

    for pool in tqdm(reader.iter_chunks(), desc="Chunks", unit="chunk"):
        synthetic = build_synthetic_genotypes(pool, segment_plans)
        writer.write_chunk(pool, synthetic)
        total_variants += pool.n_variants

    logger.info("Processed %d variants total", total_variants)
    logger.info("Finalising output (bgzip + tabix) ...")

    final_paths = writer.finalize()
    logger.info("Done. Output files:")
    for p in final_paths[:5]:
        logger.info("  %s", p)
    if len(final_paths) > 5:
        logger.info("  ... and %d more", len(final_paths) - 5)


# ---------------------------------------------------------------------------
# validate subcommand
# ---------------------------------------------------------------------------

@main.command()
@click.option("--input", "-i", "input_spec", required=True,
              help="Glob or @filelist pointing to synthetic output VCFs.")
@click.option("--reference-vcf", "-r", required=True, type=click.Path(exists=True),
              help="A merged multi-sample input VCF to compare allele frequencies against.")
@click.option("--chromosome", "-c", required=True,
              help="Chromosome to validate.")
@click.option("--verbose", is_flag=True)
def validate(
    input_spec: str,
    reference_vcf: str,
    chromosome: str,
    verbose: bool,
) -> None:
    """
    Sanity-check synthetic VCFs against the original input.

    Checks:
    1. Output VCFs are valid and parseable.
    2. Allele frequency correlation between input and output is high (r > 0.99).
    3. No output sample is identical to any input sample.
    """
    _setup_logging(verbose)

    from v_shuffler.validate import run_validate
    synth_paths = _resolve_inputs(input_spec)
    run_validate(synth_paths, Path(reference_vcf), chromosome)
