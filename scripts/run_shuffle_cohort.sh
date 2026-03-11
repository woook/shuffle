#!/bin/bash
# run_shuffle_cohort.sh — shuffle a cohort of per-sample VCFs across multiple chromosomes.
#
# Usage:
#   run_shuffle_cohort.sh [options]
#
# Options:
#   -i DIR      Directory of per-sample VCFs (required)
#   -o DIR      Output directory for final synthetic VCFs (required)
#   -m DIR      Directory of per-chromosome genetic maps, e.g. chr1.b38.gmap.gz (required)
#   -s FILE     Sex file: two columns (path sex), used to restrict chrX to female donors (optional)
#   -c CHROMS   Space-separated chromosome list, e.g. "1 2 3 22 X" (default: 1-22 X)
#   -d INT      Minimum FORMAT/DP to retain a variant (default: 20)
#   -f FIELDS   Comma-separated FORMAT fields to carry into synthetic output, e.g. "AF,DP,AD"
#   -p PATTERN  Glob pattern to select input files (default: "*.vcf.gz"); use e.g.
#               "*tnhaplotyper2*" to exclude other callers in a mixed directory
#   -n INT      Number of parallel chromosomes to process (default: nproc / 2, min 1)
#   -e PATH     Path to v-shuffler venv (default: /tmp/vshuffler-venv)
#   -k          Keep per-chromosome intermediate files after combining (default: delete)
#   -r          Resume: skip chromosomes that already have output
#
# Example:
#   bash run_shuffle_cohort.sh \
#       -i /home/wook/Downloads/v-s/h \
#       -o /home/wook/Downloads/v-s/h/s \
#       -m /tmp/shuffle_maps \
#       -s /tmp/vshuffle_sex_h.txt \
#       -d 100 \
#       -f "AF,DP,AD" \
#       -n 4

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
INPUT_DIR=""
OUTPUT_BASE=""
MAP_DIR=""
SEX_FILE_ORIG=""
CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
MIN_DP=20
FORMAT_FIELDS=""
INPUT_PATTERN="*.vcf.gz"
VENV="/tmp/vshuffler-venv"
KEEP_INTERMEDIATES=0
RESUME=0

# Default parallelism: half the CPUs (leaves headroom for bcftools I/O workers)
DEFAULT_PARALLEL=$(( $(nproc 2>/dev/null || echo 2) / 2 ))
[ "$DEFAULT_PARALLEL" -lt 1 ] && DEFAULT_PARALLEL=1
MAX_PARALLEL="$DEFAULT_PARALLEL"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while getopts "i:o:m:s:c:d:f:p:n:e:kr" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_BASE="$OPTARG" ;;
        m) MAP_DIR="$OPTARG" ;;
        s) SEX_FILE_ORIG="$OPTARG" ;;
        c) CHROMS="$OPTARG" ;;
        d) MIN_DP="$OPTARG" ;;
        f) FORMAT_FIELDS="$OPTARG" ;;
        p) INPUT_PATTERN="$OPTARG" ;;
        n) MAX_PARALLEL="$OPTARG" ;;
        e) VENV="$OPTARG" ;;
        k) KEEP_INTERMEDIATES=1 ;;
        r) RESUME=1 ;;
        *) echo "Unknown option: $opt"; exit 1 ;;
    esac
done

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_BASE" ] || [ -z "$MAP_DIR" ]; then
    echo "Error: -i, -o, and -m are required."
    echo "Run with no arguments to see usage."
    exit 1
fi

PER_CHROM_DIR="$OUTPUT_BASE/per_chrom"
LOG="$OUTPUT_BASE/shuffle.log"
WORK_BASE="$(mktemp -d -t shuffle_work.XXXXXX)"
FILTERED_DIR="$WORK_BASE/filtered"

# Clean up work directory on exit
trap 'rm -rf "$WORK_BASE"' EXIT

mkdir -p "$OUTPUT_BASE" "$PER_CHROM_DIR" "$WORK_BASE" "$FILTERED_DIR"

echo "$(date '+%H:%M:%S') ============================================" | tee -a "$LOG"
echo "$(date '+%H:%M:%S') shuffle cohort pipeline" | tee -a "$LOG"
echo "$(date '+%H:%M:%S')   input:         $INPUT_DIR" | tee -a "$LOG"
echo "$(date '+%H:%M:%S')   output:        $OUTPUT_BASE" | tee -a "$LOG"
echo "$(date '+%H:%M:%S')   maps:          $MAP_DIR" | tee -a "$LOG"
echo "$(date '+%H:%M:%S')   sex file:      ${SEX_FILE_ORIG:-none}" | tee -a "$LOG"
echo "$(date '+%H:%M:%S')   chromosomes:   $CHROMS" | tee -a "$LOG"
echo "$(date '+%H:%M:%S')   min DP:        $MIN_DP" | tee -a "$LOG"
echo "$(date '+%H:%M:%S')   format fields: ${FORMAT_FIELDS:-GT only}" | tee -a "$LOG"
echo "$(date '+%H:%M:%S')   file pattern:  $INPUT_PATTERN" | tee -a "$LOG"
echo "$(date '+%H:%M:%S')   parallel:      $MAX_PARALLEL" | tee -a "$LOG"
echo "$(date '+%H:%M:%S') ============================================" | tee -a "$LOG"

# ---------------------------------------------------------------------------
# Step 1: Apply DP filter and index
# ---------------------------------------------------------------------------
DONOR_LIST="$WORK_BASE/donor_list.txt"

echo "$(date '+%H:%M:%S') Applying DP>=$MIN_DP filter to input VCFs (pattern: $INPUT_PATTERN)..." | tee -a "$LOG"
for f in "$INPUT_DIR"/$INPUT_PATTERN; do
    base=$(basename "$f")
    out="$FILTERED_DIR/$base"
    if [ ! -f "$out" ]; then
        bcftools view -i "FORMAT/DP>=$MIN_DP" "$f" -Oz -o "$out" \
            && tabix -p vcf "$out" &
    fi
done
wait

python3 - <<PYEOF
import pathlib
files = sorted(pathlib.Path("$FILTERED_DIR").glob("*.vcf.gz"))
pathlib.Path("$DONOR_LIST").write_text("\n".join(str(f) for f in files) + "\n")
print(f"$(date '+%H:%M:%S') Filtered donor VCFs ready: {len(files)}")
PYEOF

# ---------------------------------------------------------------------------
# Per-chromosome function (runs in a background subshell)
# ---------------------------------------------------------------------------
process_chrom() {
    local CHROM="$1"
    local WORK_DIR="$WORK_BASE/${CHROM}"
    local SPLIT_DIR="$WORK_DIR/split"
    local CHROM_OUT="$PER_CHROM_DIR/${CHROM}"
    local MAP_FILE="$MAP_DIR/chr${CHROM}.b38.gmap.gz"
    local CHROM_LOG="$OUTPUT_BASE/chr${CHROM}.log"

    # Clear old log for this chromosome
    > "$CHROM_LOG"

    echo "$(date '+%H:%M:%S') [chr${CHROM}] starting" | tee -a "$LOG"
    mkdir -p "$SPLIT_DIR" "$CHROM_OUT"

    # --- Merge with reference-fill ---
    echo "$(date '+%H:%M:%S') [chr${CHROM}] merging" | tee -a "$CHROM_LOG"
    bcftools merge --missing-to-ref \
        --file-list "$DONOR_LIST" \
        -r "$CHROM" \
        -Oz -o "$WORK_DIR/merged.vcf.gz" >> "$CHROM_LOG" 2>&1
    tabix -p vcf "$WORK_DIR/merged.vcf.gz"

    # --- Normalise ---
    echo "$(date '+%H:%M:%S') [chr${CHROM}] normalising" | tee -a "$CHROM_LOG"
    bcftools norm -m -any --keep-sum AD \
        "$WORK_DIR/merged.vcf.gz" \
        -Oz -o "$WORK_DIR/merged_norm.vcf.gz" >> "$CHROM_LOG" 2>&1
    tabix -p vcf "$WORK_DIR/merged_norm.vcf.gz"
    rm "$WORK_DIR/merged.vcf.gz" "$WORK_DIR/merged.vcf.gz.tbi"

    # --- Split to per-sample ---
    echo "$(date '+%H:%M:%S') [chr${CHROM}] splitting" | tee -a "$CHROM_LOG"
    bcftools +split -Oz "$WORK_DIR/merged_norm.vcf.gz" -o "$SPLIT_DIR/" >> "$CHROM_LOG" 2>&1
    for f in "$SPLIT_DIR"/*.vcf.gz; do tabix -p vcf "$f" & done
    wait
    rm "$WORK_DIR/merged_norm.vcf.gz" "$WORK_DIR/merged_norm.vcf.gz.tbi"

    # Build split_list.txt and a chromosome-local sex file from the split paths.
    # Using Python avoids shell glob failures and ensures the sex file basenames
    # match the split VCF names (bcftools +split names files by SM tag, not by
    # the original filename).
    python3 - <<PYEOF >> "$CHROM_LOG" 2>&1
import pathlib, re
split_dir = pathlib.Path("$SPLIT_DIR")
files = sorted(split_dir.glob("*.vcf.gz"))

# split_list.txt
pathlib.Path("$WORK_DIR/split_list.txt").write_text(
    "\n".join(str(f) for f in files) + "\n"
)

# sex file: extract F/M from the sample name.
# Pattern: look for -F- or -M- between two fields that are numbers / IDs,
# e.g. -5877-F-92197814 or -9526-M-103698. Matches any panel code.
sex_lines = []
for f in files:
    m = re.search(r"-([FM])-\d{5,}", f.name)
    if m:
        sex_lines.append(f"{f}  {m.group(1)}")
pathlib.Path("$WORK_DIR/sex_file.txt").write_text(
    "\n".join(sex_lines) + "\n"
)
print(f"chr$CHROM: {len(files)} split VCFs, {len(sex_lines)} sex-mapped")
PYEOF

    # --- Shuffle ---
    echo "$(date '+%H:%M:%S') [chr${CHROM}] shuffling" | tee -a "$LOG" "$CHROM_LOG"
    SHUFFLE_ARGS=(
        --input "@$WORK_DIR/split_list.txt"
        --output-dir "$CHROM_OUT"
        --genetic-map "$MAP_FILE"
        --chromosome "$CHROM"
        --seed 42
    )
    # Only pass --sex-file if the file actually has entries; an empty sex file
    # would cause v-shuffler to report "no female donors found" on chrX.
    if [ -s "$WORK_DIR/sex_file.txt" ]; then
        SHUFFLE_ARGS+=(--sex-file "$WORK_DIR/sex_file.txt")
    fi
    [ -n "$FORMAT_FIELDS" ] && SHUFFLE_ARGS+=(--carry-format-fields "$FORMAT_FIELDS")
    "$VENV/bin/v-shuffler" shuffle "${SHUFFLE_ARGS[@]}" >> "$CHROM_LOG" 2>&1

    rm -rf "$WORK_DIR"
    echo "$(date '+%H:%M:%S') [chr${CHROM}] done" | tee -a "$LOG"
}

# ---------------------------------------------------------------------------
# Launch chromosomes with concurrency limit
# ---------------------------------------------------------------------------
declare -A CHROM_PIDS

for CHROM in $CHROMS; do
    # Resume: skip chromosomes that already have output
    if [ "$RESUME" -eq 1 ]; then
        n=$(find "$PER_CHROM_DIR/$CHROM" -name "synthetic_*.vcf.gz" 2>/dev/null | wc -l || echo 0)
        if [ "$n" -gt 0 ]; then
            echo "$(date '+%H:%M:%S') [chr${CHROM}] already done, skipping" | tee -a "$LOG"
            continue
        fi
    fi

    # Wait if at the concurrency limit
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_PARALLEL" ]; do
        sleep 2
    done

    process_chrom "$CHROM" &
    CHROM_PIDS["$CHROM"]=$!
done

# Wait for all remaining jobs and collect exit codes
FAILED_CHROMS=()
for CHROM in "${!CHROM_PIDS[@]}"; do
    if ! wait "${CHROM_PIDS[$CHROM]}"; then
        FAILED_CHROMS+=("$CHROM")
        echo "$(date '+%H:%M:%S') [chr${CHROM}] FAILED — see $OUTPUT_BASE/chr${CHROM}.log" | tee -a "$LOG"
    fi
done

if [ "${#FAILED_CHROMS[@]}" -gt 0 ]; then
    echo "$(date '+%H:%M:%S') Pipeline failed for: ${FAILED_CHROMS[*]}" | tee -a "$LOG"
    exit 1
fi

echo "$(date '+%H:%M:%S') All chromosomes complete. Combining..." | tee -a "$LOG"

# ---------------------------------------------------------------------------
# Combine per-chromosome VCFs into one file per synthetic individual
# ---------------------------------------------------------------------------
# Count synthetic outputs from the current chromosome set (not a stale directory).
N_SYNTH=0
for chrom in $CHROMS; do
    chrom_dir="$PER_CHROM_DIR/$chrom"
    if [[ -d "$chrom_dir" ]]; then
        count=$(find "$chrom_dir" -name "synthetic_*.vcf.gz" 2>/dev/null | wc -l)
        if [[ $count -gt 0 ]]; then
            N_SYNTH=$count
            break
        fi
    fi
done

if [[ $N_SYNTH -eq 0 ]]; then
    echo "$(date '+%H:%M:%S') ERROR: No synthetic outputs found for chromosomes: $CHROMS" | tee -a "$LOG"
    exit 1
fi

echo "$(date '+%H:%M:%S') Combining $N_SYNTH synthetics across $(echo $CHROMS | wc -w) chromosomes..." | tee -a "$LOG"

# Run combines in parallel too
combine_synthetic() {
    local i="$1"
    python3 - <<PYEOF
import pathlib, subprocess, sys

chroms = "$CHROMS".split()
per_chrom = pathlib.Path("$PER_CHROM_DIR")
files = []
for c in chroms:
    f = per_chrom / c / f"synthetic_$i.vcf.gz"
    if f.exists():
        files.append(str(f))

if not files:
    print(f"WARNING: no files found for synthetic_$i", file=sys.stderr)
    sys.exit(1)

out = pathlib.Path("$OUTPUT_BASE") / f"synthetic_$i.vcf.gz"
subprocess.run(["bcftools", "concat"] + files + ["-Oz", "-o", str(out)], check=True)
subprocess.run(["tabix", "-p", "vcf", str(out)], check=True)
PYEOF
}

for i in $(seq 0 $((N_SYNTH - 1))); do
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_PARALLEL" ]; do
        sleep 1
    done
    combine_synthetic "$i" &
done
wait

# ---------------------------------------------------------------------------
# Clean up and report
# ---------------------------------------------------------------------------
if [ "$KEEP_INTERMEDIATES" -eq 0 ]; then
    echo "$(date '+%H:%M:%S') Cleaning up per-chromosome intermediates..." | tee -a "$LOG"
    rm -rf "$PER_CHROM_DIR"
fi

N_OUT=$(find "$OUTPUT_BASE" -maxdepth 1 -name "synthetic_*.vcf.gz" | wc -l)
echo "$(date '+%H:%M:%S') ============================================" | tee -a "$LOG"
echo "$(date '+%H:%M:%S') Complete: $N_OUT synthetic VCFs in $OUTPUT_BASE" | tee -a "$LOG"
echo "$(date '+%H:%M:%S') ============================================" | tee -a "$LOG"
