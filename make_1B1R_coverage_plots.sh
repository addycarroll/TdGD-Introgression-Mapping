#!/bin/bash
#SBATCH --job-name=Jag1R_plots
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=Jag1R_plots_%j.out
#SBATCH --error=Jag1R_plots_%j.err

set -euo pipefail

# =============================================================================
# For each sample and MAPQ threshold, this script:
#   1. divides Jagger chromosome 1B and Weining chromosome 1R into consecutive,
#      nonoverlapping 1-Mb bins;
#   2. uses samtools bedcov to sum aligned-base depth in each bin;
#   3. records the number of alignments overlapping each bin in the same pass;
#   4. divides the depth sum by the exact bin length to obtain mean read depth;
#   5. creates a separate raw-depth PNG dot plot for each sample, with Jagger
#      chromosome 1B and Weining chromosome 1R shown side by side; and
#   6. writes the underlying tab-delimited data used for every point.
#
# Submit:
#   sbatch make_parker_1B_1R_coverage_plots.sh
# =============================================================================

# ------------------------------ EDIT HERE -----------------------------------
ANALYSIS_ROOT="/bulk/akf/Addy/Variant_Calling/1RS_1BL_analysis"
BIN_SIZE=1000000
MAPQ_THRESHOLDS=(20)

# Maximum displayed y-axis value for each sample
Y_AXIS_MAXIMA=(
    "Larry=30"
    "KS061278M-4=30"
    "KS061406LN-26=20"
    "PI487264=5"
    "PI428088=5"
    "PI471021=5"
    "PI466984=20"
    "PI471815=15"
    "PI470736=15"
    "PI414718=15"
)
# ---------------------------------------------------------------------------

ALIGNMENT_ROOT="$ANALYSIS_ROOT/alignments"
OUTPUT_ROOT="$ANALYSIS_ROOT/coverage_1Mb"
DATA_DIR="$OUTPUT_ROOT/data"
PLOT_DIR="$OUTPUT_ROOT/plots"
BED_DIR="$OUTPUT_ROOT/bed"

JAGGER_1B_ID="LR862509.1"
JAGGER_1B_LENGTH=705338699
RYE_1R_ID="Weining_1R"
RYE_1R_LENGTH=940967580

SAMPLES=(
    "Larry"
    "KS061278M-4"
    "KS061406LN-26"
    "PI487264"
    "PI428088"
    "PI471021"
    "PI466984"
    "PI471815"
    "PI470736"
    "PI414718"
)
EXPECTED_TOTAL_BINS=$(
    awk -v a="$JAGGER_1B_LENGTH" -v b="$RYE_1R_LENGTH" -v w="$BIN_SIZE" \
        'BEGIN {print int((a+w-1)/w) + int((b+w-1)/w)}'
)

module purge
module load SAMtools

for program in samtools awk; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: required program is unavailable: $program" >&2
        exit 1
    fi
done

SAMTOOLS_VERSION="$(samtools --version | awk 'NR == 1 {print $2}')"

mkdir -p "$DATA_DIR" "$PLOT_DIR" "$BED_DIR"

TARGET_BED="$BED_DIR/Jagger_1B_Weining_1R.${BIN_SIZE}bp_bins.bed"
TARGET_BED_PARTIAL="${TARGET_BED}.partial.${SLURM_JOB_ID:-manual}"

make_bins() {
    local chromosome="$1"
    local length="$2"
    local start=0
    local end
    local bin_number=1

    while (( start < length )); do
        end=$((start + BIN_SIZE))
        if (( end > length )); then
            end="$length"
        fi
        printf '%s\t%d\t%d\t%s_bin_%04d\n' \
            "$chromosome" "$start" "$end" "$chromosome" "$bin_number"
        start="$end"
        bin_number=$((bin_number + 1))
    done
}

bed_is_valid() {
    [[ -s "$TARGET_BED" ]] || return 1
    [[ "$(wc -l < "$TARGET_BED")" -eq "$EXPECTED_TOTAL_BINS" ]] || return 1
    awk -v chr1b="$JAGGER_1B_ID" -v len1b="$JAGGER_1B_LENGTH" \
        -v chr1r="$RYE_1R_ID" -v len1r="$RYE_1R_LENGTH" \
        -v width="$BIN_SIZE" '
        BEGIN {ok=1; prev_chr=""; prev_end=0}
        NF != 4 {ok=0}
        $1 != chr1b && $1 != chr1r {ok=0}
        $2 < 0 || $3 <= $2 || ($3-$2) > width {ok=0}
        $1 == prev_chr && $2 != prev_end {ok=0}
        $1 != prev_chr && $2 != 0 {ok=0}
        {prev_chr=$1; prev_end=$3; last[$1]=$3; count[$1]++}
        END {
            if (last[chr1b] != len1b || last[chr1r] != len1r) ok=0
            if (count[chr1b] != int((len1b+width-1)/width)) ok=0
            if (count[chr1r] != int((len1r+width-1)/width)) ok=0
            exit !ok
        }' "$TARGET_BED"
}

if ! bed_is_valid; then
    echo "Creating nonoverlapping ${BIN_SIZE}-bp bins."
    {
        make_bins "$JAGGER_1B_ID" "$JAGGER_1B_LENGTH"
        make_bins "$RYE_1R_ID" "$RYE_1R_LENGTH"
    } > "$TARGET_BED_PARTIAL"
    mv "$TARGET_BED_PARTIAL" "$TARGET_BED"

    if ! bed_is_valid; then
        echo "ERROR: generated BED file failed validation: $TARGET_BED" >&2
        exit 1
    fi
else
    echo "Using validated existing BED file: $TARGET_BED"
fi

table_is_valid() {
    local table="$1"
    local sample="$2"
    local mapq="$3"

    [[ -s "$table" ]] || return 1
    [[ "$(awk 'NR > 1 {n++} END {print n+0}' "$table")" \
        -eq "$EXPECTED_TOTAL_BINS" ]] || return 1
    awk -F '\t' -v sample="$sample" -v mapq="$mapq" \
        -v chr1b="$JAGGER_1B_ID" -v chr1r="$RYE_1R_ID" '
        NR == 1 {
            if ($0 != "Sample\tMAPQ_min\tChromosome\tStart_bp\tEnd_bp\tMidpoint_Mb\tBin_length_bp\tDepth_sum\tRead_count\tMean_depth") exit 1
            next
        }
        $1 != sample || $2 != mapq {exit 1}
        $3 != chr1b && $3 != chr1r {exit 1}
        $4 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ || $7 !~ /^[0-9]+$/ || $8 !~ /^[0-9]+$/ || $9 !~ /^[0-9]+$/ {exit 1}
        $10 < 0 {exit 1}
        END {if (NR <= 1) exit 1}
        ' "$table"
}

echo "Coverage calculation started: $(date)"
echo "Samtools version: $SAMTOOLS_VERSION"
echo "Bin size: $BIN_SIZE bp"
echo "MAPQ thresholds: ${MAPQ_THRESHOLDS[*]}"
echo

for sample in "${SAMPLES[@]}"; do
    sample_dir="$ALIGNMENT_ROOT/$sample"
    bam="$sample_dir/${sample}.Jagger_plus_Weining_1R.mapped.sorted.bam"

    if [[ ! -s "$bam" ]]; then
        echo "ERROR: BAM is missing or empty: $bam" >&2
        exit 1
    fi
    samtools quickcheck -v "$bam"

    for chromosome in "$JAGGER_1B_ID" "$RYE_1R_ID"; do
        if ! samtools view -c "$bam" "$chromosome" >/dev/null 2>&1; then
            echo "ERROR: indexed access failed for $sample, $chromosome." >&2
            echo "Expected a readable CSI index at ${bam}.csi" >&2
            exit 1
        fi
    done

    for mapq in "${MAPQ_THRESHOLDS[@]}"; do
        table="$DATA_DIR/${sample}.MAPQ${mapq}.mean_depth_1Mb.tsv"
        partial="${table}.partial.${SLURM_JOB_ID:-manual}"

        if table_is_valid "$table" "$sample" "$mapq"; then
            echo "Reusing validated table: $table"
            continue
        fi

        echo "Calculating $sample at MAPQ >= $mapq: $(date)"
        {
            printf 'Sample\tMAPQ_min\tChromosome\tStart_bp\tEnd_bp\tMidpoint_Mb\tBin_length_bp\tDepth_sum\tRead_count\tMean_depth\n'
            samtools bedcov -j -c -Q "$mapq" "$TARGET_BED" "$bam" |
                awk -v OFS='\t' -v sample="$sample" -v mapq="$mapq" '
                    {
                        length_bp=$3-$2
                        midpoint_mb=(($2+$3)/2)/1000000
                        mean_depth=$5/length_bp
                        print sample, mapq, $1, $2, $3,
                              sprintf("%.6f", midpoint_mb), length_bp, $5, $6,
                              sprintf("%.8f", mean_depth)
                    }'
        } > "$partial"

        mv "$partial" "$table"
        if ! table_is_valid "$table" "$sample" "$mapq"; then
            echo "ERROR: coverage table failed validation: $table" >&2
            exit 1
        fi
    done
done

COMBINED_TSV="$DATA_DIR/all_samples.mean_depth_1Mb.tsv"
COMBINED_PARTIAL="${COMBINED_TSV}.partial.${SLURM_JOB_ID:-manual}"

{
    printf 'Sample\tMAPQ_min\tChromosome\tStart_bp\tEnd_bp\tMidpoint_Mb\tBin_length_bp\tDepth_sum\tRead_count\tMean_depth\n'
    for sample in "${SAMPLES[@]}"; do
        for mapq in "${MAPQ_THRESHOLDS[@]}"; do
            tail -n +2 "$DATA_DIR/${sample}.MAPQ${mapq}.mean_depth_1Mb.tsv"
        done
    done
} > "$COMBINED_PARTIAL"
mv "$COMBINED_PARTIAL" "$COMBINED_TSV"

module purge
module load R

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript is unavailable after loading the R module." >&2
    echo "Run 'module spider R' and load an available R version." >&2
    exit 1
fi

echo "Creating Parker-style plots: $(date)"

Rscript --vanilla - "$COMBINED_TSV" "$PLOT_DIR" "${Y_AXIS_MAXIMA[@]}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
plot_dir <- args[2]
y_axis_specs <- args[-c(1, 2)]

d <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)

required <- c("Sample", "MAPQ_min", "Chromosome", "End_bp",
              "Midpoint_Mb", "Mean_depth")
if (!all(required %in% names(d))) {
  stop("Combined table is missing one or more required columns.")
}

panel_sample_order <- c(
  "Larry", "KS061278M-4", "KS061406LN-26",
  "PI487264", "PI428088", "PI471021", "PI466984", "PI471815",
  "PI470736", "PI414718"
)

if (length(y_axis_specs) != length(panel_sample_order)) {
  stop("Y_AXIS_MAXIMA must contain exactly one entry for every sample.")
}

y_axis_parts <- strsplit(y_axis_specs, "=", fixed = TRUE)
if (any(lengths(y_axis_parts) != 2L)) {
  stop("Every Y_AXIS_MAXIMA entry must use the form sample=maximum.")
}

y_axis_names <- vapply(y_axis_parts, `[`, character(1), 1)
y_axis_values <- suppressWarnings(as.numeric(
  vapply(y_axis_parts, `[`, character(1), 2)
))

if (anyDuplicated(y_axis_names) ||
    !setequal(y_axis_names, panel_sample_order) ||
    any(!is.finite(y_axis_values)) || any(y_axis_values <= 0)) {
  stop(paste(
    "Y_AXIS_MAXIMA must name every sample exactly once and use",
    "positive numeric maxima."
  ))
}

y_axis_max <- setNames(y_axis_values, y_axis_names)
chromosome_order <- c("LR862509.1", "Weining_1R")
chromosome_labels <- c(
  "LR862509.1" = "Jagger chromosome 1B",
  "Weining_1R" = "Weining chromosome 1R"
)
sample_labels <- c(
  "Larry" = "Larry",
  "KS061278M-4" = "KS061278M-4",
  "KS061406LN-26" = "KS061406LN~26",
  "PI487264" = "PI487264",
  "PI428088" = "PI428088",
  "PI471021" = "PI471021",
  "PI466984" = "PI466984",
  "PI471815" = "PI471815",
  "PI470736" = "PI470736",
  "PI414718" = "PI414718"
)
sample_colors <- c(
  "Larry" = "#2C7FB8",
  "KS061278M-4" = "#D95F0E",
  "KS061406LN-26" = "#E6AB02",
  "PI487264" = "#7570B3",
  "PI428088" = "#1B9E77",
  "PI471021" = "#E7298A",
  "PI466984" = "#66A61E",
  "PI471815" = "#A6761D",
  "PI470736" = "#1F78B4",
  "PI414718" = "#666666"
)

d$Sample <- factor(d$Sample, levels = panel_sample_order)
d$Chromosome <- factor(d$Chromosome, levels = chromosome_order)
d <- d[order(d$MAPQ_min, d$Sample, d$Chromosome, d$Midpoint_Mb), ]

draw_grid <- function() {
  abline(h = axTicks(2), col = "grey90", lwd = 0.7)
}

draw_one_sample <- function(sample, q) {
  z <- d[d$MAPQ_min == q & d$Sample == sample, ]
  ylim <- c(0, unname(y_axis_max[sample]))

  if (nrow(z) == 0L) {
    stop(paste("No plotting data found for", sample, "at MAPQ", q))
  }

  par(mfrow = c(1, 2), mar = c(4.1, 4.4, 3.0, 1.0),
      oma = c(1.0, 1.2, 3.0, 0.5), las = 1)

  for (chromosome in chromosome_order) {
    p <- z[z$Chromosome == chromosome, ]
    if (nrow(p) == 0L) {
      stop(paste("No plotting data found for", sample, chromosome,
                 "at MAPQ", q))
    }
    plot(p$Midpoint_Mb, p$Mean_depth,
         type = "n", xlim = c(0, max(p$End_bp) / 1e6), ylim = ylim,
         xlab = "Chromosome position (Mb)",
         ylab = "Mean read depth per 1-Mb bin",
         main = chromosome_labels[chromosome], cex.main = 0.96)
    draw_grid()
    points(p$Midpoint_Mb, p$Mean_depth, pch = 16, cex = 0.55,
           col = sample_colors[sample])
    box()
  }

  mtext(paste0(sample_labels[sample],
               " — coverage profiles (MAPQ ≥ ", q,
               "; y-axis maximum = ", format(y_axis_max[sample]), ")"),
        side = 3, outer = TRUE, line = 1.0, font = 2)
}

write_png <- function(stem, sample, q) {
  png(file.path(plot_dir, paste0(stem, ".png")),
      width = 3300, height = 1650, res = 300, type = "cairo")
  draw_one_sample(sample, q)
  dev.off()
}

for (q in sort(unique(d$MAPQ_min))) {
  for (sample in panel_sample_order) {
    write_png(paste0("Parker_style_", sample, "_MAPQ", q), sample, q)
  }
}
RSCRIPT

for sample in "${SAMPLES[@]}"; do
    for mapq in "${MAPQ_THRESHOLDS[@]}"; do
        expected_plot="$PLOT_DIR/Parker_style_${sample}_MAPQ${mapq}.png"
        if [[ ! -s "$expected_plot" ]]; then
            echo "ERROR: expected plot is missing or empty: $expected_plot" >&2
            exit 1
        fi
    done
done

RUN_INFO="$OUTPUT_ROOT/coverage_plot_run_${SLURM_JOB_ID:-manual}.txt"
{
    echo "completed=$(date --iso-8601=seconds)"
    echo "analysis_root=$ANALYSIS_ROOT"
    echo "bin_size_bp=$BIN_SIZE"
    echo "mapq_thresholds=${MAPQ_THRESHOLDS[*]}"
    echo "samtools_version=$SAMTOOLS_VERSION"
    echo "R_version=$(Rscript --version 2>&1)"
    echo "statistic=depth_sum_divided_by_exact_bin_length"
    echo "bedcov_options=-j; -c; -Q threshold"
    echo "normalization=none"
    echo "plot_layout=one sample per PNG; Jagger 1B and Weining 1R side by side"
    printf 'y_axis_maximum=%s\n' "${Y_AXIS_MAXIMA[@]}"
    echo "combined_data=$COMBINED_TSV"
    echo "plot_directory=$PLOT_DIR"
} > "$RUN_INFO"

echo
echo "Coverage plotting completed successfully: $(date)"
echo "Combined plotting data: $COMBINED_TSV"
echo "Plot directory: $PLOT_DIR"
echo "Run information: $RUN_INFO"
