#!/opt/homebrew/bin/bash

############################################################################
## Project: SNUseq project
## Script: Compute and plot SNU-seq metagenes on ATAC peaks and FANTOM5 
##        enhancers sorted by various annotations/ratios/signals
## Author: Umut Gerlevik
## Date: July 9, 2025
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# ──────────────── Configuration ────────────────
threads=8

# Base Paths
prefix="/MellorLab/SNUseqProject"
input_base="${prefix}/4_Anna/101_input_files"
peak_base_dir="${prefix}/4_Anna/1_fromAnna/2_ATACpeaks"
fantom_base_dir="${prefix}/3_publicData_HEK293T/FANTOM5"

# Output
base_outdir="${prefix}/4_Anna/104_SNUseq_heatmaps_sortedPeaks"
mkdir -p "$base_outdir"

# ──────────────── Inputs ────────────────
# Directories
log2_k27ac_dir="${input_base}/log2_K27ac"
log2_k4me3_dir="${input_base}/log2_K4me3"
log2_ratio_dir="${input_base}/log2_K27ac_to_K4me3_ratios"
snuseq_dir="${input_base}/log2_SNUseq_mergedStrands"

# SNUseq Files (The signal to plot)
snuseq_files=( "$snuseq_dir/SNUseq_0_concatenatedSum_log2.bw" )
snuseq_labels=( "concat" )

# Reference Signals (The signals to sort by)
k27ac_refs=( "K27ac_0_concatenatedSum_log2" )
k4me3_refs=( "K4me3_0_concatenatedSum_log2" )
ratio_refs=( "K27ac_0_concatenatedSum_log2_div_K4me3_0_concatenatedSum_log2_log2ratio" )

# Peak Sets (Regions)
peak_sets=(
  "$peak_base_dir/ATAC_peaks_filtered.bed"
  "$peak_base_dir/ATAC_peaks_filtered_fwd.bed"
  "$peak_base_dir/ATAC_peaks_filtered_rev.bed"
  "$peak_base_dir/ATAC_peaks_filtered_unstranded.bed"
  "$fantom_base_dir/FANTOM5_filtered.bed"
)

# ──────────────── Functions ────────────────

# 1. Sort BED file based on a BigWig signal
sort_bed_by_signal() {
  local bedfile="$1"
  local signal_bw="$2"
  local output_dir="$3"

  local base=$(basename "$bedfile" .bed)
  local signal_name=$(basename "$signal_bw" .bw)

  local npz="${output_dir}/${base}_${signal_name}.npz"
  local tab="${output_dir}/${base}_${signal_name}.tab"
  local sorted="${output_dir}/${base}_sortedBy_${signal_name}.bed"

  echo "   🧩 Sorting $base by signal..." >&2

  multiBigwigSummary BED-file \
    --BED "$bedfile" \
    -b "$signal_bw" \
    -out "$npz" \
    --outRawCounts "$tab" \
    -p $threads

  # Use R to sort the tab output and create a new sorted BED
  Rscript --vanilla - "$tab" "$sorted" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
df <- read.table(args[1], header=TRUE, sep="\t", comment.char="")
colnames(df)[1:3] <- c("chr", "start", "end")
# Handle NAs by treating them as 0
df[is.na(df[,4]), 4] <- 0
# Sort descending
df_sorted <- df[order(-df[,4]), c(1,2,3)]
write.table(df_sorted, file=args[2], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF

  echo "$sorted"
}

# 2. Main processing pipeline for a specific reference signal
process_reference_signal() {
  local peak_file="$1"
  local signal_path="$2"
  local signal_label="$3" # e.g., "sortedBy_K27ac..."
  
  local peak_subset=$(basename "$peak_file" .bed)
  local outdir="${base_outdir}/${peak_subset}/${signal_label}"
  mkdir -p "$outdir"

  echo "🔄 Processing: $peak_subset | $signal_label"

  # 1. Generate Sorted BED
  sorted_bed=$(sort_bed_by_signal "$peak_file" "$signal_path" "$outdir")

  # 2. Iterate SNUseq files to plot
  for i in "${!snuseq_files[@]}"; do
    local bw="${snuseq_files[$i]}"
    local label="${snuseq_labels[$i]}"
    local matrix="${outdir}/matrix_SNUseq_${label}.gz"
    local pdf="${outdir}/heatmap_SNUseq_${label}.pdf"

    # Compute Matrix
    if [[ ! -f "$matrix" ]]; then
      computeMatrix reference-point \
        -S "$bw" \
        -R "$sorted_bed" \
        --referencePoint center \
        --beforeRegionStartLength 4000 \
        --afterRegionStartLength 4000 \
        --binSize 10 \
        --sortRegions keep \
        --missingDataAsZero \
        --averageTypeBins mean \
        -o "$matrix" \
        -p $threads
    fi

    # Determine Plot Settings based on Peak Set
    local z_args=""
    
    # logic: If filename implies FANTOM5, apply specific scaling
    if [[ "$peak_subset" == *"FANTOM5"* ]]; then
      z_args="--zMin 0 --zMax 0.001"
      echo "   🎨 Plotting (FANTOM5 scaling active: 0 - 0.001)..."
    else
      echo "   🎨 Plotting (Auto-scaling)..."
    fi

    plotHeatmap \
      -m "$matrix" \
      -out "$pdf" \
      --colorMap "coolwarm" \
      --missingDataColor white \
      --yAxisLabel "log2(SNU-seq $label)" \
      --refPointLabel "center" \
      --regionsLabel "$peak_subset" \
      --sortRegions keep \
      --legendLocation none \
      --heatmapHeight 18 \
      --heatmapWidth 7 $z_args
      
  done
}

# ──────────────── Main Execution Loop ────────────────

for peak_file in "${peak_sets[@]}"; do
  
  # Group 1: K27ac
  for ref in "${k27ac_refs[@]}"; do
    process_reference_signal "$peak_file" "${log2_k27ac_dir}/${ref}.bw" "sortedBy_${ref}"
  done

  # Group 2: K4me3
  for ref in "${k4me3_refs[@]}"; do
    process_reference_signal "$peak_file" "${log2_k4me3_dir}/${ref}.bw" "sortedBy_${ref}"
  done

  # Group 3: Ratios
  for ref in "${ratio_refs[@]}"; do
    process_reference_signal "$peak_file" "${log2_ratio_dir}/${ref}.bw" "sortedBy_${ref}"
  done

done

echo -e "\n✅ All heatmaps completed successfully."
