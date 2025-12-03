#!/opt/homebrew/bin/bash

############################################################################
## Project: SNUseq project
## Script: Sort regions by Ratio -> Compute Matrix -> Plot Heatmap
##        (Multi-sample & Multi-annotation support)
## Author: Umut Gerlevik
## Date:   July 26, 2025
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# ──────────────── Configuration ────────────────
threads=8

# Base Paths
prefix="/MellorLab/SNUseqProject"

snuseq_dir="${prefix}/4_Anna/101_input_files/log2_SNUseq_mergedStrands"
annot_dir_atac="${prefix}/4_Anna/1_fromAnna/2_ATACpeaks"
annot_dir_fantom="${prefix}/3_publicData_HEK293T/FANTOM5"

ratio_bw="${prefix}/4_Anna/101_input_files/log2_K27ac_to_K4me3_ratios/K27ac_0_concatenatedSum_log2_div_K4me3_0_concatenatedSum_log2_log2ratio.bw"

# Output Directories
out_dir="${prefix}/4_Anna/105_SNUseq_multiSample_multiAnnotation_heatmaps"
sorted_bed_base="${out_dir}/sorted_annotations_by_K27ac_K4me3_ratio"
mkdir -p "$out_dir" "$sorted_bed_base"

# ──────────────── Inputs ────────────────

# 1. BigWigs to Plot (SNUseq)
# Add rep1, rep2, etc. here if needed
snuseq_bw=(
  "$snuseq_dir/SNUseq_0_concatenatedSum_log2.bw"
)
snuseq_labels=( "concat" )

# 2. Annotation Groups
# We use a string of paths separated by space for the values
declare -A annotation_groups
annotation_groups["ATAC_FANTOM5"]="$annot_dir_atac/ATAC_peaks_filtered.bed $annot_dir_fantom/FANTOM5_filtered.bed"
annotation_groups["ATAC_strands"]="$annot_dir_atac/ATAC_peaks_filtered_fwd.bed $annot_dir_atac/ATAC_peaks_filtered_rev.bed $annot_dir_atac/ATAC_peaks_filtered_unstranded.bed"
annotation_groups["intergenic"]="$annot_dir_atac/ATAC_peaks_filtered_unstranded.bed $annot_dir_fantom/FANTOM5_filtered.bed"

# ──────────────── Functions ────────────────

# Sort a BED file based on the Ratio BigWig
sort_annotations_by_ratio() {
  local bedfile="$1"
  local label=$(basename "$bedfile" .bed)
  local npz="${sorted_bed_base}/${label}_ratio_summary.npz"
  local tab="${sorted_bed_base}/${label}_ratio_signal.tab"
  local sorted="${sorted_bed_base}/${label}_sortedBy_K27acK4me3_ratio.bed"

  # Only sort if not already done (optional optimization)
  if [[ ! -f "$sorted" ]]; then
    >&2 echo "   🔢 Sorting $label by K27ac/K4me3 ratio..."

    multiBigwigSummary BED-file \
      --BED "$bedfile" \
      -b "$ratio_bw" \
      -out "$npz" \
      --outRawCounts "$tab" \
      -p "$threads"

    # Use R to sort descenting; handle NAs as 0
    Rscript --vanilla - "$tab" "$sorted" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
df <- read.table(args[1], header=TRUE, sep="\t", comment.char="")
colnames(df)[1:3] <- c("chr", "start", "end")
df[is.na(df[,4]), 4] <- 0
df_sorted <- df[order(-df[,4]), c(1,2,3)]
write.table(df_sorted, file=args[2], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF
  else
    >&2 echo "   ⏩ Sorted BED exists for $label, skipping sort..."
  fi

  echo "$sorted"
}

# ──────────────── Main Loop ────────────────

for group in "${!annotation_groups[@]}"; do
  echo -e "\n📚 Processing group: $group"

  # 1. Prepare Sorted Regions for this group
  # Split string into array
  raw_beds=(${annotation_groups[$group]})
  sorted_regions=()
  region_labels=()

  for bed in "${raw_beds[@]}"; do
    # Call sorting function
    sorted_bed_path=$(sort_annotations_by_ratio "$bed")
    sorted_regions+=("$sorted_bed_path")
    region_labels+=("$(basename "$bed" .bed)")
  done

  # 2. Define Output Names
  matrix_file="${out_dir}/matrix_SNUseq_multi_${group}.gz"
  heatmap_pdf="${out_dir}/heatmap_SNUseq_multi_${group}.pdf"

  # 3. Compute Matrix
  if [[ ! -f "$matrix_file" ]]; then
    echo "   🧬 Running computeMatrix..."
    computeMatrix reference-point \
      -S "${snuseq_bw[@]}" \
      -R "${sorted_regions[@]}" \
      --samplesLabel "${snuseq_labels[@]}" \
      --referencePoint center \
      --beforeRegionStartLength 4000 \
      --afterRegionStartLength 4000 \
      --binSize 10 \
      --missingDataAsZero \
      --sortRegions no \
      --averageTypeBins mean \
      -o "$matrix_file" \
      -p "$threads"
  else
    echo "   ⏩ Matrix exists, skipping computation..."
  fi

  # 4. Determine Plotting Parameters (Z-score scaling)
  # Apply specific scaling for 'intergenic' or groups containing FANTOM data
  local z_args=""
  if [[ "$group" == "intergenic" ]] || [[ "$group" == *"FANTOM"* ]]; then
    z_args="--zMin 0 --zMax 0.001"
    echo "   🎨 Plotting heatmap (Scaled: 0 - 0.001)..."
  else
    echo "   🎨 Plotting heatmap (Auto-scaled)..."
  fi

  # 5. Plot Heatmap
  plotHeatmap \
    -m "$matrix_file" \
    -out "$heatmap_pdf" \
    --colorMap "coolwarm" \
    --missingDataColor white \
    --yAxisLabel "Signal (log2)" \
    --refPointLabel "center" \
    --sortRegions no \
    --regionsLabel "${region_labels[@]}" \
    --samplesLabel "${snuseq_labels[@]}" \
    --legendLocation none \
    --heatmapHeight 25 \
    --heatmapWidth 10 $z_args

  echo "   ✅ Done: $heatmap_pdf"
done

echo -e "\n🎉 All multi-sample SNU-seq heatmaps plotted."
