#!/opt/homebrew/bin/bash

############################################################################
## Project: SNUseq project
## Script: Compute and plot SNU-seq metagenes on ATAC peaks and FANTOM5 
##        enhancers sorted by the SNU-seq signal
## Author: Umut Gerlevik
## Date: Aug 27, 2025
############################################################################

source ~/.zshrc
mamba activate deeptools_env

threads=8

# ─────────────── Inputs ───────────────
snuseq_dir="/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_SNUseq_mergedStrands"
snuseq_files=(
  "$snuseq_dir/SNUseq_0_concatenatedSum_log2.bw"
)
snuseq_labels=("concat")

peak_base_dir="/MellorLab/SNUseqProject/4_Anna/1_fromAnna/2_ATACpeaks"
fantom_base_dir="/MellorLab/SNUseqProject/3_publicData_HEK293T/FANTOM5"

peak_sets=(
  "$peak_base_dir/ATAC_peaks_filtered.bed"
  "$peak_base_dir/ATAC_peaks_filtered_fwd.bed"
  "$peak_base_dir/ATAC_peaks_filtered_rev.bed"
  "$peak_base_dir/ATAC_peaks_filtered_unstranded.bed"
  "$fantom_base_dir/FANTOM5_filtered.bed"
)

# ─────────────── Outputs ───────────────
out_base="/MellorLab/SNUseqProject/4_Anna/109_SNUseq_heatmaps_sortedBy_SNUseq"
mkdir -p "$out_base"

# ─────────────── Main Loop: Compute Matrix ───────────────
for peak_bed in "${peak_sets[@]}"; do
  subset=$(basename "$peak_bed" .bed)
  outdir="${out_base}/${subset}_sortedBy_SNUseq"
  mkdir -p "$outdir"

  echo "🔄 Processing subset: $subset"

  for i in "${!snuseq_files[@]}"; do
    bw="${snuseq_files[$i]}"
    label="${snuseq_labels[$i]}"
    matrix="${outdir}/matrix_${label}.gz"
    
    # Check if matrix already exists to avoid re-computation (Optional check)
    # Remove if overwrite is always desired
    if [[ ! -f "$matrix" ]]; then
      computeMatrix reference-point \
        -S "$bw" \
        -R "$peak_bed" \
        --referencePoint center \
        --beforeRegionStartLength 4000 \
        --afterRegionStartLength 4000 \
        --binSize 10 \
        --sortRegions keep \
        --missingDataAsZero \
        --averageTypeBins mean \
        -o "$matrix" \
        -p "$threads"
      
      echo "✅ Computed: $matrix"
    else
      echo "⏩ Skipped (exists): $matrix"
    fi
  done
done

echo -e "\n🎉 All matrices computed."

# ─────────────── Plotting Loop ───────────────
for peak_bed in "${peak_sets[@]}"; do
  # Recalculate subset name to ensure paths match the previous loop
  subset=$(basename "$peak_bed" .bed)
  
  # Define input matrix and output PDF
  matrix="${out_base}/${subset}_sortedBy_SNUseq/matrix_concat.gz"
  outpdf="${out_base}/${subset}_sortedBy_SNUseq/heatmap_${subset}_concat_noLegend.pdf"

  echo -e "\n🎨 Replotting ${subset} (SNU concat, no legend)..."

  if [[ -f "$matrix" ]]; then
    # Set default Z-max parameters
    z_min=""
    z_max=""
    
    # Logic: Apply specific scaling only for FANTOM5
    if [[ "$subset" == "FANTOM5_filtered" ]]; then
      z_min="--zMin 0"
      z_max="--zMax 0.001"
    fi

    plotHeatmap \
      -m "$matrix" \
      -out "$outpdf" \
      --colorMap "coolwarm" \
      --missingDataColor white \
      --yAxisLabel "log2(SNU-seq concat)" \
      --refPointLabel "center" \
      --regionsLabel "${subset}" \
      --sortRegions descend \
      --legendLocation none \
      --heatmapHeight 18 \
      --heatmapWidth 7 $z_min $z_max

    echo "✅ Done: $outpdf"
  else
    echo "⚠️ Matrix not found: $matrix"
  fi
done

echo -e "\n🚀 Script execution finished."