#!/opt/homebrew/bin/bash

############################################################################
## Project: SNUseq Project
## Script: Compute and Plot log2 signal on BEDs sorted by the same signal
##        (Self-Validation)
## Author: Umut Gerlevik
## Date:   Oct 8, 2025
############################################################################

# Load environment
source ~/.zshrc
mamba activate deeptools_env

# ──────────────── Configuration ────────────────
threads=8
binSize=10

# Base Directories
base_dir="/MellorLab/SNUseqProject/4_Anna"

input_files="${base_dir}/101_input_files"
log2_k27ac_dir="${input_files}/log2_K27ac"
log2_k4me3_dir="${input_files}/log2_K4me3"
log2_ratio_dir="${input_files}/log2_K27ac_to_K4me3_ratios"
sorted_bed_root="${base_dir}/104_SNUseq_heatmaps_sortedPeaks"

# Output Directory
outdir="${base_dir}/106_validation_heatmaps_sortedBySelf"
mkdir -p "$outdir"

# ──────────────── Data Definitions ────────────────

# 1. Peak types (Regions)
peak_types=(
  "ATAC_peaks_filtered"
  "ATAC_peaks_filtered_fwd"
  "ATAC_peaks_filtered_rev"
  "ATAC_peaks_filtered_unstranded"
  "FANTOM5_filtered"
)

# 2. Signals to validate (filename prefix -> folder path)
declare -A signal_map=(
  ["K27ac_0_concatenatedSum_log2"]=$log2_k27ac_dir
  ["K4me3_0_concatenatedSum_log2"]=$log2_k4me3_dir
  ["K27ac_0_concatenatedSum_log2_div_K4me3_0_concatenatedSum_log2_log2ratio"]=$log2_ratio_dir
)

# ──────────────── Main Loop ────────────────

for signal in "${!signal_map[@]}"; do
  # Construct BigWig path
  bw="${signal_map[$signal]}/${signal}.bw"

  for peak in "${peak_types[@]}"; do
    # Define IO paths
    sorted_bed="${sorted_bed_root}/${peak}/sortedBy_${signal}/${peak}_sortedBy_${signal}.bed"
    outmat="${outdir}/matrix_${signal}_on_${peak}.gz"
    outpdf="${outdir}/heatmap_${signal}_on_${peak}.pdf"

    echo -e "\n🎯 Processing: $signal on $peak"

    # 1. Compute Matrix
    # (Check if exists to save time, remove 'if' block to force re-run)
    if [[ ! -f "$outmat" ]]; then
      echo "   🧬 Computing matrix..."
      computeMatrix reference-point \
        -S "$bw" \
        -R "$sorted_bed" \
        --referencePoint center \
        --beforeRegionStartLength 4000 \
        --afterRegionStartLength 4000 \
        --binSize $binSize \
        --missingDataAsZero \
        --averageTypeBins mean \
        --sortRegions keep \
        -p $threads \
        -o "$outmat"
    else
      echo "   ⏩ Matrix exists, skipping computation..."
    fi

    # 2. Determine Plotting Parameters
    # Logic: If FANTOM5, restrict Z-range. Else, auto-scale.
    z_args=""
    if [[ "$peak" == *"FANTOM5"* ]]; then
      z_args="--zMin 0 --zMax 0.001"
      echo "   🎨 Plotting (Scaled 0-0.001)..."
    else
      echo "   🎨 Plotting (Auto-scaled)..."
    fi

    # 3. Plot Heatmap
    plotHeatmap \
      -m "$outmat" \
      -out "$outpdf" \
      --colorMap "coolwarm" \
      --missingDataColor "white" \
      --yAxisLabel "log2 signal" \
      --refPointLabel "center" \
      --regionsLabel "${peak}" \
      --sortRegions keep \
      --legendLocation "none" \
      --heatmapHeight 16 \
      --heatmapWidth 6 $z_args

    echo "   ✅ Output: $outpdf"
  done
done

echo -e "\n🎉 All validation heatmaps completed."
