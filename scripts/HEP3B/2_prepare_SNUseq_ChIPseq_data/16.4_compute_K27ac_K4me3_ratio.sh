#!/bin/bash

############################################################################
## Project: SNUseq project
## Script: log2(K27ac / K4me3) ratio computation
## Author: Umut Gerlevik
## Date: July 9, 2025
############################################################################

source ~/.zshrc
mamba activate deeptools_env

threads=8

# Directories
base_dir="/MellorLab/SNUseqProject/4_Anna/101_input_files"
k27ac_dir="${base_dir}/log2_K27ac"
k4me3_dir="${base_dir}/log2_K4me3"
out_dir="${base_dir}/log2_K27ac_to_K4me3_ratios"

mkdir -p "$out_dir"

# Genome size
genome_file="/MellorLab/SNUseqProject/0_commonFiles/genome/2_STARgenome/chrNameLength.txt"

# Pairs to compare: Format "K27ac_label|K4me3_label"
declare -a pairs=(
  "K27ac_0_concatenatedSum_log2|K4me3_0_concatenatedSum_log2"
)

# Loop through pairs
for pair in "${pairs[@]}"; do
  k27=$(cut -d'|' -f1 <<< "$pair")
  k4=$(cut -d'|' -f2 <<< "$pair")

  k27_file="${k27ac_dir}/${k27}.bw"
  k4_file="${k4me3_dir}/${k4}.bw"

  out_name="${k27}_div_${k4}_log2ratio.bw"
  out_path="${out_dir}/${out_name}"

  echo "🔁 Computing log2 ratio: log2($k27) - log2($k4)"

  bigwigCompare \
    -b1 "$k27_file" \
    -b2 "$k4_file" \
    --operation subtract \
    --binSize 1 \
    --numberOfProcessors $threads \
    -o "$out_path"

  echo "✅ Done → $out_path"
done

echo "🎉 All log2(K27ac/K4me3) ratio BigWigs generated."
