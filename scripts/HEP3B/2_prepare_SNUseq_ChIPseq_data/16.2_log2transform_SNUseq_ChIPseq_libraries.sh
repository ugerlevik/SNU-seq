#!/opt/homebrew/bin/bash

############################################################################
## Project: SNUseq project
## Purpose: Apply log2(x + pseudocount) to SNU-seq and ChIP-seq BigWig files
## Author: Umut Gerlevik
## Date: July 8, 2025
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Constants
pseudocount=1
genome_file="/MellorLab/SNUseqProject/0_commonFiles/genome/2_STARgenome/chrNameLength.txt"

# ─────────────────────────────────────────────────────────────────────────────
# Function: Transform all BigWig files in a folder
# ─────────────────────────────────────────────────────────────────────────────
transform_bigwig_folder() {
  local dataset_name="$1"
  local input_folder="$2"
  local output_folder="$3"
  local filename_pattern="${4:-}"  # Optional 4th argument; empty = all files

  echo "📂 Scanning $input_folder for BigWig files for: $dataset_name"
  mkdir -p "$output_folder"
  local tmp_dir="${output_folder}/temp_bedgraphs"
  mkdir -p "$tmp_dir"

  if [[ -n "$filename_pattern" ]]; then
    echo "🔍 Filtering files with pattern: $filename_pattern"
    bw_files=($(find "$input_folder" -maxdepth 1 -type f \( -iname "*.bw" -o -iname "*.bigwig" \) -name "*$filename_pattern*"))
  else
    echo "🔍 Including ALL *.bw / *.bigwig files"
    bw_files=($(find "$input_folder" -maxdepth 1 -type f \( -iname "*.bw" -o -iname "*.bigwig" \)))
  fi

  if [[ ${#bw_files[@]} -eq 0 ]]; then
    echo "⚠️  No BigWig files found in $input_folder. Skipping."
    return
  fi

  echo "🔢 Found ${#bw_files[@]} files."

  for bw in "${bw_files[@]}"; do
    base=$(basename "$bw" .bw)
    base=$(basename "$base" .bigwig)

    bedgraph="${tmp_dir}/${base}.bedgraph"
    sorted_bedgraph="${tmp_dir}/${base}_sorted.bedgraph"
    log2_bedgraph="${tmp_dir}/${base}_log2.bedgraph"
    sorted_log2_bedgraph="${tmp_dir}/${base}_log2_sorted.bedgraph"
    output_bw="${output_folder}/${base}_log2.bw"

    echo "🔁 Processing: $base"

    bigWigToBedGraph "$bw" "$bedgraph"

    awk 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' "$bedgraph" \
      | LC_COLLATE=C sort -k1,1 -k2,2n > "$sorted_bedgraph"

    awk -v pc="$pseudocount" 'BEGIN {OFS="\t"} 
      NF == 4 {
        val = log($4 + pc) / log(2);
        print $1, $2, $3, val;
      }' "$sorted_bedgraph" > "$log2_bedgraph"

    LC_COLLATE=C sort -k1,1 -k2,2n "$log2_bedgraph" > "$sorted_log2_bedgraph"
    bedGraphToBigWig "$sorted_log2_bedgraph" "$genome_file" "$output_bw"

    echo "✅ Written: $output_bw"
  done

  rm -rf "$tmp_dir"
  echo "🎉 Completed: $dataset_name"
}

# ─────────────────────────────────────────────────────────────────────────────
# Run transformations for each dataset
# ─────────────────────────────────────────────────────────────────────────────
transform_bigwig_folder "SNU-seq merged strands" \
  "/MellorLab/SNUseqProject/4_Anna/2_GSE172053_RAW/SNUseq/merged_bigwigs" \
  "/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_SNUseq_mergedStrands" \
  "_0_"

transform_bigwig_folder "K27ac split strands" \
  "/MellorLab/SNUseqProject/4_Anna/2_GSE172053_RAW/K27ac" \
  "/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_K27ac" \
  "_0_"

transform_bigwig_folder "K4me3 split strands" \
  "/MellorLab/SNUseqProject/4_Anna/2_GSE172053_RAW/K4me3" \
  "/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_K4me3" \
  "_0_"

transform_bigwig_folder "SNU-seq concatenated sum" \
  "/MellorLab/SNUseqProject/4_Anna/101_input_files/pooled_concatSum_bigwigs/SNUseq" \
  "/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_SNUseq_mergedStrands" \
  "_0_"

transform_bigwig_folder "K27ac concatenated sum" \
  "/MellorLab/SNUseqProject/4_Anna/101_input_files/pooled_concatSum_bigwigs/K27ac" \
  "/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_K27ac" \
  "_0_"

transform_bigwig_folder "K4me3 concatenated sum" \
  "/MellorLab/SNUseqProject/4_Anna/101_input_files/pooled_concatSum_bigwigs/K4me3" \
  "/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_K4me3" \
  "_0_"

echo "🎉 All transformations completed successfully!"

