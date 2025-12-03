#!/opt/homebrew/bin/bash

############################################################################
## Project: SNUseq project
## Purpose: Merge replicates from the same type of libraries to enhance the
##         signal by concatenating via summing the signal across replicates.
## Author: Umut Gerlevik
## Date: July 8, 2025
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Genome size file
genome_file="/MellorLab/SNUseqProject/0_commonFiles/genome/2_STARgenome/chrNameLength.txt"

# Base output folder
out_base="/MellorLab/SNUseqProject/4_Anna/101_input_files/pooled_concatSum_bigwigs"
mkdir -p "$out_base"

# Temporary folder
tmp_dir="${out_base}/tmp_bedgraphs"
mkdir -p "$tmp_dir"

# Define assays and their input files
declare -A assays
assays["SNUseq"]="/MellorLab/SNUseqProject/4_Anna/2_GSE172053_RAW/SNUseq/merged_bigwigs/GSM5240706_0_SNU_rep1_merged.bw /MellorLab/SNUseqProject/4_Anna/2_GSE172053_RAW/SNUseq/merged_bigwigs/GSM5240707_0_SNU_rep2_merged.bw /MellorLab/SNUseqProject/4_Anna/2_GSE172053_RAW/SNUseq/merged_bigwigs/GSM5240708_0_SNU_rep3_merged.bw"
assays["K4me3"]="/MellorLab/SNUseqProject/4_Anna/2_GSE172053_RAW/K4me3/GSM5240747_0_K4me3_rep1.bigwig /MellorLab/SNUseqProject/4_Anna/2_GSE172053_RAW/K4me3/GSM5240748_0_K4me3_rep2.bigwig"
assays["K27ac"]="/MellorLab/SNUseqProject/4_Anna/2_GSE172053_RAW/K27ac/GSM5240739_0_K27ac_rep1.bigwig /MellorLab/SNUseqProject/4_Anna/2_GSE172053_RAW/K27ac/GSM5240740_0_K27ac_rep2.bigwig"

# Function to pool and write to assay-specific subfolder
pool_replicates_concat_sum() {
    local assay_name="$1"
    local bw_files_string="$2"

    IFS=' ' read -r -a bw_files <<< "$bw_files_string"

    echo "🔁 Pooling assay: $assay_name"
    
    # Extract common metadata pattern (e.g. "_0_")
    local pattern=""
    for f in "${bw_files[@]}"; do
        fname=$(basename "$f")
        [[ "$fname" =~ (_[0-9]+_) ]] && pattern="${BASH_REMATCH[1]}" && break
    done

    if [[ -z "$pattern" ]]; then
        echo "⚠️ No metadata pattern (_0_, etc.) found — defaulting to no pattern"
        pattern=""
    fi

    # Strip underscores if desired, or keep literal (your choice)
    clean_pattern="${pattern//_/}"  # turns "_0_" → "0"

    # Define output folder and file
    local assay_out_dir="${out_base}/${assay_name}"
    mkdir -p "$assay_out_dir"
    local output_bw="${assay_out_dir}/${assay_name}_${clean_pattern}_concatenatedSum.bw"

    # Intermediate files
    local merged_bedgraph="${tmp_dir}/${assay_name}_merged.bedgraph"
    local chr_filtered="${tmp_dir}/${assay_name}_merged_chr1-22.bedgraph"
    local sorted_bedgraph="${tmp_dir}/${assay_name}_merged_chr1-22_sorted.bedgraph"

    echo "  ➤ Merging with bigWigMerge..."
    bigWigMerge "${bw_files[@]}" "$merged_bedgraph"

    echo "  ➤ Filtering to chr1–22 only..."
    awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' "$merged_bedgraph" > "$chr_filtered"

    echo "  ➤ Sorting..."
    LC_COLLATE=C sort -k1,1 -k2,2n "$chr_filtered" > "$sorted_bedgraph"

    echo "  💾 Converting to BigWig: $output_bw"
    bedGraphToBigWig "$sorted_bedgraph" "$genome_file" "$output_bw"

    echo "✅ Done: $assay_name → $output_bw"
}

# Run for all assays
for assay in "${!assays[@]}"; do
    pool_replicates_concat_sum "$assay" "${assays[$assay]}"
done

echo "🗑️ Cleaning up temp files..."
rm -rf "$tmp_dir"

echo "🎉 All pooled replicate BigWigs written to: $out_base"
