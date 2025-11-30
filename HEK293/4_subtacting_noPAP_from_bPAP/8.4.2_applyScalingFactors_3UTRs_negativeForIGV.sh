#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Create bedgraphs for IGV (Reverse strand negative)
## Date: Mar 19, 2025
## Author: Umut Gerlevik
############################################################################


# Set common directories
prefix=/MellorLab/SNUseqProject/

# Directory paths
input_dir=$prefix/1_Umut/8.2_spikeNorm_3UTRnormalised_bedgraph
output_dir=$prefix/1_Umut/8.2_spikeNorm_3UTRnormalised_bedgraph/_forIGV_reverseNegative

# Create output directories if they don't exist
mkdir -p "$output_dir"

echo "[INFO] Processing reverse strand files..."
# Process reverse strand files (negate 4th column)
for file in "$input_dir"/*_rev_3UTRnormalised.bedgraph; do
    if [[ -f "$file" ]]; then
        new_filename=$(basename "$file" | sed 's/_3UTRnormalised//')
        output_file="$output_dir/$new_filename"

        # Negate the 4th column
        awk '{OFS="\t"; $4 = -$4; print}' "$file" > "$output_file"

        echo "[INFO] Processed (reverse strand): $file -> $output_file"
    fi
done

echo "[INFO] Processing forward strand files..."
# Copy forward strand files (no change in values, only filename change)
for file in "$input_dir"/*_fwd_3UTRnormalised.bedgraph; do
    if [[ -f "$file" ]]; then
        new_filename=$(basename "$file" | sed 's/_3UTRnormalised//')
        output_file="$output_dir/$new_filename"

        # Copy file without modification
        cp "$file" "$output_file"

        echo "[INFO] Copied (forward strand): $file -> $output_file"
    fi
done

echo "[INFO] All files processed and saved to $output_dir."

