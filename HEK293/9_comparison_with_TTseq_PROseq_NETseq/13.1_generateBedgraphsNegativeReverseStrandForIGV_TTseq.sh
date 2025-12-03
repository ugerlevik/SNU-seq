#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Create bedgraphs for IGV (reverse strand negative)
##                and remove lines with 0 signal (if any) (GSM5452296 and GSM5452297)
## Date: Apr 3, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Set common directories
prefix=/MellorLab/SNUseqProject/

# Directory paths
input_bigwig_folder=$prefix/2_TTseq_Phil/1.1_bigwigs
output_bedgraph_folder=$prefix/2_TTseq_Phil/1.2_bedgraphs
negative_bedgraph_folder=$prefix/2_TTseq_Phil/1.2_bedgraphs/reverseNegative

# Create output directory if it doesn't exist
mkdir -p "$output_bedgraph_folder"
mkdir -p "$negative_bedgraph_folder"

# Loop through each BigWig file in the input folder
for bigwig_file in "$input_bigwig_folder"/*.bw; do
    
    bigwig_basename=$(basename "$bigwig_file" .bw)

    output_bedgraph_file="$output_bedgraph_folder/${bigwig_basename}.bedgraph"

    bigWigToBedGraph "$bigwig_file" "$output_bedgraph_file"

done

for bigwig_file in "$input_bigwig_folder"/*.bw; do
    bigwig_basename=$(basename "$bigwig_file" .bw)
    
    # Check if the file is negative strand based on its name (_Minus)
    if [[ "$bigwig_basename" == *"_Minus"* ]]; then
        echo "Processing negative strand BigWig: $bigwig_file"
        output_bedgraph_file="$negative_bedgraph_folder/${bigwig_basename}.bedgraph"
        
        # Convert BigWig to BedGraph, invert the score, and remove 0 value lines
        bigWigToBedGraph "$bigwig_file" /dev/stdout | \
        awk 'BEGIN {OFS="\t"} {if ($4 != 0) print $1, $2, $3, -$4}' > "$output_bedgraph_file"
        echo "Negative strand BedGraph saved to: $output_bedgraph_file"
    else
        echo "Processing positive strand BigWig: $bigwig_file"
        output_bedgraph_file="$negative_bedgraph_folder/${bigwig_basename}.bedgraph"
        
        # Convert BigWig to BedGraph and remove 0 value lines
        bigWigToBedGraph "$bigwig_file" /dev/stdout | \
        awk 'BEGIN {OFS="\t"} {if ($4 != 0) print $1, $2, $3, $4}' > "$output_bedgraph_file"
        echo "Positive strand BedGraph saved to: $output_bedgraph_file"
    fi
done

echo "BigWig to BedGraph conversion complete."
