#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Create bedgraphs for IGV (reverse strand negative)
##                and remove lines with 0 signal (if any)
## Date: Apr 4, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Set common directories
prefix=/MellorLab/SNUseqProject/

# Directory paths
input_bigwig_folder=$prefix/2_TTseq_Phil/1.4_genebodyNormalised_bedgraph
output_bedgraph_folder=$prefix/2_TTseq_Phil/1.4_genebodyNormalised_bedgraph/reverseNegative

# Create output directory if it doesn't exist
mkdir -p "$output_bedgraph_folder"

# Loop through each BigWig file in the input folder
for bigwig_file in "$input_bigwig_folder"/*_genebodyNormalised.bedgraph; do
    bigwig_basename=$(basename "$bigwig_file" _genebodyNormalised.bedgraph)
    
    # Check if the file is negative strand based on its name (_rev)
    if [[ "$bigwig_basename" == *"_rev"* ]]; then
        echo "Processing negative strand BigWig: $bigwig_file"
        output_bedgraph_file="$output_bedgraph_folder/${bigwig_basename}.bedgraph"
        
        # Convert BigWig to BedGraph, invert the score, and remove 0 value lines
        awk 'BEGIN {OFS="\t"} {print $1, $2, $3, -$4}' "$bigwig_file" > "$output_bedgraph_file"
        echo "Negative strand BedGraph saved to: $output_bedgraph_file"
    else
        echo "Processing positive strand BigWig: $bigwig_file"
        output_bedgraph_file="$output_bedgraph_folder/${bigwig_basename}.bedgraph"
        
        # Convert BigWig to BedGraph and remove 0 value lines
        awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}' "$bigwig_file" > "$output_bedgraph_file"
        echo "Positive strand BedGraph saved to: $output_bedgraph_file"
    fi
done
