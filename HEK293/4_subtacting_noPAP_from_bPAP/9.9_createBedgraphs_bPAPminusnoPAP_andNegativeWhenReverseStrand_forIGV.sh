#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Convert BigWig files to BedGraph with negative strand
##               signal values inverted (multiplied by -1) and remove lines
##               with 0 signal.
## Date: Apr 7, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Set common directories
prefix=/MellorLab/SNUseqProject/

# Directory paths
input_bigwig_folder=$prefix/1_Umut/9.1_bPAPminusnoPAP_bigwig
output_bedgraph_folder=$prefix/1_Umut/9.1_bPAPminusnoPAP_bigwig/reverseNegative
negative_bedgraph_folder=$prefix/1_Umut/9.1_bPAPminusnoPAP_bigwig/reverseNegative

# Create output directory if it doesn't exist
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
    if [[ "$bigwig_basename" == *"_rev"* ]]; then
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
