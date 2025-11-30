#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Remove negative and zero values from bPAP minus noPAP
## Date: Mar 12, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/MellorLab/SNUseqProject

# Directories
inputDir="$prefix/1_Umut/9.1_bPAPminusnoPAP_bigwig/negativesNotRemoved"
outputDir="$prefix/1_Umut/9.1_bPAPminusnoPAP_bigwig"
filteredBedgraphDir="$outputDir/filtered_bedgraphs"

# Temporary directory for intermediate BedGraph files
tempDir="$outputDir/temp_bedgraphs"
mkdir -p "$tempDir"
mkdir -p "$filteredBedgraphDir" 

# Function to filter out negative and zero values from BedGraph and convert back to BigWig
filter_bedgraph_and_convert() {
  bigwig_file=$1
  bedgraph_file="$tempDir/$(basename "$bigwig_file" .bw).bedgraph"
  filtered_bedgraph_file="$filteredBedgraphDir/$(basename "$bigwig_file" _wNeg.bw).bedgraph"
  
  # Remove "_wNeg" from the output filename
  output_file="$outputDir/$(basename "$bigwig_file" _wNeg.bw).bw"

  # Convert BigWig to BedGraph
  bigWigToBedGraph "$bigwig_file" "$bedgraph_file"

  # Filter out negative and zero values
  awk '$4 > 0' "$bedgraph_file" > "$filtered_bedgraph_file"

  # Convert filtered BedGraph back to BigWig
  bedGraphToBigWig "$filtered_bedgraph_file" "$common_prefix/genome/2_STARgenome/chrNameLength.txt" "$output_file"

  # Clean up temporary BedGraph file
  rm "$bedgraph_file"
}

# Loop through the BigWig files and filter them
for bigwig_file in "$inputDir"/*.bw; do
  filter_bedgraph_and_convert "$bigwig_file"
  echo "Filtered $bigwig_file"
done

# Remove the temporary directory
rmdir "$tempDir"

echo "All endogenous-subtracted BigWig files have been filtered and saved."

