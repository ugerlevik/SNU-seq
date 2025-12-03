#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Calculate total size and signal before and after filtering
##                out the negative and zero values from bPAP minus noPAP samples
## Date: Mar 12, 2025
## Author: Umut Gerlevik
############################################################################

# Load environment
source ~/.zshrc
mamba activate deeptools_env

# Set directories and parameters
common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/MellorLab/SNUseqProject

inputDir="$prefix/1_Umut/8.2_spikeNorm_3UTRnormalised_bedgraph"
filteredDir="$prefix/1_Umut/9.1_bPAPminusnoPAP_bigwig/filtered_bedgraphs"

outputFile="$prefix/1_Umut/9.1_bPAPminusnoPAP_bigwig/size_and_signal_stats.txt"

# Function to calculate total size and signal
calculate_size_and_signal() {
  file=$1
  awk '{
    size += ($3 - $2 + 1);
    signal += $4 * ($3 - $2 + 1);
  } END {
    print size "\t" signal;
  }' "$file"
}

# Write header to the output file
echo -e "File\tTotal_Size_Before\tTotal_Signal_Before\tTotal_Size_After\tTotal_Signal_After" > "$outputFile"

# Loop through the original files and corresponding filtered files
for original_file in "$inputDir"/*.bedgraph; do
    
  filename=$(basename "$original_file" "_3UTRnormalised.bedgraph")
  filtered_file="$filteredDir/${filename}_3UTRnormalised_minus_noPAP.bedgraph"

  # Skip files containing "noPAP" or "unlabelled" in their names
  if [[ "$filename" != *"bPAP"* || "$filename" == *"unlabelled"* ]]; then
    continue
  fi
  
  # Calculate size and signal for both files
  before_stats=$(calculate_size_and_signal "$original_file")
  after_stats=$(calculate_size_and_signal "$filtered_file")

  # Output the results
  echo -e "${filename}\t${before_stats}\t${after_stats}" >> "$outputFile"
done

echo "Size and signal statistics calculated and saved to $outputFile."


