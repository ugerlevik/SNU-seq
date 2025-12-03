#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Subtract noPAP samples from bPAP samples 
##                (at single nucleotide resolution)
## Date: Mar 11, 2025
## Author: Umut Gerlevik
############################################################################


# Load environment
source ~/.zshrc
mamba activate deeptools_env

threads=8

prefix=/MellorLab/SNUseqProject

# Directories
inputDir="$prefix/1_Umut/8.3_spikeNorm_3UTRnormalised_bigwig"
outputDir="$prefix/1_Umut/9.1_bPAPminusnoPAP_bigwig/negativesNotRemoved"

# Create output directory if it doesn't exist
mkdir -p "$outputDir"

# List of corresponding files with strand information for noPAP and bPAP
file_pairs=(
  "labelled_bPAP_fwd_3UTRnormalised labelled_noPAP_fwd_3UTRnormalised"
  "labelled_bPAP_rev_3UTRnormalised labelled_noPAP_rev_3UTRnormalised"
  "labelled_bPAP_rRNA_fwd_3UTRnormalised labelled_noPAP_rRNA_fwd_3UTRnormalised"
  "labelled_bPAP_rRNA_rev_3UTRnormalised labelled_noPAP_rRNA_rev_3UTRnormalised"
  "total_bPAP_fwd_3UTRnormalised total_noPAP_fwd_3UTRnormalised"
  "total_bPAP_rev_3UTRnormalised total_noPAP_rev_3UTRnormalised"
  "total_bPAP_rRNA_fwd_3UTRnormalised total_noPAP_rRNA_fwd_3UTRnormalised"
  "total_bPAP_rRNA_rev_3UTRnormalised total_noPAP_rRNA_rev_3UTRnormalised"
  "unlabelled_bPAP_fwd_3UTRnormalised unlabelled_noPAP_fwd_3UTRnormalised"
  "unlabelled_bPAP_rev_3UTRnormalised unlabelled_noPAP_rev_3UTRnormalised"
)

# Apply bigwigCompare subtraction
for pair in "${file_pairs[@]}"; do
  # Read the filenames and library type from the pair
  read -r bPAP_file noPAP_file <<< "$pair"
  
  # Define file paths
  bPAP_file_path="$inputDir/${bPAP_file}.bw"
  noPAP_file_path="$inputDir/${noPAP_file}.bw"
  
  # Define output file paths
  output_file="$outputDir/${bPAP_file}_minus_noPAP_wNeg.bw"
  
  # Subtract noPAP from bPAP
  bigwigCompare -b1 "$bPAP_file_path" -b2 "$noPAP_file_path" --operation subtract --binSize 1 -p $threads -o "$output_file"
  
  echo "Subtraction of noPAP from bPAP completed for $bPAP_file"
done

echo "All noPAP subtractions completed. "

