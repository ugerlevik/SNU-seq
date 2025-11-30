#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Compute matrix for PCA (raw)
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

# Load environment
source ~/.zshrc
mamba activate deeptools_env

threads=8

# Set common directories
prefix=/MellorLab/SNUseqProject

# Set directory paths
bwFolder=$prefix/1_Umut/4.4_raw_bigwigFiles_bothStrand_genomecov3end
matrixFolder=$prefix/1_Umut/4.6_PCA_matrices_raw

# Create directories if they don't exist
mkdir -p "$matrixFolder"

# All bw files in 8_allReadyForFurtherAnalysis
allbwFiles=$(ls "$bwFolder"/*labelled*.bw)

multiBigwigSummary bins --bwfiles $allbwFiles \
                        --binSize 10000 -p $threads \
                        --smartLabels \
                        --outFileName $matrixFolder/multiBigwigSummary_results_10kb_onlyLabelled.npz \
                        --outRawCounts $matrixFolder/multiBigwigSummary_results_10kb_counts_onlyLabelled.tsv

