#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Plot correlations (on raw data)
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
matrixFolder=$prefix/1_Umut/4.6_PCA_matrices_raw
outputDir=$prefix/outputs/4.2_PCAplots_raw

mkdir -p $outputDir

plotPCA -in $matrixFolder/multiBigwigSummary_results_10kb_onlyLabelled.npz \
  -o "$outputDir/raw_PCAplot_onlyLabelled.pdf" \
  --plotFileFormat pdf

