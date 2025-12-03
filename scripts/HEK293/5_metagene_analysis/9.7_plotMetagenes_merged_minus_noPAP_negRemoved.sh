#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Plot metagenes (with the bPAP minus noPAP, negatives removed)
## Date: Mar 12, 2025
## Author: Umut Gerlevik
############################################################################

# Load environment
source ~/.zshrc
mamba activate deeptools_env

threads=8

# Set common directories
prefix=/MellorLab/SNUseqProject

# Set directory paths
inputDir=$prefix/1_Umut/9.3_matrices_bPAPminusnoPAP_negRemoved
outputDir=$prefix/outputs/9.3_metagenePlots_bPAPminusnoPAP_negRemoved

mkdir -p $outputDir

################################################################################################
### All together
################################################################################################
plotHeatmap -m $inputDir/matrix_fwd_rev_combined_TSS_to_TES.gz \
  -out "$outputDir/allTogether_TSStoTES_minus_noPAP_negRemoved.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "upper-left" \
  --heatmapHeight 20 \
  --heatmapWidth 8 


################################################################################################
### Labelled only (unlabelled included)
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAP\"" "\"labelled_noPAP\"" "\"labelled_bPAP_minus_noPAP\""
  "\"labelled_bPAP_rRNA\"" "\"labelled_noPAP_rRNA\"" "\"labelled_bPAP_rRNA_minus_noPAP\""
  "\"unlabelled_bPAP\"" "\"unlabelled_noPAP\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/labelled_TSStoTES_minus_noPAP_negRemoved.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "upper-left" \
  --heatmapHeight 20 \
  --heatmapWidth 8


################################################################################################
### Total only
################################################################################################
sampleOrderTotal=(
  "\"total_bPAP\"" "\"total_noPAP\"" "\"total_bPAP_minus_noPAP\""
  "\"total_bPAP_rRNA\"" "\"total_noPAP_rRNA\"" "\"total_bPAP_rRNA_minus_noPAP\""
)

subsetMatrixTotal="$outputDir/TSS_to_TES_subset_total.gz"

computeMatrixOperations subset \
  -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES.gz" \
  --samples "${sampleOrderTotal[@]}" \
  -o "$subsetMatrixTotal"

plotHeatmap \
  -m "$subsetMatrixTotal" \
  -out "$outputDir/total_TSStoTES_minus_noPAP_negRemoved.pdf" \
  --colorMap "Blues" \
  --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "upper-left" \
  --heatmapHeight 20 \
  --heatmapWidth 8

echo "All heatmap plots have been generated."

