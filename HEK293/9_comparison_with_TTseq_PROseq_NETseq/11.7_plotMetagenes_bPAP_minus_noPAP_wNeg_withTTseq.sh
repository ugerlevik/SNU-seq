#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Plot metagenes (merged) with TTseq
## Date: Apr 4, 2025
## Author: Umut Gerlevik
############################################################################

# Load environment
source ~/.zshrc
mamba activate deeptools_env

threads=8

# Set common directories
prefix=/MellorLab/SNUseqProject

# Set directory paths
inputDir=$prefix/1_Umut/11.1_matrices_bPAPminusnoPAP_wNegatives_withTTseq
outputDir=$prefix/outputs/11.1_metagenePlots_bPAPminusnoPAP_wNegatives_withTTseq

mkdir -p $outputDir


################################################################################################
### Labelled only (unlabelled included)
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAP\"" "\"unlabelled_bPAP\"" "\"labelled_noPAP\"" "\"labelled_bPAP_minus_noPAP_wNeg\"" "\"TTseq_A\"" "\"TTseq_B\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/labelled_TSStoTES_minus_noPAP_wNegatives_withTTseq.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8


################################################################################################
### Labelled only (unlabelled included) -rRNA
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAP_rRNA\"" "\"unlabelled_bPAP\"" "\"labelled_noPAP_rRNA\"" "\"labelled_bPAP_rRNA_minus_noPAP_wNeg\"" "\"TTseq_A\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/labelledminusrRNA_TSStoTES_minus_noPAP_wNegatives_withTTseq.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8

################################################################################################
### Labelled only (unlabelled included) -rRNA without TTseq
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAP_rRNA\"" "\"unlabelled_bPAP\"" "\"labelled_noPAP_rRNA\"" "\"labelled_bPAP_rRNA_minus_noPAP_wNeg\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/labelledminusrRNA_TSStoTES_minus_noPAP_wNegatives.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8


################################################################################################
### TTseq only
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"TTseq_A\"" "\"TTseq_B\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/only_TTseq.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8


################################################################################################
### bpapminusnoPAP and TTseq only
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAP_minus_noPAP_wNeg\"" "\"TTseq_A\"" 
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/bPAPminusNoPAP_TTseq.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8


################################################################################################
### bpapminusnoPAP -rRNA and TTseq only
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAP_rRNA_minus_noPAP_wNeg\"" "\"TTseq_A\"" 
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/bPAPminusNoPAPminusrRNA_TTseq.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8


echo "All heatmap plots have been generated."

