#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Plot metagenes of SNU-seq (bPAP minus no bPAP negatives 
##                 removed), sf4sU-seq, PRO-seq, TT-seq, mNET-seq in HEK293 cells
## Date: Sep 27, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Directories
prefix="/MellorLab/SNUseqProject"
inputDir="$prefix/1_Umut/14.2_matrices_bPAPminusnoPAP_negRemoved_withTTseqPROseqNETseq"
outputDir="$prefix/outputs/14.1_metagenePlots_bPAPminusnoPAP_negRemoved_withTTseqPROseqNETseq"
mkdir -p "$outputDir"


################################################################################################
### Fig_labelledrRNAminusMinusNoPAPnegativesRemoved_sf4sUseq_TTseq_PROseq_NETseq_TSStoTES
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_rRNA_bPAPminusNoPAP\"" "\"GSM5452295_HEK_FII_5\"" "\"GSM5452295_HEK_FII_3\"" 
  "\"TTseq\"" "\"GSM4730176_TTseq\"" "\"GSM4730174_PROseq\"" "\"GSM7990390_mNETseq\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/Fig_labelledrRNAminusMinusNoPAPnegativesRemoved_sf4sUseq_TTseq_PROseq_NETseq_TSStoTES.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8

################################################################################################
### Fig_labelledMinusNoPAPnegativesRemoved_sf4sUseq_TTseq_PROseq_NETseq_TSStoTES
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAPminusNoPAP\"" "\"GSM5452295_HEK_FII_5\"" "\"GSM5452295_HEK_FII_3\"" 
  "\"TTseq\"" "\"GSM4730176_TTseq\"" "\"GSM4730174_PROseq\"" "\"GSM7990390_mNETseq\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/Fig_labelledMinusNoPAPnegativesRemoved_sf4sUseq_TTseq_PROseq_NETseq_TSStoTES.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8

