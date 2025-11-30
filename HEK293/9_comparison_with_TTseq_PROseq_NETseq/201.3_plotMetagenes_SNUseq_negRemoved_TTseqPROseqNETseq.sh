#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Plot metagenes of SNU-seq (neg removed),
##                 PRO-seq, TT-seq, mNET-seq in HEK cells
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

: '
################################################################################################
### Fig_labelled_minusNoPAP_negativesRemoved_TSStoTES
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAP\"" "\"unlabelled_bPAP\"" "\"labelled_noPAP\"" "\"labelled_bPAPminusNoPAP\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/Fig1B_labelled_minusNoPAP_negativesRemoved_TSStoTES.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8


################################################################################################
### Fig_alternative_labelled_rRNAminus_minusNoPAP_negativesRemoved_TSStoTES
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAP_rRNA\"" "\"unlabelled_bPAP\"" "\"labelled_noPAP_rRNA\"" "\"labelled_rRNA_bPAPminusNoPAP\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/Fig1B_alternative_labelled_rRNAminus_minusNoPAP_negativesRemoved_TSStoTES.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8



################################################################################################
### Fig_labelled_minusNoPAP_negativesRemoved_TSStoTES
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAP_rRNA\"" "\"labelled_noPAP_rRNA\"" "\"labelled_rRNA_bPAPminusNoPAP\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/FigS2B_labelled_minusNoPAP_negativesRemoved_TSStoTES.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8


################################################################################################
### Fig_labelled_minusNoPAP_negativesRemoved_TSStoTES
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAP\"" "\"labelled_noPAP\"" "\"labelled_bPAPminusNoPAP\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/FigS2B_labelled_minusNoPAP_negativesRemoved_TSStoTES.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8



################################################################################################
### Fig_alternative_labelled_rRNAminus_minusNoPAP_negativesRemoved_withTTseqPROseqNETseq_TSStoTES
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_rRNA_bPAPminusNoPAP\"" "\"TTseq\"" "\"GSM4730176_TTseq\"" "\"GSM4730174_PROseq\"" "\"GSM6738185_mNETseq\""
   "\"GSM4964193_NETseq\"" "\"GSM7990390_mNETseq\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/Fig2C_alternative_labelled_rRNAminus_minusNoPAP_negativesRemoved_withTTseqPROseqNETseq_TSStoTES.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8

################################################################################################
### Fig_alternative_labelled_minusNoPAP_negativesRemoved_withTTseqPROseqNETseq_TSStoTES
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAPminusNoPAP\"" "\"TTseq\"" "\"GSM4730176_TTseq\"" "\"GSM4730174_PROseq\"" "\"GSM6738185_mNETseq\""
  "\"GSM4964193_NETseq\"" "\"GSM7990390_mNETseq\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/Fig2C_alternative_labelled_minusNoPAP_negativesRemoved_withTTseqPROseqNETseq_TSStoTES.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8


################################################################################################
### only_SNUseq_TTseq
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_bPAPminusNoPAP\"" "\"TTseq\"" "\"GSM4730176_TTseq\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/only_SNUseq_TTseq.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8

################################################################################################
### only_SNUseq_rRNAminus_TTseq
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"labelled_rRNA_bPAPminusNoPAP\"" "\"TTseq\"" "\"GSM4730176_TTseq\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/only_SNUseq_rRNAminus_TTseq.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8
'


################################################################################################
### Fig_total_negativesRemoved_TSStoTES
################################################################################################
# List of samples to include and their desired order
sampleOrderLabelled=(
  "\"total_bPAP\"" "\"total_noPAP\"" "\"total_bPAPminusNoPAP\""
  "\"total_bPAP_rRNA\"" "\"total_noPAP_rRNA\"" "\"total_rRNA_bPAPminusNoPAP\""
)

# Create a temporary matrix with only the selected samples in the specified order
subsetMatrixLabelled="$outputDir/TSS_to_TES_subset_allLabelled.gz"

# Subset and reorder the matrix
computeMatrixOperations subset -m "$inputDir/matrix_fwd_rev_combined_TSS_to_TES_scaled.gz" \
    --samples "${sampleOrderLabelled[@]}" \
    -o "$subsetMatrixLabelled"

# Plot the heatmap for labelled samples
plotHeatmap -m "$subsetMatrixLabelled" \
  -out "$outputDir/Fig_total_negativesRemoved_TSStoTES.pdf" \
  --colorMap "Blues" --missingDataColor white \
  --yAxisLabel "Average Signal" \
  --legendLocation "none" \
  --heatmapHeight 20 \
  --heatmapWidth 8
