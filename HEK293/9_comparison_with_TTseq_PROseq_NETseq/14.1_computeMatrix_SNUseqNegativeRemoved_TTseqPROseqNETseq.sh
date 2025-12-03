#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Compute matrices for metagenes of SNU-seq (neg removed),
##                 PRO-seq, TT-seq, mNET-seq in HEK293 cells
## Date: Sep 27, 2025
## Author: Umut Gerlevik
############################################################################

# Load environment
source ~/.zshrc
mamba activate deeptools_env

threads=8

# Set common directories
common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/MellorLab/SNUseqProject

# Set directory paths
# SNU from 1_Umut/8.3_spikeNorm_3UTRnormalised_bigwig and 9.1_bPAPminusnoPAP_bigwig
# TT-seq from 2_TTseq_Phil/1.5_genebodyNormalised_bigwig (GSM5452296)
## others from 3_publicData_HEK293T/TTseq_NETseq_PROseq: GSM4730174, GSM4730176 and GSM7990390
bwFolder1=$prefix/1_Umut/14.1_comparisonsWithOtherMetods

matrixFolder=$prefix/1_Umut/14.2_matrices_bPAPminusnoPAP_negRemoved_withTTseqPROseqNETseq

# Annotation files for forward and reverse strands
annotationFileFwd="$common_prefix/annotations/finalAnnotationSubset_fwd.gtf"
annotationFileRev="$common_prefix/annotations/finalAnnotationSubset_rev.gtf"

# Create directories if they don't exist
mkdir -p "$matrixFolder"

# Parameters for scale-regions
regionBodyLengthScaleRegions=5000
beforeRegionStartLengthScaleRegions=2000
afterRegionStartLengthScaleRegions=2000
binSizeScaleRegions=20

# All bw files
allbwFiles=("$bwFolder1"/*.bw)

# Separate the forward and reverse strand files
fwdbwFiles=()
revbwFiles=()
for file in "${allbwFiles[@]}"; do
  if [[ "$file" == *"_fwd.bw" ]]; then
    fwdbwFiles+=("$file")
  elif [[ "$file" == *"_rev.bw" ]]; then
    revbwFiles+=("$file")
  fi
done

# Forward strand labels
samplesLabelFwd=()
for file in "${fwdbwFiles[@]}"; do
    sample=$(basename "$file" | sed -E 's/_fwd.bw//')
    samplesLabelFwd+=("$sample")
done
samplesLabelFwdString=$(printf "\"%s\" " "${samplesLabelFwd[@]}")

# Reverse strand labels
samplesLabelRev=()
for file in "${revbwFiles[@]}"; do
    sample=$(basename "$file" | sed -E 's/_rev.bw//')
    samplesLabelRev+=("$sample")
done
samplesLabelRevString=$(printf "\"%s\" " "${samplesLabelRev[@]}")

# Function to compute matrices for scale-regions
compute_scale_regions_matrix() {
  local matrixName=$1

  # ComputeMatrix for forward strand
  computeMatrix scale-regions \
    -R "$annotationFileFwd" \
    -S "${fwdbwFiles[@]}" \
    --regionBodyLength "$regionBodyLengthScaleRegions" \
    --beforeRegionStartLength "$beforeRegionStartLengthScaleRegions" \
    --afterRegionStartLength "$afterRegionStartLengthScaleRegions" \
    --binSize "$binSizeScaleRegions" \
    --missingDataAsZero --skipZeros \
    --averageTypeBins "mean" \
    --samplesLabel $samplesLabelFwdString \
    -o "$matrixFolder/matrix_fwd_$matrixName.gz" \
    --outFileNameMatrix "$matrixFolder/matrix_fwd_$matrixName.tab" \
    --outFileSortedRegions "$matrixFolder/matrix_fwd_$matrixName.bed" \
    --transcriptID "sequence_feature" \
    --transcript_id_designator "gene_id"\
    -p $threads

  echo "Matrix computation for the forward strand ($matrixName) completed."

  # ComputeMatrix for reverse strand
  computeMatrix scale-regions \
    -R "$annotationFileRev" \
    -S "${revbwFiles[@]}" \
    --regionBodyLength "$regionBodyLengthScaleRegions" \
    --beforeRegionStartLength "$beforeRegionStartLengthScaleRegions" \
    --afterRegionStartLength "$afterRegionStartLengthScaleRegions" \
    --binSize "$binSizeScaleRegions" \
    --missingDataAsZero --skipZeros \
    --averageTypeBins "mean" \
    --samplesLabel $samplesLabelRevString \
    -o "$matrixFolder/matrix_rev_$matrixName.gz" \
    --outFileNameMatrix "$matrixFolder/matrix_rev_$matrixName.tab" \
    --outFileSortedRegions "$matrixFolder/matrix_rev_$matrixName.bed" \
    --transcriptID "sequence_feature" \
    --transcript_id_designator "gene_id"\
    -p $threads

  echo "Matrix computation for the reverse strand ($matrixName) completed."

  # Combine forward and reverse matrices
  computeMatrixOperations rbind -m "$matrixFolder/matrix_fwd_$matrixName.gz" "$matrixFolder/matrix_rev_$matrixName.gz" -o "$matrixFolder/matrix_fwd_rev_combined_$matrixName.gz"

  echo "Forward and reverse matrices combined for $matrixName."
}

# Compute matrices for different scenarios
echo "Now TSS to TES"
compute_scale_regions_matrix "TSS_to_TES"

