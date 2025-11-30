#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Compute matrices for metagene profiles (bPAP minus noPAP, removed negatives)
## Date: Mar 12, 2025
## Author: Umut Gerlevik
############################################################################

# Load environment
source ~/.zshrc
mamba activate deeptools_env

threads=8

# Set common directories
common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/Documents/MellorLab/SNUseqProject

# Set directory paths
bwFolder1=$prefix/1_Umut/8.3_spikeNorm_3UTRnormalised_bigwig
bwFolder2=$prefix/1_Umut/9.1_bPAPminusnoPAP_bigwig

matrixFolder=$prefix/1_Umut/9.3_matrices_bPAPminusnoPAP_negRemoved

# Annotation files for forward and reverse strands
annotationFileFwd="$common_prefix/annotations/finalAnnotationSubset_fwd.gtf"
annotationFileRev="$common_prefix/annotations/finalAnnotationSubset_rev.gtf"

# Create directories if they don't exist
mkdir -p "$matrixFolder"

# Parameters for scale-regions
regionBodyLengthScaleRegions=5000
beforeRegionStartLengthScaleRegions=2000
afterRegionStartLengthScaleRegions=2000
binSizeScaleRegions=100

# All bw files in 8_allReadyForFurtherAnalysis
allbwFiles=("$bwFolder1"/*_3UTRnormalised.bw "$bwFolder2"/*_3UTRnormalised_minus_noPAP.bw)

# Separate the forward and reverse strand files
fwdbwFiles=()
revbwFiles=()
for file in "${allbwFiles[@]}"; do
  if [[ "$file" == *"fwd_3UTRnormalised.bw" || "$file" == *"fwd_3UTRnormalised_minus_noPAP.bw" ]]; then
    fwdbwFiles+=("$file")
  elif [[ "$file" == *"rev_3UTRnormalised.bw" || "$file" == *"rev_3UTRnormalised_minus_noPAP.bw" ]]; then
    revbwFiles+=("$file")
  fi
done

# Automatically generate samplesLabel for forward strand
samplesLabelFwd=()
for file in "${fwdbwFiles[@]}"; do
  if [[ "$file" == *"minus_noPAP.bw" ]]; then
    sample=$(basename "$file" | sed -E 's/_fwd_3UTRnormalised_minus_noPAP.bw/_minus_noPAP/')
  else
    sample=$(basename "$file" | sed -E 's/_fwd_3UTRnormalised.bw//')
  fi
  if [[ ! " ${samplesLabelFwd[*]} " =~ " ${sample} " ]]; then
    samplesLabelFwd+=("$sample")
  fi
done
samplesLabelFwdString=$(printf "\"%s\" " "${samplesLabelFwd[@]}")

# Automatically generate samplesLabel for reverse strand
samplesLabelRev=()
for file in "${revbwFiles[@]}"; do
  if [[ "$file" == *"minus_noPAP.bw" ]]; then
    sample=$(basename "$file" | sed -E 's/_rev_3UTRnormalised_minus_noPAP.bw/_minus_noPAP/')
  else
    sample=$(basename "$file" | sed -E 's/_rev_3UTRnormalised.bw//')
  fi
  if [[ ! " ${samplesLabelRev[*]} " =~ " ${sample} " ]]; then
    samplesLabelRev+=("$sample")
  fi
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

