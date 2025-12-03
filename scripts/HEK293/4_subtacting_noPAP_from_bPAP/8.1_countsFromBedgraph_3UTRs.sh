#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Count 3' UTRs in a stranded way, for scaling between bPAP
##                and noPAP samples as they should have similar 3' UTR counts
## Date: Mar 11, 2025
## Author: Umut Gerlevik
############################################################################

# "sum" is used to count all bedgraph signal in the region

source ~/.zshrc
mamba activate deeptools_env

# Set directories
common_prefix="/MellorLab/SNUseqProject/0_commonFiles"
prefix="/MellorLab/SNUseqProject"

inputDir="$prefix/1_Umut/7.1_merged_bedgraph_strandSpecific"

annotationDir="$common_prefix/annotations"
sortedAnnoDir="$annotationDir/sorted"

strandedOutputDir="$prefix/1_Umut/8.1_countsRegions_fromBedgraph/individualResults/strandsSeparated"
mergedOutputDir="$prefix/1_Umut/8.1_countsRegions_fromBedgraph/individualResults"

mkdir -p "$sortedAnnoDir"
mkdir -p "$strandedOutputDir"

# Define annotation feature types
features=("filtered_3UTRs")

# Define conditions
conditions=("labelled_bPAP" "labelled_bPAP_rRNA" "labelled_noPAP" "labelled_noPAP_rRNA" \
            "total_bPAP" "total_bPAP_rRNA" "total_noPAP" "total_noPAP_rRNA" \
            "unlabelled_bPAP" "unlabelled_noPAP")

# Sort annotation files
for strand in "fwd" "rev"; do
    for feature in "${features[@]}"; do
        annotation_file="$annotationDir/${feature}_${strand}.bed"
        sorted_file="$sortedAnnoDir/${feature}_${strand}_sorted.bed"

        if [[ -f "$annotation_file" ]]; then
            sort -k1,1 -k2,2n "$annotation_file" > "$sorted_file"
            echo "[INFO] Sorted: $sorted_file"
        else
            echo "[WARNING] Missing annotation file: $annotation_file"
        fi
    done
done

echo "[INFO] Annotation files sorted and stored in $sortedAnnoDir."

# Loop through each bedGraph file
for bedgraph_file in "$inputDir"/*.bedgraph; do
    fileName=$(basename "$bedgraph_file" _merged.bedgraph)

    # Detect strand based on filename pattern
    if echo "$fileName" | grep -q "_fwd"; then
        strand="fwd"
    elif echo "$fileName" | grep -q "_rev"; then
        strand="rev"
    else
        echo "[WARNING] Skipping unrecognized file: $fileName"
        continue
    fi

    echo "[INFO] Processing $fileName (strand: $strand)"

    for feature in "${features[@]}"; do
        annotation_file="$sortedAnnoDir/${feature}_${strand}_sorted.bed"
        output_file="$strandedOutputDir/${fileName}_${feature}.bedgraph"

        if [[ -f "$annotation_file" ]]; then
            bedtools map -a "$annotation_file" -b "$bedgraph_file" -c 4 -o sum > "$output_file"
            echo "[INFO] Coverage calculated for $fileName -> $feature"
        else
            echo "[WARNING] Annotation file missing: $annotation_file"
        fi
    done
done

echo "[INFO] All stranded coverage calculations completed. Proceeding to merge forward and reverse strands."

# Merge forward and reverse strand results for each feature
for feature in "${features[@]}"; do
    for condition in "${conditions[@]}"; do
        fwd_file="$strandedOutputDir/${condition}_fwd_${feature}.bedgraph"
        rev_file="$strandedOutputDir/${condition}_rev_${feature}.bedgraph"
        merged_file="$mergedOutputDir/${condition}_merged_${feature}.bedgraph"

        if [[ -f "$fwd_file" && -f "$rev_file" ]]; then
            echo "[INFO] Merging $fwd_file and $rev_file into $merged_file"

            # Use awk to concatenate the files row-wise (`rbind` manner)
            awk 'FNR==1 && NR!=1 {next} {print}' "$fwd_file" "$rev_file" > "$merged_file"

            echo "[INFO] Merged file created: $merged_file"
        else
            echo "[WARNING] Missing files for merging: $fwd_file or $rev_file"
        fi
    done
done

echo "[INFO] All merging completed. Merged results saved to $mergedOutputDir."
