#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Apply DESeq2 size factors to the BedGraph files
## Date: Mar 5, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Set common directories
common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/MellorLab/SNUseqProject

# Directory paths
input_bedgraph_folder=$prefix/1_Umut/5.5_filtered_bedgraph_strandSpecific

output_bedgraph_folder=$prefix/1_Umut/6.2_spikeinNormalised_bedgraph
output_bigwig_folder=$prefix/1_Umut/6.3_spikeinNormalised_bigwig
output_bedgraph_bothStrand_folder=$prefix/1_Umut/6.4_spikeinNormalised_bedgraph_bothStrand
output_bigwig_bothStrand_folder=$prefix/1_Umut/6.5_spikeinNormalised_bigwig_bothStrand

# Create output directories if they don't exist
mkdir -p "$output_bedgraph_folder"
mkdir -p "$output_bigwig_folder"
mkdir -p "$output_bedgraph_bothStrand_folder"
mkdir -p "$output_bigwig_bothStrand_folder"

scaling_factors_file=$prefix/1_Umut/6.1_featureCounts_spikeins/DESeq2_scalingFactors_spikeins.txt
genome_file="$common_prefix/genome/2_STARgenome/chrNameLength.txt"

# Loop through each BedGraph file and apply scaling factors
for bedgraph_file in "$input_bedgraph_folder"/*.bedgraph; do
    bedgraph_basename=$(basename "$bedgraph_file" _sorted.bedgraph)
    
    # Extract the sample name (remove _fwd or _rev)
    sample_name=$(echo "$bedgraph_basename" | sed -E 's/_fwd$|_rev$//')

    # Get the scaling factor for this sample by grepping the scaling_factors_file
    scaling_factor=$(grep "^$sample_name\t" "$scaling_factors_file" | cut -f2)

    if [ -z "$scaling_factor" ]; then
        echo "Warning: No scaling factor found for $sample_name. Skipping file."
        continue
    fi

    # Output normalized BedGraph and BigWig files
    scaled_bedgraph_file="$output_bedgraph_folder/${bedgraph_basename}_spikeinNormalised.bedgraph"
    scaled_bigwig_file="$output_bigwig_folder/${bedgraph_basename}_spikeinNormalised.bw"

    echo "Processing $bedgraph_file with scaling factor $scaling_factor"

    # Apply the scaling factor to the BedGraph file
    awk -v scale="$scaling_factor" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4 / scale}' "$bedgraph_file" > "$scaled_bedgraph_file"

    # Convert scaled BedGraph to BigWig
    bedGraphToBigWig "$scaled_bedgraph_file" "$genome_file" "$scaled_bigwig_file"

    echo "Scaled BedGraph saved to $scaled_bedgraph_file"
    echo "BigWig saved to $scaled_bigwig_file"
done

# Step 2: Combine forward and reverse BedGraphs for bothStrand using bedtools unionbedg
for bedgraph_file_fwd in "$output_bedgraph_folder"/*_fwd_spikeinNormalised.bedgraph; do
    
    base_name=$(basename "$bedgraph_file_fwd" "_fwd_spikeinNormalised.bedgraph")
    bedgraph_file_rev="$output_bedgraph_folder/${base_name}_rev_spikeinNormalised.bedgraph"
    
    if [ ! -f "$bedgraph_file_rev" ]; then
        echo "Warning: No reverse strand file found for $base_name. Skipping."
        continue
    fi
    
    # Output for merged BedGraph and BigWig
    merged_bedgraph_file="$output_bedgraph_bothStrand_folder/${base_name}_bothStrand_spikeinNormalised.bedgraph"
    merged_bigwig_file="$output_bigwig_bothStrand_folder/${base_name}_bothStrand_spikeinNormalised.bw"

    echo "Merging forward and reverse BedGraph files for $base_name"

    # Step 2.1: Use bedtools unionbedg to merge forward and reverse strands
    bedtools unionbedg -i "$bedgraph_file_fwd" "$bedgraph_file_rev" | \
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4 + $5)}' > "$merged_bedgraph_file.tmp"

    # Step 2.2: Sort the merged BedGraph file
    LC_COLLATE=C sort -k1,1 -k2,2n "$merged_bedgraph_file.tmp" > "$merged_bedgraph_file"

    # Step 2.3: Convert to BigWig for bothStrand 1_Umut
    bedGraphToBigWig "$merged_bedgraph_file" "$genome_file" "$merged_bigwig_file"
    
    # Clean up intermediate files
    rm "$merged_bedgraph_file.tmp"

    echo "Processed $bedgraph_file_fwd: Merged BedGraph and BigWig saved to $merged_bedgraph_file and $merged_bigwig_file"
done

echo "Normalisation and conversion complete."

