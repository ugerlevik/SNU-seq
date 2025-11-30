#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Remove Top and Bottom 0.5% of Signal
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Set common directories
common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/MellorLab/SNUseqProject

# Directories for input and output
input_bedgraph_strandSpecific=$prefix/1_Umut/5.1_halfFiltered_bedgraph_strandSpecific

output_bedgraph_strandSpecific=$prefix/1_Umut/5.5_filtered_bedgraph_strandSpecific
output_bigwig_strandSpecific=$prefix/1_Umut/5.6_filtered_bigwig_strandSpecific
output_bedgraph_bothStrand=$prefix/1_Umut/5.7_filtered_bedgraph_bothStrand
output_bigwig_bothStrand=$prefix/1_Umut/5.8_filtered_bigwig_bothStrands

# Create directories for output
mkdir -p "$output_bedgraph_strandSpecific"
mkdir -p "$output_bigwig_strandSpecific"
mkdir -p "$output_bedgraph_bothStrand"
mkdir -p "$output_bigwig_bothStrand"

# Set genome size file for conversion to BigWig
genome_file="$common_prefix/genome/2_STARgenome/chrNameLength.txt"

# Function to calculate quantiles and filter bedGraph files
process_bedgraph_file () {
    local bedgraphFile=$1
    local outputTrimmedBedGraph=$2
    local sortedTrimmedBedGraph=$3
    local outputTrimmedBigWig=$4

    # Step 1: Identify top and bottom 0.5% quantiles
    total_lines=$(awk 'END {print NR}' "$bedgraphFile")
    bottomQuantile=$(awk '{print $4}' "$bedgraphFile" | sort -n | awk -v total_lines="$total_lines" 'NR == int(0.005 * total_lines) {print; exit}')
    topQuantile=$(awk '{print $4}' "$bedgraphFile" | sort -n | awk -v total_lines="$total_lines" 'NR == int(0.995 * total_lines) {print; exit}')

    echo "Processing $bedgraphFile: Bottom 0.5% quantile = $bottomQuantile, Top 0.5% quantile = $topQuantile"

    # Step 2: Filter out values outside the 0.5% range
    awk -v min="$bottomQuantile" -v max="$topQuantile" \
    '{ if ($4 >= min && $4 <= max) {print $0} }' \
    "$bedgraphFile" > "$outputTrimmedBedGraph"

    # Step 3: Sort the filtered BedGraph
    LC_COLLATE=C sort -k1,1 -k2,2n "$outputTrimmedBedGraph" > "$sortedTrimmedBedGraph"

    # Calculate the number of removed rows and the percentage of rows removed
    filtered_lines=$(awk 'END {print NR}' "$outputTrimmedBedGraph")
    removed_lines=$((total_lines - filtered_lines))
    removed_percentage=$(awk -v removed="$removed_lines" -v total="$total_lines" 'BEGIN {printf "%.2f", (removed/total) * 100}')

    echo "Rows removed: $removed_lines / $total_lines (${removed_percentage}%)"

    # Step 4: Convert to BigWig
    bedGraphToBigWig "$sortedTrimmedBedGraph" "$genome_file" "$outputTrimmedBigWig"

    echo "Processed $bedgraphFile: Trimmed BedGraph and BigWig saved to $sortedTrimmedBedGraph and $outputTrimmedBigWig"

    # Remove intermediate files
    rm "$outputTrimmedBedGraph"
}

# Step 1: Loop through each forward and reverse BedGraph file to apply top and bottom 0.5% filtering
for strand in fwd rev; do
    
    if [ "$strand" == "fwd" ]; then
        bedgraphList=$(ls $input_bedgraph_strandSpecific/*_fwd.bedgraph)
    elif [ "$strand" == "rev" ]; then
        bedgraphList=$(ls $input_bedgraph_strandSpecific/*_rev.bedgraph)
    fi

    for bedgraphFile in $bedgraphList; do
        bedgraphFileName=$(basename "$bedgraphFile")
        baseName=$(echo "$bedgraphFileName" | sed 's/.bedgraph//')

        # Define output BedGraph and BigWig files
        outputTrimmedBedGraph="$output_bedgraph_strandSpecific/${baseName}.bedgraph"
        sortedTrimmedBedGraph="$output_bedgraph_strandSpecific/${baseName}_sorted.bedgraph"
        outputTrimmedBigWig="$output_bigwig_strandSpecific/${baseName}.bw"

        process_bedgraph_file "$bedgraphFile" "$outputTrimmedBedGraph" "$sortedTrimmedBedGraph" "$outputTrimmedBigWig"
    done
done

# Step 2: Combine forward and reverse BedGraphs for bothStrand using bedtools unionbedg
for bedgraphFile in $(ls $output_bedgraph_strandSpecific/*_fwd_sorted.bedgraph); do
    
    baseName=$(basename "$bedgraphFile" | sed 's/_fwd_sorted.bedgraph//')
    revBedGraph="$output_bedgraph_strandSpecific/${baseName}_rev_sorted.bedgraph"
    outputBedGraphbothStrand="$output_bedgraph_bothStrand/${baseName}_bothStrand.bedgraph"
    outputBigWigbothStrand="$output_bigwig_bothStrand/${baseName}_bothStrand.bw"

    # Step 2.1: Use bedtools unionbedg to merge forward and reverse strands
    bedtools unionbedg -i "$bedgraphFile" "$revBedGraph" | \
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5)}' > "$outputBedGraphbothStrand.tmp"

    # Step 2.2: Sort the merged BedGraph file
    LC_COLLATE=C sort -k1,1 -k2,2n "$outputBedGraphbothStrand.tmp" > "$outputBedGraphbothStrand"

    # Step 2.3: Convert to BigWig for bothStrand data
    bedGraphToBigWig "$outputBedGraphbothStrand" "$genome_file" "$outputBigWigbothStrand"
    
    # Clean up intermediate files
    rm "$outputBedGraphbothStrand.tmp"

    echo "Processed $bedgraphFile: Merged BedGraph and BigWig saved to $outputBedGraphbothStrand and $outputBigWigbothStrand"
done


