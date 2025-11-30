#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Filter blacklist, A-stretches and non chr1 to chr22 (or spike-in)
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Set common directories
common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/MellorLab/SNUseqProject

# Set directory paths
bedgraphFolder=$prefix/1_Umut/4.1_raw_bedgraphFiles_genomecov3end

# Directories for filtered output
filteredFolder=$prefix/1_Umut
filtered_bedgraph_strandSpecific=$filteredFolder/5.1_halfFiltered_bedgraph_strandSpecific
filtered_bigwig_strandSpecific=$filteredFolder/5.2_halfFiltered_bigwig_strandSpecific
filtered_bedgraph_bothStrand=$filteredFolder/5.3_halfFiltered_bedgraph_bothStrand
filtered_bigwig_bothStrand=$filteredFolder/5.4_halfFiltered_bigwig_bothStrand

# Create directories for filtered output
mkdir -p "$filtered_bedgraph_strandSpecific"
mkdir -p "$filtered_bigwig_strandSpecific"
mkdir -p "$filtered_bedgraph_bothStrand"
mkdir -p "$filtered_bigwig_bothStrand"

# Bed file for blacklist and A-stretch regions
filterBedFile="$common_prefix/filterRegions/merged_blacklist_KevinRoyA_sorted_merged.bed"

# Set genome size file for conversion to BigWig
genome_file="$common_prefix/genome/2_STARgenome/chrNameLength.txt"

# Step 1: Filter forward and reverse strand BedGraphs
for strand in fwd rev; do
    
    if [ "$strand" == "fwd" ]; then
        bedgraphList=$(ls $bedgraphFolder/*_fwd_sorted.bedgraph)
    elif [ "$strand" == "rev" ]; then
        bedgraphList=$(ls $bedgraphFolder/*_rev_sorted.bedgraph)
    fi

    # Loop through each BedGraph file in the directory
    for bedgraphFile in $bedgraphList; do
        
        bedgraphFileName=$(basename "$bedgraphFile")
        baseName=$(echo "$bedgraphFileName" | sed 's/_sorted.bedgraph//')

        # Define output BedGraph and BigWig files for filtered 1_Umut (without "filtered" or "sorted" in the name)
        outputFilteredBedGraph="$filtered_bedgraph_strandSpecific/${baseName}.bedgraph"
        outputBigWigStrandSpecific="$filtered_bigwig_strandSpecific/${baseName}.bw"

        # Step 1.1: Subtract blacklist and A-stretches using bedtools subtract
        bedtools subtract -a "$bedgraphFile" -b "$filterBedFile" > "$outputFilteredBedGraph.tmp"

        # Step 1.2: Remove non-chr1 to chr22 regions and keep UGSPIKE
        awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ || $1 == "UGSPIKE"' "$outputFilteredBedGraph.tmp" > "$outputFilteredBedGraph"

        # Step 1.3: Convert sorted BedGraph to BigWig for strand-specific 1_Umut
        bedGraphToBigWig "$outputFilteredBedGraph" "$genome_file" "$outputBigWigStrandSpecific"
        
        # Clean up intermediate files
        rm "$outputFilteredBedGraph.tmp"

        echo "Processed $bedgraphFile: Filtered BedGraph and BigWig saved to $outputFilteredBedGraph and $outputBigWigStrandSpecific"
    done
done

# Step 2: Combine forward and reverse BedGraphs for bothStrand using bedtools unionbedg
for bedgraphFile in $(ls $filtered_bedgraph_strandSpecific/*_fwd.bedgraph); do
    
    baseName=$(basename "$bedgraphFile" | sed 's/_fwd.bedgraph//')
    revBedGraph="$filtered_bedgraph_strandSpecific/${baseName}_rev.bedgraph"
    outputBedGraphbothStrand="$filtered_bedgraph_bothStrand/${baseName}_bothStrand.bedgraph"
    outputBigWigbothStrand="$filtered_bigwig_bothStrand/${baseName}_bothStrand.bw"

    # Step 2.1: Use bedtools unionbedg to merge forward and reverse strands
    bedtools unionbedg -i "$bedgraphFile" "$revBedGraph" | \
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5)}' > "$outputBedGraphbothStrand.tmp"

    # Step 2.2: Sort the merged BedGraph file
    LC_COLLATE=C sort -k1,1 -k2,2n "$outputBedGraphbothStrand.tmp" > "$outputBedGraphbothStrand"

    # Step 2.3: Convert to BigWig for bothStrand 1_Umut
    bedGraphToBigWig "$outputBedGraphbothStrand" "$genome_file" "$outputBigWigbothStrand"
    
    # Clean up intermediate files
    rm "$outputBedGraphbothStrand.tmp"

    echo "Processed $bedgraphFile: Merged BedGraph and BigWig saved to $outputBedGraphbothStrand and $outputBigWigbothStrand"
done

