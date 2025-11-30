#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: featureCounts spikeins
## Date: Mar 5, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

threads=8

# Set common directories
common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/MellorLab/SNUseqProject

# Set directory paths
mappedFile_dir=$prefix/1_Umut/2_mapped_bamFiles
output_dir=$prefix/1_Umut/6.1_featureCounts_spikeins

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Set 3' UTR GTF annotation file
gtf=$common_prefix/annotations/onlySpikein.gtf

# Get list of BAM files in mapped directory
bamlist=$(ls $mappedFile_dir/*.bam)

# Loop through each BAM file and run featureCounts
for bamfile in $bamlist; do
    bamBaseName=$(basename "$bamfile" | sed 's/_trimmed_Aligned.sortedByCoord.out.bam//')

    echo "Processing $bamfile..."
    
    # Run featureCounts for each BAM file
    featureCounts -T $threads -a $gtf -g "gene_id" -t "exon" -s 1 -o $output_dir/${bamBaseName}_featureCounts_spikeins.txt "$bamfile"
    
    echo "Finished counting for $bamfile"
done

# Combine featureCounts outputs
echo "Merging featureCounts results..."


featureCountsList=$(ls $output_dir/*_featureCounts_spikeins.txt | grep -v "summary")

# Merge all featureCounts outputs into a single file and clean up
csvtk join -t -f "Geneid,Start,Length,End,Chr,Strand" $featureCountsList | \
csvtk cut -t -f "-Start,-Chr,-End,-Length,-Strand" | \
sed 's/_trimmed_Aligned.sortedByCoord.out.bam//g' | \
sed 's|/MellorLab/SNUseqProject/1_Umut/2_mapped_bamFiles/||g' > $output_dir/merged_spikeins_featureCounts.txt

echo "Merged counts saved to $output_dir/merged_spikeins_featureCounts.txt"

