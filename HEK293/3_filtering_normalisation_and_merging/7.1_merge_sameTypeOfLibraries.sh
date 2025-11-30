#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Merge same type of libraries by using 
##                  spike-in normalized BedGraph files
## Date: Mar 5, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Set directories and parameters
common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/MellorLab/SNUseqProject

# Define output directories for strand-specific and bothStrand files
input_bedgraph_folder=$prefix/1_Umut/6.2_spikeinNormalised_bedgraph
input_bigwig_folder=$prefix/1_Umut/6.3_spikeinNormalised_bigwig
input_bedgraph_bothStrand=$prefix/1_Umut/6.4_spikeinNormalised_bedgraph_bothStrand
input_bigwig_bothStrand=$prefix/1_Umut/6.5_spikeinNormalised_bigwig_bothStrand

output_bedgraph_strandSpecific=$prefix/1_Umut/7.1_merged_bedgraph_strandSpecific
output_bigwig_strandSpecific=$prefix/1_Umut/7.2_merged_bigwig_strandSpecific
output_bedgraph_bothStrand=$prefix/1_Umut/7.3_merged_bedgraph_bothStrand
output_bigwig_bothStrand=$prefix/1_Umut/7.4_merged_bigwig_bothStrand

genome_file="$common_prefix/genome/2_STARgenome/chrNameLength.txt"

# Create output directories if they don't exist
mkdir -p "$output_bedgraph_strandSpecific"
mkdir -p "$output_bigwig_strandSpecific"
mkdir -p "$output_bedgraph_bothStrand"
mkdir -p "$output_bigwig_bothStrand"

# Function to merge files for each group
merge_group_files () {
    local group="$1"
    local sample_list="$2"

    # Define merged output files for each strand
    local merged_bedgraph_fwd="$output_bedgraph_strandSpecific/${group}_fwd_merged.bedgraph"
    local merged_bedgraph_rev="$output_bedgraph_strandSpecific/${group}_rev_merged.bedgraph"
    local merged_bigwig_fwd="$output_bigwig_strandSpecific/${group}_fwd_merged.bw"
    local merged_bigwig_rev="$output_bigwig_strandSpecific/${group}_rev_merged.bw"

    echo "Merging forward and reverse strands for group: $group"

    # Gather forward and reverse strand files
    local file_list_fwd=""
    local file_list_rev=""
    for sample in $sample_list; do
        file_fwd="$input_bedgraph_folder/${sample}_fwd_spikeinNormalised.bedgraph"
        file_rev="$input_bedgraph_folder/${sample}_rev_spikeinNormalised.bedgraph"

        # Append files if they exist and print info
        if [ -f "$file_fwd" ]; then
            file_list_fwd+="$file_fwd "
            echo "Adding $file_fwd to forward strand merge for group $group"
        else
            echo "Warning: Forward strand file not found for sample $sample, skipping."
        fi

        if [ -f "$file_rev" ]; then
            file_list_rev+="$file_rev "
            echo "Adding $file_rev to reverse strand merge for group $group"
        else
            echo "Warning: Reverse strand file not found for sample $sample, skipping."
        fi
    done

    # Check if there's only one file in file_list_fwd and file_list_rev
    if [ $(echo "$file_list_fwd" | wc -w) -eq 1 ]; then
        cp "$file_list_fwd" "$merged_bedgraph_fwd"
        cp "$input_bigwig_folder/${group}_fwd_spikeinNormalised.bw" "$merged_bigwig_fwd"
    elif [ -n "$file_list_fwd" ]; then
        bedtools unionbedg -i $file_list_fwd | \
        awk 'BEGIN {OFS="\t"} {sum=0; for (i=4; i<=NF; i++) sum+=$i; print $1, $2, $3, sum}' > "$merged_bedgraph_fwd.tmp"
        LC_COLLATE=C sort -k1,1 -k2,2n "$merged_bedgraph_fwd.tmp" > "$merged_bedgraph_fwd"
        rm "$merged_bedgraph_fwd.tmp"
        bedGraphToBigWig "$merged_bedgraph_fwd" "$genome_file" "$merged_bigwig_fwd"
    fi

    if [ $(echo "$file_list_rev" | wc -w) -eq 1 ]; then
        cp "$file_list_rev" "$merged_bedgraph_rev"
        cp "$input_bigwig_folder/${group}_rev_spikeinNormalised.bw" "$merged_bigwig_rev"
    elif [ -n "$file_list_rev" ]; then
        bedtools unionbedg -i $file_list_rev | \
        awk 'BEGIN {OFS="\t"} {sum=0; for (i=4; i<=NF; i++) sum+=$i; print $1, $2, $3, sum}' > "$merged_bedgraph_rev.tmp"
        LC_COLLATE=C sort -k1,1 -k2,2n "$merged_bedgraph_rev.tmp" > "$merged_bedgraph_rev"
        rm "$merged_bedgraph_rev.tmp"
        bedGraphToBigWig "$merged_bedgraph_rev" "$genome_file" "$merged_bigwig_rev"
    fi

    echo "Merged forward and reverse strands for $group saved as BigWig."

    # Combine forward and reverse merged BedGraphs into bothStrand files
    local both_strand_bedgraph="$output_bedgraph_bothStrand/${group}_bothStrand_merged.bedgraph"
    local both_strand_bigwig="$output_bigwig_bothStrand/${group}_bothStrand_merged.bw"

    # Check if there’s only one file in both forward and reverse lists
    if [ $(echo "$file_list_fwd" | wc -w) -eq 1 ] && [ $(echo "$file_list_rev" | wc -w) -eq 1 ]; then
        cp "$input_bedgraph_bothStrand/${group}_bothStrand_spikeinNormalised.bedgraph" "$both_strand_bedgraph"
        cp "$input_bigwig_bothStrand/${group}_bothStrand_spikeinNormalised.bw" "$both_strand_bigwig"
    else
        # Use bedtools unionbedg to merge the forward and reverse merged files
        bedtools unionbedg -i "$merged_bedgraph_fwd" "$merged_bedgraph_rev" | \
        awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4 + $5)}' > "$both_strand_bedgraph.tmp"
        LC_COLLATE=C sort -k1,1 -k2,2n "$both_strand_bedgraph.tmp" > "$both_strand_bedgraph"
        rm "$both_strand_bedgraph.tmp"

        # Convert the combined BedGraph to BigWig
        bedGraphToBigWig "$both_strand_bedgraph" "$genome_file" "$both_strand_bigwig"
    fi

    echo "BothStrand Merged BigWig for $group saved to $both_strand_bigwig"
}

# Define sample lists for each group

merge_group_files "labelled_bPAP" "labelled_H1_bPAP labelled_H2_bPAP labelled_S1_bPAP"
merge_group_files "labelled_noPAP" "labelled_H1_noPAP labelled_H2_noPAP labelled_S1_noPAP"
merge_group_files "labelled_bPAP_rRNA" "labelled_H3_bPAP labelled_H4_bPAP labelled_S3_bPAP"
merge_group_files "labelled_noPAP_rRNA" "labelled_H3_noPAP labelled_H4_noPAP labelled_S3_noPAP"
merge_group_files "total_bPAP" "total_H1_bPAP total_H2_bPAP total_S1_bPAP total_S4_bPAP"
merge_group_files "total_noPAP" "total_H1_noPAP total_H2_noPAP total_S1_noPAP total_S4_noPAP"
merge_group_files "total_bPAP_rRNA" "total_H3_bPAP total_H4_bPAP total_S3_bPAP"
merge_group_files "total_noPAP_rRNA" "total_H3_noPAP total_H4_noPAP total_S3_noPAP"

: '
# These are put manually because only one sample exists in these groups
merge_group_files "unlabelled_bPAP" "unlabelled_S4_bPAP"
merge_group_files "unlabelled_noPAP" "unlabelled_S4_noPAP"
'

echo "Merging and conversion for all groups complete."

