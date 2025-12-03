#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Splicing efficiency (SE) calculation using SPLICE-q
## Date: Mar 11, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate spliceQ_env

# Set directories
common_prefix="/MellorLab/SNUseqProject/0_commonFiles"
prefix="/MellorLab/SNUseqProject"

input_dir="$prefix/1_Umut/2_mapped_bamFiles"
gtf_file="$common_prefix/genome/1_rawGenomeFiles/gencode46spikes_seqidFixed.gtf"

output_dir="$prefix/1_Umut/3.1_spliceQ_splicingEfficiency"
mkdir -p $output_dir

# Function to run SPLICE-q on one sample
run_spliceq() {
    bam_file="$1"
    gtf_file="$2"
    output_dir="$3"

    sample_name=$(basename "$bam_file" _trimmed_Aligned.sortedByCoord.out.bam)
    echo "Processing $sample_name"
    
    # Run SPLICE-q
    SPLICE-q.py -b "$bam_file" -g "$gtf_file" -o "$output_dir/${sample_name}_spliceQ.tsv"

    # Check if SPLICE-q completed successfully
    if [ $? -ne 0 ]; then
        echo "Error: SPLICE-q failed for $sample_name"
        return 1
    fi
    
    # Verify that the expected output file was created
    if [ ! -f "$output_dir/${sample_name}_spliceQ.tsv" ]; then
        echo "Error: Output file for $sample_name not found"
        return 1
    fi
    
    echo "Completed $sample_name successfully"
    return 0
}

export -f run_spliceq  # Export the function for parallel

# Use GNU Parallel to process all BAM files in parallel
ls $input_dir/*.bam | parallel --halt now,fail=1 -j 4 run_spliceq {} "$gtf_file" "$output_dir"

