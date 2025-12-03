#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Adapter and quality trimming
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/MellorLab/SNUseqProject

rawFileList=($prefix/1_Umut/0_raw_fastqFiles/*.fastq.gz)
adapterFile=$common_prefix/sequencingAdapters/QuantSeqAdaptersIndex.fa

trimmedFile_dir=$prefix/1_Umut/1_trimmed_fastqFiles
trimmingStats_dir=$prefix/outputs/1_trimmingStats

mkdir -p $trimmedFile_dir
mkdir -p $trimmingStats_dir

for file in ${rawFileList[@]}
do
	
	echo "Now working on $file"
	
	## Run fastqc on raw data
	fastqc $file -o $trimmingStats_dir

	## Adapters, polyA and quality trimming
	bbduk.sh in=$file out="$trimmedFile_dir/`basename $file .fastq.gz`_trimmed.fastq.gz" \
	ref=$adapterFile k=12 ktrim=r useshortkmers=t mink=5 threads=10 literal=AAAAAAAAAAAA \
	qtrim=w trimq=20 minlength=0 stats=$trimmingStats_dir/`basename $file .fastq.gz`_bbduk_stats.out \
	2>&1 | tee $trimmingStats_dir/`basename $file .fastq.gz`.out

	## Run fastqc to check after trimming
	fastqc $trimmedFile_dir/`basename $file .fastq.gz`_trimmed.fastq.gz -o $trimmingStats_dir
	
done
