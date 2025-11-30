#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Map fastq files to the genome
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

common_prefix=/MellorLab/SNUseqProject/0_commonFiles
prefix=/MellorLab/SNUseqProject

threads=8

trimmedFile_dir=$prefix/1_Umut/1_trimmed_fastqFiles
trimmedFileList=($trimmedFile_dir/*_trimmed.fastq.gz)

genome_dir=$common_prefix/genome/2_STARgenome
mappedFile_dir=$prefix/1_Umut/2_mapped_bamFiles

mappingStats_dir=$prefix/outputs/2_mappingStats

mkdir -p $mappedFile_dir
mkdir -p $mappingStats_dir

# Mapping to the genome w/ spike-ins
for file in ${trimmedFileList[@]}
do
	echo "Now working on $file"

	STAR \
	--runMode alignReads \
	--runThreadN $threads \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes All \
	--genomeDir $genome_dir \
	--readFilesIn $file \
	--readFilesCommand gunzip -c \
	--outFileNamePrefix $mappedFile_dir/`basename $file _trimmed.fastq.gz`_trimmed_ \
	--outFilterMultimapNmax 10 \
	--winAnchorMultimapNmax 50 \
	--alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 10 \
	--outFilterMismatchNoverReadLmax 1 \
	--outMultimapperOrder Random \
	--alignEndsType EndToEnd \
    --alignIntronMin 11 \
	--alignIntronMax 0 \
	--alignMatesGapMax 0

	samtools index $mappedFile_dir/`basename $file _trimmed.fastq.gz`_trimmed_Aligned.sortedByCoord.out.bam

done

mv $mappedFile_dir/*.out $mappingStats_dir/.
mv $mappedFile_dir/*.tab $mappingStats_dir/.

