#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Build STAR genome with hg38
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

threads=8

prefix=/MellorLab/SNUseqProject/0_commonFiles

genome_dir=$prefix/genome
genome_fasta=$genome_dir/1_rawGenomeFiles/GRCh38spikes.fa

# Merge annotations with spikes and fix seqid based on the fasta headers in the genome
genome_gtf=$genome_dir/1_rawGenomeFiles/GRCh38.p14_GENCODE/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz
spike_gtf=$genome_dir/1_rawGenomeFiles/SpikeIns/allSpikeIns.gtf
gunzip -c $genome_gtf | grep -v ^## | cat $spike_gtf - > tmp.gtf
grep ">" $genome_fasta | sort > tmp_fasta_headers.txt
cut -f1 tmp.gtf |  uniq | sort > tmp_gtf_headers.txt
RScript $prefix/scripts/1.2_correctGTFandFASTAheaders.R
gxf2bed -i $genome_dir/1_rawGenomeFiles/gencode46spikes_seqidFixed.gtf -o $genome_dir/1_rawGenomeFiles/gencode46spikes_seqidFixed.bed
rm tmp_*_headers.txt tmp.gtf

# Set the corrected GTF file path
genome_gtf=$genome_dir/1_rawGenomeFiles/gencode46spikes_seqidFixed.gtf

# Build STAR genome
star_genome_dir=$genome_dir/2_STARgenome
mkdir -p $star_genome_dir

STAR \
    --runThreadN $threads \
    --runMode genomeGenerate \
    --genomeDir $star_genome_dir \
    --genomeFastaFiles $genome_fasta \
    --sjdbGTFfile $genome_gtf \
    --sjdbOverhang 100
