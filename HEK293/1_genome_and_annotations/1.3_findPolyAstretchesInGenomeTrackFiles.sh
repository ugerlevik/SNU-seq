#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Find the genomic polyA stretches
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

: '
# Mask reads primed at genomic A-stretches
* adapted https://github.com/manschmi/MexNab_3seq of Schmid et al. (2018), who adapted from Roy et al. (2016).

Notes from https://github.com/manschmi/MexNab_3seq of Schmid et al. (2018), quoted as it is:

"Next step is to mask genomic pA reads from the genome. 
This could possible be included in the above step. But the approach I used is to create a “masking file”. 
i.e., a bed file of all genomic positions to be filtered. 
(This should allow for a more flexible approach for comparing different filters.) 
The main script to create the masking file is again implemented in python 2, 
which takes various parameters (how many nucleotides X in a window of size Y) to control what is being masked. 
The criteria are taken from Kevin Roy, Chanfreau lab. Which he assembled based on yeast Lexogen data. They are:  

1)	all regions with >= 4A within 6 nucleotides but no C or Ts (downstream of position).  
2)	all regions with >= 12A within 18 nucleotides (no CT restrictions, downstream of position).  
3)	All regions with >=15A within 18 nucleotides (upstream of position).  

I would say number 3 is optional, and the technical reason for this filtering is somewhat obscure. 
Apparently QuantSeq protocol can produce artifact reads containing polyT stretches, 
and these will map to A-rich sequences in the genome, which can be filtered away using criterion 3. 
On our first run this was not obvious and I did not include that.

Scans the genome for genomic A stretches and prints them in bed format.
Usage:
    flag_genomic_As_KevinRoy_like.py -i genome.fa -l [word_len] -A [min As within word] \
        -maxCT [max C+Ts within word] -el [extend output left] -er [extend output right] -o [output_file]
-i ... genome in fasta format
-l ... length of motif, ie word length (default = 6)
-A ... minimum number of A for word to be flagged (default = 4)
-maxCT ... maximum number of C and Ts for word to be flagged (default = 2)
-el ... extend each single nucleotide hit by this number upstream (default = 0)
-er ... extend each single nucleotide hit by this number downstream (default = 0)
-o ... output file

Intended for combination as input for bedtools subtract of A-tail based 3p end sequencing."
'

source ~/.zshrc
mamba activate deeptools_env

prefix="/MellorLab/SNUseqProject/0_commonFiles"

rawGenome_dir="$prefix/genome/1_rawGenomeFiles"
kevinroy_dir="$prefix/filterRegions"

flag_genomic_As_KevinRoy_like="$prefix/scripts/flag_genomic_As_KevinRoy_like.py"

## make all uppercase in the genome
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' $rawGenome_dir/GRCh38.p14.fa > $rawGenome_dir/GRCh38.p14_uppercase.fa
genomeUpper="$rawGenome_dir/GRCh38.p14_uppercase.fa"

## get A6 regions
mamba activate py27 # environment with python2
python $flag_genomic_As_KevinRoy_like $genomeUpper -l 6 -A 4 -maxCT 0 -o $kevinroy_dir/GRCh38.genome_flagged_A4in6.bed
mamba activate deeptools_env

awk '{if($6=="+") print $0}' $kevinroy_dir/GRCh38.genome_flagged_A4in6.bed | \
awk '$2 >= 0' | \
sort -k1,1 -k2,2n -o $kevinroy_dir/GRCh38.genome_flagged_A4in6_plus.bed
bedtools merge -d 2 -i $kevinroy_dir/GRCh38.genome_flagged_A4in6_plus.bed > $kevinroy_dir/GRCh38.genome_flagged_A4in6_plus_merged.bed

awk '{if($6=="-") print $0}' $kevinroy_dir/GRCh38.genome_flagged_A4in6.bed | \
awk '$2 >= 0' | \
sort -k1,1 -k2,2n -o $kevinroy_dir/GRCh38.genome_flagged_A4in6_minus.bed
bedtools merge -d 2 -i $kevinroy_dir/GRCh38.genome_flagged_A4in6_minus.bed > $kevinroy_dir/GRCh38.genome_flagged_A4in6_minus_merged.bed


## get A18 regions
mamba activate py27
python $flag_genomic_As_KevinRoy_like $genomeUpper -l 18 -A 12 -maxCT 6 -o $kevinroy_dir/GRCh38.genome_flagged_A12in18.bed
mamba activate deeptools_env

awk '{if($6=="+") print $0}' $kevinroy_dir/GRCh38.genome_flagged_A12in18.bed | \
awk '$2 >= 0' | \
sort -k1,1 -k2,2n -o $kevinroy_dir/GRCh38.genome_flagged_A12in18_plus.bed
bedtools merge -d 2 -i $kevinroy_dir/GRCh38.genome_flagged_A12in18_plus.bed > $kevinroy_dir/GRCh38.genome_flagged_A12in18_plus_merged.bed

awk '{if($6=="-") print $0}' $kevinroy_dir/GRCh38.genome_flagged_A12in18.bed | \
awk '$2 >= 0' | \
sort -k1,1 -k2,2n -o $kevinroy_dir/GRCh38.genome_flagged_A12in18_minus.bed
bedtools merge -d 2 -i $kevinroy_dir/GRCh38.genome_flagged_A12in18_minus.bed > $kevinroy_dir/GRCh38.genome_flagged_A12in18_minus_merged.bed


## get A15 out of A18 regions
mamba activate py27
python $flag_genomic_As_KevinRoy_like $genomeUpper -l 18 -A 15 -maxCT 3 -o $kevinroy_dir/GRCh38.genome_flagged_A15in18.bed
mamba activate deeptools_env

#these are for upstream filtering ... need to shift
awk '{OFS="\t"}{if($6=="+"){$2+=19; $3+=19; print $0}}' $kevinroy_dir/GRCh38.genome_flagged_A15in18.bed | \
awk '$2 >= 0' | \
sort -k1,1 -k2,2n -o $kevinroy_dir/GRCh38.genome_flagged_A15in18_plus.bed
bedtools merge -d 2 -i $kevinroy_dir/GRCh38.genome_flagged_A15in18_plus.bed > $kevinroy_dir/GRCh38.genome_flagged_A15in18_plus_merged.bed

awk '{OFS="\t"}{if($6=="-"){ $2-=19; $3-=19; print $0}}' $kevinroy_dir/GRCh38.genome_flagged_A15in18.bed | \
awk '$2 >= 0' | \
sort -k1,1 -k2,2n -o $kevinroy_dir/GRCh38.genome_flagged_A15in18_minus.bed
bedtools merge -d 2 -i $kevinroy_dir/GRCh38.genome_flagged_A15in18_minus.bed > $kevinroy_dir/GRCh38.genome_flagged_A15in18_minus_merged.bed

## merge the above
cat $kevinroy_dir/GRCh38.genome_flagged_A12in18_plus_merged.bed \
    $kevinroy_dir/GRCh38.genome_flagged_A4in6_plus_merged.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i stdin > $kevinroy_dir/GRCh38.genome_flagged_KevinRoy_plus.bed
cat $kevinroy_dir/GRCh38.genome_flagged_A12in18_minus_merged.bed $kevinroy_dir/GRCh38.genome_flagged_A4in6_minus_merged.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i stdin > $kevinroy_dir/GRCh38.genome_flagged_KevinRoy_minus.bed

## single file
awk '{if(NR == FNR){strand="+"}else{strand="-"};print $0"\t.\t0\t"strand}' \
    $kevinroy_dir/GRCh38.genome_flagged_KevinRoy_plus.bed \
    $kevinroy_dir/GRCh38.genome_flagged_KevinRoy_minus.bed | \
    sort -k1,1 -k2,2n > $kevinroy_dir/GRCh38.p14_KevinRoyAregions.bed

## single mask human only
#grep -v -E 'SPIKEIN' $kevinroy_dir/GRCh38.p14_KevinRoyAregions.bed > $kevinroy_dir/genomicAmask_noSpike.bed

## count positions
awk '{sum+=($3-$2)}END{print sum}' $kevinroy_dir/GRCh38.p14_KevinRoyAregions.bed
# 186108971

## genome size
head $rawGenome_dir/../2_STARgenome/chrNameLength.txt
: '
chr1    248956422
chr1_KI270706v1_random  175055
chr1_KI270707v1_random  32032
chr1_KI270708v1_random  127682
chr1_KI270709v1_random  66860
chr1_KI270710v1_random  40176
chr1_KI270711v1_random  42210
chr1_KI270712v1_random  176043
chr1_KI270713v1_random  40745
chr1_KI270714v1_random  41717
'

awk '{sum+=$2}END{print sum*2}' $rawGenome_dir/../2_STARgenome/chrNameLength.txt
# 6596861272 (total nt in both strands!!)
## 186108971/6596861272 = 0.02821174545 --> 2.8% of genome are masked


# Clean everything
rm $genomeUpper
rm $kevinroy_dir/GRCh38.genome_flagged_A4in6*
rm $kevinroy_dir/GRCh38.genome_flagged_A12in18*
rm $kevinroy_dir/GRCh38.genome_flagged_A15in18*
rm $kevinroy_dir/GRCh38.genome_flagged_KevinRoy_*us.bed

