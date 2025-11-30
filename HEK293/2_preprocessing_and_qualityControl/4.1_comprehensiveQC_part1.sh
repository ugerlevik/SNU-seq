#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Comprehensive QC - Part 1
## Date: Sep 29, 2024
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

common_prefix=/MellorLab/SNUseqProject/commonFiles
prefix=/MellorLab/SNUseqProject

mappedFile_dir=$prefix/1_Umut/2_mapped_bamFiles
mappingStats_dir=$prefix/outputs/2_mappingStats

bamFileList=$(ls $mappedFile_dir/*.bam)

gtf=$common_prefix/genome/1_rawGenomeFiles/gencode46spikes_seqidFixed.gtf
bed=$common_prefix/genome/1_rawGenomeFiles/gencode46spikes_seqidFixed.bed

hkgenes=$common_prefix/genome/3_annotations/sources/hg38.HouseKeepingGenes.bed

biotypes_header=$common_prefix/multiqcHeaders/biotypes_header.txt
alias mqc_features_stat.py=$prefix/scripts/4.1.1_mqc_features_stat.py

for file in ${bamFileList[@]}
do
    bamBaseName=`basename $file _trimmed_Aligned.sortedByCoord.out.bam`
    echo " Now working on $bamBaseName "

    bam_stat.py -i $file -q 20 > $mappingStats_dir/${bamBaseName}_bam_stat.txt
    clipping_profile.py -i $file -q 20 -s "SE" -o $mappingStats_dir/${bamBaseName}
    # deletion_profile.py -i $file -q 20 -l 101 -o $mappingStats_dir/${bamBaseName}
    infer_experiment.py -i $file  -q 20 -r $bed> $mappingStats_dir/${bamBaseName}_infer_experiment.txt
    insertion_profile.py -i $file -q 20 -s "SE" -o $mappingStats_dir/${bamBaseName}
    junction_annotation.py -i $file -q 20 -r $bed -o $mappingStats_dir/${bamBaseName}
    junction_saturation.py -i $file -q 20 -r $bed -o $mappingStats_dir/${bamBaseName}
    mismatch_profile.py -i $file -q 20 -l 101 -o $mappingStats_dir/${bamBaseName}
    read_distribution.py -i $file -r $hkgenes > $mappingStats_dir/${bamBaseName}_housekeepingGenes_read_distribution.txt
    read_duplication.py -i $file -q 20 -o $mappingStats_dir/${bamBaseName}
    RPKM_saturation.py -i $file -q 20 -r $bed -d "++,--" -o $mappingStats_dir/${bamBaseName}
    
    featureCounts -a $gtf -g "gene_id" -s 1 -o $mappingStats_dir/${bamBaseName}_featureCounts_geneid $file
    featureCounts -a $gtf -g "gene_type" -s 1 -o $mappingStats_dir/${bamBaseName}_featureCounts_genetype $file
    cut -f 1,7 $mappingStats_dir/${bamBaseName}_featureCounts_genetype | tail -n +3 | cat $biotypes_header - > $mappingStats_dir/${bamBaseName}_biotype_counts_mqc.txt
    mqc_features_stat.py $mappingStats_dir/${bamBaseName}_biotype_counts_mqc.txt -s ${bamBaseName} -f rRNA -o $mappingStats_dir/${bamBaseName}_biotype_counts_gs_mqc.tsv

    stringtie $file -v -e --fr -G $gtf -o $mappingStats_dir/${bamBaseName}_transcripts.gtf \
    -A $mappingStats_dir/${bamBaseName}.gene_abund.txt -C $mappingStats_dir/${bamBaseName}.cov_refs.gtf \
    -b $mappingStats_dir/${bamBaseName}_ballgown
done

bamList=$(ls $mappedFile_dir/*.bam | sed -e '$ ! s/$/,/g' | tr -d '\n')
# housekeeping might be useful for geneBody_coverage.py, for the others, whole genome annotations in bed12 format should be good
geneBody_coverage.py -i $bamList  -r $hkgenes -o $mappingStats_dir/output_geneBody_coverage_housekeeping
#geneBody_coverage.py -i $bamList  -r $bed -o $mappingStats_dir/output_geneBody_coverage_all

tin.py -i $bamList -r $hkgenes
sed 's/_trimmed_Aligned.sortedByCoord.out.bam//g' *.summary.txt | awk 'NR==1; NR > 1 && !/^Bam_file/' > $mappingStats_dir/TIN.summary.txt
mv *.summary.txt $mappingStats_dir

#featureCountsList=($mappingStats_dir/*_featureCounts_geneid)
featureCountsList=$(ls $mappingStats_dir/*_featureCounts_geneid)
csvtk join -t -f "Geneid,Start,Length,End,Chr,Strand" $featureCountsList | csvtk cut -t -f "-Start,-Chr,-End,-Length,-Strand" | sed 's/_trimmed_Aligned.sortedByCoord.out.bam//g' > $mappingStats_dir/merged_gene_featureCounts.txt

# cleanup
rm tmp_summary