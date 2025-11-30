#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Comprehensive QC - Part 2
## Date: Sep 29, 2024
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba deactivate

common_prefix=/MellorLab/SNUseqProject/commonFiles
prefix=/MellorLab/SNUseqProject

mappingStats_dir=$prefix/outputs/2_mappingStats
cd $mappingStats_dir
countFileList=$(ls *_featureCounts_geneid)

mdsplot_header=$common_prefix/multiqcHeaders/mdsplot_header.txt
heatmap_header=$common_prefix/multiqcHeaders/heatmap_header.txt

$prefix/scripts/4.2.1_edgeR_heatmap_MDS.R $countFileList
cat $mdsplot_header edgeR_MDS_Aplot_coordinates_mqc.csv | sed -e 's/X.MellorLab.SNUseqProject.1_Umut.2_mapped_bamFiles.//g' > tmp_file
mv tmp_file edgeR_MDS_Aplot_coordinates_mqc.csv
cat $heatmap_header log2CPM_sample_distances_mqc.csv | sed -e 's/X.MellorLab.SNUseqProject.1_Umut.2_mapped_bamFiles.//g' > tmp_file
mv tmp_file log2CPM_sample_distances_mqc.csv
