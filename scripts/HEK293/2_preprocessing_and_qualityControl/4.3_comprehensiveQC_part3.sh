#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Comprehensive QC - Part 3
## Date: Sep 29, 2024
## Author: Umut Gerlevik
############################################################################

SECONDS=0

source ~/.zshrc
mamba activate multiqc_env

prefix=/MellorLab/SNUseqProject

logFiles_dir=$prefix/scripts/logs
trimmingStats_dir=$prefix/outputs/1_trimmingStats
mappingStats_dir=$prefix/outputs/2_mappingStats
multiqcSummary_dir=$prefix/outputs/3_multiqcSummary

mkdir -p $multiqcSummary_dir

multiqc $logFiles_dir $trimmingStats_dir $mappingStats_dir -o $multiqcSummary_dir
