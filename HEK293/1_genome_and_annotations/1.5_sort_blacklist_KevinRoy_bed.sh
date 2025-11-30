#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Sort and merge the regions to be masked
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

source ~/.zshrc
mamba activate deeptools_env

cat /MellorLab/SNUseqProject/0_commonFiles/filterRegions/merged_blacklist_KevinRoyA.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i stdin > /MellorLab/SNUseqProject/0_commonFiles/filterRegions/merged_blacklist_KevinRoyA_sorted_merged.bed
