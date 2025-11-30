#!/bin/bash

############################################################################
## Project: SNUseq project
## Script purpose: Rename fastq files
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

## mapping file: new names on the first column, old names on the second column
source_folder=/MellorLab/_fromLexogen_allSNUseqRawFiles_Umut_2022_2024
target_folder=/MellorLab/SNUseqProject/1_Umut/0_raw_fastqFiles

mapping_file=$target_folder/renameFastqMapping.txt
# These names and renames work on my local machine. However, in the GEO upload,
# these fastq files were named differently, so please follow the descriptions
# in the GEO, to match and rename correctly if needed.

tr -d '\r' < $mapping_file |
while read new_name current_name seq_name; do
    echo "Copying $current_name as $new_name" 
    cp ${source_folder}/$current_name ${target_folder}/$new_name
done
