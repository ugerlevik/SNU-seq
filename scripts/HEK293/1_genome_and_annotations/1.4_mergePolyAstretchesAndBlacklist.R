############################################################################
## Project: SNUseq project
## Script purpose: Merge genomic-A and blacklist regions to mask together
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

library(data.table)

prefix <- "/MellorLab/SNUseqProject/0_commonFiles"

kevin <- fread(paste0(prefix, "/filterRegions/GRCh38.p14_KevinRoyAregions.bed"))
kevin <- kevin[, c(1:3, 6)]
colnames(kevin)[4] <- "V4"

black <- fread(paste0(prefix, "/filterRegions/hg38-blacklist.v2.bed"))
black <- black[, c(1:3)]
black$V4 <- "."

merged <- rbind(black, kevin)
write.table(merged, paste0(prefix, "/filterRegions/merged_blacklist_KevinRoyA.bed"), 
            quote = F, row.names = F, col.names = F, sep = "\t")
