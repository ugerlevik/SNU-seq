############################################################################
## Project: SNUseq project
## Script purpose: Prepare the genome annotations
## Date: Mar 13, 2025
## Author: Umut Gerlevik
############################################################################

prefix <- "/MellorLab/SNUseqProject/0_commonFiles"

setwd(prefix)

#############################################
# 1) Load libraries and data
#############################################
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
suppressWarnings(suppressPackageStartupMessages(library(plyranges)))

# Set output directory
output_dir <- "annotations/"

# Import your housekeeping genes BED file
# https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38.HouseKeepingGenes.bed.gz/download
hkg <- import("annotations/sources/hg38.HouseKeepingGenes.bed")
length(unique(hkg$name)) # 3798

#############################################
# 2) Collapse transcripts into gene-level
#############################################
genes <- hkg %>%
  filter(seqnames %in% paste0("chr", 1:22)) 
length(unique(genes$name)) # 3684

#############################################
# 3) Filter gene-level intervals
#   (length >= 1 kb and discard genes < 3.5 kb from each other)
#############################################

# a) Remove genes that are within 3.5 kb of each other
remove_close_regions_3.5kb <- function(gr) {
  # Sort ignoring strand so that close regions across strands are compared
  gr <- sort(gr, ignore.strand = TRUE)
  
  if (length(gr) == 0) return(gr)  # Return if empty
  
  # Find all regions that are within 3500 bp of any other region
  overlaps <- findOverlaps(gr, gr, maxgap = 3499, ignore.strand = TRUE)
  
  # Get unique indices of regions that overlap with another region (excluding self-hits)
  regions_to_remove <- unique(queryHits(overlaps)[queryHits(overlaps) != subjectHits(overlaps)])
  
  # Remove all overlapping regions
  filtered_gr <- gr[-regions_to_remove]
  
  return(filtered_gr)
}
genes_filtered <- remove_close_regions_3.5kb(genes)
# genes_filtered <- sort(genes_filtered, ignore.strand = TRUE)
length(unique(genes_filtered$name)) # 2783

# b) Keep genes >= 1 kb
genes_filtered <- genes_filtered[width(genes_filtered) >= 1000]
length(unique(genes_filtered$name)) # 2774

genes_filtered[genes_filtered$name %in% genes_filtered$name[which(duplicated(genes_filtered$name))]]

genes_filtered <- as.data.frame(genes_filtered)
# Now 'genes_filtered' are our final set of genes: 
# length >= 1 kb, and at least 3.5 kb from each other.

#############################################
# 4) Save
#############################################
write.table(genes_filtered %>% select(seqnames, start, end, name, strand), 
            paste0(output_dir, "housekeeping_genes.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(genes_filtered %>% filter(strand == "+") %>% select(seqnames, start, end, name), 
            paste0(output_dir, "housekeeping_genes_fwd.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(genes_filtered %>% filter(strand == "-") %>% select(seqnames, start, end, name), 
            paste0(output_dir, "housekeeping_genes_rev.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
