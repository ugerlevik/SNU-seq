############################################################################
## Project: SNUseq project
## Script purpose: Prepare filtered intergenic FANTOM5 enhancer annotations
## Date: July 13, 2025
## Author: Umut Gerlevik
############################################################################

suppressPackageStartupMessages({
  library(rtracklayer)
  library(plyranges)
})

prefix <- "/MellorLab/SNUseqProject"
output_dir <- file.path(prefix, "3_publicData_HEK293T/FANTOM5")

# -------------------------
# 1. Load FANTOM5 Enhancers
# -------------------------
fan5 <- import(file.path(output_dir, "F5.hg38.enhancers.bed"))
length(unique(fan5$name)) # 63285
fan5 <- fan5 %>% filter(seqnames %in% paste0("chr", 1:22))
fan5_original <- fan5
length(unique(fan5_original$name)) # 61825 in chr1-22

# -------------------------------
# 2. Remove enhancers near genes
# -------------------------------
genes <- import.gff(file.path(prefix, "0_commonFiles/genome/1_rawGenomeFiles/gencode46spikes_seqidFixed.gtf"))
genes <- genes[genes$type == "gene" & seqnames(genes) %in% paste0("chr", 1:22)]

# Expand gene body ±5 kb and exclude overlapping enhancers
genes_5kb <- resize(genes, width = width(genes) + 10000, fix = "center")
fan5 <- fan5[!fan5 %over% genes_5kb]
length(unique(fan5$name)) # 12093

# -------------------------------
# 3. Remove enhancer pairs < 5 kb
# -------------------------------
remove_close_regions_5kb <- function(gr) {
  gr <- sort(gr, ignore.strand = TRUE)
  if (length(gr) == 0) return(gr)
  overlaps <- findOverlaps(gr, gr, maxgap = 4999, ignore.strand = TRUE)
  remove_idx <- unique(queryHits(overlaps)[queryHits(overlaps) != subjectHits(overlaps)])
  gr[-remove_idx]
}

fan5_filtered <- remove_close_regions_5kb(fan5)
length(unique(fan5_filtered$name)) # 6911

# -------------------------------
# 4. Sort and Export
# -------------------------------
seqlevels(fan5_filtered) <- sortSeqlevels(seqlevels(fan5_filtered))
fan5_filtered <- sort(fan5_filtered, ignore.strand = TRUE)

write.table(as.data.frame(fan5_filtered)[, 1:3],
            file.path(output_dir, "FANTOM5_filtered.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

bad_cols <- vapply(mcols(fan5_filtered), function(col) is.list(col) || is(col, "List"), logical(1))
if (any(bad_cols)) {
  message("Removing non-atomic metadata columns: ", paste(names(which(bad_cols)), collapse = ", "))
  mcols(fan5_filtered)[, bad_cols] <- NULL
}
export.gff(fan5_filtered, file.path(output_dir, "FANTOM5_filtered.gtf"))

# -------------------------------
# 5. Summary
# -------------------------------
cat("✅ Final Enhancers:", length(fan5_filtered), "\n")
cat("Filtered from:", length(fan5_original), "\n")
