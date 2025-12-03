############################################################################
## Project: SNUseq project
## Script purpose: Prepare ATAC-seq reference point annotations
##                See GSE172053 for Anna's data used here
## Date: July 13, 2025
## Author: Umut Gerlevik
############################################################################

suppressPackageStartupMessages({
  library(rtracklayer)
  library(plyranges)
})

# ───────────────────────────────────────────────
# 1. Paths
# ───────────────────────────────────────────────
prefix <- "/MellorLab/SNUseqProject"
output_dir <- file.path(prefix, "4_Anna/1_fromAnna/2_ATACpeaks")
peaks_file <- file.path(output_dir, "ATAC_MergedPeaks.gff")
gtf_file <- file.path(prefix, "0_commonFiles/genome/1_rawGenomeFiles/gencode46spikes_seqidFixed.gtf")

# ───────────────────────────────────────────────
# 2. Load ATAC peaks and filter to chr1–22
# ───────────────────────────────────────────────
peaks <- import.gff(peaks_file)
length(unique(peaks$gene_id)) # 64536
peaks <- peaks[seqnames(peaks) %in% paste0("chr", 1:22)]
length(unique(peaks$gene_id)) # 63618

# ───────────────────────────────────────────────
# 3. Load genes and keep chr1–22 only
# ───────────────────────────────────────────────
genes <- import(gtf_file)
genes <- genes[genes$type == "gene" & seqnames(genes) %in% paste0("chr", 1:22)]

# ───────────────────────────────────────────────
# 4. Join nearest gene and assign strand (full set)
# ───────────────────────────────────────────────
peaks_joined <- peaks %>%
  join_nearest(genes, suffix = c("", ".gene"), distance = TRUE)

peaks_joined <- peaks_joined[peaks_joined$distance <= 1000]
peaks_joined$assigned_strand <- genes@strand[match(peaks_joined$gene_id.gene, genes$gene_id)]

# ───────────────────────────────────────────────
# 5. Split into fwd/rev and remove 5kb-close peaks
# ───────────────────────────────────────────────
remove_close_regions_5kb <- function(gr) {
  gr <- sort(gr, ignore.strand = TRUE)
  hits <- findOverlaps(gr, gr, maxgap = 4999, ignore.strand = TRUE)
  to_remove <- unique(queryHits(hits)[queryHits(hits) != subjectHits(hits)])
  gr[-to_remove]
}

peaks_fwd <- peaks_joined[peaks_joined$assigned_strand == "+"]
peaks_rev <- peaks_joined[peaks_joined$assigned_strand == "-"]

peaks_fwd <- remove_close_regions_5kb(peaks_fwd)
length(unique(peaks_fwd$gene_id)) # 17367
peaks_rev <- remove_close_regions_5kb(peaks_rev)
length(unique(peaks_rev$gene_id)) # 16570

# ───────────────────────────────────────────────
# 6. Prepare unstranded peaks: >1kb from genes, exclude already assigned gene_ids
# ───────────────────────────────────────────────
gene_ids_stranded <- unique(peaks_joined$gene_id)

# Remove peaks linked to <1kb annotations
peaks_unstranded <- peaks[!peaks$gene_id %in% gene_ids_stranded]

# Also exclude peaks overlapping within ±1kb window (additional check)
genes_1kb <- resize(genes, width = width(genes) + 2000, fix = "center")
peaks_unstranded <- peaks_unstranded[!(peaks_unstranded %over% genes_1kb)]

# Enforce ≥5kb spacing
peaks_unstranded <- remove_close_regions_5kb(peaks_unstranded)
length(unique(peaks_unstranded$gene_id)) # 13434

# ───────────────────────────────────────────────
# 7. All peaks (chr1-22 + 5kb spacing)
# ───────────────────────────────────────────────
peaks_all <- remove_close_regions_5kb(peaks)
length(unique(peaks_all$gene_id)) # 45289

# ───────────────────────────────────────────────
# 7.5 Make a Venn Diagram to check
# ───────────────────────────────────────────────
library(ggVennDiagram)

# Extract unique gene IDs from each group
venn_sets <- list(
  FWD         = unique(peaks_fwd$gene_id),
  REV         = unique(peaks_rev$gene_id),
  Unstranded  = unique(peaks_unstranded$gene_id)#,
  # All         = unique(peaks_all$gene_id)
)

# Plot
ggVennDiagram(venn_sets, label_alpha = 0, category.names = names(venn_sets)) +
  ggplot2::theme_void() +
  ggplot2::ggtitle("Venn diagram of ATAC peak gene_id assignments")

# ───────────────────────────────────────────────
# 8. Export BED and GTF
# ───────────────────────────────────────────────
export_bed <- function(gr, path, extra_cols = NULL) {
  df <- as.data.frame(gr)[, c("seqnames", "start", "end")]
  if (!is.null(extra_cols)) {
    df <- cbind(df, mcols(gr)[, extra_cols, drop = FALSE])
  }
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

export_bed(peaks_all,        file.path(output_dir, "ATAC_peaks_filtered.bed"))
export_bed(peaks_fwd,        file.path(output_dir, "ATAC_peaks_filtered_fwd.bed"), extra_cols = "assigned_strand")
export_bed(peaks_rev,        file.path(output_dir, "ATAC_peaks_filtered_rev.bed"), extra_cols = "assigned_strand")
export_bed(peaks_unstranded, file.path(output_dir, "ATAC_peaks_filtered_unstranded.bed"))

export.gff(peaks_all,        file.path(output_dir, "ATAC_peaks_filtered.gtf"))
export.gff(peaks_fwd,        file.path(output_dir, "ATAC_peaks_filtered_fwd.gtf"))
export.gff(peaks_rev,        file.path(output_dir, "ATAC_peaks_filtered_rev.gtf"))
export.gff(peaks_unstranded, file.path(output_dir, "ATAC_peaks_filtered_unstranded.gtf"))

message("✅ ATAC-seq annotationsgenerated and exported.")
