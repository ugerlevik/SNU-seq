############################################################################
## Project: SNUseq project
## Script purpose: Prepare the genome annotations
## Date: Mar 4, 2025
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
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Import your GENCODE v46 annotation GTF (with spike, seqid fixed, etc.)
gencode <- import("genome/1_rawGenomeFiles/gencode46spikes_seqidFixed.gtf")
length(unique(gencode$gene_id)) # 70611 + 4 spike-ins

#############################################
# 2) Collapse transcripts into gene-level
#############################################
# GENCODE has lines with type == "gene" that already define each gene’s span
# We can use them directly, or collapse all transcripts for each gene_id manually.
# Here, we'll rely on the gene entries:
genes <- gencode %>%
  filter(type == "gene",
         # gene_type == "protein_coding",
         # grepl("^ENSG", gene_id),
         seqnames %in% paste0("chr", 1:22)
         ) 
length(unique(genes$gene_id)) # 59950

# Now 'genes' should be one line per gene in GENCODE. 
# # If GENCODE has multiple lines for the same gene, we can reduce them:
# genes <- genes %>%
#   group_by(gene_id) %>%
#   reduce_ranges_directed() %>%
#   ungroup()
# length(unique(genes$gene_id))  # 59950
# Note: reduce_ranges() from plyranges merges overlapping intervals for each group.

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
length(unique(genes_filtered$gene_id)) # 12224

# b) Keep genes >= 1 kb
genes_filtered <- genes_filtered[width(genes_filtered) >= 1000]
length(unique(genes_filtered$gene_id)) # 6962

# c) Keep "protein_coding" genes
genes_filtered <- genes_filtered %>%
  filter(gene_id %in% gencode$gene_id[gencode$gene_type == "protein_coding"])
length(unique(genes_filtered$gene_id)) # 2807

# Now 'genes_filtered' are our final set of genes: protein-coding, 
# length >= 1 kb, and at least 3.5 kb from each other.


################################################################################################
# 4) Find 3'UTRs within ±100 bp of gene end UTRs are not annotated as 5' or 3' in GENCODE GTF, 
# so we needed to infer based on strand and gene end as below. However, later disciovered that 
# GENCODE GFF3 formatted files have explicitly categorised 3' UTR annotations. So those can be
# used directly in the future analysis instead of this code block determining 3' UTRs.
################################################################################################
# 3' UTRs: focus on protein-coding genes and transcripts
utr3 <- gencode %>%
  filter(type == "UTR", 
         seqnames %in% paste0("chr", 1:22),
         gene_type == "protein_coding",
         strand != "*",  # Only valid stranded regions
         transcript_support_level %in% 1:5) 
# %>%
#   group_by(gene_id, transcript_id) %>%
#   reduce_ranges_directed() %>%
#   ungroup()

# Annotate 3' ends of genes based on strand
genes_ends <- genes_filtered %>%
  as.data.frame() %>% 
  mutate(
    gene_tes = ifelse(strand == "+", end, start)  # TES: end for "+" strand, start for "-"
  ) %>%
  select(gene_id, strand, gene_tes)

# Merge 3' UTRs with gene ends
utr3_with_gene_ends <- merge(as.data.frame(utr3), as.data.frame(genes_ends), by = "gene_id")
# any(utr3_with_gene_ends$strand.x != utr3_with_gene_ends$strand.y)

# Filter UTRs based on their proximity to gene TES
utr3_filtered <- utr3_with_gene_ends %>%
  filter((strand.x == "+" & abs(end - gene_tes) <= 100) | 
           (strand.x == "-" & abs(start - gene_tes) <= 100)) %>%
  mutate(strand = strand.x) %>% 
  as_granges()
length(utr3_filtered$gene_id) # 4153

# Collapse overlapping UTRs for each gene (if multiple transcripts per gene)
utr3_filtered_unique <- utr3_filtered %>%
  group_by(gene_id) %>%
  reduce_ranges_directed() %>%
  ungroup()

# Manually select or merge any duplicated gene_id for separate
utr3_filtered_unique$gene_id[duplicated(utr3_filtered_unique$gene_id)]

# utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000167104.11"]
# utr3_filtered_unique <- utr3_filtered_unique %>%
#   filter(start != 33044045)
# end(utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000167104.11"]) <- 33044047
# utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000167104.11"]
# 
# utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000103266.11"]
# utr3_filtered_unique <- utr3_filtered_unique %>%
#   filter(start != 682807)
# end(utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000103266.11"]) <- 682870
# utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000103266.11"]
# 
# utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000204348.10"]
# utr3_filtered_unique <- utr3_filtered_unique %>%
#   filter(start != 31969907)
# end(utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000204348.10"]) <- 31970024
# utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000204348.10"]
# 
# utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000171811.14"]
# utr3_filtered_unique <- utr3_filtered_unique %>%
#   filter(start != 132808443)
# end(utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000171811.14"]) <- 132808794
# utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000171811.14"]
# 
# utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000187726.9"]
# utr3_filtered_unique <- utr3_filtered_unique %>%
#   filter(start != 73970283)
# end(utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000187726.9"]) <- 73970366
# utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000187726.9"]
# 
utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000213339.9"]
utr3_filtered_unique <- utr3_filtered_unique %>%
  filter(start != 10713369)
end(utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000213339.9"]) <- 10713437
utr3_filtered_unique[utr3_filtered_unique$gene_id == "ENSG00000213339.9"]

utr3_filtered_unique$gene_id[duplicated(utr3_filtered_unique$gene_id)]

length(utr3_filtered_unique$gene_id) # 2394

# Export BED files for 3' UTRs
write.table(as.data.frame(utr3_filtered_unique) %>% select(seqnames, start, end, gene_id, strand),
            paste0(output_dir, "filtered_3UTRs.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(as.data.frame(utr3_filtered_unique) %>% filter(strand == "+") %>% select(seqnames, start, end, gene_id), 
            paste0(output_dir, "filtered_3UTRs_fwd.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(as.data.frame(utr3_filtered_unique) %>% filter(strand == "-") %>% select(seqnames, start, end, gene_id), 
            paste0(output_dir, "filtered_3UTRs_rev.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

export.gff(utr3_filtered_unique, paste0(output_dir, "filtered_3UTRs.gtf"))
export.gff(utr3_filtered_unique %>% filter(strand == "+"), paste0(output_dir, "filtered_3UTRs_fwd.gtf"))
export.gff(utr3_filtered_unique %>% filter(strand == "-"), paste0(output_dir, "filtered_3UTRs_rev.gtf"))

# # Check how many of the final gene subset have 3' UTR or not
# genes_filtered_ids <- genes_filtered$gene_id
# utr3_gene_ids <- unique(utr3_filtered_unique$gene_id)
# 
# genes_with_utr <- intersect(genes_filtered_ids, utr3_gene_ids)
# num_genes_with_utr <- length(genes_with_utr)
# 
# genes_without_utr <- setdiff(genes_filtered_ids, utr3_gene_ids)
# num_genes_without_utr <- length(genes_without_utr)
# 
# cat("Number of genes in genes_filtered WITH UTR annotations:", num_genes_with_utr, "\n")
# # Number of genes in genes_filtered WITH UTR annotations: 1958
# cat("Number of genes in genes_filtered WITHOUT UTR annotations:", num_genes_without_utr, "\n")
# # Number of genes in genes_filtered WITHOUT UTR annotations: 336


#############################################
# 5) Prepare data frames for genes and 3'UTRs
#############################################
library(dplyr)
utr3_df <- as.data.frame(utr3_filtered_unique) %>%
  rowwise() %>%
  mutate(
    # For + strand, the 3'UTR starts at 'start'
    # For - strand, the 3'UTR "start" in a transcription sense is actually the max genomic coordinate
    # but we’ll call it "utr3_start_coord" to unify usage.
    utr3_start_coord = ifelse(strand == "+", start, end)
  ) %>%
  ungroup() %>%
  select(gene_id, strand, utr3_start_coord)

genes_df <- as.data.frame(genes_filtered) %>%
  rowwise() %>%
  mutate(
    gene_start = start, 
    gene_end   = end,
    # TSS coordinate:
    gene_tss   = ifelse(strand == "+", start, end),
    # TES coordinate:
    gene_tes   = ifelse(strand == "+", end, start)
  ) %>%
  ungroup() %>%
  select(seqnames, start, end, strand, gene_id, gene_tss, gene_tes, width)

genes_with_utr <- left_join(genes_df, utr3_df, by = c("gene_id", "strand")) %>% 
  filter(!is.na(utr3_start_coord))
length(unique(genes_with_utr$gene_id)) # 2394


write.table(genes_with_utr %>% select(seqnames, start, end, gene_id, strand), 
            paste0(output_dir, "finalAnnotationSubset.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(genes_with_utr %>% filter(strand == "+") %>% select(seqnames, start, end, gene_id), 
            paste0(output_dir, "finalAnnotationSubset_fwd.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(genes_with_utr %>% filter(strand == "-") %>% select(seqnames, start, end, gene_id), 
            paste0(output_dir, "finalAnnotationSubset_rev.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

export.gff(genes_with_utr, paste0(output_dir, "finalAnnotationSubset.gtf"))
export.gff(genes_with_utr %>% filter(strand == "+"), paste0(output_dir, "finalAnnotationSubset_fwd.gtf"))
export.gff(genes_with_utr %>% filter(strand == "-"), paste0(output_dir, "finalAnnotationSubset_rev.gtf"))

#############################################
# 6) Define regions
#############################################
# TSS region: [TSS-50, TSS+200]
# GeneBody region: [TSS+201, 3'UTR_start-1]
# GeneEnd region: [3'UTR_start, TES+50]
# Readthrough region: [TES+51, TES+3000]

define_subregions_for_gene <- function(chr, strand, gene_id,
                                       gene_tss, gene_tes, utr3_start_coord) {
  # We assume gene_tss < gene_tes if strand == "+"
  # and gene_tss > gene_tes if strand == "-".
  # But we want to produce intervals in ascending genomic coordinates.
  
  # We'll define numeric starts/ends for each region, then reorder with min(), max().
  
  # TSS region
  if (strand == "+") {
    tss_start <- gene_tss - 50
    tss_end   <- gene_tss + 200
  } else {
    # If strand == "-", TSS is the higher coordinate
    # so TSS-50 means +50 in genomic space
    tss_start <- gene_tss + 50
    tss_end   <- gene_tss - 200
  }
  
  # GeneBody region
  # from (TSS+201) to (3'UTR_start-1) for + strand
  # or from (TSS-201) to (3'UTR_start-1) for - strand,
  # but we unify logic with the same approach:
  if (strand == "+") {
    body_start <- gene_tss + 201
    body_end   <- utr3_start_coord - 1
  } else {
    body_start <- utr3_start_coord + 1
    body_end   <- gene_tss - 201
  }
  
  # GeneEnd region
  # from 3'UTR_start to (TES+50) for +,
  # from (TES-50) to 3'UTR_start for -.
  if (strand == "+") {
    end_start <- utr3_start_coord
    end_end   <- gene_tes + 50
  } else {
    end_start <- gene_tes - 50
    end_end   <- utr3_start_coord
  }
  
  # Readthrough region
  # from (TES+51) to (TES+3000) for +,
  # from (TES-3000) to (TES-51) for -.
  if (strand == "+") {
    rt_start <- gene_tes + 51
    rt_end   <- gene_tes + 3000
  } else {
    rt_start <- gene_tes - 3000
    rt_end   <- gene_tes - 51
  }
  
  # Build intervals with ascending genomic coords
  build_region <- function(s, e, region_type) {
    # reorder so start < end
    rstart <- min(s, e)
    rend   <- max(s, e)
    if (rend - rstart < 1) return(NULL)  # skip zero/negative-length intervals
    GRanges(seqnames = chr,
            ranges   = IRanges(start = rstart, end = rend),
            strand   = strand,
            gene_id  = gene_id,
            region_type = region_type)
  }
  
  gr_list <- GRangesList()
  gr_list[[1]] <- build_region(tss_start, tss_end, "TSS")
  gr_list[[2]] <- build_region(body_start, body_end, "GeneBody")
  gr_list[[3]] <- build_region(end_start, end_end, "GeneEnd")
  gr_list[[4]] <- build_region(rt_start, rt_end, "Readthrough")
  
  # Return a flattened GRanges
  unlist(gr_list)
}

# Single core version
# all_regions_list <- vector("list", nrow(genes_with_utr))
# 
# for (i in seq_len(nrow(genes_with_utr))) {
#   rowdat <- genes_with_utr[i, ]
#   
#   # Extract needed columns
#   chr  <- rowdat$seqnames
#   strd <- as.character(rowdat$strand)
#   gid  <- rowdat$gene_id
#   tss  <- as.numeric(rowdat$gene_tss)
#   tes  <- as.numeric(rowdat$gene_tes)
#   utr3 <- as.numeric(rowdat$utr3_start_coord)
#   
#   # Call our helper
#   subregions <- define_subregions_for_gene(chr, strd, gid, tss, tes, utr3)
#   
#   # Store
#   all_regions_list[[i]] <- subregions
# }

# Parallel version
library(doParallel)
nCores <- parallel::detectCores(logical = TRUE) - 2
cl <- makeCluster(nCores)
registerDoParallel(cl)

all_regions_list <- foreach(i = seq_len(nrow(genes_with_utr)), 
                            .packages = c("GenomicRanges", "plyranges")) %dopar% {
                              rowdat <- genes_with_utr[i, ]
                              
                              # Extract columns
                              chr  <- rowdat$seqnames
                              strd <- as.character(rowdat$strand)
                              gid  <- rowdat$gene_id
                              tss  <- as.numeric(rowdat$gene_tss)
                              tes  <- as.numeric(rowdat$gene_tes)
                              utr3 <- as.numeric(rowdat$utr3_start_coord)
                              
                              # Call your subregion function
                              subregions <- define_subregions_for_gene(chr, strd, gid, tss, tes, utr3)
                              
                              # Return result
                              subregions
                            }
stopCluster(cl)

# Combine all GRanges
all_regions <- do.call(c, all_regions_list)

# Filter out any NULL or zero-length intervals
all_regions <- all_regions[width(all_regions) > 0]
length(all_regions)

# Check if any overlaps
ov <- findOverlaps(all_regions, all_regions, ignore.strand = FALSE)
ov_df <- as.data.frame(ov)
ov_df_self <- subset(ov_df, queryHits != subjectHits)
if (nrow(ov_df_self) == 0) {
  cat("No overlaps detected among the intervals in all_regions.\n")
} else {
  cat("Overlaps detected! Number of overlapping pairs:", nrow(ov_df_self), "\n")
  head(ov_df_self)  # Preview the first few overlapping pairs
}

# Separate each region
genes_with_utr_gr <- as_granges(genes_with_utr)
TSS_gr <- all_regions[all_regions$region_type == "TSS"]
GeneBody_gr <- all_regions[all_regions$region_type == "GeneBody"]
GeneEnd_gr <- all_regions[all_regions$region_type == "GeneEnd"]
Readthrough_gr <- all_regions[all_regions$region_type == "Readthrough"]

# Export each subset as BED and GTF
TSS_gr %>% as.data.frame() %>% select(seqnames, start, end, gene_id, strand) %>% 
  write.table(., paste0(output_dir, "/TSS.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
TSS_gr %>% as.data.frame() %>% filter(strand == "+") %>% select(seqnames, start, end, gene_id) %>% 
  write.table(., paste0(output_dir, "/TSS_fwd.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
TSS_gr %>% as.data.frame() %>% filter(strand == "-") %>% select(seqnames, start, end, gene_id) %>% 
  write.table(., paste0(output_dir, "/TSS_rev.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
export.gff(TSS_gr, paste0(output_dir, "/TSS.gtf"))
export.gff(TSS_gr %>% filter(strand == "+"), paste0(output_dir, "/TSS_fwd.gtf"))
export.gff(TSS_gr %>% filter(strand == "-"), paste0(output_dir, "/TSS_rev.gtf"))

GeneBody_gr %>% as.data.frame() %>% select(seqnames, start, end, gene_id, strand) %>% 
  write.table(., paste0(output_dir, "/GeneBody.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
GeneBody_gr %>% as.data.frame() %>% filter(strand == "+") %>% select(seqnames, start, end, gene_id) %>% 
  write.table(., paste0(output_dir, "/GeneBody_fwd.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
GeneBody_gr %>% as.data.frame() %>% filter(strand == "-") %>% select(seqnames, start, end, gene_id) %>% 
  write.table(., paste0(output_dir, "/GeneBody_rev.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
export.gff(GeneBody_gr, paste0(output_dir, "/GeneBody.gtf"))
export.gff(GeneBody_gr %>% filter(strand == "+"), paste0(output_dir, "/GeneBody_fwd.gtf"))
export.gff(GeneBody_gr %>% filter(strand == "-"), paste0(output_dir, "/GeneBody_rev.gtf"))

GeneEnd_gr %>% as.data.frame() %>% select(seqnames, start, end, gene_id, strand) %>% 
  write.table(., paste0(output_dir, "/GeneEnd.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
GeneEnd_gr %>% as.data.frame() %>% filter(strand == "+") %>% select(seqnames, start, end, gene_id) %>% 
  write.table(., paste0(output_dir, "/GeneEnd_fwd.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
GeneEnd_gr %>% as.data.frame() %>% filter(strand == "-") %>% select(seqnames, start, end, gene_id) %>% 
  write.table(., paste0(output_dir, "/GeneEnd_rev.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
export.gff(GeneEnd_gr, paste0(output_dir, "/GeneEnd.gtf"))
export.gff(GeneEnd_gr %>% filter(strand == "+"), paste0(output_dir, "/GeneEnd_fwd.gtf"))
export.gff(GeneEnd_gr %>% filter(strand == "-"), paste0(output_dir, "/GeneEnd_rev.gtf"))

Readthrough_gr %>% as.data.frame() %>% select(seqnames, start, end, gene_id, strand) %>% 
  write.table(., paste0(output_dir, "/Readthrough.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
Readthrough_gr %>% as.data.frame() %>% filter(strand == "+") %>% select(seqnames, start, end, gene_id) %>% 
  write.table(., paste0(output_dir, "/Readthrough_fwd.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
Readthrough_gr %>% as.data.frame() %>% filter(strand == "-") %>% select(seqnames, start, end, gene_id) %>% 
  write.table(., paste0(output_dir, "/Readthrough_rev.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
export.gff(Readthrough_gr, paste0(output_dir, "/Readthrough.gtf"))
export.gff(Readthrough_gr %>% filter(strand == "+"), paste0(output_dir, "/Readthrough_fwd.gtf"))
export.gff(Readthrough_gr %>% filter(strand == "-"), paste0(output_dir, "/Readthrough_rev.gtf"))

