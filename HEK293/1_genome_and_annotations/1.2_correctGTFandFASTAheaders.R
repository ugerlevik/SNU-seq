############################################################################
## Project: SNUseq project
## Script purpose: Correct GTF and FASTA headers for genome
## Date: Mar 4, 2025
## Author: Umut Gerlevik
############################################################################

setwd("/MellorLab/SNUseqProject/0_commonFiles")

gtf_head <- readLines("tmp_gtf_headers.txt")
fa_head <- readLines("tmp_fasta_headers.txt")

fa_head <- gsub(">", "", fa_head)

any(gtf_head %in% fa_head)

common <- intersect(fa_head, gtf_head)

gtf_head <- gsub("\\.", "v", gtf_head)

noncommon <- setdiff(gtf_head, common)

tmp <- data.frame(original = c(), found = c())
za <- c()
zu <- c()
for(i in noncommon) {
  if(any(grepl(i, fa_head))) {
    tmp <- rbind(tmp, data.frame(original = i,
                                 found = fa_head[grepl(i, fa_head)]))
  } else if (length(fa_head[grepl(i, fa_head)]) > 1) {
    za <- c(za, i)
  } else {
    zu <- c(zu, i)
  }
}

library(rtracklayer)
gencode <- readGFF("tmp.gtf")

gencode$seqid <- as.character(gencode$seqid)
gencode$seqid <- gsub("\\.", "v", gencode$seqid)

all(unique(gencode$seqid) %in% unique(tmp$original)) # only noncommon
unique(gencode$seqid)[!(unique(gencode$seqid) %in% unique(tmp$original))]
all(unique(tmp$original) %in% unique(gencode$seqid))

for(i in tmp$original) {
  gencode$seqid[gencode$seqid == i] <- tmp$found[tmp$original == i]
}

export.gff(gencode, "/MellorLab/SNUseqProject/commonFiles/genome/1_rawGenomeFiles/gencode46spikes_seqidFixed.gtf", format = "gtf")
