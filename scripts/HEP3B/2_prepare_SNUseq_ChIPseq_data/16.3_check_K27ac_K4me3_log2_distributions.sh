#!/bin/bash

############################################################################
## Project: SNUseq project
## Script: Plot log2(K27ac) vs log2(K4me3) signal distributions to check for 
##        consistency across replicates and concatenated data.
## Author: Umut Gerlevik
## Date: July 9, 2025
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Input BigWigs
k27ac_dir="/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_K27ac"
k4me3_dir="/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_K4me3"

k27ac_rep1="${k27ac_dir}/GSM5240739_0_K27ac_rep1_log2.bw"
k27ac_rep2="${k27ac_dir}/GSM5240740_0_K27ac_rep2_log2.bw"
k27ac_concat="${k27ac_dir}/K27ac_0_concatenatedSum_log2.bw"

k4me3_rep1="${k4me3_dir}/GSM5240747_0_K4me3_rep1_log2.bw"
k4me3_rep2="${k4me3_dir}/GSM5240748_0_K4me3_rep2_log2.bw"
k4me3_concat="${k4me3_dir}/K4me3_0_concatenatedSum_log2.bw"

# Region BED (ATAC peaks as example)
bed_file="/MellorLab/SNUseqProject/4_Anna/1_fromAnna/2_ATACpeaks/ATAC_peaks_filtered.bed"

# Output
outdir="/MellorLab/SNUseqProject/4_Anna/102_distribution_validation_K27ac_K4me3"
mkdir -p "$outdir"

npz="${outdir}/log2_K27ac_K4me3_signal_summary.npz"
tab="${outdir}/log2_K27ac_K4me3_signal_summary.tab"
plot="${outdir}/violin_log2_K27ac_K4me3.pdf"

echo "🔁 Running multiBigwigSummary..."
multiBigwigSummary BED-file \
  --BED "$bed_file" \
  -b "$k27ac_rep1" "$k27ac_rep2" "$k27ac_concat" \
     "$k4me3_rep1" "$k4me3_rep2" "$k4me3_concat" \
  -out "$npz" \
  --outRawCounts "$tab" \
  -p 8

echo "📊 Generating violin plot in R..."

Rscript --vanilla - "$tab" "$plot" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
tabfile <- args[1]
plotfile <- args[2]

df <- read.table(tabfile, header=TRUE, sep="\t", comment.char="")
colnames(df)[1:3] <- c("chr", "start", "end")

df_clean <- df[complete.cases(df), ]

# Convert to long format
library(reshape2)
library(ggplot2)

df_melt <- melt(df_clean[, -c(1:3)], variable.name="Sample", value.name="Log2Signal")

ggplot(df_melt, aes(x=Sample, y=Log2Signal, fill=Sample)) +
  geom_violin(trim=FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, outlier.size=0.3, alpha=0.3) +
  theme_bw() +
  labs(title="Distribution of log2(K27ac) and log2(K4me3) signals",
       y="log2(signal + pseudocount)", x=NULL) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="none")

ggsave(plotfile, width=8, height=6)
EOF

echo "✅ Done. Plot saved to: $plot"
