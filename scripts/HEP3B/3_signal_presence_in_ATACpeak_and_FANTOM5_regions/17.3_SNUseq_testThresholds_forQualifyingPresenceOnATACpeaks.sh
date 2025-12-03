#!/usr/bin/env bash

############################################################################
## Project: SNUseq project
## Script:  Test SNU-seq thresholds for qualifying presence on ATAC peaks
## Author:  Umut Gerlevik
## Date:    Aug 27, 2025
############################################################################

source ~/.zshrc
mamba activate deeptools_env

threads=8

# ─────────── Inputs ───────────
ATAC_BED="/MellorLab/SNUseqProject/4_Anna/1_fromAnna/2_ATACpeaks/ATAC_peaks_filtered.bed"
K27AC_BW="/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_K27ac/K27ac_0_concatenatedSum_log2.bw"
SNUSEQ_BW="/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_SNUseq_mergedStrands/SNUseq_0_concatenatedSum_log2.bw"
K4ME3_BW="/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_K4me3/K4me3_0_concatenatedSum_log2.bw"

# Fixed thresholds (log2(count+1))
K27_THR=1.25
K4_THR=1.25

# Sweep these for SNU-seq (log2(count+1))
SNU_THRESH_LIST=("0" "0.5" "1" "1.25" "1.5" "2")

# ─────────── Outputs ───────────
OUT_DIR="/MellorLab/SNUseqProject/4_Anna/111b_ATAC_SNUseq_threshold_sweep"
mkdir -p "$OUT_DIR"

TAB="${OUT_DIR}/ATAC_K27ac_SNUseq_K4me3_means.tab"
NPZ="${OUT_DIR}/ATAC_K27ac_SNUseq_K4me3_means.npz"
SUMMARY="${OUT_DIR}/summary_SNUseq_threshold_sweep.tsv"
HIST_SNU="${OUT_DIR}/SNUseq_on_ATAC_histogram.pdf"
QFILE="${OUT_DIR}/quantiles_SNUseq.txt"

# ─────────── 1) Mean signal per peak (K27ac, SNU-seq, K4me3) ───────────
echo "▶ multiBigwigSummary on ATAC peaks (K27ac + SNU-seq + K4me3)…"
multiBigwigSummary BED-file \
  --BED "$ATAC_BED" \
  -b "$K27AC_BW" "$SNUSEQ_BW" "$K4ME3_BW" \
  -out "$NPZ" \
  --outRawCounts "$TAB" \
  -p $threads

TOTAL=$(wc -l < "$ATAC_BED")
echo "Number of ATAC peaks: $TOTAL"

# ─────────── 1.5) Quick quantiles for SNU-seq ───────────
awk 'NR>1 {print $5}' "$TAB" | \
Rscript -e 'x<-scan("stdin",quiet=TRUE); x<-x[is.finite(x)]; print(quantile(x,c(.5,.8,.9,.95,.99),na.rm=TRUE))' \
| tee "$QFILE"

# ─────────── 2) Sweep SNU thresholds (K27 & K4 fixed) ───────────
echo -e "snuseq_threshold\tclass\tcount\tpercent" > "$SUMMARY"

pct() { awk -v c="$1" -v t="$TOTAL" 'BEGIN{printf("%.2f", (t?100*c/t:0))}'; }

for STHR in "${SNU_THRESH_LIST[@]}"; do
  echo "▶ K27ac ≥ $K27_THR ; K4me3 ≥ $K4_THR ; sweeping SNU-seq > $STHR"

  BOTH=$(awk -v k27="$K27_THR" -v k4thr="$K4_THR" -v sthr="$STHR" '
    NR==1{next}
    {ka=($4=="nan"||$4=="NaN"?0:$4); su=($5=="nan"||$5=="NaN"?0:$5); k4=($6=="nan"||$6=="NaN"?0:$6)}
    ka>=k27 && k4>=k4thr && su>sthr {c++} END{print c+0}' "$TAB")

  K4pos_SNUneg=$(awk -v k27="$K27_THR" -v k4thr="$K4_THR" -v sthr="$STHR" '
    NR==1{next}
    {ka=($4=="nan"||$4=="NaN"?0:$4); su=($5=="nan"||$5=="NaN"?0:$5); k4=($6=="nan"||$6=="NaN"?0:$6)}
    ka>=k27 && k4>=k4thr && su<=sthr {c++} END{print c+0}' "$TAB")

  K4neg_SNUpos=$(awk -v k27="$K27_THR" -v k4thr="$K4_THR" -v sthr="$STHR" '
    NR==1{next}
    {ka=($4=="nan"||$4=="NaN"?0:$4); su=($5=="nan"||$5=="NaN"?0:$5); k4=($6=="nan"||$6=="NaN"?0:$6)}
    ka>=k27 && k4<k4thr  && su>sthr {c++} END{print c+0}' "$TAB")

  Neither=$(awk -v k27="$K27_THR" -v k4thr="$K4_THR" -v sthr="$STHR" '
    NR==1{next}
    {ka=($4=="nan"||$4=="NaN"?0:$4); su=($5=="nan"||$5=="NaN"?0:$5); k4=($6=="nan"||$6=="NaN"?0:$6)}
    ka>=k27 && k4<k4thr  && su<=sthr {c++} END{print c+0}' "$TAB")

  printf "%s\tSNU+_K4+\t%d\t%s\n" "$STHR" "$BOTH"           "$(pct "$BOTH")"           >> "$SUMMARY"
  printf "%s\tSNU-_K4+\t%d\t%s\n" "$STHR" "$K4pos_SNUneg"  "$(pct "$K4pos_SNUneg")"   >> "$SUMMARY"
  printf "%s\tSNU+_K4-\t%d\t%s\n" "$STHR" "$K4neg_SNUpos"  "$(pct "$K4neg_SNUpos")"   >> "$SUMMARY"
  printf "%s\tSNU-_K4-\t%d\t%s\n" "$STHR" "$Neither"       "$(pct "$Neither")"        >> "$SUMMARY"
done

echo "✅ Wrote summary → $SUMMARY"
echo "✅ Quantiles → $QFILE"

# ─────────── 3) Histogram of SNU-seq ───────────
Rscript --vanilla - "$TAB" "$HIST_SNU" "$(IFS=,; echo "${SNU_THRESH_LIST[*]}")" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
tab    <- args[1]
outS   <- args[2]
thrCsv <- args[3]
ths <- as.numeric(strsplit(thrCsv, ",")[[1]])

df <- read.table(tab, header=TRUE, sep="\t", comment.char="")
s <- df[,5]; s[!is.finite(s)] <- NA; s <- s[!is.na(s)]

suppressPackageStartupMessages(library(ggplot2))

p <- ggplot(data.frame(val=s), aes(x=val)) +
  geom_histogram(bins=100) +
  geom_density(linewidth=0.7, alpha=0.2) +
  labs(title="SNU-seq on ATAC peaks",
       subtitle=paste("Thresholds:", paste(ths, collapse=", ")),
       x="mean log2(count+1)", y="Count") +
  theme_minimal(base_size=12)

for (t in ths) p <- p + geom_vline(xintercept=t, linetype="dashed")

ggsave(outS, p, width=7, height=4.5)
message("Histogram written to: ", outS)
EOF

echo "✅ Histogram → $HIST_SNU"