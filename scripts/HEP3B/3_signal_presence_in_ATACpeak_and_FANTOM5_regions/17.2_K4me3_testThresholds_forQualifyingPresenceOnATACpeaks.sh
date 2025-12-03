#!/usr/bin/env bash

############################################################################
## Project: SNUseq project
## Script:  Test K4m3 thresholds for qualifying presence on ATAC peaks
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

# Fixed thresholds
K27_THR=1.25   # chosen cutoff for K27ac
SNU_THR=0      # SNU-seq already filtered

# Sweep these for K4me3
K4_THRESH_LIST=("0.5" "1.0" "1.25" "1.5" "2.0")

# ─────────── Outputs ───────────
OUT_DIR="/MellorLab/SNUseqProject/4_Anna/111_ATAC_K4me3_threshold_sweep"
mkdir -p "$OUT_DIR"

TAB="${OUT_DIR}/ATAC_K27ac_SNUseq_K4me3_means.tab"
NPZ="${OUT_DIR}/ATAC_K27ac_SNUseq_K4me3_means.npz"
SUMMARY="${OUT_DIR}/summary_K4me3_threshold_sweep.tsv"
HIST4="${OUT_DIR}/K4me3_on_ATAC_histogram.pdf"
QFILE="${OUT_DIR}/quantiles_K4me3.txt"

# ─────────── 1) Mean signal per peak ───────────
echo "▶ multiBigwigSummary on ATAC peaks (K27ac + SNU-seq + K4me3)…"
multiBigwigSummary BED-file \
  --BED "$ATAC_BED" \
  -b "$K27AC_BW" "$SNUSEQ_BW" "$K4ME3_BW" \
  -out "$NPZ" \
  --outRawCounts "$TAB" \
  -p $threads

TOTAL=$(wc -l < "$ATAC_BED")
echo "Number of ATAC peaks: $TOTAL"

# ─────────── 1.5) Quick quantiles for K4me3 ───────────
awk 'NR>1 && $6!="nan" && $6!="NaN" {print $6}' "$TAB" | \
Rscript -e 'x<-scan("stdin",quiet=TRUE); x<-x[is.finite(x)]; print(quantile(x,c(.5,.8,.9,.95,.99),na.rm=TRUE))' | tee "$QFILE"

# ─────────── 2) Sweep K4me3 thresholds ───────────
echo -e "k4me3_threshold\tclass\tcount\tpercent" > "$SUMMARY"

pct() { awk -v c="$1" -v t="$TOTAL" 'BEGIN{printf("%.2f", (t?100*c/t:0))}'; }

for K4THR in "${K4_THRESH_LIST[@]}"; do
  echo "▶ K27ac ≥ $K27_THR ; SNU > $SNU_THR ; sweeping K4me3 ≥ $K4THR"

  # Classes: (a) K27ac present is fixed, so we look at SNU vs K4 combos
  BOTH=$(awk -v k27="$K27_THR" -v sthr="$SNU_THR" -v k4thr="$K4THR" '
    NR==1{next}
    {ka=($4=="nan"?0:$4); su=($5=="nan"?0:$5); k4=($6=="nan"?0:$6)}
    ka>=k27 && su>sthr && k4>=k4thr {c++} END{print c+0}' "$TAB")

  SNUpos_K4neg=$(awk -v k27="$K27_THR" -v sthr="$SNU_THR" -v k4thr="$K4THR" '
    NR==1{next}
    {ka=($4=="nan"?0:$4); su=($5=="nan"?0:$5); k4=($6=="nan"?0:$6)}
    ka>=k27 && su>sthr && k4<k4thr {c++} END{print c+0}' "$TAB")

  SNUneg_K4pos=$(awk -v k27="$K27_THR" -v sthr="$SNU_THR" -v k4thr="$K4THR" '
    NR==1{next}
    {ka=($4=="nan"?0:$4); su=($5=="nan"?0:$5); k4=($6=="nan"?0:$6)}
    ka>=k27 && su<=sthr && k4>=k4thr {c++} END{print c+0}' "$TAB")

  Neither=$(awk -v k27="$K27_THR" -v sthr="$SNU_THR" -v k4thr="$K4THR" '
    NR==1{next}
    {ka=($4=="nan"?0:$4); su=($5=="nan"?0:$5); k4=($6=="nan"?0:$6)}
    ka>=k27 && su<=sthr && k4<k4thr {c++} END{print c+0}' "$TAB")

  printf "%s\tBoth_SNU+K4+\t%d\t%s\n" "$K4THR" "$BOTH"        "$(pct "$BOTH")"        >> "$SUMMARY"
  printf "%s\tSNU+_K4-\t%d\t%s\n"      "$K4THR" "$SNUpos_K4neg" "$(pct "$SNUpos_K4neg")" >> "$SUMMARY"
  printf "%s\tSNU-_K4+\t%d\t%s\n"      "$K4THR" "$SNUneg_K4pos" "$(pct "$SNUneg_K4pos")" >> "$SUMMARY"
  printf "%s\tNeither\t%d\t%s\n"       "$K4THR" "$Neither"      "$(pct "$Neither")"      >> "$SUMMARY"
done

echo "✅ Wrote summary → $SUMMARY"
echo "✅ Quantiles → $QFILE"

# ─────────── 3) Histogram of K4me3 ───────────
Rscript --vanilla - "$TAB" "$HIST4" "$(IFS=,; echo "${K4_THRESH_LIST[*]}")" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
tab    <- args[1]
out4   <- args[2]
thrCsv <- args[3]
ths <- as.numeric(strsplit(thrCsv, ",")[[1]])

df <- read.table(tab, header=TRUE, sep="\t", comment.char="")
k4 <- df[,6]; k4[!is.finite(k4)] <- NA; k4 <- k4[!is.na(k4)]

suppressPackageStartupMessages(library(ggplot2))

p <- ggplot(data.frame(val=k4), aes(x=val)) +
  geom_histogram(bins=100) +
  geom_density(linewidth=0.7, alpha=0.2) +
  labs(title="K4me3 on ATAC peaks",
       subtitle=paste("Thresholds:", paste(ths, collapse=", ")),
       x="mean log2(count+1)", y="Count") +
  theme_minimal(base_size=12)

for (t in ths) p <- p + geom_vline(xintercept=t, linetype="dashed")

ggsave(out4, p, width=7, height=4.5)
message("Histogram written to: ", out4)
EOF

echo "✅ Histogram → $HIST4"
