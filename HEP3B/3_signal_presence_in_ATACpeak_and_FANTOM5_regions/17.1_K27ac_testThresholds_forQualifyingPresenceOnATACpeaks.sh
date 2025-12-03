#!/usr/bin/env bash

############################################################################
## Project: SNUseq project
## Script:  Test K27ac thresholds for qualifying presence on ATAC peaks
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

# Thresholds to try for K27ac presence (in log2(count+1) units)
K27_THRESH_LIST=("0.5" "1" "1.25" "1.5" "2")

# SNU‑seq presence threshold (default 0 → >0 means ≥1 read)
SNU_THR="${SNU_THR:-0}"

# ─────────── Outputs ───────────
OUT_DIR="/MellorLab/SNUseqProject/4_Anna/110_ATAC_K27ac_vs_SNUseq_presence_threshSweep"
mkdir -p "$OUT_DIR"

TAB="${OUT_DIR}/ATAC_K27ac_SNUseq_means.tab"
NPZ="${OUT_DIR}/ATAC_K27ac_SNUseq_means.npz"
SUMMARY="${OUT_DIR}/summary_threshold_sweep.tsv"
HIST_PDF="${OUT_DIR}/K27ac_on_ATAC_histogram.pdf"

# ─────────── 1) Mean signal per peak for both tracks ───────────
# TAB columns: chr start end K27ac_mean SNUseq_mean   (header on line 1)
echo "▶ Running multiBigwigSummary on ATAC peaks (K27ac + SNU-seq)…"
multiBigwigSummary BED-file \
  --BED "$ATAC_BED" \
  -b "$K27AC_BW" "$SNUSEQ_BW" \
  -out "$NPZ" \
  --outRawCounts "$TAB" \
  -p $threads

awk 'NR>1{print $4}' "$TAB" | \
Rscript -e 'x<-scan("stdin",quiet=TRUE); print(quantile(x,c(.5,.8,.9,.95,.99)))'

TOTAL=$(wc -l < "$ATAC_BED")
echo "Number of ATAC peaks: $TOTAL"

# ─────────── 2) Sweep thresholds & classify ───────────
echo -e "k27ac_threshold\tclass\tcount\tpercent" > "$SUMMARY"

for THR in "${K27_THRESH_LIST[@]}"; do
  echo "▶ Threshold K27ac ≥ ${THR} (SNU present > ${SNU_THR})"

  K27_ONLY=$(awk -v thr="$THR" -v sthr="$SNU_THR" 'NR>1 && $4>=thr && $5<=sthr {c++} END{print c+0}' "$TAB")
  SNU_ONLY=$(awk -v thr="$THR" -v sthr="$SNU_THR" 'NR>1 && $4<thr  && $5>sthr  {c++} END{print c+0}' "$TAB")
  BOTH=$(     awk -v thr="$THR" -v sthr="$SNU_THR" 'NR>1 && $4>=thr && $5>sthr  {c++} END{print c+0}' "$TAB")
  NEITHER=$(  awk -v thr="$THR" -v sthr="$SNU_THR" 'NR>1 && $4<thr  && $5<=sthr {c++} END{print c+0}' "$TAB")

  pct() { awk -v c="$1" -v t="$TOTAL" 'BEGIN{printf("%.2f", (t?100*c/t:0))}'; }

  printf "%s\tK27ac_only\t%d\t%s\n" "$THR" "$K27_ONLY" "$(pct "$K27_ONLY")"   >> "$SUMMARY"
  printf "%s\tSNUseq_only\t%d\t%s\n" "$THR" "$SNU_ONLY" "$(pct "$SNU_ONLY")"   >> "$SUMMARY"
  printf "%s\tBoth\t%d\t%s\n"        "$THR" "$BOTH"      "$(pct "$BOTH")"      >> "$SUMMARY"
  printf "%s\tNeither\t%d\t%s\n"     "$THR" "$NEITHER"   "$(pct "$NEITHER")"   >> "$SUMMARY"

  # Optional: write BEDs per class per threshold (uncomment if you want files)
  # thr_tag="thr${THR}"
  # mkdir -p "${OUT_DIR}/beds_${thr_tag}"
  # awk -v thr="$THR" -v sthr="$SNU_THR" 'BEGIN{FS=OFS="\t"} NR==FNR{a[++i]=$0; next}
  #   NR>1 && $4>=thr && $5<=sthr {print a[NR-1]}' "$ATAC_BED" "$TAB" > "${OUT_DIR}/beds_${thr_tag}/K27ac_only.bed"
  # awk -v thr="$THR" -v sthr="$SNU_THR" 'BEGIN{FS=OFS="\t"} NR==FNR{a[++i]=$0; next}
  #   NR>1 && $4<thr && $5>sthr {print a[NR-1]}' "$ATAC_BED" "$TAB" > "${OUT_DIR}/beds_${thr_tag}/SNUseq_only.bed"
  # awk -v thr="$THR" -v sthr="$SNU_THR" 'BEGIN{FS=OFS="\t"} NR==FNR{a[++i]=$0; next}
  #   NR>1 && $4>=thr && $5>sthr {print a[NR-1]}' "$ATAC_BED" "$TAB" > "${OUT_DIR}/beds_${thr_tag}/Both.bed"
  # awk -v thr="$THR" -v sthr="$SNU_THR" 'BEGIN{FS=OFS="\t"} NR==FNR{a[++i]=$0; next}
  #   NR>1 && $4<thr && $5<=sthr {print a[NR-1]}' "$ATAC_BED" "$TAB" > "${OUT_DIR}/beds_${thr_tag}/Neither.bed"
done

echo "✅ Wrote summary with threshold sweep → $SUMMARY"

# ─────────── 3) Plot histogram of K27ac means with vertical lines at thresholds ───────────
Rscript --vanilla - "$TAB" "$HIST_PDF" "$(IFS=,; echo "${K27_THRESH_LIST[*]}")" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
tabfile <- args[1]
outpdf  <- args[2]
thr_csv <- args[3]

df <- read.table(tabfile, header=TRUE, sep="\t", comment.char="")
# Columns: chr start end K27ac SNUseq
k <- df[,4]

ths <- as.numeric(strsplit(thr_csv, ",")[[1]])

suppressPackageStartupMessages(library(ggplot2))

p <- ggplot(data.frame(k27 = k), aes(x = k27)) +
  geom_histogram(bins = 100) +
  geom_density(linewidth = 0.7, alpha = 0.2) +
  labs(
    title = "K27ac log2(count+1) mean on ATAC peaks",
    subtitle = if (length(ths) > 0) paste("Thresholds:", paste(ths, collapse = ", ")) else NULL,
    x = "K27ac mean (log2(count+1))",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

# Add vertical lines one by one (can't add a list directly)
if (length(ths) > 0) {
  for (t in ths) {
    p <- p + geom_vline(xintercept = t, linetype = "dashed")
  }
}

ggsave(outpdf, p, width = 7, height = 4.5)
message("Histogram written to: ", outpdf)
EOF

echo "✅ Histogram saved → $HIST_PDF"

echo -e "\nNext: open the histogram and compare the summary rows across thresholds.\nYou can set SNU_THR (default 0) like: SNU_THR=1 ./102.2_....sh if you also want a stricter SNU-seq presence cutoff."