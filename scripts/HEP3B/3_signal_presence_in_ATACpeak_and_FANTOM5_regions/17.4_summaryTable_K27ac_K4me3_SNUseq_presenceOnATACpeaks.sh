#!/usr/bin/env bash

############################################################################
## Project: SNUseq visualization
## Script:  Final three-way summary (K4me3, K27ac, SNU-seq) with thresholds
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

# Thresholds (fixed)
K27_THR=1.25
K4_THR=1.25
SNU_THR=0

# ─────────── Outputs ───────────
OUT_DIR="/MellorLab/SNUseqProject/4_Anna/112_ATAC_threeway_summary"
mkdir -p "$OUT_DIR"

TAB="${OUT_DIR}/ATAC_means.tab"
NPZ="${OUT_DIR}/ATAC_means.npz"
SUMMARY="${OUT_DIR}/threeway_summary.tsv"

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

# ─────────── 2) Classify into 8 groups ───────────
echo -e "K4me3\tK27ac\tSNU-seq\tcount\tpercent" > "$SUMMARY"

awk -v k27thr="$K27_THR" -v k4thr="$K4_THR" -v snuthr="$SNU_THR" -v total="$TOTAL" '
  BEGIN{ FS=OFS="\t" }
  NR==1{ next }
  {
    # Treat non-finite as 0
    k27 = ($4=="nan"||$4=="NaN"||$4=="inf"||$4=="-inf") ? 0 : $4+0
    snu = ($5=="nan"||$5=="NaN"||$5=="inf"||$5=="-inf") ? 0 : $5+0
    k4  = ($6=="nan"||$6=="NaN"||$6=="inf"||$6=="-inf") ? 0 : $6+0

    k27p = (k27 >= k27thr) ? "YES" : "NO"
    snup = (snu >  snuthr) ? "YES" : "NO"
    k4p  = (k4  >= k4thr ) ? "YES" : "NO"

    key = k4p OFS k27p OFS snup
    cnt[key]++
  }
  END{
    vals[1] = "YES"; vals[0] = "NO"
    for (i1=1; i1>=0; i1--)      # K4me3 YES/NO
      for (i2=1; i2>=0; i2--)    # K27ac YES/NO
        for (i3=1; i3>=0; i3--){ # SNU-seq YES/NO
          k = vals[i1] OFS vals[i2] OFS vals[i3]
          c = (k in cnt) ? cnt[k] : 0
          pct = (total ? 100*c/total : 0)
          printf "%s\t%d\t%.2f\n", k, c, pct
        }
  }' "$TAB" >> "$SUMMARY"

echo "✅ Wrote summary → $SUMMARY"

# ─────────── 3) Export per-class regions, sorted by SNU-seq ───────────
echo "▶ Writing per-class files sorted by SNU-seq…"
CLASS_DIR_TSV="${OUT_DIR}/classes_tsv"
CLASS_DIR_BED="${OUT_DIR}/classes_bed"
TMP_DIR="${OUT_DIR}/tmp_classes"
mkdir -p "$CLASS_DIR_TSV" "$CLASS_DIR_BED" "$TMP_DIR"

# We’ll create 8 temp files keyed by class tag, then sort each by SNU (col4)
# Class tags: K4pos_K27pos_SNUpos, K4pos_K27pos_SNUneg, ..., K4neg_K27neg_SNUneg

# 3.1 Emit per-class temp TSVs: chr  start  end  SNU  K27  K4
awk -v k27thr="$K27_THR" -v k4thr="$K4_THR" -v snuthr="$SNU_THR" \
    -v tmpdir="$TMP_DIR" 'BEGIN{FS=OFS="\t"}
  # First file: ATAC_BED → store the original BED lines by rank
  FNR==NR {
    bed[++n] = $0
    next
  }
  # Second file: ATAC_means.tab  (chr start end K27 SNU K4)
  NR==1 { next }  # skip header
  {
    # Row index matches bed[] (1-based): current record # among data lines
    idx = NR-1

    chr=$1; st=$2; en=$3
    k27 = ($4=="nan"||$4=="NaN"||$4=="inf"||$4=="-inf")?0:$4+0
    snu = ($5=="nan"||$5=="NaN"||$5=="inf"||$5=="-inf")?0:$5+0
    k4  = ($6=="nan"||$6=="NaN"||$6=="inf"||$6=="-inf")?0:$6+0

    k27p = (k27 >= k27thr) ? "pos" : "neg"
    snup = (snu >  snuthr) ? "pos" : "neg"
    k4p  = (k4  >= k4thr ) ? "pos" : "neg"

    tag = "K4" k4p "_K27" k27p "_SNU" snup

    # TSV line: chr  start  end  SNU  K27  K4  (keep numeric signals for sorting/inspection)
    print chr, st, en, snu, k27, k4 > (tmpdir "/" tag ".tsv")
  }
' "$ATAC_BED" "$TAB"

# 3.2 For each class: sort by SNU (col4) descending and write TSV + BED
for f in "${TMP_DIR}"/K4*_K27*_SNU*.tsv; do
  [ -s "$f" ] || continue
  base=$(basename "$f" .tsv)

  # (a) Sorted TSV with header
  {
    echo -e "chr\tstart\tend\tSNU\tK27ac\tK4me3"
    sort -k4,4nr "$f"
  } > "${CLASS_DIR_TSV}/${base}_sortedBySNU.tsv"

  # (b) BED4 for IGV; name field carries short annotation + SNU value
  #     IGV will load BED3/4 fine; score is optional. We keep it simple.
  awk 'BEGIN{FS=OFS="\t"}
       NR==1{next}
       {printf "%s\t%s\t%s\tSNU=%.3f;K27=%.3f;K4=%.3f\n", $1,$2,$3,$4,$5,$6}' \
      "${CLASS_DIR_TSV}/${base}_sortedBySNU.tsv" > "${CLASS_DIR_BED}/${base}_sortedBySNU.bed"

done

# 3.3 (Optional) Also write an “all peaks” file sorted by SNU
ALL_TSV="${OUT_DIR}/all_peaks_sortedBySNU.tsv"
ALL_BED="${OUT_DIR}/all_peaks_sortedBySNU.bed"
{
  echo -e "chr\tstart\tend\tSNU\tK27ac\tK4me3"
  awk 'NR>1{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6}' "$TAB" | sort -k4,4nr
} > "$ALL_TSV"
awk 'BEGIN{FS=OFS="\t"} NR==1{next}
     {printf "%s\t%s\t%s\tSNU=%.3f;K27=%.3f;K4=%.3f\n", $1,$2,$3,$4,$5,$6}' \
    "$ALL_TSV" > "$ALL_BED"

echo "✅ Per-class TSVs → $CLASS_DIR_TSV"
echo "✅ Per-class BEDs → $CLASS_DIR_BED"
echo "✅ All-peaks (SNU-sorted) → $ALL_TSV ; $ALL_BED"