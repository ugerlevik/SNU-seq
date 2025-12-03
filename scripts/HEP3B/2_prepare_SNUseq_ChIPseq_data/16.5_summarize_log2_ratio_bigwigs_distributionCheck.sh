#!/bin/bash

############################################################################
## Project: SNUseq project
## Purpose: Summarize log2(K27ac/K4me3) ratio BigWigs into a combined TSV
## Author: Umut Gerlevik
## Date: July 10, 2025
############################################################################

source ~/.zshrc
mamba activate deeptools_env

# Input & output setup
ratio_dir="/MellorLab/SNUseqProject/4_Anna/101_input_files/log2_K27ac_to_K4me3_ratios"
genome_file="/MellorLab/SNUseqProject/0_commonFiles/genome/2_STARgenome/chrNameLength.txt"
tmp_dir="${ratio_dir}/tmp_stats"
output_tsv="${ratio_dir}/log2_K27ac_to_K4me3_ratio_stats.tsv"

mkdir -p "$tmp_dir"
echo -e "File\tMin\tMax\tMean\tMedian" > "$output_tsv"

echo "📊 Computing distribution stats and saving to TSV..."

# Loop through BigWigs
for bw in "$ratio_dir"/*.bw; do
  [[ -e "$bw" ]] || continue
  base=$(basename "$bw" .bw)
  bg="$tmp_dir/${base}.bedgraph"

  echo "🔁 $base"

  # Convert BigWig to BedGraph
  bigWigToBedGraph "$bw" "$bg"

  # Compute stats in R and append to TSV
  Rscript --vanilla - "$bg" "$base" "$output_tsv" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
bgfile <- args[1]
basename <- args[2]
outfile <- args[3]

vals <- read.table(bgfile)[,4]
vals <- vals[is.finite(vals)]

if (length(vals) > 0) {
  stats <- c(min(vals), max(vals), mean(vals), median(vals))
  write.table(
    data.frame(File=basename, t(stats)),
    file=outfile,
    append=TRUE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE
  )
}
EOF

done

rm -rf "$tmp_dir"

echo "✅ Done. Summary saved to:"
echo "$output_tsv"
