############################################################################
## Project: SNUseq project
## Script purpose: Rescale data in deepTools computeMatrix output matrix file
##                to adjust specific samples (e.g. PRO-seq, TT-seq) while
##                keeping others unchanged, for better visualisation in
##                metagene heatmaps/plots
## Date: Sep 27, 2025
## Author: Umut Gerlevik
############################################################################

library(data.table)
library(jsonlite)
library(R.utils)   # for gzip()

# === INPUTS ===
matrix_file <- "/MellorLab/SNUseqProject/1_Umut/14.2_matrices_bPAPminusnoPAP_negRemoved_withTTseqPROseqNETseq/matrix_fwd_rev_combined_TSS_to_TES.gz"
out_file <- sub(".gz$", "_scaled", matrix_file)

# === READ HEADER ===
h <- scan(matrix_file, n = 1, sep = "\n", what = character(), quiet = TRUE)
params <- fromJSON(txt = gsub("^@", "", h))   # strip leading "@"

# Define your scaling factors here (names must match sample_labels)
scale_factors <- setNames(rep(1, length(params$sample_labels)), params$sample_labels)

scale_factors["\"GSM4730174_PROseq\""] <- 0.025
scale_factors["\"GSM4730176_TTseq\""]  <- 0.038
scale_factors["\"GSM7990390_mNETseq\""] <- 6
scale_factors["\"GSM5452295_HEK_FII_5\""] <- 0.02
scale_factors["\"GSM5452295_HEK_FII_3\""] <- 0.02

# === READ MATRIX ===
mat <- fread(matrix_file, skip = 1, header = FALSE)

# How many samples and bins?
n_samples <- length(params$sample_labels)
n_bins <- (ncol(mat) - 6) / n_samples

# Label columns: first 6 metadata cols + per-sample bins
bin_cols <- unlist(lapply(params$sample_labels, function(s) paste0(s, "_bin", seq_len(n_bins))))
colnames(mat) <- c("chr","start","end","feature","dot","strand", bin_cols)

# === APPLY SCALING ===
for (s in params$sample_labels) {
  factor <- scale_factors[[s]]
  if (is.null(factor)) factor <- 1   # default: no scaling
  sel <- grep(paste0("^", s, "_bin"), colnames(mat))
  mat[, (sel) := lapply(.SD, function(x) x * factor), .SDcols = sel]
}

# === WRITE OUTPUT ===
# 1. Write back the original JSON header line exactly
writeLines(h, out_file)

# 2. Append scaled matrix
fwrite(mat, file = out_file, sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE, row.names = FALSE)

# 3. Compress again
gzip(out_file, overwrite = TRUE)
