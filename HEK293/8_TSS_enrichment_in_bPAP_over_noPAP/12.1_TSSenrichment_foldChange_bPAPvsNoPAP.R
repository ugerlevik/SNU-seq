############################################################################
## Project: SNUseq project
## Script purpose: Quantify the fold enrichment of labelled_bPAP vs 
##                labelled_noPAP in the TSS region
## Date: Nov 25, 2025
## Author: Umut Gerlevik
############################################################################

library(dplyr)
library(readr)
library(ggplot2)

# 1. Setup Directories
# --------------------
input_dir  <- "/MellorLab/SNUseqProject/1_Umut/10.1_countsRegions_fromBedgraph"
output_dir <- "/MellorLab/SNUseqProject/outputs/10.3_TSS_enrichment_bPAP"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. Load Data
# --------------------
counts_file <- file.path(input_dir, "final_merged_counts_ordered_normalized.csv")
if (!file.exists(counts_file)) {
  stop(paste("[ERROR] Input file not found:", counts_file))
}

df <- read.csv(counts_file)
print(paste("[INFO] Loaded data with", nrow(df), "genes."))

# 3. Calculate Enrichment
# --------------------
# Target columns: "labelled_bPAP_TSS" vs "labelled_noPAP_TSS"

# Define a small pseudocount to handle division by zero or log(0)
# Since data is normalized by length, values can be small. 
# We use a small epsilon relative to the data scale.
pseudocount <- 0.001 

enrichment_df <- df %>%
  dplyr::select(GeneID, labelled_bPAP_TSS, labelled_noPAP_TSS) %>%
  # Filter out genes with absolutely 0 signal in BOTH conditions to reduce noise
  filter(labelled_bPAP_TSS > 0 | labelled_noPAP_TSS > 0) %>%
  mutate(
    # Calculate simple Ratio
    Ratio_bPAP_vs_noPAP = (labelled_bPAP_TSS + pseudocount) / (labelled_noPAP_TSS + pseudocount),
    
    # Calculate Log2 Fold Change (Standard metric for enrichment)
    # Log2FC > 0 indicates enrichment in bPAP
    Log2FC_bPAP_vs_noPAP = log2(Ratio_bPAP_vs_noPAP)
  ) %>%
  arrange(desc(Log2FC_bPAP_vs_noPAP))

# 4. Global Statistics
# --------------------
mean_fc   <- mean(enrichment_df$Log2FC_bPAP_vs_noPAP, na.rm = TRUE)
median_fc <- median(enrichment_df$Log2FC_bPAP_vs_noPAP, na.rm = TRUE)
up_genes  <- sum(enrichment_df$Log2FC_bPAP_vs_noPAP > 1) # > 2-fold enrichment
down_genes <- sum(enrichment_df$Log2FC_bPAP_vs_noPAP < -1)

summary_text <- paste0(
  "Global Enrichment Statistics (TSS Region - bPAP vs noPAP):\n",
  "--------------------------------------------------------\n",
  "Total Genes Analyzed: ", nrow(enrichment_df), "\n",
  "Mean Log2FC:          ", round(mean_fc, 3), "\n",
  "Median Log2FC:        ", round(median_fc, 3), "\n",
  "Genes Enriched (>2x): ", up_genes, "\n",
  "Genes Depleted (<0.5x): ", down_genes, "\n"
)

cat(summary_text)

# 5. Save Results
# --------------------
# Save the gene-level table
write.csv(enrichment_df, file.path(output_dir, "TSS_enrichment_bPAP_vs_noPAP.csv"), row.names = FALSE)

# Save the summary stats
writeLines(summary_text, file.path(output_dir, "global_statistics_summary.txt"))

# 6. Visualization (Histogram of Enrichment)
# --------------------
p <- ggplot(enrichment_df, aes(x = Log2FC_bPAP_vs_noPAP)) +
  geom_histogram(binwidth = 0.2, fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
  geom_vline(aes(xintercept = median_fc), color="red", linetype="dashed", size=1) +
  theme_minimal() +
  labs(
    title = "Distribution of TSS Signal Enrichment (labelled bPAP vs noPAP)",
    subtitle = paste("Median Log2FC =", round(median_fc, 3)),
    x = "Log2 Fold Change (bPAP / noPAP)",
    y = "Count of Genes"
  )

ggsave(file.path(output_dir, "plot_TSS_enrichment_distribution.pdf"), plot = p, width = 8, height = 6)

print(paste("[INFO] Analysis complete. Results saved to:", output_dir))