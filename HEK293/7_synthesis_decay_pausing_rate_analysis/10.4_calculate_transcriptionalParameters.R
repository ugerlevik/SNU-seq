############################################################################
## Project: SNUseq project
## Script purpose: Calculate synthesis, decay, and pausing index
## Date: Mar 22, 2025
## Author: Umut Gerlevik
############################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)

# Set working directories
output_dir <- "/MellorLab/SNUseqProject/1_Umut/10.1_countsRegions_fromBedgraph"

# Load the merged counts file
merged_file <- file.path(output_dir, "final_merged_counts_ordered_normalized.csv")
data <- read.csv(merged_file)

# Set labeling duration
labeling_time <- 10  # 10 minutes

# Identify relevant conditions
labelled_conditions <- grep("^labelled_bPAP", colnames(data), value = TRUE)
total_conditions <- grep("^total_noPAP", colnames(data), value = TRUE)

# Ensure correct matches (rRNA vs non-rRNA)
labelled_geneBody_cols <- grep("^labelled_bPAP.*_GeneBody$", labelled_conditions, value = TRUE)
total_whole_cols <- grep("^total_noPAP.*finalAnnotationSubset$", total_conditions, value = TRUE)

# Compute Decay & Synthesis Rates for bPAP
for (labelled_col in labelled_geneBody_cols) {
  if (grepl("_rRNA_", labelled_col)) {
    total_col <- total_whole_cols[grepl("_rRNA_", total_whole_cols)]
  } else {
    total_col <- total_whole_cols[!grepl("_rRNA_", total_whole_cols)]
  }
  
  if (length(total_col) == 1) {
    gene_prefix <- gsub("_GeneBody", "", labelled_col)
    
    # Compute decay rate
    decay_col <- paste0(gene_prefix, "_DecayRate")
    data[[decay_col]] <- ifelse(
      data[[total_col]] > 0 & data[[labelled_col]] > 0, 
      (-1 / labeling_time) * log10(1 - (data[[labelled_col]] / data[[total_col]])),
      NA
    )
    
    # Compute synthesis rate
    synthesis_col <- paste0(gene_prefix, "_SynthesisRate")
    data[[synthesis_col]] <- data[[decay_col]] * data[[total_col]]
  }
}

summary(is.nan(data$labelled_bPAP_DecayRate) | is.na(data$labelled_bPAP_DecayRate))
summary(is.nan(data$labelled_bPAP_rRNA_DecayRate) | is.na(data$labelled_bPAP_rRNA_DecayRate))

synth_bPAP <- data %>% 
  filter(!(is.nan(data$labelled_bPAP_DecayRate) | is.na(data$labelled_bPAP_DecayRate))) %>% 
  select(c(GeneID, colnames(data)[grepl("^labelled_bPAP", colnames(data)) & !grepl("rRNA", colnames(data))]))

synth_bPAP_rRNA <- data %>% 
  filter(!(is.nan(data$labelled_bPAP_rRNA_DecayRate) | is.na(data$labelled_bPAP_rRNA_DecayRate))) %>% 
  select(c(GeneID, colnames(data)[grepl("^labelled_bPAP_rRNA", colnames(data))]))


# Compute Pausing Index (TSS/GeneBody) for bPAP and noPAP
for (cond in gsub("_GeneBody", "", grep("^labelled_bPAP.*_GeneBody$", labelled_conditions, value = TRUE))) {
  tss_col <- paste0(cond, "_TSS")
  geneBody_col <- paste0(cond, "_GeneBody")
  
  pausing_col <- paste0(cond, "_PausingIndex")
  data[[pausing_col]] <- ifelse(
    data[[geneBody_col]] > 0, 
    data[[tss_col]] / data[[geneBody_col]], 
    NA
  )
}

for (cond in gsub("_GeneBody", "", grep("^labelled_noPAP.*_GeneBody$", colnames(data), value = TRUE))) {
  tss_col <- paste0(cond, "_TSS")
  geneBody_col <- paste0(cond, "_GeneBody")
  
  pausing_col <- paste0(cond, "_PausingIndex")
  data[[pausing_col]] <- ifelse(
    data[[geneBody_col]] > 0, 
    data[[tss_col]] / data[[geneBody_col]], 
    NA
  )
}

summary(is.nan(data$labelled_bPAP_PausingIndex) | is.na(data$labelled_bPAP_PausingIndex))
summary(is.nan(data$labelled_bPAP_rRNA_PausingIndex) | is.na(data$labelled_bPAP_rRNA_PausingIndex))
summary(is.nan(data$labelled_noPAP_PausingIndex) | is.na(data$labelled_noPAP_PausingIndex))
summary(is.nan(data$labelled_noPAP_rRNA_PausingIndex) | is.na(data$labelled_noPAP_rRNA_PausingIndex))

pausing_bPAP <- data %>% 
  filter(!(is.nan(data$labelled_bPAP_PausingIndex) | is.na(data$labelled_bPAP_PausingIndex))) %>% 
  select(c(GeneID, colnames(data)[grepl("^labelled_bPAP", colnames(data)) & !grepl("rRNA", colnames(data))]))

pausing_bPAP_rRNA <- data %>% 
  filter(!(is.nan(data$labelled_bPAP_rRNA_PausingIndex) | is.na(data$labelled_bPAP_rRNA_PausingIndex))) %>% 
  select(c(GeneID, colnames(data)[grepl("^labelled_bPAP_rRNA", colnames(data))]))

pausing_noPAP <- data %>% 
  filter(!(is.nan(data$labelled_noPAP_PausingIndex) | is.na(data$labelled_noPAP_PausingIndex))) %>% 
  select(c(GeneID, colnames(data)[grepl("^labelled_noPAP", colnames(data)) & !grepl("rRNA", colnames(data))]))

pausing_noPAP_rRNA <- data %>% 
  filter(!(is.nan(data$labelled_noPAP_rRNA_PausingIndex) | is.na(data$labelled_noPAP_rRNA_PausingIndex))) %>% 
  select(c(GeneID, colnames(data)[grepl("^labelled_noPAP_rRNA", colnames(data))]))


# Compute Termination Efficiency (GeneEnd / Readthrough) for labelled bPAP, bPAP_rRNA, noPAP, noPAP_rRNA
for (cond in c("labelled_bPAP", "labelled_bPAP_rRNA", "labelled_noPAP", "labelled_noPAP_rRNA")) {
  geneEnd_col <- paste0(cond, "_GeneEnd")
  readthrough_col <- paste0(cond, "_Readthrough")
  term_eff_col <- paste0(cond, "_TerminationEfficiency")
  
  data[[term_eff_col]] <- ifelse(
    data[[readthrough_col]] > 0,
    data[[geneEnd_col]] / data[[readthrough_col]],
    NA
  )
}


# Save the final computed table
output_file <- file.path(output_dir, "final_merged_transcriptional_rates.csv")
write.csv(data, output_file, row.names = FALSE)

print("[INFO] Transcriptional rates (synthesis, decay, pausing index, termination efficiency) computed and saved as final_merged_transcriptional_rates.csv")