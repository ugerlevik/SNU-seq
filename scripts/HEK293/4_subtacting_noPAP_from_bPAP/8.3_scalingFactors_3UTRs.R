############################################################################
## Project: SNUseq project
## Script purpose: Calculate size factors for 3'UTR normalisation (bPAP vs noPAP)
## Date: Mar 11, 2025
## Author: Umut Gerlevik
############################################################################

# Load necessary libraries
library(DESeq2)
library(dplyr)

# Set file paths
prefix <- "/MellorLab/SNUseqProject/1_Umut/8.1_countsRegions_fromBedgraph"
counts_file <- paste0(prefix, "/final_merged_counts_ordered_normalized.csv")

# Load counts data
counts_data <- read.csv(counts_file, header = TRUE)

# Set rownames as the GeneID column and remove it from the data frame
rownames(counts_data) <- counts_data$GeneID
counts_data <- counts_data %>% select(-GeneID)

# DESeq2 requires integer counts
counts_data <- counts_data %>% mutate(across(everything(), ~ as.integer(round(.))))

# Define groups for scaling
bpap_samples <- c("labelled_bPAP", "labelled_noPAP")
bpap_rRNA_samples <- c("labelled_bPAP_rRNA", "labelled_noPAP_rRNA")
total_samples <- c("total_bPAP", "total_noPAP")
total_rRNA_samples <- c("total_bPAP_rRNA", "total_noPAP_rRNA")

# Extract only the relevant columns for normalization
counts_bpap <- counts_data %>% select(all_of(bpap_samples))
counts_bpap_rRNA <- counts_data %>% select(all_of(bpap_rRNA_samples))
counts_total <- counts_data %>% select(all_of(total_samples))
counts_total_rRNA <- counts_data %>% select(all_of(total_rRNA_samples))

# Prepare metadata
colData_bpap <- data.frame(row.names = colnames(counts_bpap), condition = bpap_samples)
colData_bpap_rRNA <- data.frame(row.names = colnames(counts_bpap_rRNA), condition = bpap_rRNA_samples)
colData_total <- data.frame(row.names = colnames(counts_total), condition = total_samples)
colData_total_rRNA <- data.frame(row.names = colnames(counts_total_rRNA), condition = total_rRNA_samples)

# Create DESeq2 datasets
dds_bpap <- DESeqDataSetFromMatrix(countData = counts_bpap, colData = colData_bpap, design = ~ condition)
dds_bpap_rRNA <- DESeqDataSetFromMatrix(countData = counts_bpap_rRNA, colData = colData_bpap_rRNA, design = ~ condition)
dds_total <- DESeqDataSetFromMatrix(countData = counts_total, colData = colData_total, design = ~ condition)
dds_total_rRNA <- DESeqDataSetFromMatrix(countData = counts_total_rRNA, colData = colData_total_rRNA, design = ~ condition)

# Estimate size factors
dds_bpap <- estimateSizeFactors(dds_bpap)
dds_bpap_rRNA <- estimateSizeFactors(dds_bpap_rRNA)
dds_total <- estimateSizeFactors(dds_total)
dds_total_rRNA <- estimateSizeFactors(dds_total_rRNA)

# Extract scaling factors
sf_bpap <- sizeFactors(dds_bpap)
sf_bpap_rRNA <- sizeFactors(dds_bpap_rRNA)
sf_total <- sizeFactors(dds_total)
sf_total_rRNA <- sizeFactors(dds_total_rRNA)

# Set bPAP scaling factor to 1 and adjust noPAP accordingly
sf_bpap["labelled_noPAP"] <- sf_bpap["labelled_noPAP"] / sf_bpap["labelled_bPAP"]
sf_bpap["labelled_bPAP"] <- 1

sf_bpap_rRNA["labelled_noPAP_rRNA"] <- sf_bpap_rRNA["labelled_noPAP_rRNA"] / sf_bpap_rRNA["labelled_bPAP_rRNA"]
sf_bpap_rRNA["labelled_bPAP_rRNA"] <- 1

sf_total["total_noPAP"] <- sf_total["total_noPAP"] / sf_total["total_bPAP"]
sf_total["total_bPAP"] <- 1

sf_total_rRNA["total_noPAP_rRNA"] <- sf_total_rRNA["total_noPAP_rRNA"] / sf_total_rRNA["total_bPAP_rRNA"]
sf_total_rRNA["total_bPAP_rRNA"] <- 1

# Combine scaling factors into a single dataframe
sf_combined <- data.frame(Sample = c(names(sf_bpap), names(sf_bpap_rRNA), names(sf_total), names(sf_total_rRNA)),
                          ScalingFactor = c(sf_bpap, sf_bpap_rRNA, sf_total, sf_total_rRNA))

# Add unchanged conditions with a scaling factor of 1
unchanged_conditions <- setdiff(colnames(counts_data), c(bpap_samples, bpap_rRNA_samples, total_samples, total_rRNA_samples))
unchanged_sf <- data.frame(Sample = unchanged_conditions, ScalingFactor = 1)

# Merge all scaling factors
final_sf <- rbind(sf_combined, unchanged_sf)

# Print scaling factors
print(final_sf)

# Save scaling factors to a file
write.table(final_sf, file = paste0(prefix, "/DESeq2_scalingFactors_3UTRs.txt"), 
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

