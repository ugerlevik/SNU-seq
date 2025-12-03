############################################################################
## Project: SNUseq project
## Script purpose: Calculate size factors for spike-in normalisation
## Date: Mar 5, 2025
## Author: Umut Gerlevik
############################################################################

# Load necessary libraries
library(DESeq2)
library(dplyr)
# library(ggplot2)
# library(ggrepel)
# library(reshape2)

# Set file paths
prefix <- "/MellorLab/SNUseqProject/1_Umut/6.1_featureCounts_spikeins"
counts_file <- paste0(prefix, "/merged_spikeins_featureCounts.txt")

# Load counts data
counts_data <- read.delim(counts_file, header = TRUE)

# Set rownames as the Geneid column and remove it from the data frame
rownames(counts_data) <- counts_data$Geneid
counts_data <- counts_data %>% select(-Geneid)

# Prepare metadata for samples (you can adapt this part based on your experiment design)
sample_conditions <- c("labelled_bPAP", "labelled_noPAP", "labelled_bPAP", 
                       "labelled_noPAP", "labelled_bPAP_rRNA", "labelled_noPAP_rRNA", 
                       "labelled_bPAP_rRNA", "labelled_noPAP_rRNA", "labelled_bPAP", 
                       "labelled_noPAP", "labelled_bPAP", "labelled_noPAP", 
                       "labelled_bPAP_rRNA", "labelled_noPAP_rRNA", "total_bPAP", "total_noPAP", 
                       "total_bPAP", "total_noPAP", "total_bPAP_rRNA", "total_noPAP_rRNA", 
                       "total_bPAP_rRNA", "total_noPAP_rRNA", "total_bPAP", "total_noPAP", 
                       "total_bPAP", "total_noPAP", "total_bPAP_rRNA", "total_noPAP_rRNA", 
                       "total_bPAP", "total_noPAP", "unlabelled_bPAP", "unlabelled_noPAP")

# Create DESeq2 data object
colData <- data.frame(row.names = colnames(counts_data), condition = sample_conditions)
dds_spikein <- DESeqDataSetFromMatrix(countData = counts_data,
                                    colData = colData, 
                                    design = ~ condition)
dds_sf_spikein <- sizeFactors(estimateSizeFactors(dds_spikein))
print(dds_sf_spikein)


# Save scaling factors to a file
write.table(dds_sf_spikein, file = paste0(prefix, "/DESeq2_scalingFactors_spikeins.txt"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)


