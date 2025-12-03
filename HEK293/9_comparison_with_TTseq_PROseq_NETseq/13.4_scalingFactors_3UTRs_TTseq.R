############################################################################
## Project: SNUseq project
## Script purpose: Calculate size factors for 3'UTR normalisation (bPAP vs TTseq)
## Date: Apr 4, 2025
## Author: Umut Gerlevik
############################################################################

# Load necessary libraries
library(DESeq2)
library(dplyr)

# Set file paths
prefix <- "/MellorLab/SNUseqProject/2_TTseq_Phil/1.3_countsRegions_fromBedgraph"
counts_file <- paste0(prefix, "/final_merged_counts_ordered_normalized.csv")

# Load counts data
counts_data <- read.csv(counts_file, header = TRUE)

# Set rownames as the GeneID column and remove it from the data frame
rownames(counts_data) <- counts_data$GeneID
counts_data <- counts_data %>% select(-GeneID)

# DESeq2 requires integer counts
counts_data <- counts_data %>% mutate(across(everything(), ~ as.integer(round(.))))

# Define groups for scaling
samples <- c("labelled_bPAP", "TTseq_A", "TTseq_B")

# Extract only the relevant columns for normalization
counts <- counts_data %>% select(all_of(samples))

# Prepare metadata
colData <- data.frame(row.names = colnames(counts), condition = samples)

# Create DESeq2 datasets
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)

# Estimate size factors
dds <- estimateSizeFactors(dds)

# Extract scaling factors
sf <- sizeFactors(dds)

# Set bPAP scaling factor to 1 and adjust noPAP accordingly
sf["TTseq_A"] <- sf["TTseq_A"] / sf["labelled_bPAP"]
sf["TTseq_B"] <- sf["TTseq_B"] / sf["labelled_bPAP"]
sf["labelled_bPAP"] <- 1

# Combine scaling factors into a single dataframe
sf_combined <- data.frame(Sample = c(names(sf)),
                          ScalingFactor = c(sf))


# Merge all scaling factors
final_sf <- rbind(sf_combined)

# Print scaling factors
print(final_sf)

# Save scaling factors to a file
write.table(final_sf, file = paste0(prefix, "/DESeq2_scalingFactors_3UTRs.txt"), 
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

