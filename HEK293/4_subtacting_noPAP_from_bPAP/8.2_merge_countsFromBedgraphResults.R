############################################################################
## Project: SNUseq project
## Script purpose: Merge counts from 3' UTR regions (no need normalization by region length)
## Date: Mar 11, 2025
## Author: Umut Gerlevik
############################################################################

library(dplyr)
library(tidyr)
library(readr)

# Set working directory to merged files
merged_dir <- "/MellorLab/SNUseqProject/1_Umut/8.1_countsRegions_fromBedgraph/individualResults"
output_dir <- "/MellorLab/SNUseqProject/1_Umut/8.1_countsRegions_fromBedgraph"

# Get all merged txt files
merged_files <- list.files(merged_dir, pattern = "_merged_.*\\.bedgraph$", full.names = TRUE)

# Check file names
print(merged_files)

# Initialize an empty list to store data frames
data_list <- list()

# Read each file and store it in the list
for (file in merged_files) {
  # Extract condition name from filename
  condition <- gsub(".*?/|_merged_.*.bedgraph", "", file)  # Extracts condition name
  feature <- gsub(".*_merged_|.bedgraph", "", file)  # Extracts feature name
  
  # Read the file
  df <- read.delim(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Rename columns (assuming the format: chr, start, end, ..., value)
  count_col <- paste0(condition, "_", feature)
  colnames(df) <- c("Chr", "Start", "End", "GeneID", count_col)
  
  # Compute region length
  df <- df %>% mutate(RegionLength = as.numeric(End) - as.numeric(Start))
  
  # Convert score column to numeric
  df[[count_col]][df[[count_col]] == "."] <- 0
  df[[count_col]] <- as.numeric(df[[count_col]])
  
  # Normalize count by region length
  # df <- df %>% mutate(!!count_col := get(count_col) / RegionLength)
  
  # Keep only relevant columns
  df <- df %>% select(c("GeneID", count_col))
  
  # Store in list
  data_list[[paste0(condition, "_", feature)]] <- df
}

# Merge all data frames by "GeneID"
final_merged_table <- Reduce(function(x, y) full_join(x, y, by = c("GeneID")), data_list)

# Extract **full** condition names dynamically
conditions <- unique(gsub("_(finalAnnotationSubset|TSS|GeneBody|GeneEnd|Readthrough|filtered_3UTRs)", "", colnames(final_merged_table)[-1]))

# Define feature order
feature_order <- c("filtered_3UTRs")

# Generate ordered column names
ordered_columns <- c("GeneID", unlist(sapply(conditions, function(cond) {
  paste0(cond, "_", feature_order)
})))

# Reorder columns based on ordered_columns
final_merged_table <- final_merged_table %>% select(all_of(ordered_columns))

# Save the final merged table
write.csv(final_merged_table, file.path(output_dir, "/final_merged_counts_ordered_normalized.csv"), row.names = FALSE)

print("[INFO] All conditions merged, normalized, and ordered into final_merged_counts_ordered_normalized.csv")

