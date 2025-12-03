############################################################################
## Project: SNUseq project
## Script purpose: Count distributions and percentages in TSS, GeneBody, 
##                GeneEnd, Readthrough regions of the features
## Date: Mar 12, 2025
## Author: Umut Gerlevik
############################################################################

# Load necessary libraries
library(dplyr)
library(readr)
library(reshape2)

# Define input and output directories
input_dir <- "/MellorLab/SNUseqProject/1_Umut/10.1_countsRegions_fromBedgraph"
output_dir <- "/MellorLab/SNUseqProject/outputs/10.1_countsRegions_summary"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load the count matrix
counts_file <- file.path(input_dir, "final_merged_counts_ordered_normalized.csv")
final_counts <- read.csv(counts_file)

# Extract condition names and regions (excluding "finalAnnotationSubset")
columns <- colnames(final_counts)[-1]  # Exclude GeneID
conditions <- unique(gsub("_(finalAnnotationSubset|TSS|GeneBody|GeneEnd|Readthrough)", "", columns))
regions <- c("TSS", "GeneBody", "GeneEnd", "Readthrough")

# Initialize summary dataframe
summary_data <- data.frame()

# Calculate total counts and percentage distribution for each condition and region
for (condition in conditions) {
  total_counts <- sum(final_counts[[paste0(condition, "_TSS")]], na.rm = TRUE) +
    sum(final_counts[[paste0(condition, "_GeneBody")]], na.rm = TRUE) +
    sum(final_counts[[paste0(condition, "_GeneEnd")]], na.rm = TRUE) +
    sum(final_counts[[paste0(condition, "_Readthrough")]], na.rm = TRUE)
  
  for (region in regions) {
    region_counts <- sum(final_counts[[paste0(condition, "_", region)]], na.rm = TRUE)
    percentage <- ifelse(total_counts > 0, (region_counts / total_counts) * 100, 0)
    
    summary_data <- rbind(summary_data, data.frame(
      Condition = condition,
      Region = region,
      Total_Counts = region_counts,
      Percentage_of_Total = percentage
    ))
  }
}

# Reshape the table
tmp1 <- dcast(summary_data, Region ~ Condition, value.var = "Total_Counts")
colnames(tmp1) <- c("Region", paste0(unique(summary_data$Condition), "_TotalCounts"))
tmp2 <- dcast(summary_data, Region ~ Condition, value.var = "Percentage_of_Total")
colnames(tmp2) <- c("Region", paste0(unique(summary_data$Condition), "_Percentage"))

summary_data <- merge(tmp1, tmp2)

# Define the desired row order
row_order <- c("TSS", "GeneBody", "GeneEnd", "Readthrough")

# Reorder rows based on the defined order
summary_data <- summary_data[match(row_order, summary_data$Region), ]

# Extract column names
count_cols <- grep("_TotalCounts$", colnames(summary_data), value = TRUE)
percent_cols <- grep("_Percentage$", colnames(summary_data), value = TRUE)

# Ensure counts and percentages are in the same condition order
conditions <- gsub("_TotalCounts", "", count_cols)

# Construct new column order
new_col_order <- c("Region", unlist(sapply(conditions, function(cond) {
  c(paste0(cond, "_TotalCounts"), paste0(cond, "_Percentage"))
})))

# Reorder columns
summary_data <- summary_data[, new_col_order]

# Save the summary table
summary_file <- file.path(output_dir, "region_counts_summary.csv")
write.csv(summary_data, summary_file, row.names = FALSE)

# Print completion message
print("[INFO] Summary table saved to region_counts_summary.csv")
