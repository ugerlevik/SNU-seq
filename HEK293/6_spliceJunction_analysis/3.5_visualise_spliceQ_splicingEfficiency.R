############################################################################
## Project: SNUseq project
## Script purpose: Visualise spliceQ results
## Date: Mar 11, 2025
## Author: Umut Gerlevik
############################################################################

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# Set the directories
input_dir <- "/MellorLab/SNUseqProject/1_Umut/3.1_spliceQ_splicingEfficiency"
output_dir <- "/MellorLab/SNUseqProject/outputs/3.1_spliceQ_splicingEfficiency"
if(!dir.exists(output_dir)) dir.create(output_dir)

# List all Splice-Q output files
file_list <- list.files(input_dir, pattern = "*_spliceQ.tsv", full.names = TRUE)

# Define a function to load each sample's data and add a column for the sample name
load_spliceq_data <- function(file) {
  sample_name <- gsub("_spliceQ.tsv", "", basename(file))  # Extract the sample name
  data <- read.csv(file, sep = "\t", header = TRUE)  # Load the Splice-Q output file
  
  # Convert 'chr' column to character if it exists
  if ("chr" %in% colnames(data)) {
    data$chr <- as.character(data$chr)
  }
  
  data$Sample <- sample_name  # Add sample name as a column
  return(data)
}

# Load all files into a list of data frames
spliceq_data <- lapply(file_list, load_spliceq_data)


##############################################
## Data point-level comparison
##############################################

# Combine all data frames into one large data frame
combined_data <- bind_rows(spliceq_data)

# Select relevant columns for analysis (gene_ID, intron_ID, and SE which is the 'score' column)
combined_data <- combined_data %>%
  select(gene_ID, intron_ID, Sample, score) %>%
  rename(SE = score)  # Rename the score column to SE (Splicing Efficiency)

# Separate sample types
unique(combined_data$Sample)

# combined_data <- combined_data %>%
#   mutate(Sample_Type = case_when(
#     grepl("labelled", Sample) & grepl("bPAP", Sample) & !grepl("S4", Sample) & !grepl("H3|H4|S3", Sample) ~ "labelled_bPAP",
#     grepl("labelled", Sample) & grepl("bPAP", Sample) & !grepl("S4", Sample) & grepl("H3|H4|S3", Sample) ~ "labelled_bPAP_rRNA",
#     grepl("total", Sample) & grepl("bPAP", Sample) & !grepl("H3|H4|S3", Sample) ~ "total_bPAP",
#     grepl("total", Sample) & grepl("bPAP", Sample) & grepl("H3|H4|S3", Sample) ~ "total_bPAP_rRNA",
#     grepl("labelled", Sample) & grepl("noPAP", Sample) & !grepl("S4", Sample) & !grepl("H3|H4|S3", Sample) ~ "labelled_noPAP",
#     grepl("labelled", Sample) & grepl("noPAP", Sample) & !grepl("S4", Sample) & grepl("H3|H4|S3", Sample) ~ "labelled_noPAP_rRNA",
#     grepl("total", Sample) & grepl("noPAP", Sample) & !grepl("H3|H4|S3", Sample) ~ "total_noPAP",
#     grepl("total", Sample) & grepl("noPAP", Sample) & grepl("H3|H4|S3", Sample) ~ "total_noPAP_rRNA",
#     grepl("unlabelled", Sample) & grepl("bPAP", Sample) ~ "unlabelled_bPAP",
#     grepl("unlabelled", Sample) & grepl("noPAP", Sample) ~ "unlabelled_noPAP",
#     TRUE ~ "Other"  # Catch other types if necessary
#   ))
combined_data <- combined_data %>%
  mutate(Sample_Type = case_when(
    grepl("labelled", Sample) & grepl("bPAP", Sample) & !grepl("S4", Sample) & !grepl("H3|H4|S3", Sample) ~ "labelled_bPAP",
    grepl("labelled", Sample) & grepl("bPAP", Sample) & !grepl("S4", Sample) & grepl("H3|H4|S3", Sample) ~ "labelled_bPAP",
    grepl("total", Sample) & grepl("bPAP", Sample) & !grepl("H3|H4|S3", Sample) ~ "total_bPAP",
    grepl("total", Sample) & grepl("bPAP", Sample) & grepl("H3|H4|S3", Sample) ~ "total_bPAP",
    grepl("labelled", Sample) & grepl("noPAP", Sample) & !grepl("S4", Sample) & !grepl("H3|H4|S3", Sample) ~ "labelled_noPAP",
    grepl("labelled", Sample) & grepl("noPAP", Sample) & !grepl("S4", Sample) & grepl("H3|H4|S3", Sample) ~ "labelled_noPAP",
    grepl("total", Sample) & grepl("noPAP", Sample) & !grepl("H3|H4|S3", Sample) ~ "total_noPAP",
    grepl("total", Sample) & grepl("noPAP", Sample) & grepl("H3|H4|S3", Sample) ~ "total_noPAP",
    grepl("unlabelled", Sample) & grepl("bPAP", Sample) ~ "unlabelled_bPAP",
    grepl("unlabelled", Sample) & grepl("noPAP", Sample) ~ "unlabelled_noPAP",
    TRUE ~ "Other"  # Catch other types if necessary
  ))

summary(as.factor(combined_data$Sample_Type))
tmp_tbl <- as.data.frame(table(combined_data$Sample, combined_data$Sample_Type))
tmp_tbl <- tmp_tbl %>% filter(Freq > 0)
tmp_tbl

# Count the number of data points for each sample type
sample_counts <- combined_data %>%
  group_by(Sample_Type) %>%
  summarise(n = n())

# Extract the number of data points
n_labelled_bPAP <- sample_counts %>% filter(Sample_Type == "labelled_bPAP") %>% pull(n)
# n_labelled_bPAP_rRNA <- sample_counts %>% filter(Sample_Type == "labelled_bPAP_rRNA") %>% pull(n)
n_labelled_noPAP <- sample_counts %>% filter(Sample_Type == "labelled_noPAP") %>% pull(n)
# n_labelled_noPAP_rRNA <- sample_counts %>% filter(Sample_Type == "labelled_noPAP_rRNA") %>% pull(n)
n_total_bPAP <- sample_counts %>% filter(Sample_Type == "total_bPAP") %>% pull(n)
# n_total_bPAP_rRNA <- sample_counts %>% filter(Sample_Type == "total_bPAP_rRNA") %>% pull(n)
n_total_noPAP <- sample_counts %>% filter(Sample_Type == "total_noPAP") %>% pull(n)
# n_total_noPAP_rRNA <- sample_counts %>% filter(Sample_Type == "total_noPAP_rRNA") %>% pull(n)
n_unlabelled_bPAP <- sample_counts %>% filter(Sample_Type == "unlabelled_bPAP") %>% pull(n)
n_unlabelled_noPAP <- sample_counts %>% filter(Sample_Type == "unlabelled_noPAP") %>% pull(n)

# Create explicit labels with sample counts
labelled_bPAP_label <- paste0("Labelled + bPAP\n(n = ", n_labelled_bPAP, ")")
# labelled_bPAP_rRNA_label <- paste0("labelled_bPAP_rRNA\n(n = ", n_labelled_bPAP_rRNA, ")")
labelled_noPAP_label <- paste0("labelled_noPAP\n(n = ", n_labelled_noPAP, ")")
# labelled_noPAP_rRNA_label <- paste0("labelled_noPAP_rRNA\n(n = ", n_labelled_noPAP_rRNA, ")")
total_bPAP_label <- paste0("total_bPAP\n(n = ", n_total_bPAP, ")")
# total_bPAP_rRNA_label <- paste0("total_bPAP_rRNA\n(n = ", n_total_bPAP_rRNA, ")")
total_noPAP_label <- paste0("Total no bPAP\n(n = ", n_total_noPAP, ")")
# total_noPAP_rRNA_label <- paste0("total_noPAP_rRNA\n(n = ", n_total_noPAP_rRNA, ")")
unlabelled_bPAP_label <- paste0("unlabelled_bPAP\n(n = ", n_unlabelled_bPAP, ")")
unlabelled_noPAP_label <- paste0("unlabelled_noPAP\n(n = ", n_unlabelled_noPAP, ")")

# Re-assign the labels to the Sample_Type factor with the explicit labels
# combined_data$Sample_Type <- factor(combined_data$Sample_Type, 
#                                     levels = c("labelled_bPAP", "labelled_bPAP_rRNA", "labelled_noPAP", "labelled_noPAP_rRNA", 
#                                                "total_bPAP", "total_bPAP_rRNA", "total_noPAP", "total_noPAP_rRNA", 
#                                                "unlabelled_bPAP", "unlabelled_noPAP"), 
#                                     labels = c(labelled_bPAP_label, labelled_bPAP_rRNA_label, labelled_noPAP_label, labelled_noPAP_rRNA_label,
#                                                total_bPAP_label, total_bPAP_rRNA_label, total_noPAP_label, total_noPAP_rRNA_label,
#                                                unlabelled_bPAP_label, unlabelled_noPAP_label))
combined_data$Sample_Type <- factor(combined_data$Sample_Type,
                                    levels = c("labelled_bPAP", "labelled_noPAP",
                                               "total_bPAP", "total_noPAP",
                                               "unlabelled_bPAP", "unlabelled_noPAP"),
                                    labels = c(labelled_bPAP_label, labelled_noPAP_label,
                                               total_bPAP_label, total_noPAP_label,
                                               unlabelled_bPAP_label, unlabelled_noPAP_label))

# Filter the data for sample types
combined_data <- combined_data %>% 
  filter(Sample_Type %in% c(labelled_bPAP_label, total_noPAP_label))
combined_data$Sample_Type <- droplevels(combined_data$Sample_Type)

# Add statistical comparisons (e.g., Wilcoxon test) between groups using the new labels
stat_test <- compare_means(SE ~ Sample_Type, data = combined_data, 
                           method = "wilcox.test", paired = FALSE)

# Define color-blind friendly palettes for boxplots and jitter points using the new labels
boxplot_colors <- c(setNames(c("#56B4E9", "#E69F00"), c(labelled_bPAP_label, total_noPAP_label))) 

# Plot SE for all data points across SNU-seq vs RNA-seq with updated x-axis labels and color palette
ggplot(combined_data, aes(x = Sample_Type, y = SE, fill = Sample_Type)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  # Boxplot with no outliers
  geom_jitter(aes(color = Sample), width = 0.25, size = 1, alpha = 0.6) +  # Jitter points for data visibility
  labs(
    title = "Splicing Efficiency: SNU-seq vs RNA-seq",
    x = "Sample Type", y = "Splicing Efficiency"
  ) +
  theme_minimal() +
  scale_fill_manual(values = boxplot_colors) +  # Apply custom boxplot colors
  scale_color_manual(values = rainbow(length(unique(combined_data$Sample)))) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 36, face = "bold"),
    axis.text.y = element_text(size = 24),
    axis.text.x = element_text(size = 26, face = "bold", hjust = 0.5),
    legend.position = "none"  # Remove legend to simplify plot
  ) +
  stat_compare_means(
    method = "wilcox.test", label = "p.signif", size = 12, tip.length = 0.01,
    comparisons = list(
      c(labelled_bPAP_label, total_noPAP_label)
    ))

# Save the plot with the original filename
ggsave(paste0(output_dir, "/labelled_bPAP_vs_total_noPAP_splicing_efficiency_comparison_all_points.pdf"), width = 8, height = 11)

rm(list = setdiff(ls(), c("boxplot_colors", "output_dir", "spliceq_data")))


##############################################
## Library-level comparison
##############################################

# Combine all data frames into one large data frame
combined_data <- bind_rows(spliceq_data)

# Select relevant columns for analysis (gene_ID, intron_ID, and SE which is the 'score' column)
combined_data <- combined_data %>%
  select(gene_ID, intron_ID, Sample, score) %>%
  rename(SE = score)  # Rename the score column to SE (Splicing Efficiency)

# Separate sample types
unique(combined_data$Sample)

combined_data <- combined_data %>%
  mutate(Sample_Type = case_when(
    grepl("labelled", Sample) & grepl("bPAP", Sample) & !grepl("S4", Sample) & !grepl("H3|H4|S3", Sample) ~ "labelled_bPAP",
    grepl("labelled", Sample) & grepl("bPAP", Sample) & !grepl("S4", Sample) & grepl("H3|H4|S3", Sample) ~ "labelled_bPAP",
    grepl("total", Sample) & grepl("bPAP", Sample) & !grepl("H3|H4|S3", Sample) ~ "total_bPAP",
    grepl("total", Sample) & grepl("bPAP", Sample) & grepl("H3|H4|S3", Sample) ~ "total_bPAP",
    grepl("labelled", Sample) & grepl("noPAP", Sample) & !grepl("S4", Sample) & !grepl("H3|H4|S3", Sample) ~ "labelled_noPAP",
    grepl("labelled", Sample) & grepl("noPAP", Sample) & !grepl("S4", Sample) & grepl("H3|H4|S3", Sample) ~ "labelled_noPAP",
    grepl("total", Sample) & grepl("noPAP", Sample) & !grepl("H3|H4|S3", Sample) ~ "total_noPAP",
    grepl("total", Sample) & grepl("noPAP", Sample) & grepl("H3|H4|S3", Sample) ~ "total_noPAP",
    grepl("unlabelled", Sample) & grepl("bPAP", Sample) ~ "unlabelled_bPAP",
    grepl("unlabelled", Sample) & grepl("noPAP", Sample) ~ "unlabelled_noPAP",
    TRUE ~ "Other"  # Catch other types if necessary
  ))

summary(as.factor(combined_data$Sample_Type))
tmp_tbl <- as.data.frame(table(combined_data$Sample, combined_data$Sample_Type))
tmp_tbl <- tmp_tbl %>% filter(Freq > 0)
tmp_tbl

# Filter to include only SNU-seq and RNA-seq data points for visualization
summary_stats <- combined_data %>%
  filter(Sample_Type %in% c("labelled_bPAP", "total_noPAP")) %>%
  group_by(Sample_Type, Sample) %>%
  summarise(mean_SE = mean(SE, na.rm = TRUE))

# Count the number of samples for RNA-seq and SNU-seq
sample_counts <- summary_stats %>%
  group_by(Sample_Type) %>%
  summarise(n = n())

# Extract counts for RNA-seq and SNU-seq
n_total_noPAP <- sample_counts %>% filter(Sample_Type == "total_noPAP") %>% pull(n)
n_labelled_bPAP <- sample_counts %>% filter(Sample_Type == "labelled_bPAP") %>% pull(n)

# Create explicit labels with sample counts for RNA-seq and SNU-seq
total_noPAP_label <- paste0("Total no bPAP\n(n = ", n_total_noPAP, ")")
labelled_bPAP_label <- paste0("Labelled + bPAP\n(n = ", n_labelled_bPAP, ")")

# Update the factor levels to include the sample counts in the labels
summary_stats$Sample_Type <- factor(summary_stats$Sample_Type, 
                                    levels = c("labelled_bPAP", "total_noPAP"), 
                                    labels = c(labelled_bPAP_label, total_noPAP_label))

# Add statistical comparisons (e.g., Wilcoxon test) between groups using the new labels
stat_test <- compare_means(mean_SE ~ Sample_Type, data = summary_stats, 
                           method = "wilcox.test", paired = FALSE)

# Define color-blind friendly palettes for boxplots and jitter points using the new labels
jitter_colors <- c(setNames(c("#0072B2", "#D55E00"), c(labelled_bPAP_label, total_noPAP_label)))  # Dark Blue for SNU-seq, Dark Orange for RNA-seq
boxplot_colors <- c(setNames(c("#56B4E9", "#E69F00"), c(labelled_bPAP_label, total_noPAP_label)))  # Light Blue for SNU-seq, Light Orange for RNA-seq

# Plot SE across SNU-seq vs RNA-seq
ggplot(summary_stats, aes(x = Sample_Type, y = mean_SE, fill = Sample_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Remove outliers to avoid visual clutter
  geom_jitter(aes(color = Sample_Type), width = 0.2, size = 4, alpha = 0.8) +  # Add jitter for individual points
  labs(title = "Splicing Efficiency: SNU-seq vs RNA-seq",
       x = "Sample Type", y = "Mean Splicing Efficiency") +
  scale_y_continuous(limits = c(0.9, 1)) +
  theme_minimal() +
  scale_fill_manual(values = boxplot_colors) +  # Use custom boxplot colors
  scale_color_manual(values = jitter_colors) +  # Use custom jitter colors
  theme(
    plot.title = element_blank(),  # Remove the plot title
    axis.title.x = element_blank(),  # Remove the x-axis title
    axis.title.y = element_text(size = 28, face = "bold"),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 22, face = "bold"),
    legend.position = "none"  # Remove the legend to simplify the plot
  ) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", size = 12, tip.length = 0.01,
                     comparisons = list(c(labelled_bPAP_label, total_noPAP_label)))  # Add statistical comparisons

# Save the plot
ggsave(paste0(output_dir, "/labelled_bPAP_vs_total_noPAP_splicing_efficiency_comparison.pdf"), width = 6, height = 9)

