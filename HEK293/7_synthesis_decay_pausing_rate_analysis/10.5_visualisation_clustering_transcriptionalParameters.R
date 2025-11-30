############################################################################
## Project: SNUseq project
## Script purpose: Visualise counts, synthesis, decay, pausing
## Date: Mar 22, 2025
## Author: Umut Gerlevik
############################################################################

# Libraries
library(dplyr)
library(ggpubr)
library(ggplot2)

input_dir <- "/MellorLab/SNUseqProject/1_Umut/10.1_countsRegions_fromBedgraph"
output_dir <- "/MellorLab/SNUseqProject/outputs/10.2_transcriptionalModelling_plots"
if(!dir.exists(output_dir)) dir.create(output_dir)

all_datasets <- read.csv(file.path(input_dir, "/final_merged_transcriptional_rates.csv"))

############################
## PAUSE bPAP vs noPAP
############################

# Merge "labelled bPAP" with "labelled bPAP -rRNA" (mean of both)
# Merge "labelled noPAP" with "labelled noPAP -rRNA" (mean of both)
data <- data.frame(
  GeneID = all_datasets$GeneID,
  PauseIndex = rowMeans(cbind(all_datasets$labelled_bPAP_PausingIndex, all_datasets$labelled_bPAP_rRNA_PausingIndex), na.rm = TRUE),
  Condition = "labelled + bPAP"
)

data_noPAP <- data.frame(
  GeneID = all_datasets$GeneID,
  PauseIndex = rowMeans(cbind(all_datasets$labelled_noPAP_PausingIndex, all_datasets$labelled_noPAP_rRNA_PausingIndex), na.rm = TRUE),
  Condition = "labelled no bPAP"
)

# Combine both groups
data <- rbind(data, data_noPAP)

# Remove NA and NaN values
data <- data %>% filter(!is.na(PauseIndex) & !is.nan(PauseIndex))

# Remove outliers using 1.5 * IQR rule
Q1 <- quantile(data$PauseIndex, 0.25, na.rm = TRUE)
Q3 <- quantile(data$PauseIndex, 0.75, na.rm = TRUE)
IQR_value <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

data <- data %>% filter(PauseIndex >= lower_bound & PauseIndex <= upper_bound)

# Compute sample sizes (n)
sample_sizes <- data %>% group_by(Condition) %>% summarise(n = n())

# Update condition names with (n = ...) in labels
data$Condition <- recode_factor(data$Condition,
                                "labelled + bPAP" = paste0("labelled + bPAP\n(n = ", sample_sizes$n[sample_sizes$Condition == "labelled + bPAP"], ")"),
                                "labelled no bPAP" = paste0("labelled no bPAP\n(n = ", sample_sizes$n[sample_sizes$Condition == "labelled no bPAP"], ")")
)

# Define color palette manually **AFTER** renaming conditions
color_palette <- setNames(
  c("#4E79A7", "#E15759"),  # Define colors
  levels(data$Condition)  # Assign colors to condition labels
)

# Define pairwise comparisons **AFTER** renaming conditions
# comparisons <- list(c(levels(data$Condition)[1], levels(data$Condition)[2]))

# Plotting the violin + boxplot
ggplot(data, aes(x = Condition, y = PauseIndex, fill = Condition)) +
  # geom_violin(trim = TRUE, color = "black", alpha = 0.7, stat = T) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.25, alpha = 0.1, color = "darkgrey") +  # Add jittered points for data visibility
  scale_fill_manual(values = color_palette) +  # Use custom colors for each condition
  scale_y_continuous(limits = c(-0.1, 7), breaks = seq(0, 8, 1), minor_breaks = seq(0, 8, 0.5)) +
  labs(
    title = "Pausing Index Across Conditions (Merged)", 
    y = "Pausing Index", 
    x = NULL  # Remove x-axis label
  ) +  
  # stat_compare_means(comparisons = comparisons, label = "p.signif", method = "wilcox.test",
  #                    vjust = 0.65, size = 6, step.increase = 0.1, bracket.size = 1.0) +  # Add significance annotation
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.position = "none",  # Remove legend for boxplots
    plot.title = element_blank()
  )

# Save the plot
ggsave(filename = "pausing_index_boxplot_bPAPvsnoPAP_merged.pdf", plot = last_plot(), path = output_dir,
       width = 8, height = 5, units = "in")


rm(list = setdiff(ls(), c("all_datasets", "output_dir")))


############################
## TERMINATION EFFICIENCY: bPAP vs noPAP
############################

# Merge labelled_bPAP and labelled_bPAP_rRNA
term_data <- data.frame(
  GeneID = all_datasets$GeneID,
  TerminationEfficiency = rowMeans(cbind(
    all_datasets$labelled_bPAP_TerminationEfficiency,
    all_datasets$labelled_bPAP_rRNA_TerminationEfficiency
  ), na.rm = TRUE),
  Condition = "labelled + bPAP"
)

# Merge labelled_noPAP and labelled_noPAP_rRNA
term_data_noPAP <- data.frame(
  GeneID = all_datasets$GeneID,
  TerminationEfficiency = rowMeans(cbind(
    all_datasets$labelled_noPAP_TerminationEfficiency,
    all_datasets$labelled_noPAP_rRNA_TerminationEfficiency
  ), na.rm = TRUE),
  Condition = "labelled no bPAP"
)

# Combine both
term_data <- rbind(term_data, term_data_noPAP)

# Remove NA/NaN
term_data <- term_data %>% filter(!is.na(TerminationEfficiency) & !is.nan(TerminationEfficiency))

# Remove outliers using 1.5 * IQR
Q1 <- quantile(term_data$TerminationEfficiency, 0.25, na.rm = TRUE)
Q3 <- quantile(term_data$TerminationEfficiency, 0.75, na.rm = TRUE)
IQR_value <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value
term_data <- term_data %>% filter(TerminationEfficiency >= lower_bound & TerminationEfficiency <= upper_bound)

# Add (n =) to condition labels
term_sample_sizes <- term_data %>% group_by(Condition) %>% summarise(n = n())
term_data$Condition <- recode_factor(
  term_data$Condition,
  "labelled + bPAP" = paste0("labelled + bPAP\n(n = ", term_sample_sizes$n[term_sample_sizes$Condition == "labelled + bPAP"], ")"),
  "labelled no bPAP" = paste0("labelled no bPAP\n(n = ", term_sample_sizes$n[term_sample_sizes$Condition == "labelled no bPAP"], ")")
)

# Color palette
term_color_palette <- setNames(
  c("#4E79A7", "#E15759"),
  levels(term_data$Condition)
)

# Plot
ggplot(term_data, aes(x = Condition, y = TerminationEfficiency, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.25, alpha = 0.1, color = "darkgrey") +
  scale_fill_manual(values = term_color_palette) +
  labs(
    title = "Termination Efficiency Across Conditions (Merged)",
    y = "Termination Efficiency",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.position = "none",
    plot.title = element_blank()
  )

# Save
ggsave(
  filename = "termination_efficiency_boxplot_bPAPvsnoPAP_merged.pdf",
  plot = last_plot(),
  path = output_dir,
  width = 8,
  height = 5,
  units = "in"
)


##########################################
## Synthesis rate & counts - bPAP
##########################################
data <- data.frame(
  GeneID = all_datasets$GeneID,
  GeneBodyCounts = all_datasets$labelled_bPAP_GeneBody,
  SynthesisRate = all_datasets$labelled_bPAP_SynthesisRate,
  DecayRate = all_datasets$labelled_bPAP_DecayRate,
  Condition = "labelled + bPAP"
)

# Remove NA and NaN values
data <- data %>% filter(!is.na(SynthesisRate) & !is.nan(SynthesisRate))
data <- data %>% filter(!is.na(DecayRate) & !is.nan(DecayRate))

# Remove outliers using 1.5 * IQR rule
Q1 <- quantile(data$SynthesisRate, 0.25, na.rm = TRUE)
Q3 <- quantile(data$SynthesisRate, 0.75, na.rm = TRUE)
IQR_value <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value
data <- data %>% filter(SynthesisRate >= lower_bound & SynthesisRate <= upper_bound)

Q1 <- quantile(data$DecayRate, 0.25, na.rm = TRUE)
Q3 <- quantile(data$DecayRate, 0.75, na.rm = TRUE)
IQR_value <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value
data <- data %>% filter(DecayRate >= lower_bound & DecayRate <= upper_bound)

# Compute sample sizes (n)
sample_sizes <- data %>% group_by(Condition) %>% summarise(n = n())

# Update condition names with (n = ...) in labels
data$Condition <- recode_factor(data$Condition,
                                "labelled + bPAP" = paste0("labelled + bPAP\n(n = ", sample_sizes$n[sample_sizes$Condition == "labelled + bPAP"], ")"))

# Define a function to calculate the total within-cluster sum of squares
wss <- function(x, k) {
  kmeans(x, k, nstart = 10, iter.max = 50)$tot.withinss
}

# Set the range of k values to test
k.values <- 1:15

# Calculate WSS for different k values for H1 synthesis rates
wss_values_bPAP <- sapply(k.values, function(k) wss(data$SynthesisRate, k))
elbow_data_bPAP <- data.frame(k = k.values, wss = wss_values_bPAP)

pdf(file.path(output_dir, "synthesisRate_numberOfClusters_elbow.pdf"))
ggplot(elbow_data_bPAP, aes(x = k, y = wss)) +
  geom_point() +
  geom_line() +
  labs(title = "Elbow Method for Optimal Number of Clusters - bPAP",
       x = "Number of Clusters (k)",
       y = "Total Within-Cluster Sum of Squares (WSS)") +
  theme_minimal()
dev.off()

# Use the optimal number of clusters
k <- 3 
set.seed(123)
kmeans_bPAP <- kmeans(data$SynthesisRate, centers = k)

# Merge clusters with the data and name the clusters
data$synthesisRate_cluster <- factor(kmeans_bPAP$cluster)
ggplot(data, aes(x = synthesisRate_cluster, y = log2(SynthesisRate))) + 
  geom_boxplot() +
  labs(title = "Synthesis Rate Clusters", x = "Cluster", y = "Synthesis Rate")
data$synthesisRate_cluster <- ifelse(data$synthesisRate_cluster == 1, "Middle Synthesis", 
                                                ifelse(data$synthesisRate_cluster == 2, "Slow Synthesis", "Fast Synthesis"))
data$synthesisRate_cluster <- factor(data$synthesisRate_cluster, levels = c("Slow Synthesis", "Middle Synthesis", "Fast Synthesis"))
summary(data$synthesisRate_cluster)
ggplot(data, aes(x = synthesisRate_cluster, y = log2(SynthesisRate))) + 
  geom_boxplot() +
  labs(title = "Synthesis Rate Clusters", x = "Cluster", y = "Synthesis Rate")

# Function to create violin and box plots for a feature
violin_colors <- c("#a6dba0", "#ffbf00", "#1b9e77")  # Colors for the violins (light green, yellow, green)
boxplot_colors <- "#FFFFFF"  # White fill for boxplots inside the violins

create_violin_box_plot <- function(data, feature, feature_label, suffix, violin_colors, boxplot_colors) {

  data <- data %>%
    mutate(feature_value = .[[paste0(feature, suffix)]])
  
  # Filter for finite values
  data <- data[is.finite(data$feature_value), ]
  
  # Define comparisons for significance testing
  comparisons <- list(c(paste0("Slow\nSynthesis\n(n = ", n_slow, ")"), paste0("Middle\nSynthesis\n(n = ", n_middle, ")")),
                      c(paste0("Slow\nSynthesis\n(n = ", n_slow, ")"), paste0("Fast\nSynthesis\n(n = ", n_fast,")")),
                      c(paste0("Middle\nSynthesis\n(n = ", n_middle, ")"), paste0("Fast\nSynthesis\n(n = ", n_fast,")")))
  
  # Plot
  p <- ggplot(data, aes_string(x = "synthesisRate_cluster", y = "feature_value", fill = "synthesisRate_cluster")) +  
    geom_violin(trim = TRUE, color = "black", alpha = 0.9) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = boxplot_colors, color = "black", alpha = 0.6) +
    scale_fill_manual(values = violin_colors) +
    labs(y = feature_label, x = "") +
    theme_pubr() +
    theme(
      legend.position = "none",
      # axis.title.y = element_text(face = "bold", size = 26),
      axis.title.y = element_text(face = "bold", size = 20),
      axis.text.x = element_text(face = "bold", size = 19),
      # axis.text.y = element_text(size = 19)
      axis.text.y = element_text(size = 12)
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    stat_compare_means(comparisons = comparisons, label = "p.signif", method = "wilcox.test", 
                       vjust = 0.65, size = 14, step.increase = 0.1, bracket.size = 1.0)  # Add significance annotation
  
  return(p)
}

# Data prep
bPAP <- data

# Define n numbers in groups
summary(bPAP$synthesisRate_cluster)
n_fast <- summary(bPAP$synthesisRate_cluster)["Fast Synthesis"]
n_middle <- summary(bPAP$synthesisRate_cluster)["Middle Synthesis"]
n_slow <- summary(bPAP$synthesisRate_cluster)["Slow Synthesis"]
bPAP$synthesisRate_cluster <- as.character(bPAP$synthesisRate_cluster)

bPAP$synthesisRate_cluster <- ifelse(bPAP$synthesisRate_cluster == "Fast Synthesis", paste0("Fast\nSynthesis\n(n = ", n_fast,")"), 
                                          ifelse(bPAP$synthesisRate_cluster == "Slow Synthesis", paste0("Slow\nSynthesis\n(n = ", n_slow, ")"), paste0("Middle\nSynthesis\n(n = ", n_middle, ")")))

bPAP$synthesisRate_cluster <- factor(bPAP$synthesisRate_cluster, levels = c(paste0("Slow\nSynthesis\n(n = ", n_slow, ")"),
                                                                                      paste0("Middle\nSynthesis\n(n = ", n_middle, ")"),
                                                                                      paste0("Fast\nSynthesis\n(n = ", n_fast,")")))
summary(bPAP$synthesisRate_cluster)

# Create individual plots for the three features using the H1 dataset
plot_counts <- create_violin_box_plot(bPAP, "GeneBodyCounts", "Normalised Counts", "", violin_colors, boxplot_colors)
plot_synthesis <- create_violin_box_plot(bPAP, "SynthesisRate", "Synthesis Rate (A.U.)", "", violin_colors, boxplot_colors)
plot_decay <- create_violin_box_plot(bPAP, "DecayRate", "Decay Rate (A.U.)", "", violin_colors, boxplot_colors)

# Correlation plot
correlation_data <- bPAP %>%
  dplyr::select(GeneBodyCounts, SynthesisRate) %>%
  filter(is.finite(GeneBodyCounts) & is.finite(SynthesisRate))

correlation_plot <- ggplot(correlation_data, aes(x = SynthesisRate, y = GeneBodyCounts)) +
  geom_point(color = "#1b9e77", size = 3, alpha = 0.7) +  # Scatter points
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # Trend line with confidence interval
  labs(x = "Synthesis Rate (A.U.)", y = "Normalised Counts") +  # Axis labels
  theme_pubr() +
  theme(
    # axis.title.x = element_text(face = "bold", size = 22),
    # axis.title.y = element_text(face = "bold", size = 22),
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    # axis.text.x = element_text(face = "bold", size = 18),
    axis.text.x = element_text( size = 12),
    # axis.text.y = element_text(face = "bold", size = 18)
    axis.text.y = element_text(size = 12)
  ) +
  stat_cor(method = "kendall", size = 10, p.accuracy = 0.0001, cor.coef.name = "tau", label.sep = "\n")  # Add correlation coefficient

# Arrange the plots in 1x4 format
final_plot <- ggarrange(plot_synthesis, plot_decay, plot_counts, correlation_plot, ncol = 4, nrow = 1)

# Save the plot to a PDF
pdf(file.path(output_dir, "synthesisRate_ViolinBoxPlots_correlation_bPAP.pdf"), height = 6, width = 24)
print(final_plot)
dev.off()


# Arrange the plots in 2x2 format
final_plot <- ggarrange(plot_synthesis, plot_decay, plot_counts, correlation_plot, ncol = 2, nrow = 2, 
                        labels = c("A", "B", "C", "D"), font.label = list(size = 44))

# Save the plot to a PDF
pdf(file.path(output_dir, "synthesisRate_ViolinBoxPlots_correlation_bPAP_forThesis.pdf"), height = 8, width = 12)
print(final_plot)
dev.off()
