############################################################################
## Project: SNUseq project
## Script purpose: Visualize statistics on region and signal deletion after subtraction
## Date: Mar 12, 2025
## Author: Umut Gerlevik
############################################################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# Set file paths
prefix <- "/MellorLab/SNUseqProject"

# Load the data
stats <- read.table(paste0(prefix, "/1_Umut/9.1_bPAPminusnoPAP_bigwig/size_and_signal_stats.txt"), header = TRUE, sep = "\t")

output_folder <- paste0(prefix, "/outputs/9.2_quantificationPlots_negRemoved_bPAPminusnoPAP")

if(!dir.exists(output_folder)) dir.create(output_folder)

# Clean up file names to remove strand information
stats$File <- gsub("_rev", "", stats$File)
stats$File <- gsub("_fwd", "", stats$File)

# Summarize and combine the statistics for each file
data_combined <- stats %>%
  group_by(File) %>%
  summarize(
    Combined_Signal_Before = sum(Total_Signal_Before),
    Combined_Signal_After = sum(Total_Signal_After),
    Combined_Deleted_Signal = Combined_Signal_Before - Combined_Signal_After,
    Combined_Size_Before = sum(Total_Size_Before),
    Combined_Size_After = sum(Total_Size_After),
    Combined_Deleted_Size = Combined_Size_Before - Combined_Size_After,
    Percentage_Deleted_Size = (Combined_Deleted_Size / Combined_Size_Before) * 100,
    Percentage_Deleted_Signal = (Combined_Deleted_Signal / Combined_Signal_Before) * 100
  )

# Create individual plots with adjusted text size and bar width
plot1 <- ggplot(data_combined, aes(x = File)) +
  geom_bar(aes(y = Combined_Signal_Before, fill = "Signal Before"), stat = "identity", position = "dodge", width = 0.6) +
  geom_bar(aes(y = Combined_Signal_After, fill = "Signal After"), stat = "identity", position = "dodge", width = 0.6) +
  labs(title = "Total Signal Before and After Subtraction", y = "Total Signal", x = "File") +
  scale_fill_manual(values = c("Signal Before" = "steelblue", "Signal After" = "tomato")) +
  theme_minimal() +
  theme(axis.text = element_text(angle = 45, hjust = 1, size = 12), legend.title = element_blank(), legend.position = c(0.1, 0.9),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))

plot2 <- ggplot(data_combined, aes(x = File)) +
  geom_bar(aes(y = Combined_Size_Before, fill = "Size Before"), stat = "identity", position = "dodge", width = 0.6) +
  geom_bar(aes(y = Combined_Size_After, fill = "Size After"), stat = "identity", position = "dodge", width = 0.6) +
  labs(title = "Total Size Before and After Subtraction", y = "Total Size", x = "File") +
  scale_fill_manual(values = c("Size Before" = "lightgreen", "Size After" = "orange")) +
  theme_minimal() +
  theme(axis.text = element_text(angle = 45, hjust = 1, size = 12), legend.title = element_blank(), legend.position = c(0.1, 0.9),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))

plot3 <- ggplot(data_combined, aes(x = File, y = Percentage_Deleted_Signal, fill = "Deleted Signal")) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "Percentage of Deleted Signal After Subtraction", y = "Percentage Deleted Signal (%)", x = "File") +
  scale_fill_manual(values = c("Deleted Signal" = "purple")) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text = element_text(angle = 45, hjust = 1, size = 12), legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))

plot4 <- ggplot(data_combined, aes(x = File, y = Percentage_Deleted_Size, fill = "Deleted Size")) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "Percentage of Deleted Size After Subtraction", y = "Percentage Deleted Size (%)", x = "File") +
  scale_fill_manual(values = c("Deleted Size" = "darkred")) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text = element_text(angle = 45, hjust = 1, size = 12), legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))

# Combine all plots into one figure
combined_plot <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)

# Save the combined plot as a PDF
ggsave(file.path(output_folder, "/deletedSignalStatistics_bPAPminusNoPAP.pdf"), combined_plot, width = 14, height = 12)



############## ONLY LABELLED #################
# Subset the data
data_combined <- data_combined %>% 
  filter(grepl("labelled", File))

# Create individual plots with adjusted text size and bar width
plot1 <- ggplot(data_combined, aes(x = File)) +
  geom_bar(aes(y = Combined_Signal_Before, fill = "Signal Before"), stat = "identity", position = "dodge", width = 0.6) +
  geom_bar(aes(y = Combined_Signal_After, fill = "Signal After"), stat = "identity", position = "dodge", width = 0.6) +
  labs(title = "Total Signal Before and After Subtraction", y = "Total Signal", x = "File") +
  scale_fill_manual(values = c("Signal Before" = "steelblue", "Signal After" = "tomato")) +
  scale_y_continuous(breaks = seq(0, 20000000, 2000000)) +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 12), 
        axis.text.x = element_text(size = 16, face = "bold"),
        legend.title = element_blank(), legend.position = c(0.1, 0.9),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"))

plot2 <- ggplot(data_combined, aes(x = File)) +
  geom_bar(aes(y = Combined_Size_Before, fill = "Size Before"), stat = "identity", position = "dodge", width = 0.6) +
  geom_bar(aes(y = Combined_Size_After, fill = "Size After"), stat = "identity", position = "dodge", width = 0.6) +
  labs(title = "Total Size Before and After Subtraction", y = "Total Size", x = "File") +
  scale_fill_manual(values = c("Size Before" = "lightgreen", "Size After" = "orange")) +
  scale_y_continuous(breaks = seq(0, 20000000, 2000000)) +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 12), 
        axis.text.x = element_text(size = 16, face = "bold"),
        legend.title = element_blank(), legend.position = c(0.1, 0.9),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"))

plot3 <- ggplot(data_combined, aes(x = File, y = Percentage_Deleted_Signal, fill = "Deleted Signal")) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "Percentage of Deleted Signal After Subtraction", y = "Percentage Deleted Signal (%)", x = "File") +
  scale_fill_manual(values = c("Deleted Signal" = "purple")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 16, face = "bold"), 
        legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"))

plot4 <- ggplot(data_combined, aes(x = File, y = Percentage_Deleted_Size, fill = "Deleted Size")) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "Percentage of Deleted Size After Subtraction", y = "Percentage Deleted Size (%)", x = "File") +
  scale_fill_manual(values = c("Deleted Size" = "darkred")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 16, face = "bold"),
        legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"))

# Combine all plots into one figure
combined_plot <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2,
                           labels = c("A", "B", "C", "D"), font.label = list(size = 44))

# Save the combined plot as a PDF
ggsave(file.path(output_folder, "/deletedSignalStatistics_bPAPminusNoPAP_onlyLabelled.pdf"), combined_plot, width = 13, height = 12)
