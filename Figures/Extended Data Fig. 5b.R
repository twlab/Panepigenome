# Load required libraries
library(ggplot2)
library(dplyr)

# Read data correctly



data <- read.table("input.log", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# data$count<-data$count/1000000

# Ensure correct column names
# colnames(data) <- c("id", "pop", "count", "type", "hap","category")

# Define custom colors for populations (Nature-style palette)
population_colors <- c("EAS" = "#2ca02c",
                       "SAS" = "#ff7f0e",
                       "AFR" = "#1f77b4",
                       "EUR" = "#9467bd",
                       "AMR" = "#d62728")

# Convert 'type' to a factor with the desired order
data$type <- factor(data$type, levels = c("SNP", "Indel-insertion", "Indel-deletion", "SV-insertion", "SV-deletion"))

data<-data[data$type == "SV-insertion"  | data$type == "SV-deletion",]

mean(data$count)
# [1] 0.02250672

# Save the plot as a high-resolution PDF
pdf("Length_SV_CpG_normalized.pdf", width = 4, height = 3.5)  # Adjust size as needed

ggplot(data, aes(x = type, y = count)) +
  geom_jitter(
    aes(color = pop),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
    size = 1.8,
    alpha = 0.8,
    show.legend = c(shape = FALSE)  # remove category legend only
  ) +
  geom_hline(yintercept = 0.0105, color = "grey", linetype = "dashed", size = 1) +  # ← Add this line
  scale_color_manual(values = population_colors) +
  scale_shape_manual(values = c(16, 17, 15, 18, 3, 4)) +
  theme_classic(base_size = 14) +
  labs(
    y = "CpG density across SVs",
    color = "Population"
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 1)
    # no guide for shape needed
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "top",
    legend.box = "vertical",
    legend.spacing.y = unit(1, "pt"),
    legend.margin = margin(t = 0, b = 0),
    plot.title = element_blank()
  )




# Close the PDF device
dev.off()
