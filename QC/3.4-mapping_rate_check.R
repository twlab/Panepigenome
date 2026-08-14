# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read data
data <- read.delim("example_input_3.2and3.3and3.4.log", header = T, sep = "\t", stringsAsFactors = FALSE)
data<-data[,c(2,5)]
colnames(data) <- c("FilePath", "Value")
data$Value <- as.numeric((data$Value))
data$Sample <- basename(data$FilePath)

# Compute boxplot statistics
q1 <- quantile(data$Value, 0.25)
q3 <- quantile(data$Value, 0.75)
iqr <- q3 - q1
lower_bound <- q1 - 1.5 * iqr
upper_bound <- q3 + 1.5 * iqr

# Mark outliers
data <- data %>%
  mutate(Outlier = ifelse(Value < lower_bound | Value > upper_bound, "Outlier", "Normal"))

# Save outliers to a text file
outliers <- data %>% filter(Outlier == "Outlier")
write.table(outliers, file = "outliers_mapping_rate.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Create the plot
p<-ggplot(data, aes(x = "", y = Value)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", color = "black") +  # Remove default outliers
  geom_jitter(aes(color = Outlier), width = 0.1, alpha = 0.8, size = 2) +  # Custom outlier coloring
  scale_color_manual(values = c("Normal" = "black", "Outlier" = "red")) +  # Define colors
  # geom_hline(yintercept = 0.99, linetype = "dashed", color = "blue", linewidth = 1) +
  labs(
    # title = "Sequence",
    y = "Map rate",
    x = "",
    color = "Point Type"
  ) +
  theme_minimal()#+coord_cartesian(ylim = c(99, 100))
ggsave(filename = "Mappingrate.pdf", plot = p, width = 4, height = 3)

