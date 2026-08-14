library(tidyverse)

# Step 1: Load and reshape data
data <- read.delim("example_input_4.2and4.3_methylation_check.log", sep = "\t", stringsAsFactors = FALSE)
data$Sample <- basename(data$File)

# Select relevant columns (File and Bin1–Bin5)
colnames(data)[5:9] <- c("Bin1[0-20]", "Bin2[20-40]", "Bin3[40-60]", "Bin4[60-80]", "Bin5[80-100]")
# Select relevant columns (File and Bin1–Bin5)
bin_cols <- colnames(data)[5:9]
data_long <- reshape(
  data[c("Sample", bin_cols)],
  varying = bin_cols,
  v.names = "Value",
  timevar = "Bin",
  times = bin_cols,
  direction = "long"
)

# Step 2: Identify outliers per Bin
outlier_info <- data_long %>%
  group_by(Bin) %>%
  summarise(
    Q1 = quantile(Value, 0.25),
    Q3 = quantile(Value, 0.75),
    IQR = Q3 - Q1,
    Lower = Q1 - 1.5 * IQR,
    Upper = Q3 + 1.5 * IQR,
    .groups = "drop"
  )

data_long <- data_long %>%
  left_join(outlier_info, by = "Bin") %>%
  mutate(
    Outlier = ifelse(Value < Lower | Value > Upper, "Outlier", "Normal")
  )

# Step 3: Identify samples that have at least one outlier
outlier_samples <- data_long %>%
  group_by(Sample) %>%
  summarise(AnyOutlier = any(Outlier == "Outlier"), .groups = "drop")

# Step 4: Extract all unique outlier entries
outliers <- data_long %>%
  filter(Outlier == "Outlier") %>%
  distinct(Sample, Bin, Value)

# Step 5: Save to a tab-delimited text file
write.table(outliers, file = "outliers_methylation_category.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Merge info back to data
data_long <- data_long %>%
  left_join(outlier_samples, by = "Sample")

# Step 4: Plot
p<-ggplot(data_long, aes(x = Bin, y = Value, group = Sample)) +
  geom_boxplot(aes(group = Bin), outlier.shape = NA, fill = "lightblue", color = "black") +
  geom_jitter(aes(color = Outlier), width = 0.1, alpha = 0.9, size = 2) +
  geom_line(data = data_long %>% filter(AnyOutlier), aes(group = Sample), color = "red", size = 0.7, alpha = 0.7) +
  scale_color_manual(values = c("Normal" = "black", "Outlier" = "red")) +
  labs(
    title = "Methylation Distribution Across Bins",
    x = "Bin (methylation level range)",
    y = "Proportion",
    color = "Point Type"
  ) +
  theme_minimal()
ggsave(filename = "methylation_category.pdf", plot = p, width = 6, height = 4)


