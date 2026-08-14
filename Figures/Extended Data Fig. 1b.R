library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 14))

# ---------------------------
# Read data
# ---------------------------
data <- read.table(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Bed-assembly/Analysis/merged_matched_col2.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# ---------------------------
# Extract feature columns
# ---------------------------
fisher_cols <- grep("\\.fisher_results\\.tsv$", colnames(data), value = TRUE)

feature_names <- fisher_cols %>%
  gsub("\\.fisher_results\\.tsv$", "", .) %>%
  gsub("_", " ", .)

heatmap_matrix <- as.matrix(data[, fisher_cols])
colnames(heatmap_matrix) <- feature_names

# ---------------------------
# Long format
# ---------------------------
plot_df <- data.frame(
  continent = data$continent,
  heatmap_matrix,
  stringsAsFactors = FALSE
) %>%
  pivot_longer(
    cols = -continent,
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  filter(is.finite(Value))

# ---------------------------
# Mean ± SD
# ---------------------------
summary_df <- plot_df %>%
  group_by(continent, Feature) %>%
  summarise(
    MeanValue = mean(Value, na.rm = TRUE),
    n = sum(!is.na(Value)),
    SD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    xmin = pmax(MeanValue - SD, 0),
    xmax = MeanValue + SD
  )

# ---------------------------
# Order features by overall mean
# ---------------------------
feature_order <- summary_df %>%
  group_by(Feature) %>%
  summarise(
    OverallMean = mean(MeanValue, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(OverallMean) %>%
  pull(Feature)

summary_df <- summary_df %>%
  mutate(
    Feature = factor(Feature, levels = feature_order),
    continent = factor(continent, levels = c("AFR", "AMR", "EAS", "EUR", "SAS"))
  )

pd <- position_dodge(width = 0.65)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(summary_df, aes(x = MeanValue, y = Feature, color = continent)) +
  geom_linerange(
    aes(xmin = xmin, xmax = xmax),
    position = pd,
    linewidth = 0.8,
    alpha = 0.75
  ) +
  geom_point(
    position = pd,
    size = 3.0,
    alpha = 0.95
  ) +
  scale_color_manual(
    values = c(
      AFR = "#0072B2",
      AMR = "#D55E00",
      EAS = "#009E73",
      EUR = "#CC79A7",
      SAS = "#E69F00"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 0.81),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Proportion of non-reference CpGs overlapping genomic elements",
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 11),
    panel.grid = element_blank()
  )

p

ggsave(
  "Non_ref_feature_overlap_by_continent_dot_errorbar_mean_SD.pdf",
  plot = p,
  width = 7.2,
  height = 5.2,
  units = "in"
)
