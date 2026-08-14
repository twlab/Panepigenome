library(dplyr)
library(ggplot2)
library(grid)
library(scales)

# ============================================================
# 1. Read data
# ============================================================
data <- read.delim(
  "../../input",
  sep = "\t",
  stringsAsFactors = FALSE
)

goodfile <- read.table(
  "../goodfile.log",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Keep only BAM files in the self-defined filtered set
data_self <- data %>%
  filter(File %in% goodfile$V1)

# ============================================================
# 2. Summarize genomic coverage by sample
# ============================================================
summed_self <- data_self %>%
  group_by(Sample) %>%
  summarise(
    Sum_Coverage = sum(
      `genomic.coverage...5.`,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(is.finite(Sum_Coverage))

# Remove HG00272
summed_self <- summed_self %>%
  filter(Sample != "HG00272")

# Replace or add HG02080 with the specified coverage value
summed_self <- summed_self %>%
  filter(Sample != "HG02080") %>%
  bind_rows(
    tibble(
      Sample = "HG02080",
      Sum_Coverage = 22.32858976
    )
  ) %>%
  arrange(Sample)

# Confirm final sample number
cat("Final number of filtered samples:", nrow(summed_self), "\n")

# Save the final sample-level coverage table
write.table(
  summed_self,
  file = "self_defined_sample_coverage.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 3. Identify samples below the 30× threshold
# ============================================================
low_coverage_self <- summed_self %>%
  filter(Sum_Coverage < 30) %>%
  arrange(Sum_Coverage)

write.table(
  low_coverage_self,
  file = "self_defined_low_coverage_samples_below30x.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "Samples below 30×:",
  nrow(low_coverage_self),
  "\n"
)

print(low_coverage_self)

# ============================================================
# 4. Summary statistics
# ============================================================
summary_stats <- summed_self %>%
  summarise(
    n = n(),
    mean = mean(Sum_Coverage, na.rm = TRUE),
    median = median(Sum_Coverage, na.rm = TRUE),
    q1 = quantile(Sum_Coverage, 0.25, na.rm = TRUE),
    q3 = quantile(Sum_Coverage, 0.75, na.rm = TRUE),
    min = min(Sum_Coverage, na.rm = TRUE),
    max = max(Sum_Coverage, na.rm = TRUE),
    n_below_30 = sum(Sum_Coverage < 30, na.rm = TRUE)
  )

print(summary_stats)

write.table(
  summary_stats,
  file = "self_defined_coverage_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 5. Prepare plotting variables
# ============================================================
n_samples <- nrow(summed_self)
median_cov <- median(
  summed_self$Sum_Coverage,
  na.rm = TRUE
)

summary_label <- paste0(
  "Median = ",
  sprintf("%.1f×", median_cov)
)

# Classify samples relative to the threshold
plot_data <- summed_self %>%
  mutate(
    Coverage_status = if_else(
      Sum_Coverage < 30,
      "Below 30×",
      "At least 30×"
    )
  )

point_colours <- c(
  "At least 30×" = "#42771F",
  "Below 30×" = "#8B1A1A"
)

# ============================================================
# 6. plot
# ============================================================
set.seed(123)

p <- ggplot(
  plot_data,
  aes(
    x = 1,
    y = Sum_Coverage
  )
) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    fill = "#B7CDA8",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(colour = Coverage_status),
    width = 0.07,
    height = 0,
    shape = 16,
    size = 1.3,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = point_colours
  ) +
  geom_hline(
    yintercept = 30,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  annotate(
    geom = "text",
    x = 1.15,
    y = 30,
    label = "30× threshold",
    hjust = 0,
    vjust = -0.45,
    size = 3.0,
    colour = "black",
    family = "sans"
  ) +
  annotate(
    geom = "text",
    x = 0.81,
    y = max(
      plot_data$Sum_Coverage,
      na.rm = TRUE
    ),
    label = summary_label,
    hjust = 0,
    vjust = 1,
    size = 3.2,
    colour = "black",
    family = "sans"
  ) +
  scale_x_continuous(
    breaks = 1,
    labels = paste0(
      "Filtered samples\n(n = ",
      n_samples,
      ")"
    ),
    limits = c(0.75, 1.32),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(
      mult = c(0.02, 0.06)
    )
  ) +
  labs(
    x = NULL,
    y = "Genomic coverage per sample (×)"
  ) +
  theme_classic(
    base_size = 12,
    base_family = "sans"
  ) +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.75
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.75
    ),
    axis.ticks.x = element_blank(),
    axis.ticks.length = unit(1.4, "mm"),
    axis.text = element_text(
      colour = "black",
      size = 11
    ),
    axis.text.x = element_text(
      colour = "black",
      size = 11,
      angle = 0,
      hjust = 0.5,
      vjust = 1,
      lineheight = 0.95,
      margin = margin(t = 4)
    ),
    axis.title.y = element_text(
      colour = "black",
      size = 12,
      margin = margin(r = 4)
    ),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      t = 5,
      r = 8,
      b = 8,
      l = 5,
      unit = "pt"
    )
  )


# ============================================================
# 7. Save figure
# ============================================================
ggsave(
  filename = "Self_defined_genomic_coverage.pdf",
  plot = p,
  width = 1.8,
  height = 2.8
)
