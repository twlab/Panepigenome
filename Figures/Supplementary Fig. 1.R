library(ggplot2)
library(dplyr)
library(data.table)

# --------------------------------------------------
# 1. Load and clean data
# --------------------------------------------------
df_combined <- fread("rmse.csv")

df_combined <- df_combined %>%
  filter(
    is.finite(rmse),
    !is.na(pop_label),
    !is.na(split)
  ) %>%
  mutate(
    split = factor(
      split,
      levels = c("train", "test"),
      labels = c("Training", "Testing")
    )
  )

# Check the order
print(levels(df_combined$split))

# --------------------------------------------------
# 2. Calculate and save summary statistics
# --------------------------------------------------
summary_stats <- df_combined %>%
  group_by(pop_label, split) %>%
  summarise(
    n = n(),
    mean = mean(rmse, na.rm = TRUE),
    median = median(rmse, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

fwrite(
  summary_stats,
  file = "rmse_mean_median.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# --------------------------------------------------
# 3. Calculate overall whisker range
# --------------------------------------------------
y_range <- df_combined %>%
  group_by(pop_label, split) %>%
  summarise(
    lower = boxplot.stats(rmse)$stats[1],
    upper = boxplot.stats(rmse)$stats[5],
    .groups = "drop"
  ) %>%
  summarise(
    ymin = min(lower, na.rm = TRUE),
    ymax = max(upper, na.rm = TRUE)
  )

padding <- 0.05 * (y_range$ymax - y_range$ymin)

# --------------------------------------------------
# 4. Plot
# --------------------------------------------------
p <- ggplot(
  df_combined,
  aes(
    x = pop_label,
    y = rmse,
    fill = split
  )
) +
  geom_boxplot(
    width = 0.65,
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    linewidth = 0.4,
    colour = "black"
  ) +
  stat_summary(
    aes(group = split),
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.75),
    shape = 21,
    size = 2.2,
    stroke = 0.4,
    fill = "white",
    colour = "black",
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Training" = "#f8d4bf",
      "Testing" = "#d189bb"
    ),
    breaks = c("Training", "Testing"),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "RMSE"
  ) +
  coord_cartesian(
    ylim = c(
      y_range$ymin - padding,
      y_range$ymax + padding
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.35
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.3
    ),
    axis.text = element_text(
      colour = "black",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    axis.title.y = element_text(
      colour = "black",
      size = 12,
      margin = margin(r = 6)
    ),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 12),
    legend.key.width = grid::unit(0.45, "cm"),
    legend.key.height = grid::unit(0.3, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

# --------------------------------------------------
# 5. Save
# --------------------------------------------------
ggsave(
  filename = "rmse_boxplot.pdf",
  plot = p,
  width = 3.5,
  height = 2.5
)
