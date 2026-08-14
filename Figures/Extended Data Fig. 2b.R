library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 12))

df <- tribble(
  ~metric, ~tested, ~outlier, ~accumulated,
  "Median ML value",        1158, 3,  NA,
  "ML value category",      1158, 5,  5,
  "% low-quality reads",    1158, 2,  7,
  "MapQ 30",                1158, 0,  7,
  "Alignment metrics",      1158, 4,  9,
  "Median methylation",     1157, 13, 19,
  "Methylation category",   1157, 44, 47
) %>%
  mutate(
    metric = factor(metric, levels = metric),
    outlier_pct = outlier / tested * 100,
    accumulated_pct = accumulated / tested * 100
  )

p <- ggplot(df, aes(x = metric)) +
  geom_col(
    aes(y = outlier_pct),
    width = 0.58,
    fill = "grey75",
    color = "black",
    linewidth = 0.25
  ) +
  geom_line(
    aes(y = accumulated_pct, group = 1),
    linewidth = 0.55,
    color = "#7B3294",
    na.rm = TRUE
  ) +
  geom_point(
    aes(y = accumulated_pct),
    size = 2.6,
    shape = 21,
    fill = "#7B3294",
    color = "black",
    stroke = 0.25,
    na.rm = TRUE
  ) +
  geom_text(
    aes(y = outlier_pct, label = outlier),
    vjust = -0.45,
    size = 3.0,
    color = "black"
  ) +
  geom_text(
    aes(y = accumulated_pct, label = accumulated),
    vjust = -0.9,
    size = 3.0,
    color = "#7B3294",
    na.rm = TRUE
  ) +
  scale_y_continuous(
    labels = percent_format(scale = 1, accuracy = 1),
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    x = NULL,
    y = "Outlier samples (%)"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 8, 5, 5)
  )

# p

ggsave(
  "QC_outlier_accumulation.pdf",
  plot = p,
  width = 6,
  height = 3.2,
  device = "pdf"
)
