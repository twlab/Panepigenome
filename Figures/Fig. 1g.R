library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 14))

# ---------------------------
# Read input data
# ---------------------------
cpg_nodes_growth <- read_tsv(
  "cpg_nodes_panacus_growth2.tsv",
  skip = 6,
  col_names = c("count", "[0-0.01]", "[0.01-0.1]", "[0.1-0.9]", "[0.9-1.0]")
)

cpg_edges_growth <- read_tsv(
  "cpg_edges_panacus_growth2.tsv",
  skip = 6,
  col_names = c("count", "[0-0.01]", "[0.01-0.1]", "[0.1-0.9]", "[0.9-1.0]")
)

# ---------------------------
# Combine nodes and edges
# ---------------------------
cpg_growth_df <- bind_rows(cpg_nodes_growth, cpg_edges_growth) %>%
  group_by(count) %>%
  summarise(
    `[0-0.01]`  = sum(`[0-0.01]`, na.rm = TRUE),
    `[0.01-0.1]` = sum(`[0.01-0.1]`, na.rm = TRUE),
    `[0.1-0.9]`  = sum(`[0.1-0.9]`, na.rm = TRUE),
    `[0.9-1.0]`  = sum(`[0.9-1.0]`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -count,
    names_to = "Frequency",
    values_to = "cpgs"
  ) %>%
  mutate(
    Frequency = factor(
      Frequency,
      levels = c("[0.9-1.0]", "[0.1-0.9]", "[0.01-0.1]", "[0-0.01]"),
      labels = c(
        "High frequency (AF ≥ 90%)",
        "Intermediate frequency (10% ≤ AF < 90%)",
        "Low frequency (1% ≤ AF < 10%)",
        "Rare frequency (AF < 1%)"
      )
    )
  )

# ---------------------------
# Nature-style colors
# ---------------------------
freq_colors <- c(
  "High frequency (AF ≥ 90%)" = "#2166AC",
  "Intermediate frequency (10% ≤ AF < 90%)" = "#67A9CF",
  "Low frequency (1% ≤ AF < 10%)" = "#FDAE61",
  "Rare frequency (AF < 1%)" = "#B2182B"
)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(
  cpg_growth_df,
  aes(x = count, y = cpgs, fill = Frequency, color = Frequency,alpha = 0.6)
) +
  geom_area(
    position = "identity",
    alpha = 0.6,
    linewidth = 0.20
  ) +
  scale_fill_manual(values = freq_colors) +
  scale_color_manual(values = freq_colors) +
  scale_y_continuous(
    labels = label_number(scale = 1e-6, suffix = "M", accuracy = 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_x_continuous(
    breaks = c(1, 100, 200, 300, 400, 464),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Number of haplotypes",
    y = "Cumulative number of CpGs",
    fill = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 13, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.42, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

# p

ggsave(
  "CpG_saturation_HPRC2.pdf",
  plot = p,
  width = 6,
  height = 3.6
)
