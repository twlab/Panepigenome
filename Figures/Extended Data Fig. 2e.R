library(dplyr)
library(ggplot2)
library(viridis)
library(scales)
library(grid)

data <- fread(
  "/scratch/zdong/Projects/PanEpiG/Benchmark/HiFi/Pbmm2/pb-CpG-tools/NA19909_matched.txt"
)

data <- data %>%
  transmute(
    V4 = suppressWarnings(as.numeric(V4)) / 100,
    V5 = suppressWarnings(as.numeric(V5)) / 100
  ) %>%
  filter(
    is.finite(V4),
    is.finite(V5)
  )

# Confirm that the values are actually between 0 and 1
print(range(data$V4))
print(range(data$V5))
print(nrow(data))

# Keep only biologically valid methylation values
data <- data %>%
  filter(
    between(V4, 0, 1),
    between(V5, 0, 1)
  ) %>%
  mutate(
    abs_diff = abs(V4 - V5)
  )

# ============================================================
# Correlation
# ============================================================
cor_res <- cor.test(
  data$V4,
  data$V5,
  method = "pearson"
)

r_label <- paste0(
  "n = ", format(nrow(data), big.mark = ","),
  "\nPearson’s r = ",
  sprintf("%.3f", unname(cor_res$estimate)),
  "\nP < 2.2 × 10⁻¹⁶"
)

# ============================================================
# Plot
# ============================================================
p <- ggplot(
  data,
  aes(
    x = V4,
    y = V5
  )
) +
  stat_bin_2d(
    aes(fill = after_stat(count)),
    binwidth = c(0.01, 0.01),
    boundary = 0,
    drop = TRUE
  ) +
  scale_fill_viridis_c(
    option = "mako",
    direction = 1,
    trans = scales::transform_log10(),
    name = "CpG count",
    breaks = scales::breaks_log(n = 5),
    labels = scales::label_number(
      scale_cut = scales::cut_short_scale()
    ),
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3, "mm"),
      barheight = unit(25, "mm"),
      ticks = TRUE,
      frame.colour = "black",
      frame.linewidth = 0.3
    )
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "white",
    linewidth = 0.5
  ) +
  annotate(
    geom = "text",
    x = 0.03,
    y = 0.97,
    label = r_label,
    hjust = 0,
    vjust = 1,
    size = 3.2,
    lineheight = 1.05,
    colour = "black"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = 0)
  ) +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, 1),
    ratio = 1,
    clip = "on"
  ) +
  labs(
    x = "Methylation level (2023/12/07)",
    y = "Methylation level (2023/07/28)"
  ) +
  theme_classic(
    base_size = 12,
    base_family = "sans"
  ) +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.7
    ),
    axis.ticks.length = unit(1.3, "mm"),
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    axis.title.x = element_text(
      margin = margin(t = 4)
    ),
    axis.title.y = element_text(
      margin = margin(r = 4)
    ),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      size = 11,
      colour = "black"
    ),
    legend.text = element_text(
      size = 10,
      colour = "black"
    ),
    plot.margin = margin(
      4, 4, 4, 4,
      unit = "mm"
    )
  )

ggsave(
  filename = "correlation_time_NA19909.pdf",
  plot = p,
  width = 4.2,
  height = 3.0
)

