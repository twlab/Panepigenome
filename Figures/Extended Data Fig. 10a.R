library(dplyr)
library(ggplot2)
library(scales)

df_plot <- data.frame(
  model = c(
    "Self var-CpG copy number",
    "Self var-CpG copy number",
    "Nearby var-CpG copy number",
    "Nearby var-CpG copy number"
  ),
  set = c(
    "All associated var-CpGs",
    "Shared across three models",
    "All associated var-CpGs",
    "Shared across three models"
  ),
  prop_positive = c(
    0.6663199,
    0.9954899,
    0.5272111,
    0.6261134
  )
) %>%
  mutate(
    percent = prop_positive * 100,
    model = factor(
      model,
      levels = c("Self var-CpG copy number", "Nearby var-CpG copy number")
    ),
    set = factor(
      set,
      levels = c("All associated var-CpGs", "Shared across three models")
    ),
    label = sprintf("%.1f%%", percent)
  )

p <- ggplot(df_plot, aes(x = percent, y = model, color = set)) +
  geom_segment(
    aes(x = 0, xend = percent, yend = model),
    color = "grey70",
    linewidth = 0.45
  ) +
  geom_point(
    size = 3.5,
    alpha = 0.95
  ) +
  geom_text(
    aes(label = label),
    hjust = -0.25,
    size = 3.2,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "All associated var-CpGs" = "grey45",
      "Shared across three models" = "#D55E00"
    )
  ) +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 108),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Positive effect direction (%)",
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 9.5),
    panel.grid = element_blank(),
    plot.margin = margin(5, 18, 5, 5)
  )

# p

ggsave(
  "positive_effect_direction_varCpGs.pdf",
  plot = p,
  width = 5.0,
  height = 2.5,
  device = "pdf"
)
