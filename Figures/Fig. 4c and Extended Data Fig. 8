library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 12))

df <- tribble(
  ~n_groups, ~`Shared CpG methylation`, ~`CpG gain/loss`, ~`Integrated model`,
  1,16123,1302486,711354,
  2,	2055,597309,274556,
  3	,23,172942,56638,
  4	,1,31825,4174,
  5	,0,3223,301,
)

plot_df <- df %>%
  pivot_longer(
    -n_groups,
    names_to = "signature",
    values_to = "count"
  ) %>%
  mutate(
    count_plot = ifelse(count == 0, 0.5, count),
    label = comma(count),
    signature = factor(
      signature,
      levels = c(
        "Shared CpG methylation",
        "CpG gain/loss",
        "Integrated model"
      )
    )
  )

p <- ggplot(plot_df, aes(x = n_groups, y = count_plot, color = signature),alpha = 0.6) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.3) +
  geom_text(
    aes(label = label),
    size = 2.7,
    vjust = -0.8,
    show.legend = FALSE
  ) +
  scale_y_log10(
    labels = label_number(scale_cut = cut_short_scale()),
    breaks = c(0.5, 1, 10, 100, 1000, 10000, 100000, 1000000),
    limits = c(0.45, 2e6)
  ) +
  scale_x_continuous(
    breaks = 1:5,
    limits = c(1, 5)
  ) +
  scale_color_manual(
    values = c(
      "Shared CpG methylation" = "#0e7175",
      "CpG gain/loss" = "#fd7901",
      "Integrated model" = "#c35bca"
    )
  ) +
  labs(
    x = "Number of continental-group comparisons",
    y = "Number of population-specific signatures",
    color = NULL
  ) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 9.5),
    panel.grid = element_blank()
  )

# p

ggsave(
  "population_specific_signature_overlap.pdf",
  plot = p,
  width = 5.2,
  height = 3.6,
  device = "pdf"
)
