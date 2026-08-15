library(ggplot2)
library(scales)

pie_df <- data.frame(
  feature = c(
    "Self var-CpG\ncopy number",
    "Nearby var-CpG\ncopy number",
    "Nearby non-CpG-\naltering variants"
  ),
  count = c(3042, 3730, 2097)
)

pie_df$proportion <- pie_df$count / sum(pie_df$count)

pie_df$label <- sprintf(
  "%s\n%.1f%%",
  comma(pie_df$count),
  100 * pie_df$proportion
)

pie_df$feature <- factor(
  pie_df$feature,
  levels = c(
    "Self var-CpG\ncopy number",
    "Nearby var-CpG\ncopy number",
    "Nearby non-CpG-\naltering variants"
  )
)

nature_cols <- c(
  "Self var-CpG\ncopy number"        = "#3B6FB6",
  "Nearby var-CpG\ncopy number"       = "#6BAED6",
  "Nearby non-CpG-\naltering variants" = "#D98C3F"
)

p_pie <- ggplot(
  pie_df,
  aes(x = "", y = count, fill = feature)
) +
  geom_col(
    width = 1,
    colour = "white",
    linewidth = 0.6
  ) +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    colour = "black",
    size = 3.5,
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = nature_cols,
    breaks = levels(pie_df$feature)
  ) +
  labs(fill = NULL) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    legend.text = element_text(
      size = 9,
      colour = "black"
    ),
    legend.key.size = unit(0.45, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(
  "shared_varCpG_largest_effect_pie.pdf",
  p_pie,
  width = 3.7,
  height = 2.4,
  units = "in"
)
