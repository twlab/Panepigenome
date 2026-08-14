library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 12))

df <- tibble(
  event = rep(c("CGI gain", "CGI loss"), each = 3),
  category = rep(
    c("SINE/AluY-driven", "Other MEI-driven", "Non-MEI-driven"),
    2
  ),
  percent = c(
    83.1 * 93.5 / 100,   # CGI gain: SINE/AluY-driven
    83.1 * 6.5 / 100,    # CGI gain: other MEI-driven
    100 - 83.1,          # CGI gain: non-MEI-driven
    
    90.9 * 92.7 / 100,   # CGI loss: SINE/AluY-driven
    90.9 * 7.3 / 100,    # CGI loss: other MEI-driven
    100 - 90.9           # CGI loss: non-MEI-driven
  )
) %>%
  mutate(
    event = factor(
      event,
      levels = c("CGI gain", "CGI loss"),
      labels = c("CGI gains", "CGI losses")
    ),
    category = factor(
      category,
      levels = c(
        "SINE/AluY-driven",
        "Other MEI-driven",
        "Non-MEI-driven"
      )
    ),
    label = ifelse(
      percent >= 5,
      paste0(round(percent, 1), "%"),
      ""
    )
  )

bar_colors <- c(
  "SINE/AluY-driven" = "#D55E00",
  "Other MEI-driven" = "#E69F00",
  "Non-MEI-driven"   = "grey85"
)

p <- ggplot(
  df,
  aes(
    x = event,
    y = percent,
    fill = category,alpha = 0.6
  )
) +
  geom_col(
    width = 0.62,
    colour = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.4,
    colour = "black"
  ) +
  scale_fill_manual(
    values = bar_colors,
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "CGI events (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    axis.text.x = element_text(
      colour = "black",
      size = 11
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 11
    ),
    axis.title.y = element_text(
      colour = "black",
      size = 12
    ),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.35, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(4, 6, 4, 4)
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )


ggsave(
  "CGI_gain_loss_MEI_stacked_bar.pdf",
  plot = p,
  width = 2.8,
  height = 2.8
)
