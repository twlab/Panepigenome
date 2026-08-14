library(ggplot2)
library(dplyr)
library(tidyr)
library(maps)
library(scales)
library(grid)

continental_df <- tibble(
  continental_group = c("AFR", "AMR", "EAS", "EUR", "SAS"),
  longitude = c(18, -72, 118, 12, 77),
  latitude  = c(3, 8, 38, 52, 21),
  
  `Methylation difference` = c(
    6170, 148, 1417, 248, 26
  ),
  
  `Integrated model` = c(
    311627, 9430, 86946, 12267, 5939
  ),
  
  `CpG gain/loss` = c(
    613983, 25790, 195965, 27172, 16335
  )
)

plot_df <- continental_df %>%
  pivot_longer(
    cols = c(
      `Methylation difference`,
      `CpG gain/loss`,
      `Integrated model`
    ),
    names_to = "category",
    values_to = "value"
  ) %>%
  mutate(
    category = factor(
      category,
      levels = c(
        "Methylation difference",
        "CpG gain/loss",
        "Integrated model"
      )
    )
  )

bar_height_max <- 18
bar_width      <- 3.8
bar_spacing    <- 4.6

max_transformed <- max(sqrt(plot_df$value), na.rm = TRUE)

plot_df <- plot_df %>%
  mutate(
    category_number = as.numeric(category),
    
    bar_x = longitude + (category_number - 2) * bar_spacing,
    
    xmin = bar_x - bar_width / 2,
    xmax = bar_x + bar_width / 2,
    
    ymin = latitude,
    ymax = latitude + sqrt(value) / max_transformed * bar_height_max
  )

label_df <- continental_df %>%
  mutate(
    label_y = latitude + bar_height_max + 4.2
  )

world_map <- map_data("world") %>%
  filter(region != "Antarctica")

bar_colours <- c(
  "Methylation difference" = "#0e7175",
  "CpG gain/loss"          = "#fd7901",
  "Integrated model"       = "#c35bca"
)

# -----------------------------------
# Mini y-axis for each continent
# -----------------------------------
axis_breaks <- c(0, 300000, 600000)

axis_df <- continental_df %>%
  mutate(
    axis_x = longitude - 2.2 * bar_spacing   # put axis left of 3 bars
  ) %>%
  tidyr::crossing(value = axis_breaks) %>%
  mutate(
    y = latitude + sqrt(value) / max_transformed * bar_height_max,
    label = case_when(
      value == 0 ~ "0",
      # value == 200000 ~ "0.2M",
      # value == 400000 ~ "0.4M",
      value == 300000 ~ "0.3M",
      value == 600000 ~ "0.6M"
    )
  )

axis_line_df <- continental_df %>%
  mutate(
    axis_x = longitude - 2.2 * bar_spacing,
    y0 = latitude,
    y1 = latitude + bar_height_max
  )

p_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "#E6E6E6",
    colour = NA
  ) +
  
  # baseline under each continent group
  geom_segment(
    data = continental_df,
    aes(
      x = longitude - 1.5 * bar_spacing,
      xend = longitude + 1.5 * bar_spacing,
      y = latitude,
      yend = latitude
    ),
    linewidth = 0.3,
    colour = "black"
  ) +
  
  # bars
  geom_rect(
    data = plot_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = category
    ),
    colour = "black",
    alpha = 0.65,
    linewidth = 0.2
  ) +
  
  # continent labels
  geom_text(
    data = label_df,
    aes(
      x = longitude,
      y = label_y,
      label = continental_group
    ),
    size = 2.8,
    fontface = "bold",
    colour = "black"
  ) +
  
  # local vertical axis line for each continent
  geom_segment(
    data = axis_line_df,
    aes(
      x = axis_x, xend = axis_x,
      y = y0, yend = y1
    ),
    linewidth = 0.7,
    colour = "black"
  ) +
  
  # local axis tick marks
  geom_segment(
    data = axis_df,
    aes(
      x = axis_x - 1.2,
      xend = axis_x,
      y = y,
      yend = y
    ),
    linewidth = 0.7,
    colour = "black"
  ) +
  
  # local axis labels
  geom_text(
    data = axis_df,
    aes(
      x = axis_x - 1.8,
      y = y,
      label = label
    ),
    hjust = 1,
    size = 3,
    colour = "black"
  ) +
  
  scale_fill_manual(
    values = bar_colours,
    name = NULL
  ) +
  coord_quickmap(
    xlim = c(-170, 180),
    ylim = c(-58, 88),
    expand = FALSE,
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(
      size = 12,
      colour = "black"
    ),
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(2, 3, 2, 3)
  )

ggsave(
  "world_map_population_specific_CpGs_1vsother_p_axis_each_continent.pdf",
  plot = p_map,
  width = 5.8,
  height = 3.6
)
