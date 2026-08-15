library(ggplot2)
library(dplyr)
library(maps)
library(ggrepel)
library(scales)
library(grid)
library(scatterpie)

# --------------------------------------------------
# Input data
# --------------------------------------------------
af_df <- tibble(
  population  = c("AFR", "AMR", "EAS", "EUR", "SAS"),
  longitude   = c(18, -72, 118, 12, 77),
  latitude    = c(3, 8, 38, 52, 21),
  AF          = c(0.1536, 0.03977, 0.01471, 0.03333, 0.006944),
  sample_size = c(70, 44, 51, 30, 36)
)

stopifnot(
  all(is.finite(af_df$longitude)),
  all(is.finite(af_df$latitude)),
  all(is.finite(af_df$AF)),
  all(af_df$AF >= 0 & af_df$AF <= 1),
  all(is.finite(af_df$sample_size)),
  all(af_df$sample_size > 0)
)

# --------------------------------------------------
# Prepare data for pie chart
# --------------------------------------------------
af_df <- af_df %>%
  mutate(
    AluY_insertion = AF,
    Other          = 1 - AF,
    radius         = rescale(sqrt(sample_size), to = c(4.5, 8.5)),
    label_txt      = paste0(population, "\n", percent(AF, accuracy = 0.1))
  )

# --------------------------------------------------
# World map background
# --------------------------------------------------
world_map <- map_data("world") %>%
  filter(region != "Antarctica")

# --------------------------------------------------
# Plot
# --------------------------------------------------
p_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "#E6E6E6",
    colour = NA
  ) +
  geom_scatterpie(
    data = af_df,
    aes(x = longitude, y = latitude, r = radius),
    cols = c("AluY_insertion", "Other"),
    color = "black",
    linewidth = 0.25
  ) +
  geom_text_repel(
    data = af_df,
    aes(
      x = longitude,
      y = latitude,
      label = label_txt
    ),
    size = 3.6,
    colour = "black",
    lineheight = 0.95,
    box.padding = 0.35,
    point.padding = 0.45,
    min.segment.length = 0,
    segment.colour = "#737373",
    segment.linewidth = 0.22,
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_fill_manual(
    values = c(
      "AluY_insertion" = "#3B4CC0",
      "Other"          = "#F2F2F2"
    ),
    labels = c(
      "AluY_insertion" = "Insertion",
      "Other"          = "Non-insertion"
    ),
    name = NULL
  ) +
  coord_quickmap(
    xlim = c(-170, 180),
    ylim = c(-58, 82),
    expand = FALSE,
    clip = "off"
  ) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 11, colour = "black"),
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(2, 3, 2, 3)
  )

ggsave(
  "world_map_AluY_LINC01572_pie_sampleSize_labels.pdf",
  plot = p_map,
  width = 3.8,
  height = 3
)
