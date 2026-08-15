library(ggplot2)
library(maps)
library(dplyr)
library(ggrepel)

theme_set(theme_classic(base_size = 10))

# ---------------------------
# Data
# ---------------------------
data <- data.frame(
  pop = c("CDX", "CHB", "JPT", "KHV", "CHS", "BEB", "GIH", "ITU", "PJL", "STU", 
          "ASW", "ACB", "ESN", "GWD", "LWK", "MSL", "YRI", "GBR", "FIN", "IBS", 
          "TSI", "CEU", "CLM", "MXL", "PEL", "PUR", "ASL", "MKK"),
  region = c("East Asia", "East Asia", "East Asia", "East Asia", "East Asia", 
             "South Asia", "South Asia", "South Asia", "South Asia", "South Asia", 
             "Africa", "Africa", "Africa", "Africa", "Africa", "Africa", "Africa", 
             "Europe", "Europe", "Europe", "Europe", "Europe", "America", "America", 
             "America", "America", "Africa", "Africa"),
  lat = c(22, 39.916666, 35.68, 10.78, 23.13333, 23.7, 29.7589, 52.486243, 
          31.554606, 54.00, 35.483, 13.1, 9.06666, 13.454876, -1.27, 8.48, 7.4, 
          52.486243, 60.17, 40.38, 42.1, 40.767, 4.58333, 34.0544, -12.04, 
          18.4, 38.627003, 0.0236),
  lng = c(100.78, 116.383333, 139.68, 106.68, 113.266667, 90.35, -95.3677, 
          2, 74.357158, -1.890401, -97.53333, -59.62, 7.483333, -16.579032, 
          36.61, -13.23, 3.92, -1.890401, 24.93, -3.72, 12, -111.8904, 
          -74.066666, -118.2439, -77.03, -66.1, -90.199402, 37.9062),
  set = c(5, 5, 16, 13, 11, 10, 8, 3, 10, 5, 6, 12, 10, 17, 6, 9, 8, 
          10, 8, 3, 8, 0, 11, 2, 15, 16, 1, 1)
)

data <- data %>%
  filter(set > 0)

# ---------------------------
# colors
# ---------------------------
region_colors <- c(
  "Africa"     = "#1F78B4",
  "America"    = "#E31A1C",
  "East Asia"  = "#33A02C",
  "Europe"     = "#6A3D9A",
  "South Asia" = "#FF7F00"
)

world_map <- map_data("world")

# ---------------------------
# Plot
# ---------------------------
p <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    color = "white",
    linewidth = 0.08
  ) +
  geom_point(
    data = data,
    aes(x = lng, y = lat, fill = region, size = set),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = data,
    aes(x = lng, y = lat, label = pop),
    size = 2.5,
    fontface = "bold",
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    segment.linewidth = 0.2,
    segment.color = "grey50"
  ) +
  scale_fill_manual(values = region_colors) +
  scale_size_continuous(
    range = c(2.2, 6.5),
    breaks = c(1, 5, 10, 15),
    name = "Samples"
  ) +
  coord_quickmap(xlim = c(-170, 170), ylim = c(-55, 75), expand = FALSE) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Continental group"
  ) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

p

ggsave(
  "global_sample_distribution.pdf",
  plot = p,
  width = 7.2,
  height = 4.0,
  units = "in",
  device = cairo_pdf
)
