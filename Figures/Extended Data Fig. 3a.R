library(ggplot2)
library(dplyr)
library(scales)
library(grid)

library(ggplot2)
library(dplyr)
library(tidyr)

# ---- Load data ----
data <- read.table("input.log", header = FALSE)
data <- data[, c("V5", "V3", "V6")]  # freq, meth, pop
colnames(data) <- c("freq", "meth", "pop")
data$pop <- as.factor(data$pop)

# ---- Reshape to long format ----
data_long <- data %>%
  pivot_longer(cols = c(freq, meth), names_to = "Type", values_to = "Methylation") %>%
  mutate(Type = ifelse(Type == "freq", "Ref", "Nonref"))

median(data_long$Methylation[data_long$Type=='Nonref'])

median(data_long$Methylation[data_long$Type=='Ref'])


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdist)
library(grid)

# ==================================================
# 1. Load and reshape data
# ==================================================
data <- read.table(
  "input.log",
  header = FALSE,
  stringsAsFactors = FALSE
)

# V5 = reference methylation
# V3 = non-reference methylation
# V6 = continental group
data <- data[, c("V5", "V3", "V6")]

colnames(data) <- c(
  "Reference",
  "Non-reference",
  "pop"
)

data_long <- data %>%
  pivot_longer(
    cols = c("Reference", "Non-reference"),
    names_to = "Type",
    values_to = "Methylation"
  ) %>%
  filter(
    is.finite(Methylation),
    pop %in% c("AFR", "AMR", "EAS", "EUR", "SAS")
  ) %>%
  mutate(
    pop = factor(
      pop,
      levels = c("AFR", "AMR", "EAS", "EUR", "SAS")
    ),
    Type = factor(
      Type,
      levels = c("Reference", "Non-reference")
    )
  )

# Increase spacing between continental groups
population_spacing <- 1.35

data_long <- data_long %>%
  mutate(
    pop_number = as.numeric(pop) * population_spacing,
    x_position = case_when(
      Type == "Reference" ~ pop_number - 0.18,
      Type == "Non-reference" ~ pop_number + 0.18
    )
  )

# Check group sizes
print(table(data_long$pop, data_long$Type))

# Overall medians
median_summary <- data_long %>%
  group_by(Type) %>%
  summarise(
    n = n(),
    median = median(Methylation, na.rm = TRUE),
    mean = mean(Methylation, na.rm = TRUE),
    .groups = "drop"
  )

print(median_summary)

# Save summary statistics by continental group
summary_stats <- data_long %>%
  group_by(pop, Type) %>%
  summarise(
    n = n(),
    mean = mean(Methylation, na.rm = TRUE),
    median = median(Methylation, na.rm = TRUE),
    q1 = quantile(Methylation, 0.25, na.rm = TRUE),
    q3 = quantile(Methylation, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  summary_stats,
  file = "population_reference_nonreference_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==================================================
# 2. Sample points for visualization
# ==================================================
# Density and box plots use all observations.


set.seed(123)

point_data <- data_long %>%
  group_by(pop, Type) %>%
  mutate(random_order = runif(n())) %>%
  arrange(random_order, .by_group = TRUE) %>%
  slice_head(n = 300) %>%
  ungroup() %>%
  select(-random_order)

# ==================================================
# 3. Colours
# ==================================================
type_colours <- c(
  "Reference" = "#4D4D4D",
  "Non-reference" = "#7B3294"
)

# ==================================================
# 4. Raincloud plot
# ==================================================
p <- ggplot(
  data_long,
  aes(
    x = x_position,
    y = Methylation,
    fill = Type,
    colour = Type,
    group = interaction(pop, Type)
  )
) +
  
  # Half-violin density
  ggdist::stat_halfeye(
    adjust = 0.6,
    width = 0.48,
    justification = -0.28,
    slab_alpha = 0.65,
    slab_colour = "black",
    slab_linewidth = 0.25,
    point_colour = NA,
    interval_colour = NA,
    .width = 0,
    show.legend = TRUE
  ) +
  
  # Box plot
  geom_boxplot(
    width = 0.10,
    outlier.shape = NA,
    linewidth = 0.35,
    colour = "black",
    alpha = 0.90,
    show.legend = FALSE
  ) +
  
  # Sampled observations
  geom_point(
    data = point_data,
    aes(
      x = x_position + 0.075,
      y = Methylation,
      colour = Type
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.035,
      height = 0,
      seed = 123
    ),
    size = 0.55,
    alpha = 0.25,
    stroke = 0,
    show.legend = FALSE
  ) +
  
  # White point indicates median
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.5,
    stroke = 0.35,
    fill = "white",
    colour = "black",
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = type_colours,
    breaks = c("Reference", "Non-reference"),
    name = NULL
  ) +
  
  scale_colour_manual(
    values = type_colours,
    breaks = c("Reference", "Non-reference"),
    name = NULL
  ) +
  
  scale_x_continuous(
    breaks = seq_along(levels(data_long$pop)) * population_spacing,
    labels = levels(data_long$pop),
    expand = expansion(add = c(0.55, 0.65))
  ) +
  
  scale_y_continuous(
    breaks = seq(50, 100, 10),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  
  coord_cartesian(
    ylim = c(50, 100),
    clip = "on"
  ) +
  
  labs(
    x = NULL,
    y = "DNA methylation"
  ) +
  
  theme_classic(
    base_size = 12,
    base_family = "sans"
  ) +
  
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.7
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.7
    ),
    axis.ticks.length = unit(1.5, "mm"),
    axis.text = element_text(
      colour = "black",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    axis.title.y = element_text(
      colour = "black",
      size = 12,
      margin = margin(r = 5)
    ),
    legend.position = "top",
    legend.justification = "left",
    legend.direction = "horizontal",
    legend.text = element_text(
      colour = "black",
      size = 12
    ),
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.30, "cm"),
    legend.spacing.x = unit(2, "mm"),
    plot.margin = margin(
      t = 3,
      r = 3,
      b = 3,
      l = 3,
      unit = "mm"
    )
  )



# ==================================================
# 5. Save
# ==================================================
ggsave(
  filename = "pop_ref_nonref_methylation_raincloud.pdf",
  plot = p,
  width = 5,
  height = 3.2
)

