library(tidyverse)
library(scales)
library(grid)

theme_set(
  theme_classic(
    base_size = 10,
    base_family = "Arial"
  )
)

# ============================================================
# 1. Input counts
# ============================================================

before_total <- 1994
after_total  <- 1891
shared       <- 1657

before_only <- before_total - shared
after_only  <- after_total - shared

union_total <- before_only + shared + after_only

# ============================================================
# 2. Prepare plotting data
# ============================================================

plot_df <- tibble(
  Category = factor(
    c("Before only", "Shared", "After only"),
    levels = c("Before only", "Shared", "After only")
  ),
  Count = c(
    before_only,
    shared,
    after_only
  )
) %>%
  mutate(
    proportion = Count / union_total,
    label = paste0(
      comma(Count),
      "\n",
      percent(proportion, accuracy = 0.1)
    )
  )

category_colours <- c(
  "Before only" = "#BDBDBD",
  "Shared"      = "#4D4D4D",
  "After only"  = "#A6761D"
)

# Positions for annotations below the stacked bar
before_midpoint <- before_total / 2

after_start <- before_only
after_midpoint <- after_start + after_total / 2

# ============================================================
# 3. Nature-style overlap-composition plot
# ============================================================

library(tidyverse)
library(scales)
library(grid)

theme_set(theme_void(base_size = 10))

# ---------------------------
# Input counts
# ---------------------------
before_total <- 1994
after_total  <- 1891
shared       <- 1657

before_only <- before_total - shared
after_only  <- after_total - shared
union_total <- before_only + shared + after_only

# ---------------------------
# Explicit segment positions
# ---------------------------
plot_df <- tibble(
  Category = factor(
    c("Before only", "Shared", "After only"),
    levels = c("Before only", "Shared", "After only")
  ),
  Count = c(before_only, shared, after_only)
) %>%
  mutate(
    xmin = lag(cumsum(Count), default = 0),
    xmax = cumsum(Count),
    xmid = (xmin + xmax) / 2,
    Percent_of_union = Count / union_total * 100,
    label = paste0(
      comma(Count),
      "\n",
      sprintf("%.1f%%", Percent_of_union)
    )
  )

category_colours <- c(
  "Before only" = "#caeacb",
  "Shared"      = "#95a1ae",
  "After only"  = "#cbcfd6"
)

# Before set spans: before-only + shared
before_xmin <- 0
before_xmax <- before_only + shared

# After set spans: shared + after-only
after_xmin <- before_only
after_xmax <- union_total

# ---------------------------
# Nature-style overlap plot
# ---------------------------
p <- ggplot() +
  
  # Main composition bar
  geom_rect(
    data = plot_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = 0.78,
      ymax = 1.22,
      fill = Category
    ),
    colour = "white",
    linewidth = 0.4
  ) +
  
  # Segment labels
  geom_text(
    data = plot_df,
    aes(
      x = xmid,
      y = 1,
      label = label
    ),
    size = 3.1,
    lineheight = 0.9,
    colour = "white"
  ) +
  
  # Before-adjustment bracket
  annotate(
    "segment",
    x = before_xmin,
    xend = before_xmax,
    y = 0.58,
    yend = 0.58,
    linewidth = 0.45
  ) +
  annotate(
    "segment",
    x = before_xmin,
    xend = before_xmin,
    y = 0.54,
    yend = 0.62,
    linewidth = 0.45
  ) +
  annotate(
    "segment",
    x = before_xmax,
    xend = before_xmax,
    y = 0.54,
    yend = 0.62,
    linewidth = 0.45
  ) +
  annotate(
    "text",
    x = (before_xmin + before_xmax) / 2,
    y = 0.45,
    label = paste0("Before adjustment: ", comma(before_total)),
    size = 3
  ) +
  
  # After-adjustment bracket
  annotate(
    "segment",
    x = after_xmin,
    xend = after_xmax,
    y = 0.25,
    yend = 0.25,
    linewidth = 0.45
  ) +
  annotate(
    "segment",
    x = after_xmin,
    xend = after_xmin,
    y = 0.21,
    yend = 0.29,
    linewidth = 0.45
  ) +
  annotate(
    "segment",
    x = after_xmax,
    xend = after_xmax,
    y = 0.21,
    yend = 0.29,
    linewidth = 0.45
  ) +
  annotate(
    "text",
    x = (after_xmin + after_xmax) / 2,
    y = 0.12,
    label = paste0("After adjustment: ", comma(after_total)),
    size = 3
  ) +
  
  scale_fill_manual(
    values = category_colours,
    name = NULL
  ) +
  
  scale_x_continuous(
    limits = c(0, union_total),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1.42),
    expand = c(0, 0)
  ) +
  
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  
  theme_void(base_size = 10) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(
      size = 8.5,
      colour = "black"
    ),
    legend.key.size = unit(0.32, "cm"),
    plot.margin = margin(3, 5, 3, 5)
  )



ggsave(
  "before_after_adjustment_overlap_composition.pdf",
  plot = p,
  width = 4,
  height = 2.2,
)
