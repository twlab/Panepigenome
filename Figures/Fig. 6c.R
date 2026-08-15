library(tidyverse)
library(scales)
library(grid)

theme_set(theme_classic(base_size = 12))

# ============================================================
# 1. Prepare data
# ============================================================

df <- tribble(
  ~Group, ~Context, ~Count, ~Total,
  
  "Clinically associated", "Within CGIs",
  44, 613,
  
  "Clinically associated", "Within promoters",
  156, 613,
  
  "Clinically associated", "Non-CGI and non-promoter",
  445, 613,
  
  "Benign", "Within CGIs",
  6447, 49389,
  
  "Benign", "Within promoters",
  14366, 49389,
  
  "Benign", "Non-CGI and non-promoter",
  32384, 49389
) %>%
  mutate(
    Proportion = Count / Total,
    
    Group = factor(
      Group,
      levels = c(
        "Clinically associated",
        "Benign"
      )
    ),
    
    Context = factor(
      Context,
      levels = c(
        "Non-CGI and non-promoter",
        "Within promoters",
        "Within CGIs"
      )
    ),
    
    count_label = comma(Count),
    
    percentage_label = percent(
      Proportion,
      accuracy = 0.1
    )
  )

# Check calculated values
print(
  df %>%
    select(
      Group,
      Context,
      Count,
      Total,
      Proportion
    )
)

# ============================================================
# 2. Facet labels including total numbers
# ============================================================

group_labels <- c(
  "Clinically associated" =
    "Clinically associated variants\n(n = 613 var-CpGs)",
  
  "Benign" =
    "Benign variants\n(n = 49,389 var-CpGs)"
)

# ============================================================
# 3. Colours
# ============================================================

group_colours <- c(
  "Clinically associated" = "#9E1B32",
  "Benign" = "#3B6F8F"
)

# ============================================================
# 4. lollipop plot
# ============================================================

p <- ggplot(
  df,
  aes(
    x = Proportion,
    y = Context,
    colour = Group
  )
) +
  
  # Line from zero to each observed proportion
  geom_segment(
    aes(
      x = 0,
      xend = Proportion,
      y = Context,
      yend = Context
    ),
    linewidth = 1.3,
    lineend = "round",
    show.legend = FALSE
  ) +
  
  # Endpoint
  geom_point(
    size = 4.2,
    stroke = 0,
    show.legend = FALSE
  ) +
  
  # Percentage inside or immediately adjacent to point
  geom_text(
    aes(
      label = percentage_label
    ),
    hjust = -0.35,
    size = 3.4,
    colour = "black",
    show.legend = FALSE
  ) +
  
  # Original count shown in a right-side aligned column
  geom_text(
    aes(
      x = 0.79,
      label = count_label
    ),
    hjust = 1,
    size = 3.35,
    colour = "black",
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~Group,
    nrow = 2,
    labeller = as_labeller(group_labels)
  ) +
  
  scale_colour_manual(
    values = group_colours
  ) +
  
  scale_x_continuous(
    limits = c(0, 0.82),
    breaks = seq(
      0,
      0.8,
      by = 0.2
    ),
    labels = percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  labs(
    x = "Proportion of var-CpGs",
    y = NULL
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    # ------------------------
    # Axis
    # ------------------------
    axis.text.x = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "black",
      margin = margin(r = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      colour = "black",
      margin = margin(t = 7)
    ),
    
    axis.line.x = element_line(
      linewidth = 0.6,
      colour = "black"
    ),
    
    axis.line.y = element_blank(),
    
    axis.ticks.x = element_line(
      linewidth = 0.6,
      colour = "black"
    ),
    
    axis.ticks.y = element_blank(),
    
    axis.ticks.length = unit(
      0.14,
      "cm"
    ),
    
    # ------------------------
    # Facet labels
    # ------------------------
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 11.5,
      face = "bold",
      colour = "black",
      hjust = 0,
      margin = margin(
        b = 6
      )
    ),
    
    # ------------------------
    # Panels
    # ------------------------
    panel.grid.major.x = element_line(
      linewidth = 0.3,
      colour = "grey90"
    ),
    
    panel.grid.major.y = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    panel.spacing.y = unit(
      0.65,
      "lines"
    ),
    
    # ------------------------
    # Remove legend
    # ------------------------
    legend.position = "none",
    
    # ------------------------
    # Plot margin
    # ------------------------
    plot.margin = margin(
      t = 5,
      r = 26,
      b = 5,
      l = 5
    )
  )


# ============================================================
# 5. Save
# ============================================================

ggsave(
  filename =
    "clinical_benign_varCpG_genomic_context.pdf",
  plot = p,
  width = 5,
  height = 3.2
)
