library(tidyverse)
library(scales)
library(grid)

theme_set(theme_classic(base_size = 12))

# ============================================================
# 1. Prepare data
# ============================================================

df <- tribble(
  ~Group, ~Category, ~Overlap, ~Fold, ~P,
  
  "CpG gain/loss", "Pathogenic",       7,     1.1200, 0.8,
  "CpG gain/loss", "Benign",           48087, 1.8000, 2e-16,
  "CpG gain/loss", "Risk",             50,    1.9300, 0.000336,
  "CpG gain/loss", "Protective",       10,    1.7700, 0.1815,
  "CpG gain/loss", "Clinical",         7,     1.7900, 2e-16,
  "CpG gain/loss", "Pharmacogenomic",  549,   1.7440, 2e-16,
  
  "CpG gain",      "Pathogenic",       2,     0.6885, 1,
  "CpG gain",      "Benign",           19457, 1.5660, 2e-16,
  "CpG gain",      "Risk",             19,    1.5790, 0.07,
  "CpG gain",      "Protective",       2,     0.7610, 1,
  "CpG gain",      "Clinical",         2,     1.5590, 2e-16,
  "CpG gain",      "Pharmacogenomic",  233,   1.5940, 9e-10,
  
  "CpG loss",      "Pathogenic",       5,     1.4690, 0.4,
  "CpG loss",      "Benign",           29440, 2.0220, 2e-16,
  "CpG loss",      "Risk",             31,    2.1980, 0.0004,
  "CpG loss",      "Protective",       8,     2.5980, 0.04,
  "CpG loss",      "Clinical",         5,     2.0140, 2e-16,
  "CpG loss",      "Pharmacogenomic",  317,   1.8500, 2e-16
) %>%
  mutate(
    FDR = p.adjust(P, method = "BH"),
    
    Significance = if_else(
      FDR < 0.05,
      "FDR < 0.05",
      "Not significant"
    ),
    
    Group = factor(
      Group,
      levels = c(
        "CpG gain/loss",
        "CpG gain",
        "CpG loss"
      )
    ),
    
    # Explicit y positions prevent category–label mismatch
    y_position = case_when(
      Category == "Pathogenic"      ~ 6,
      Category == "Benign"          ~ 5,
      Category == "Risk"            ~ 4,
      Category == "Protective"      ~ 3,
      Category == "Clinical"        ~ 1.8,
      Category == "Pharmacogenomic" ~ 0.6
    ),
    
    Significance = factor(
      Significance,
      levels = c(
        "FDR < 0.05",
        "Not significant"
      )
    )
  )

# Check that every category received a y position
stopifnot(!any(is.na(df$y_position)))

# ============================================================
# 2. Background regions
# ============================================================

background_df <- tibble(
  xmin = -Inf,
  xmax = Inf,
  ymin = 2.5,
  ymax = 6.5
)

# ============================================================
# 3. Plot
# ============================================================

p <- ggplot(
  df,
  aes(
    x = Fold,
    y = y_position
  )
) +
  
  # Background for clinical classes
  geom_rect(
    data = background_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = "grey96",
    colour = NA
  ) +
  
  # Separators between annotation groups
  geom_hline(
    yintercept = c(2.35, 1.2),
    colour = "grey75",
    linewidth = 0.45
  ) +
  
  # No-enrichment reference
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.45
  ) +
  
  geom_point(
    aes(
      size = Overlap,
      fill = Significance
    ),
    shape = 21,
    colour = "black",
    stroke = 0.35,
    alpha = 1
  ) +
  
  facet_wrap(
    ~Group,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = c(
      "FDR < 0.05" = "#7B3294",
      "Not significant" = "white"
    ),
    drop = FALSE
  ) +
  
  scale_size_continuous(
    trans = "log10",
    range = c(2.5, 7.5),
    breaks = c(
      2,
      10,
      100,
      1000,
      10000
    ),
    labels = comma
  ) +
  
  scale_y_continuous(
    breaks = c(
      6,
      5,
      4,
      3,
      1.8,
      0.6
    ),
    labels = c(
      "Pathogenic",
      "Benign",
      "Risk",
      "Protective",
      "Clinical",
      "Pharmacogenomic"
    ),
    limits = c(0.1, 6.5),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_x_continuous(
    breaks = c(
      0.5,
      1.0,
      1.5,
      2.0,
      2.5
    ),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  
  coord_cartesian(
    xlim = c(0.55, 2.70),
    clip = "off"
  ) +
  
  labs(
    x = "Fold enrichment",
    y = NULL,
    fill = NULL,
    size = "Overlapping variants"
  ) +
  
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        size = 4,
        stroke = 0.35
      )
    ),
    
    size = guide_legend(
      order = 2,
      title.position = "top",
      override.aes = list(
        shape = 21,
        fill = "grey70",
        colour = "black",
        stroke = 0.25
      )
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "black",
      margin = margin(r = 4)
    ),
    
    axis.title.x = element_text(
      size = 12,
      colour = "black",
      margin = margin(t = 7)
    ),
    
    axis.line = element_line(
      linewidth = 0.75,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.75,
      colour = "black"
    ),
    
    axis.ticks.length = unit(
      0.14,
      "cm"
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 11,
      face = "bold",
      colour = "black",
      margin = margin(b = 5)
    ),
    
    panel.grid = element_blank(),
    
    panel.spacing.x = unit(
      0.7,
      "lines"
    ),
    
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    
    legend.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    legend.title = element_text(
      size = 11,
      colour = "black"
    ),
    
    legend.key.height = unit(
      0.47,
      "cm"
    ),
    
    legend.key.width = unit(
      0.47,
      "cm"
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 4,
      l = 5
    )
  )

# p

# ============================================================
# 4. Save
# ============================================================

ggsave(
  filename = "clinical_variant_enrichment.pdf",
  plot = p,
  width = 7,
  height = 3.2
)
