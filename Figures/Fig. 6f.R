library(tidyverse)
library(scales)
library(grid)

theme_set(theme_classic(base_size = 12))

# ============================================================
# 1. Input data
# ============================================================

df2 <- tibble(
  group = c(
    "Clinical or pharmacogenomic\nvar-CpGs",
    "Other var-CpGs"
  ),
  count = c(
    3188,
    251228 - 3188
  ),
  total = c(
    30953,
    2816359 - 30953
  )
) %>%
  mutate(
    other = total - count,
    proportion = count / total,
    percent = proportion * 100,
    
    group = factor(
      group,
      levels = c(
        "Other var-CpGs",
        "Clinical or pharmacogenomic\nvar-CpGs"
      )
    ),
    
    value_label = paste0(
      sprintf("%.1f%%", percent),
      "\n(",
      comma(count),
      "/",
      comma(total),
      ")"
    )
  )

# ============================================================
# 2. Fisher exact test and fold enrichment
# ============================================================

clinical_count <- df2$count[
  df2$group == "Clinical or pharmacogenomic\nvar-CpGs"
]

clinical_total <- df2$total[
  df2$group == "Clinical or pharmacogenomic\nvar-CpGs"
]

other_count <- df2$count[
  df2$group == "Other var-CpGs"
]

other_total <- df2$total[
  df2$group == "Other var-CpGs"
]

fisher_matrix <- matrix(
  c(
    clinical_count,
    clinical_total - clinical_count,
    other_count,
    other_total - other_count
  ),
  nrow = 2,
  byrow = TRUE
)

ft <- fisher.test(
  fisher_matrix,
  alternative = "two.sided"
)

fold_enrichment <- (
  clinical_count / clinical_total
) / (
  other_count / other_total
)

# ============================================================
# 3. Format P value
# ============================================================

format_p <- function(p) {
  
  if (is.na(p)) {
    return("P = NA")
  }
  
  if (p < 2.2e-16) {
    return("P < 2.2 × 10⁻¹⁶")
  }
  
  if (p < 0.001) {
    exponent <- floor(log10(p))
    coefficient <- p / 10^exponent
    
    superscript <- c(
      "0" = "⁰",
      "1" = "¹",
      "2" = "²",
      "3" = "³",
      "4" = "⁴",
      "5" = "⁵",
      "6" = "⁶",
      "7" = "⁷",
      "8" = "⁸",
      "9" = "⁹",
      "-" = "⁻"
    )
    
    exponent_text <- paste0(
      superscript[
        strsplit(as.character(exponent), "")[[1]]
      ],
      collapse = ""
    )
    
    return(
      paste0(
        "P = ",
        sprintf("%.2f", coefficient),
        " × 10",
        exponent_text
      )
    )
  }
  
  paste0("P = ", sprintf("%.3f", p))
}

stat_label <- paste0(
  "Fold enrichment = ",
  sprintf("%.2f", fold_enrichment),
  "; ",
  format_p(ft$p.value)
)

# ============================================================
# 4. Colours
# ============================================================

group_colours <- c(
  "Clinical or pharmacogenomic\nvar-CpGs" = "#9E1B32",
  "Other var-CpGs" = "#3B6F8F"
)

# ============================================================
# 5. Vertical proportion plot
# ============================================================

y_max <- max(df2$percent)
bracket_y <- y_max + 1.2
tick_height <- 0.35
label_y <- bracket_y + 0.65

p2 <- ggplot(
  df2,
  aes(
    x = group,
    y = percent,
    fill = group
  )
) +
  
  geom_col(
    width = 0.58,
    colour = "black",
    linewidth = 0.4,
    show.legend = FALSE
  ) +
  
  geom_text(
    aes(label = value_label),
    vjust = -0.35,
    size = 3.3,
    colour = "black",
    lineheight = 0.95
  ) +
  
  # Statistical bracket
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = bracket_y,
    yend = bracket_y,
    linewidth = 0.55,
    colour = "black"
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = bracket_y,
    yend = bracket_y - tick_height,
    linewidth = 0.55,
    colour = "black"
  ) +
  
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = bracket_y,
    yend = bracket_y - tick_height,
    linewidth = 0.55,
    colour = "black"
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = label_y,
    label = stat_label,
    size = 3.2,
    colour = "black"
  ) +
  
  scale_fill_manual(
    values = group_colours
  ) +
  
  scale_y_continuous(
    limits = c(0, label_y + 0.8),
    breaks = seq(0, 12, by = 2),
    labels = label_percent(
      scale = 1,
      accuracy = 1
    ),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = NULL,
    y = paste0(
      "Integrated-model-only population-specific\n",
      "var-CpGs (%)"
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      size = 10.5,
      colour = "black",
      lineheight = 0.95,
      margin = margin(t = 6)
    ),
    
    axis.text.y = element_text(
      size = 10.5,
      colour = "black"
    ),
    
    axis.title.y = element_text(
      size = 12,
      colour = "black",
      margin = margin(r = 7)
    ),
    
    axis.line = element_line(
      linewidth = 0.65,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.65,
      colour = "black"
    ),
    
    axis.ticks.length = unit(
      0.13,
      "cm"
    ),
    
    panel.grid = element_blank(),
    legend.position = "none",
    
    plot.margin = margin(
      t = 7,
      r = 7,
      b = 5,
      l = 5
    )
  )


# ============================================================
# 6. Save
# ============================================================

ggsave(
  filename = "Clinical_enrichment_proportion.pdf",
  plot = p2,
  width = 3,
  height = 3.2
)
