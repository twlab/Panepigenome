library(ggplot2)
library(dplyr)
library(mgcv)
library(scales)
library(grid)

# ==================================================
# 1. Load data
# ==================================================
read_data <- function(file, group_name) {
  
  df <- read.table(
    file,
    header = FALSE,
    col.names = c("presence_count", "methylation")
  )
  
  df %>%
    transmute(
      presence_count = as.integer(presence_count),
      freq = presence_count / 440,
      meth = methylation,
      group = group_name
    ) %>%
    filter(
      is.finite(freq),
      is.finite(meth),
      between(presence_count, 0, 440),
      between(meth, 0, 100)
    )
}

data_ref <- read_data(
  "1per.ref.txt",
  "Reference"
)

data_nonref <- read_data(
  "1per.nonref.txt",
  "Non-reference"
)

plot_data <- bind_rows(
  data_ref,
  data_nonref
) %>%
  mutate(
    group = factor(
      group,
      levels = c("Reference", "Non-reference")
    )
  )

# ==================================================
# 2. Summarise methylation at each presence frequency
# ==================================================
summary_data <- plot_data %>%
  group_by(
    group,
    presence_count,
    freq
  ) %>%
  summarise(
    n = n(),
    q10 = quantile(
      meth,
      probs = 0.10,
      na.rm = TRUE
    ),
    q25 = quantile(
      meth,
      probs = 0.25,
      na.rm = TRUE
    ),
    median = median(
      meth,
      na.rm = TRUE
    ),
    q75 = quantile(
      meth,
      probs = 0.75,
      na.rm = TRUE
    ),
    q90 = quantile(
      meth,
      probs = 0.90,
      na.rm = TRUE
    ),
    mean = mean(
      meth,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(n >= 20)

# ==================================================
# 3. Spearman correlations
# ==================================================
get_spearman <- function(df) {
  
  test <- suppressWarnings(
    cor.test(
      df$freq,
      df$meth,
      method = "spearman",
      exact = FALSE
    )
  )
  
  tibble(
    rho = unname(test$estimate),
    p_value = test$p.value,
    n = nrow(df)
  )
}

stats_group <- plot_data %>%
  group_by(group) %>%
  group_modify(
    ~ get_spearman(.x)
  ) %>%
  ungroup()

print(stats_group)

# ==================================================
# 4. Format P values
# ==================================================
format_p <- function(p) {
  
  if (is.na(p)) {
    return("P = NA")
  }
  
  if (p < 2.2e-16) {
    return("P < 2.2 × 10⁻¹⁶")
  }
  
  if (p < 0.001) {
    
    exponent <- floor(log10(p))
    coefficient <- p / (10^exponent)
    
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
    
    exponent_characters <- strsplit(
      as.character(exponent),
      split = ""
    )[[1]]
    
    exponent_text <- paste0(
      superscript[exponent_characters],
      collapse = ""
    )
    
    return(
      sprintf(
        "P = %.2f × 10%s",
        coefficient,
        exponent_text
      )
    )
  }
  
  sprintf(
    "P = %.3f",
    p
  )
}

# ==================================================
# 5. Annotation labels
# ==================================================
annotation_df <- stats_group %>%
  mutate(
    p_label = vapply(
      p_value,
      format_p,
      FUN.VALUE = character(1)
    ),
    label = paste0(
      "ρ = ",
      sprintf("%.3f", rho),
      "\n",
      p_label,
      "\nn = ",
      comma(n)
    ),
    x = 0.03,
    y = 97
  )

print(annotation_df)

# ==================================================
# 6. Colours
# ==================================================
group_colours <- c(
  "Reference" = "#4D4D4D",
  "Non-reference" = "#7B3294"
)

# ==================================================
# 7. Plot
# ==================================================
p <- ggplot() +
  
  # 10th–90th percentile range
  geom_ribbon(
    data = summary_data,
    aes(
      x = freq,
      ymin = q10,
      ymax = q90,
      fill = group
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # Interquartile range
  geom_ribbon(
    data = summary_data,
    aes(
      x = freq,
      ymin = q25,
      ymax = q75,
      fill = group
    ),
    alpha = 0.32,
    colour = NA
  ) +
  
  # Median methylation
  geom_line(
    data = summary_data,
    aes(
      x = freq,
      y = median,
      colour = group
    ),
    linewidth = 0.75,
    lineend = "round"
  ) +
  
  # GAM fit calculated from the plotted 1% sample
  geom_smooth(
    data = plot_data,
    aes(
      x = freq,
      y = meth,
      colour = group
    ),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    linewidth = 0.65,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  
  geom_text(
    data = annotation_df,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3,
    lineheight = 1.05,
    colour = "black"
  ) +
  
  facet_wrap(
    ~group,
    nrow = 1
  ) +
  
  scale_colour_manual(
    values = group_colours,
    guide = "none"
  ) +
  
  scale_fill_manual(
    values = group_colours,
    guide = "none"
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    labels = percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0.01)
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(
      0,
      100,
      by = 25
    ),
    labels = function(x) {
      paste0(x, "%")
    },
    expand = expansion(
      mult = c(0, 0.01)
    )
  ) +
  
  labs(
    x = "CpG presence frequency",
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
    
    axis.ticks.length = unit(
      1.5,
      "mm"
    ),
    
    axis.text = element_text(
      colour = "black",
      size = 11
    ),
    
    axis.title = element_text(
      colour = "black",
      size = 12
    ),
    
    axis.title.x = element_text(
      margin = margin(
        t = 5
      )
    ),
    
    axis.title.y = element_text(
      margin = margin(
        r = 5
      )
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 11,
      face = "bold",
      colour = "black",
      margin = margin(
        b = 4
      )
    ),
    
    panel.spacing.x = unit(
      5,
      "mm"
    ),
    
    plot.margin = margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4,
      unit = "mm"
    )
  )


# ==================================================
# 8. Save figure
# ==================================================
ggsave(
  filename = "CpG_frequency_methylation_quantile.pdf",
  plot = p,
  width = 6.8,
  height = 3.3
)
