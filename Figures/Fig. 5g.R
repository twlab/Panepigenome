library(data.table)
library(tidyverse)
library(patchwork)
library(scales)
library(ggplot2)

theme_set(
  theme_classic(
    base_size = 12
  )
)

# ---------------------------
# Read data
# ---------------------------
re <- fread(
  "../re_hierarchical_fdr_ACAT_4.txt"
)

# ---------------------------
# Keep significant interaction pairs
# ---------------------------
sig <- re %>%
  filter(
    !is.na(I_inter_p_FDR),
    I_inter_p_FDR < 0.05
  ) %>%
  mutate(
    pair_id = row_number()
  )

# ---------------------------
# Paired Wilcoxon tests
# ---------------------------
p_geno_signed <- wilcox.test(
  sig$Y_geno_beta,
  sig$I_geno_beta,
  paired = TRUE
)$p.value

p_geno_abs <- wilcox.test(
  abs(sig$Y_geno_beta),
  abs(sig$I_geno_beta),
  paired = TRUE
)$p.value

p_meth_signed <- wilcox.test(
  sig$Y_meth_beta,
  sig$I_meth_beta,
  paired = TRUE
)$p.value

p_meth_abs <- wilcox.test(
  abs(sig$Y_meth_beta),
  abs(sig$I_meth_beta),
  paired = TRUE
)$p.value

# ---------------------------
# Helper function to prepare data
# ---------------------------
prepare_pair_df <- function(
    sig,
    before_col,
    after_col,
    use_abs = FALSE
) {
  
  before <- sig[[before_col]]
  after  <- sig[[after_col]]
  
  if (use_abs) {
    before <- abs(before)
    after  <- abs(after)
  }
  
  tibble(
    pair_id = sig$pair_id,
    before = before,
    after = after
  ) %>%
    pivot_longer(
      cols = c(before, after),
      names_to = "Model",
      values_to = "Beta"
    ) %>%
    mutate(
      Model = factor(
        Model,
        levels = c(
          "before",
          "after"
        ),
        labels = c(
          "Before\nadjustment",
          "After\nadjustment"
        )
      )
    )
}

# ---------------------------
# Prepare four datasets
# ---------------------------
geno_signed_df <- prepare_pair_df(
  sig = sig,
  before_col = "Y_geno_beta",
  after_col = "I_geno_beta",
  use_abs = FALSE
)

geno_abs_df <- prepare_pair_df(
  sig = sig,
  before_col = "Y_geno_beta",
  after_col = "I_geno_beta",
  use_abs = TRUE
)

meth_signed_df <- prepare_pair_df(
  sig = sig,
  before_col = "Y_meth_beta",
  after_col = "I_meth_beta",
  use_abs = FALSE
)

meth_abs_df <- prepare_pair_df(
  sig = sig,
  before_col = "Y_meth_beta",
  after_col = "I_meth_beta",
  use_abs = TRUE
)

# ---------------------------
# P-value formatter
# ---------------------------
format_p <- function(p) {
  if (is.na(p)) {
    return("P = NA")
  }
  
  if (p < 2.2e-16) {
    return("P < 2.2 × 10\u207b\u00b9\u2076")
  }
  
  paste0(
    "P = ",
    formatC(
      p,
      format = "e",
      digits = 2
    )
  )
}

# ---------------------------
# Plot function
# ---------------------------
make_violin_paired_plot <- function(
    df,
    ylab,
    pval,
    use_hline = TRUE,
    show_x_text = TRUE
) {
  
  df <- df %>%
    filter(
      !is.na(Beta),
      !is.na(Model),
      !is.na(pair_id)
    ) %>%
    mutate(
      Model = factor(
        Model,
        levels = c(
          "Before\nadjustment",
          "After\nadjustment"
        )
      )
    )
  
  y_max <- max(
    df$Beta,
    na.rm = TRUE
  )
  
  y_min <- min(
    df$Beta,
    na.rm = TRUE
  )
  
  y_range <- y_max - y_min
  
  if (!is.finite(y_range) || y_range == 0) {
    y_range <- max(
      abs(c(y_min, y_max)),
      na.rm = TRUE
    )
  }
  
  if (!is.finite(y_range) || y_range == 0) {
    y_range <- 1
  }
  
  y_label <- y_max + 0.14 * y_range
  
  x_theme <- if (show_x_text) {
    theme(
      axis.text.x = element_text(
        size = 11,
        colour = "black"
      ),
      axis.ticks.x = element_line(
        linewidth = 0.35,
        colour = "black"
      )
    )
  } else {
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  
  ggplot(
    df,
    aes(
      x = Model,
      y = Beta,
      fill = Model
    )
  ) +
    
    # Paired lines
    geom_line(
      aes(
        group = pair_id
      ),
      colour = "grey75",
      linewidth = 0.22,
      alpha = 0.35
    ) +
    
    # Violin distributions
    geom_violin(
      width = 0.82,
      linewidth = 0.30,
      colour = "black",
      alpha = 0.95,
      trim = FALSE,
      scale = "width"
    ) +
    
    # Compact internal boxplots
    geom_boxplot(
      width = 0.10,
      outlier.shape = NA,
      fill = "white",
      colour = "black",
      linewidth = 0.28
    ) +
    
    # Individual observations
    geom_point(
      position = position_jitter(
        width = 0.045,
        height = 0,
        seed = 1
      ),
      shape = 21,
      size = 1.15,
      colour = "black",
      stroke = 0.18,
      alpha = 0.65
    ) +
    
    # Median marker
    stat_summary(
      fun = median,
      geom = "point",
      shape = 95,
      size = 5.5,
      colour = "black"
    ) +
    
    {
      if (use_hline) {
        geom_hline(
          yintercept = 0,
          linetype = "dashed",
          colour = "grey45",
          linewidth = 0.35
        )
      }
    } +
    
    # Paired Wilcoxon P value
    annotate(
      "text",
      x = 1.5,
      y = y_label,
      label = format_p(pval),
      size = 3.9,
      colour = "black"
    ) +
    
    scale_fill_manual(
      values = c(
        "Before\nadjustment" = "#CAEACB",
        "After\nadjustment" = "#CBCFD6"
      )
    ) +
    
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.08,
          0.22
        )
      )
    ) +
    
    labs(
      x = NULL,
      y = ylab,
      title = NULL
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_blank(),
      
      axis.text.y = element_text(
        size = 11,
        colour = "black"
      ),
      
      axis.title.y = element_text(
        size = 12,
        colour = "black",
        margin = margin(
          r = 6
        )
      ),
      
      axis.line = element_line(
        linewidth = 0.35,
        colour = "black"
      ),
      
      axis.ticks.y = element_line(
        linewidth = 0.35,
        colour = "black"
      ),
      
      legend.position = "none",
      
      panel.grid = element_blank(),
      
      plot.margin = margin(
        t = 3,
        r = 7,
        b = 3,
        l = 5
      )
    ) +
    
    x_theme
}

# ---------------------------
# Create four panels
# Only the fourth panel retains x-axis labels
# ---------------------------
p1 <- make_violin_paired_plot(
  df = geno_signed_df,
  ylab = "Genetic effect size",
  pval = p_geno_signed,
  use_hline = TRUE,
  show_x_text = FALSE
)

p2 <- make_violin_paired_plot(
  df = geno_abs_df,
  ylab = expression("|Genetic effect size|"),
  pval = p_geno_abs,
  use_hline = FALSE,
  show_x_text = FALSE
)

p3 <- make_violin_paired_plot(
  df = meth_signed_df,
  ylab = "Methylation effect size",
  pval = p_meth_signed,
  use_hline = TRUE,
  show_x_text = FALSE
)

p4 <- make_violin_paired_plot(
  df = meth_abs_df,
  ylab = expression("|Methylation effect size|"),
  pval = p_meth_abs,
  use_hline = FALSE,
  show_x_text = TRUE
)

# ---------------------------
# Combine vertically
# ---------------------------
p <- p1 / p2 / p3 / p4 +
  plot_layout(
    ncol = 1,
    heights = c(
      1,
      1,
      1,
      1
    )
  ) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(
        2,
        2,
        2,
        2
      )
    )
  )


# ---------------------------
# Save
# ---------------------------
ggsave(
  filename =
    "interaction_adjusted_effect_sizes_four_panels_violin_vertical.pdf",
  plot = p,
  width = 3,
  height = 8
)
