library(tidyverse)
library(scales)
library(grid)

# ============================================================
# 1. Read and merge data
# ============================================================
ml <- read.delim(
  "../../ML-coverage/1158/match.log",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  transmute(
    Sample = basename(V1),
    Coverage = suppressWarnings(as.numeric(V4))
  )

meth <- read.delim(
  "../../Align_methyaltion_files/1158/re_methyaltion_qc_1157.tsv",
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  transmute(
    Sample = basename(
      sub(
        "model\\.pbmm2\\.combined\\.bed$",
        "bam",
        File
      )
    ),
    Methylation = suppressWarnings(as.numeric(Median))
  )

data_raw <- inner_join(
  ml,
  meth,
  by = "Sample"
) %>%
  filter(
    is.finite(Coverage),
    is.finite(Methylation)
  )

cat("Coverage records:", nrow(ml), "\n")
cat("Methylation records:", nrow(meth), "\n")
cat("Matched records:", nrow(data_raw), "\n")

# Check for duplicated sample names
# duplicated_samples <- data_raw %>%
#   count(Sample) %>%
#   filter(n > 1)
# 
# if (nrow(duplicated_samples) > 0) {
#   warning(
#     nrow(duplicated_samples),
#     " sample names occur more than once after merging."
#   )
# }

# ============================================================
# 2. Read low-quality BAM list
# ============================================================
lowq <- read.table(
  "../../lowQ_1158_self.log",
  header = FALSE,
  stringsAsFactors = FALSE
)

lowq_samples <- basename(trimws(lowq$V1))

data_filtered <- data_raw %>%
  filter(!Sample %in% lowq_samples)

cat("Original samples:", nrow(data_raw), "\n")
cat("Filtered samples:", nrow(data_filtered), "\n")
cat("Removed samples:", nrow(data_raw) - nrow(data_filtered), "\n")

# ============================================================
# 3. Statistical formatting functions
# ============================================================
format_p <- function(p) {
  
  if (is.na(p)) {
    return("NA")
  }
  
  if (p < 2.2e-16) {
    return("< 2.2 × 10\u207B\u00B9\u2076")
  }
  
  if (p < 0.001) {
    
    exponent <- floor(log10(p))
    coefficient <- p / 10^exponent
    
    exponent_text <- chartr(
      "0123456789-",
      "⁰¹²³⁴⁵⁶⁷⁸⁹⁻",
      as.character(exponent)
    )
    
    return(
      paste0(
        sprintf("%.2f", coefficient),
        " × 10",
        exponent_text
      )
    )
  }
  
  sprintf("%.3f", p)
}

# ============================================================
# 4. plotting function
# ============================================================
make_methylation_coverage_plot <- function(
    df,
    outfile,
    panel_label = NULL,
    x_limits = NULL,
    y_limits = c(0, 100)) {
  
  if (nrow(df) < 3) {
    stop("At least three observations are required.")
  }
  
  cor_res <- cor.test(
    df$Coverage,
    df$Methylation,
    method = "spearman",
    exact = FALSE
  )
  
  rho_value <- unname(cor_res$estimate)
  p_value <- cor_res$p.value
  
  statistics_label <- paste0(
    "n = ", format(nrow(df), big.mark = ","),
    "\nSpearman’s \u03C1 = ", sprintf("%.2f", rho_value),
    "\nP ", ifelse(p_value < 2.2e-16, "", "= "),
    format_p(p_value)
  )
  
  # Determine plotting range when no common range is supplied
  if (is.null(x_limits)) {
    
    x_range <- range(
      df$Coverage,
      na.rm = TRUE
    )
    
    x_padding <- diff(x_range) * 0.04
    
    if (x_padding == 0) {
      x_padding <- 0.5
    }
    
    x_limits <- c(
      x_range[1] - x_padding,
      x_range[2] + x_padding
    )
  }
  
  p <- ggplot(
    df,
    aes(
      x = Coverage,
      y = Methylation
    )
  ) +
    
    geom_point(
      shape = 16,
      size = 1.35,
      colour = "#6baed6",
      alpha = 0.55
    ) +
    
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      span = 0.75,
      se = TRUE,
      linewidth = 0.55,
      colour = "black",
      fill = "grey75",
      alpha = 0.35
    ) +
    
    annotate(
      geom = "text",
      x = x_limits[1],
      y = y_limits[2],
      label = statistics_label,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      lineheight = 1.05,
      colour = "black",
      family = "sans"
    ) +
    
    scale_x_continuous(
      breaks = pretty_breaks(n = 5),
      labels = label_number(
        accuracy = 0.1,
        trim = TRUE
      ),
      expand = expansion(
        mult = c(0.015, 0.015)
      )
    ) +
    
    scale_y_continuous(
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(
        mult = c(0, 0.015)
      )
    ) +
    
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits,
      clip = "off"
    ) +
    
    labs(
      x = "Genome coverage (\u00D7)",
      y = "Median DNA methylation"
    ) +
    
    theme_classic(
      base_size = 12,
      base_family = "sans"
    ) +
    
    theme(
      axis.line = element_line(
        colour = "black",
        linewidth = 1
      ),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 1
      ),
      axis.ticks.length = unit(
        1.4,
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
        margin = margin(t = 4)
      ),
      axis.title.y = element_text(
        margin = margin(r = 4)
      ),
      panel.grid = element_blank(),
      plot.margin = margin(
        t = 5,
        r = 6,
        b = 5,
        l = 5,
        unit = "pt"
      )
    )
  
  # Optional panel label
  if (!is.null(panel_label)) {
    
    p <- p +
      annotate(
        geom = "text",
        x = -Inf,
        y = Inf,
        label = panel_label,
        hjust = -0.15,
        vjust = 1.25,
        fontface = "bold",
        size = 3.5,
        family = "sans"
      )
  }
  
  ggsave(
    filename = outfile,
    plot = p,
    width = 4,
    height = 3.2,
    device = "pdf"
  )
  
  return(p)
}

# ============================================================
# 5. Use identical axes in the original and filtered plots
# ============================================================
common_x_range <- range(
  data_raw$Coverage,
  na.rm = TRUE
)

x_padding <- diff(common_x_range) * 0.04

if (x_padding == 0) {
  x_padding <- 0.5
}

common_x_limits <- c(
  common_x_range[1] - x_padding,
  common_x_range[2] + x_padding
)

# ============================================================
# 6. Generate plots
# ============================================================


p_filtered <- make_methylation_coverage_plot(
  df = data_filtered,
  outfile = "Coverage_methylation_filtered.pdf",
  panel_label = "b",
  x_limits = common_x_limits
)
