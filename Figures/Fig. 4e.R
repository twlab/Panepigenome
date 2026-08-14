library(data.table)
library(ggplot2)
library(scales)
library(viridis)
library(patchwork)
library(grid)
library(hexbin)

# ============================================================
# 1. General settings
# ============================================================

groups <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# Absolute deviation threshold:
# delta = expected PST - observed PST
delta_threshold <- 0.05

outlier_colours <- c(
  "Lower than expected PST"  = "#756BB1",
  "Higher than expected PST" = "#31A354"
)

group_labels <- c(
  "AFR" = "AFR",
  "AMR" = "AMR",
  "EAS" = "EAS",
  "EUR" = "EUR",
  "SAS" = "SAS"
)

# ============================================================
# 2. Function to read, merge and analyse one group
# ============================================================

prepare_group_data <- function(group) {
  
  message("Processing ", group, "...")
  
  fst_file <- paste0(
    "../PCA_Nonref_geno/Pop_con_vcftools/",
    group,
    "_vs_rest.weir.fst"
  )
  
  pst_file <- paste0(
    "../Pst_Non_ae/Pst_con_beta/",
    group,
    "_sorted_filtered.txt"
  )
  
  # ----------------------------------------------------------
  # Check files
  # ----------------------------------------------------------
  
  if (!file.exists(fst_file)) {
    stop("FST file does not exist: ", fst_file)
  }
  
  if (!file.exists(pst_file)) {
    stop("PST file does not exist: ", pst_file)
  }
  
  # ----------------------------------------------------------
  # Read FST
  # ----------------------------------------------------------
  
  fst_dt <- fread(fst_file)
  
  required_fst_columns <- c(
    "CHROM",
    "POS",
    "WEIR_AND_COCKERHAM_FST"
  )
  
  missing_fst_columns <- setdiff(
    required_fst_columns,
    names(fst_dt)
  )
  
  if (length(missing_fst_columns) > 0) {
    stop(
      "Missing FST columns in ",
      fst_file,
      ": ",
      paste(missing_fst_columns, collapse = ", ")
    )
  }
  
  fst_dt[, V1 := paste(CHROM, POS, sep = ":")]
  
  fst_dt <- fst_dt[, .(
    V1,
    FST = as.numeric(WEIR_AND_COCKERHAM_FST)
  )]
  
  # Replace negative FST estimates with zero
  fst_dt[
    is.finite(FST) & FST < 0,
    FST := 0
  ]
  
  # ----------------------------------------------------------
  # Read PST
  # ----------------------------------------------------------
  
  pst_dt <- fread(pst_file)
  
  if (!all(c("V1", "V2") %in% names(pst_dt))) {
    stop(
      "The PST file must contain columns V1 and V2: ",
      pst_file
    )
  }
  
  # Retain chromosome and position from the identifier
  pst_dt[, V1 := sub(
    "^([^:]+:[^:]+).*",
    "\\1",
    V1
  )]
  
  pst_dt[, PST := as.numeric(V2)]
  
  pst_dt <- pst_dt[, .(
    V1,
    PST
  )]
  
  # ----------------------------------------------------------
  # Merge FST and PST
  # ----------------------------------------------------------
  
  dt <- merge(
    fst_dt,
    pst_dt,
    by = "V1",
    all = FALSE,
    allow.cartesian = TRUE
  )
  
  dt <- dt[
    is.finite(FST) &
      is.finite(PST) &
      FST >= 0 &
      PST >= 0
  ]
  
  if (nrow(dt) < 3) {
    warning(
      "Too few valid loci for ",
      group,
      ". This group will be omitted."
    )
    
    return(NULL)
  }
  
  # ----------------------------------------------------------
  # Fit PST ~ FST
  # ----------------------------------------------------------
  
  lm_fit <- lm(
    PST ~ FST,
    data = dt
  )
  
  lm_summary <- summary(lm_fit)
  
  intercept <- unname(
    coef(lm_fit)[["(Intercept)"]]
  )
  
  slope <- unname(
    coef(lm_fit)[["FST"]]
  )
  
  pearson_r <- cor(
    dt$FST,
    dt$PST,
    method = "pearson",
    use = "complete.obs"
  )
  
  r_squared <- lm_summary$r.squared
  
  # ----------------------------------------------------------
  # Expected PST from the fitted model
  # ----------------------------------------------------------
  
  dt[, expected_PST := as.numeric(
    predict(
      lm_fit,
      newdata = data.frame(FST = FST)
    )
  )]
  
  # Positive delta: observed PST is lower than expected
  # Negative delta: observed PST is higher than expected
  dt[, delta := expected_PST - PST]
  
  dt[, outlier_class := fifelse(
    delta > delta_threshold,
    "Lower than expected PST",
    fifelse(
      delta < -delta_threshold,
      "Higher than expected PST",
      "Non-outlier"
    )
  )]
  
  dt[, outlier_class := factor(
    outlier_class,
    levels = c(
      "Lower than expected PST",
      "Higher than expected PST",
      "Non-outlier"
    )
  )]
  
  dt[, Group := group]
  
  # ----------------------------------------------------------
  # Save all group-specific outliers
  # ----------------------------------------------------------
  
  outlier_dt <- dt[
    outlier_class != "Non-outlier"
  ]
  
  fwrite(
    outlier_dt,
    paste0(
      "outliers_",
      group,
      "_delta",
      delta_threshold,
      ".tsv"
    ),
    sep = "\t",
    na = "NA"
  )
  
  # ----------------------------------------------------------
  # Regression line
  # ----------------------------------------------------------
  
  regression_grid <- data.table(
    FST = seq(
      min(dt$FST, na.rm = TRUE),
      max(dt$FST, na.rm = TRUE),
      length.out = 400
    )
  )
  
  regression_grid[, PST_fit := as.numeric(
    predict(
      lm_fit,
      newdata = data.frame(FST = FST)
    )
  )]
  
  regression_grid[, Group := group]
  
  # ----------------------------------------------------------
  # Summary information
  # ----------------------------------------------------------
  
  n_lower <- dt[
    outlier_class == "Lower than expected PST",
    .N
  ]
  
  n_higher <- dt[
    outlier_class == "Higher than expected PST",
    .N
  ]
  
  n_non_outlier <- dt[
    outlier_class == "Non-outlier",
    .N
  ]
  
  model_summary <- data.table(
    Group = group,
    n_total = nrow(dt),
    n_non_outlier = n_non_outlier,
    n_lower_than_expected = n_lower,
    n_higher_than_expected = n_higher,
    intercept = intercept,
    slope = slope,
    pearson_r = pearson_r,
    r_squared = r_squared
  )
  
  message(
    group,
    ": ",
    format(nrow(dt), big.mark = ","),
    " total loci; ",
    format(n_non_outlier, big.mark = ","),
    " non-outliers; ",
    format(n_lower, big.mark = ","),
    " lower-than-expected; ",
    format(n_higher, big.mark = ","),
    " higher-than-expected."
  )
  
  list(
    data = dt,
    regression = regression_grid,
    summary = model_summary
  )
}

# ============================================================
# 3. Analyse all continental groups
# ============================================================

results <- lapply(
  groups,
  prepare_group_data
)

names(results) <- groups

# Remove groups that could not be analysed
results <- results[
  !vapply(results, is.null, logical(1))
]

if (length(results) == 0) {
  stop("No group contained sufficient data for plotting.")
}

# ============================================================
# 4. Combine results
# ============================================================

all_data <- rbindlist(
  lapply(results, `[[`, "data"),
  use.names = TRUE,
  fill = TRUE
)

regression_data <- rbindlist(
  lapply(results, `[[`, "regression"),
  use.names = TRUE,
  fill = TRUE
)

model_summary <- rbindlist(
  lapply(results, `[[`, "summary"),
  use.names = TRUE,
  fill = TRUE
)

all_data[, Group := factor(
  Group,
  levels = groups
)]

regression_data[, Group := factor(
  Group,
  levels = groups
)]

model_summary[, Group := factor(
  Group,
  levels = groups
)]

# ============================================================
# 5. Save model statistics
# ============================================================

fwrite(
  model_summary,
  "FST_PST_regression_summary.tsv",
  sep = "\t",
  na = "NA"
)

print(model_summary)

# ============================================================
# 6. Separate non-outliers and outliers
# ============================================================

background_data <- all_data[
  outlier_class == "Non-outlier"
]

outlier_data <- all_data[
  outlier_class != "Non-outlier"
]

# Draw the more abundant outlier category first within each group
# so that the less abundant category remains visible on top.
outlier_counts <- outlier_data[
  ,
  .N,
  by = .(
    Group,
    outlier_class
  )
]

outlier_data <- merge(
  outlier_data,
  outlier_counts,
  by = c(
    "Group",
    "outlier_class"
  ),
  all.x = TRUE
)

# Larger categories are drawn first; smaller categories are drawn last
setorder(
  outlier_data,
  Group,
  -N
)

# ============================================================
# 7. Common plotting range
# ============================================================

common_max <- max(
  c(
    all_data$FST,
    all_data$PST
  ),
  na.rm = TRUE
)

if (!is.finite(common_max) || common_max <= 0) {
  common_max <- 1
}

# Round upward to a clear axis limit
common_max <- ceiling(common_max * 10) / 10

message("Common FST/PST axis maximum: ", common_max)

# ============================================================
# 8. Panel annotations
# ============================================================

annotation_data <- copy(model_summary)

annotation_data[, annotation_label := paste0(
  "Slope = ",
  sprintf("%.2f", slope),
  "\nr = ",
  sprintf("%.2f", pearson_r)
)]

# Bottom-right position in each panel
annotation_data[, x := common_max * 0.96]
annotation_data[, y := common_max * 0.04]

# ============================================================
# 9. Nature-style multi-panel plot
# ============================================================

p <- ggplot() +
  
  # ----------------------------------------------------------
# Hexbin density: non-outliers only
# ----------------------------------------------------------

geom_hex(
  data = background_data,
  aes(
    x = FST,
    y = PST
  ),
  bins = 65,
  linewidth = 0
) +
  
  scale_fill_viridis(
    option = "mako",
    direction = 1,
    trans = "log10",
    name = "Non-outlier loci",
    limits = c(1, 10000),
    breaks = c(
      1,
      10,
      100,
      1000,
      10000
    ),
    labels = label_number(
      scale_cut = cut_short_scale()
    ),
    oob = squish
  ) +
  
  # ----------------------------------------------------------
# FST = PST reference
# ----------------------------------------------------------

geom_abline(
  slope = 1,
  intercept = 0,
  linetype = "dashed",
  linewidth = 0.40,
  colour = "grey45"
) +
  
  # ----------------------------------------------------------
# Linear regression fit
# ----------------------------------------------------------

geom_line(
  data = regression_data,
  aes(
    x = FST,
    y = PST_fit,
    group = Group
  ),
  colour = "black",
  linewidth = 0.65,
  lineend = "round"
) +
  
  # ----------------------------------------------------------
# Outliers
# ----------------------------------------------------------

geom_point(
  data = outlier_data,
  aes(
    x = FST,
    y = PST,
    colour = outlier_class
  ),
  size = 0.28,
  alpha = 0.45,
  stroke = 0,
  show.legend = TRUE
) +
  
  scale_colour_manual(
    values = outlier_colours,
    breaks = c(
      "Lower than expected PST",
      "Higher than expected PST"
    ),
    name = NULL,
    drop = FALSE
  ) +
  
  # ----------------------------------------------------------
# Per-panel statistics
# ----------------------------------------------------------

geom_text(
  data = annotation_data,
  aes(
    x = x,
    y = y,
    label = annotation_label
  ),
  inherit.aes = FALSE,
  hjust = 1,
  vjust = 0,
  size = 3.15,
  lineheight = 1.05,
  colour = "black"
) +
  
  # ----------------------------------------------------------
# Facets
# ----------------------------------------------------------

facet_wrap(
  ~Group,
  ncol = 3,
  labeller = as_labeller(group_labels)
) +
  
  # ----------------------------------------------------------
# Shared coordinate scales
# ----------------------------------------------------------

scale_x_continuous(
  limits = c(0, common_max),
  breaks = pretty_breaks(n = 4),
  expand = expansion(
    mult = c(0, 0)
  )
) +
  
  scale_y_continuous(
    limits = c(0, common_max),
    breaks = pretty_breaks(n = 4),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  coord_fixed(
    ratio = 1,
    xlim = c(0, common_max),
    ylim = c(0, common_max),
    clip = "on"
  ) +
  
  labs(
    x = expression(F[ST]),
    y = expression(P[ST])
  ) +
  
  guides(
    fill = guide_colourbar(
      order = 1,
      title.position = "top",
      title.hjust = 0,
      barwidth = unit(3, "cm"),
      barheight = unit(0.25, "cm")
    ),
    
    colour = guide_legend(
      order = 2,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        size = 2.2,
        alpha = 1
      )
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    # --------------------------------------------------------
    # Axes
    # --------------------------------------------------------
    
    axis.title = element_text(
      size = 13,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      margin = margin(t = 6)
    ),
    
    axis.title.y = element_text(
      margin = margin(r = 6)
    ),
    
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.line = element_line(
      linewidth = 1,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 1,
      colour = "black"
    ),
    
    axis.ticks.length = unit(
      0.11,
      "cm"
    ),
    
    # --------------------------------------------------------
    # Facet labels
    # --------------------------------------------------------
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 12,
      face = "bold",
      colour = "black",
      hjust = 0.5,
      margin = margin(
        b = 4
      )
    ),
    
    # --------------------------------------------------------
    # Panels
    # --------------------------------------------------------
    
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    ),
    
    panel.spacing.x = unit(
      2.4,
      "lines"
    ),
    
    panel.spacing.y = unit(
      2,
      "lines"
    ),
    
    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------
    
    legend.position = "bottom",
    
    legend.direction = "horizontal",
    
    legend.box = "horizontal",
    
    legend.box.just = "right",
    
    legend.justification = "right",
    
    legend.title = element_text(
      size = 12,
      colour = "black"
    ),
    
    legend.text = element_text(
      size = 12,
      colour = "black"
    ),
    
    legend.key.height = unit(
      0.30,
      "cm"
    ),
    
    legend.key.width = unit(
      0.45,
      "cm"
    ),
    
    legend.spacing.x = unit(
      1,
      "cm"
    ),
    
    legend.box.spacing = unit(
      0.10,
      "cm"
    ),
    
    legend.margin = margin(
      b = 3
    ),
    
    # --------------------------------------------------------
    # Plot margins
    # --------------------------------------------------------
    
    plot.margin = margin(
      t = 5,
      r = 6,
      b = 5,
      l = 5
    )
  )


# ============================================================
# 10. Save multi-panel main figure
# ============================================================

ggsave(
  filename = "FST_PST_all_continental_groups.pdf",
  plot = p,
  width = 18,
  height = 8.5
)


