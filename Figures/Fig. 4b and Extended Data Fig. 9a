library(data.table)
library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 12))

# ---------------------------
# Read data
# ---------------------------
df <- fread(
  "enrichment_results_fdr_model1.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

# ---------------------------
# Clean numeric columns
# ---------------------------
df <- df %>%
  mutate(
    shared1_n = as.numeric(shared1_n),
    backg_n = as.numeric(backg_n),
    shared1_all = as.numeric(shared1_all),
    backg_all = as.numeric(backg_all),
    odds_ratio = as.numeric(odds_ratio),
    p_value = as.numeric(p_value),
    fdr = as.numeric(fdr)
  ) %>%
  filter(
    !is.na(pair),
    !is.na(type)
  )

# ---------------------------
# Merge duplicated features and recompute enrichment
# ---------------------------
df <- df %>%
  mutate(
    feature = case_when(
      type == "1_Active_Promoter" ~ "Active promoter",
      type == "2_Weak_Promoter" ~ "Weak promoter",
      type == "3_Poised_Promoter" ~ "Poised promoter",
      type %in% c("4_Strong_Enhancer", "5_Strong_Enhancer") ~ "Strong enhancer",
      type %in% c("6_Weak_Enhancer", "7_Weak_Enhancer") ~ "Weak enhancer",
      type == "8_Insulator" ~ "Insulator",
      type == "9_Txn_Transition" ~ "Txn transition",
      type == "10_Txn_Elongation" ~ "Txn elongation",
      type == "11_Weak_Txn" ~ "Weak transcription",
      type == "12_Repressed" ~ "Repressed",
      type == "13_Heterochrom/lo" ~ "Heterochromatin",
      type %in% c("14_Repetitive/CNV", "15_Repetitive/CNV") ~ "Repetitive/CNV",
      type == "exon" ~ "Exon",
      type == "intron" ~ "Intron",
      type == "SE" ~ "Super-enhancer",
      TRUE ~ type
    )
  ) %>%
  group_by(pair, feature) %>%
  summarise(
    shared1_n = sum(shared1_n, na.rm = TRUE),
    backg_n = sum(backg_n, na.rm = TRUE),
    shared1_all = first(shared1_all),
    backg_all = first(backg_all),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(
      matrix(
        c(
          shared1_n,
          shared1_all - shared1_n,
          backg_n,
          backg_all - backg_n
        ),
        nrow = 2,
        byrow = TRUE
      )
    )$p.value,
    odds_ratio = {
      a <- matrix(
        c(
          shared1_n,
          shared1_all - shared1_n,
          backg_n,
          backg_all - backg_n
        ),
        nrow = 2,
        byrow = TRUE
      ) + 0.1
      
      (a[1, 1] * a[2, 2]) / (a[1, 2] * a[2, 1])
    }
  ) %>%
  ungroup() %>%
  group_by(pair) %>%
  mutate(
    fdr = p.adjust(p_value, method = "fdr"),
    log2OR = log2(odds_ratio),
    sig = fdr < 0.05
  ) %>%
  ungroup()

# Save merged enrichment table
fwrite(
  as.data.table(df),
  "enrichment_results_p_model3_merged_features.tsv",
  sep = "\t"
)

# ---------------------------
# Order rows and columns
# ---------------------------
feature_order <- c(
  "Active promoter",
  "Weak promoter",
  "Poised promoter",
  "Strong enhancer",
  "Weak enhancer",
  "Insulator",
  "Txn transition",
  "Txn elongation",
  "Weak transcription",
  "Repressed",
  "Heterochromatin",
  "Repetitive/CNV",
  "Exon",
  "Intron",
  "Super-enhancer"
)

pair_order <- c(
  "AFR_vs_Others",
  "AMR_vs_Others",
  "EAS_vs_Others",
  "EUR_vs_Others",
  "SAS_vs_Others"
)

df <- df %>%
  mutate(
    feature = factor(feature, levels = rev(feature_order)),
    pair = factor(pair, levels = pair_order),
    log2OR_capped = pmax(pmin(log2OR, 3), -3)
  ) %>%
  filter(
    !is.na(feature),
    !is.na(pair),
    is.finite(log2OR_capped)
  )

# ---------------------------
# Prepare data for bar plot
# ---------------------------
df_plot <- df %>%
  mutate(
    significance = if_else(
      fdr < 0.05,
      "Significant",
      "Not significant"
    ),
    significance = factor(
      significance,
      levels = c("Significant", "Not significant")
    ),
    feature = factor(
      feature,
      levels = rev(feature_order)
    ),
    pair = factor(
      pair,
      levels = pair_order,
      labels = c(
        "AFR versus others",
        "AMR versus others",
        "EAS versus others",
        "EUR versus others",
        "SAS versus others"
      )
    ),
    log2OR_capped = pmax(
      pmin(log2OR, 3),
      -3
    )
  ) %>%
  filter(
    !is.na(feature),
    !is.na(pair),
    is.finite(log2OR_capped)
  )

# ---------------------------
# Prepare grouped bar-plot data
# ---------------------------
# continent_colours <- c(
#   "AFR_vs_Others" = "#D55E00",
#   "AMR_vs_Others" = "#CC79A7",
#   "EAS_vs_Others" = "#009E73",
#   "EUR_vs_Others" = "#0072B2",
#   "SAS_vs_Others" = "#E69F00"
# )

df_plot <- df %>%
  mutate(
    feature = factor(
      feature,
      levels = rev(feature_order)
    ),
    pair = factor(
      pair,
      levels = pair_order,
      labels = c("AFR", "AMR", "EAS", "EUR", "SAS")
    ),
    log2OR_capped = pmax(pmin(log2OR, 3), -3),
    significance = if_else(
      fdr < 0.05,
      "Significant",
      "Not significant"
    )
  ) %>%
  filter(
    !is.na(feature),
    !is.na(pair),
    is.finite(log2OR_capped)
  )

# Update colour names after relabelling pair
continent_colours <- c(
  "EAS" = "#2ca02c", 
  "SAS" = "#ff7f0e", 
  "AFR" = "#1f77b4", 
  "EUR" = "#9467bd", 
  "AMR" = "#d62728"
)

significance_colours <- c(
  "Significant" = "#B2182B",
  "Not significant" = "#BDBDBD"
)

# Same dodge setting must be used for bars and significance points
dodge_position <- position_dodge2(
  width = 0.82,
  preserve = "single",
  padding = 0.08
)

# ---------------------------
# Diverging grouped bar plot
# ---------------------------
p <- ggplot(
  df_plot,
  aes(
    x = log2OR_capped,
    y = feature,
    fill = pair,
    group = pair
  )
) +
  geom_vline(
    xintercept = 0,
    colour = "black",
    linewidth = 0.35
  ) +
  
  # Continental-group bars
  geom_col(
    position = dodge_position,
    width = 0.72,
    alpha = 0.82,
    colour = NA
  ) +
  
  # Significance marker at bar end
  geom_point(
    aes(
      x = log2OR_capped,
      colour = significance
    ),
    position = dodge_position,
    shape = 16,
    size = 1.5,
    show.legend = TRUE,alpha = 0.6
  ) +
  
  scale_fill_manual(
    values = continent_colours,
    breaks = c("AFR", "AMR", "EAS", "EUR", "SAS"),
    name = "Continental group"
  ) +
  
  scale_colour_manual(
    values = significance_colours,
    breaks = c("Significant", "Not significant"),
    name = NULL
  ) +
  
  scale_x_continuous(
    limits = c(-3.15, 3.15),
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    labels = c(
      "≤−3", "−2", "−1", "0", "1", "2", "≥3"
    ),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = expression(log[2]~"(odds ratio)"),
    y = NULL
  ) +
  
  annotate(
    "text",
    x = -1.5,
    y = Inf,
    label = "Depletion",
    vjust = 1.4,
    size = 3,
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    label = "Enrichment",
    vjust = 1.4,
    size = 3,
    fontface = "bold"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.line.y = element_blank(),
    axis.line.x = element_line(
      colour = "black",
      linewidth = 0.7
    ),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(
      colour = "black",
      linewidth = 0.7
    ),
    axis.text.x = element_text(
      size = 11,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 12,
      colour = "black",
      margin = margin(t = 6)
    ),
    legend.position = "top",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.32, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 4, 4)
  ) +
  
  guides(
    fill = guide_legend(
      order = 1,
      nrow = 1,
      byrow = TRUE
    ),
    colour = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(size = 2)
    )
  )


ggsave(
  "population_specific_varCpG_feature_enrichment_grouped_barplot_fdr_shared1.pdf",
  plot = p,
  width = 6,
  height = 5.6,
)
