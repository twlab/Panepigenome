library(data.table)
library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 12))

# ---------------------------
# Read data
# ---------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

continent_order <- c("AFR", "AMR", "EAS", "EUR", "SAS")

df <- fread(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Figures/Model_IND/CREs/Pair/p/enrichment_shared1_pairwise_all.txt",
  sep = "\t",
  header = TRUE
) %>%
  separate(
    pair,
    into = c("group1", "group2"),
    sep = "_vs_",
    remove = FALSE
  ) %>%
  mutate(
    rank1 = match(group1, continent_order),
    rank2 = match(group2, continent_order),
    
    pair = if_else(
      rank1 < rank2,
      paste(group1, "vs", group2),
      paste(group2, "vs", group1)
    ),
    
    pair = factor(
      pair,
      levels = c(
        "AFR vs AMR",
        "AFR vs EAS",
        "AFR vs EUR",
        "AFR vs SAS",
        "AMR vs EAS",
        "AMR vs EUR",
        "AMR vs SAS",
        "EAS vs EUR",
        "EAS vs SAS",
        "EUR vs SAS"
      )
    )
  ) %>%
  select(-group1, -group2, -rank1, -rank2)

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

pair_order <- unique(df$pair)

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
# heatmap
# ---------------------------
p <- ggplot(df, aes(x = pair, y = feature, fill = log2OR_capped)) +
  geom_tile(
    color = "white",
    linewidth = 0.25
  ) +
  geom_point(
    data = df %>% filter(sig),
    aes(x = pair, y = feature),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.3,
    fill = "black",
    color = "black",
    stroke = 0
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-3, 3),
    breaks = c(-3, -1.5, 0, 1.5, 3),
    labels = c("≤−3", "−1.5", "0", "1.5", "≥3"),
    name = expression(log[2]~"(odds ratio)")
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(
      size = 11,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank()
  )

ggsave(
  "population_specific_varCpG_feature_enrichment_heatmap_p_model1_merged_features.pdf",
  plot = p,
  width = 7.0,
  height = 4.6
)
