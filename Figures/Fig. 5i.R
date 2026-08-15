library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)
library(grid)

theme_set(theme_classic(base_size = 12))

# ---------------------------
# 1. Read and prepare data
# ---------------------------
data <- fread(
  "/scratch/zdong/Projects/PanEpiG/V1-9/eQTL/Juan_eqtl_gene/pop/eQTL_pop_counts_by_comparison.tsv",
  sep = "\t",
  header = TRUE
) %>%
  mutate(
    total = CpG_Shared + CpG_Model1_only + CpG_Model3_only,
    
    comparison_label = str_replace_all(
      Comparison ,
      "_vs_",
      " vs "
    ),
    
    comparison_label = factor(
      comparison_label,
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
  )

# ---------------------------
# 2. Convert to long format
# ---------------------------
plot_data <- data %>%
  select(
    comparison_label,
    total,
    CpG_Shared,
    CpG_Model1_only,
    CpG_Model3_only
  ) %>%
  pivot_longer(
    cols = c(CpG_Model3_only, CpG_Model1_only, CpG_Shared),
    names_to = "Category",
    values_to = "Count"
  ) %>%
  mutate(
    Category = factor(
      Category,
      levels = c(
        "CpG_Model3_only",
        "CpG_Model1_only",
        "CpG_Shared"
      ),
      labels = c(
        "Integrated model only",
        "CpG gain/loss only",
        "Shared between models"
      )
    )
  ) %>%
  group_by(comparison_label) %>%
  mutate(
    proportion = Count / sum(Count),
    
    count_label = case_when(
      proportion >= 0.08 & Count >= 1000 ~
        label_number(
          scale = 1e-3,
          accuracy = 0.1,
          suffix = "k"
        )(Count),
      
      proportion >= 0.08 ~
        comma(Count),
      
      TRUE ~ ""
    )
  ) %>%
  ungroup()

# Total labels at the end of each bar
total_labels <- data %>%
  transmute(
    comparison_label,
    total,
    label = label_number(
      scale = 1e-3,
      accuracy = 1,
      suffix = "k"
    )(total)
  )

# ---------------------------
# 3. colours
# ---------------------------
category_colours <- c(
  "Integrated model only" = "#c35bca",
  "CpG gain/loss only"    = "#fd7901",
  "Shared between models" = "#999999"
)

# ---------------------------
# 4. Plot
# ---------------------------
nature_plot <- ggplot(
  plot_data,
  aes(
    y = comparison_label,
    x = Count,
    fill = Category
  )
) +
  geom_col(
    width = 0.66,
    colour = "white",
    linewidth = 0.25
  ) +
  
  # Labels inside sufficiently large segments
  geom_text(
    aes(label = count_label),
    position = position_stack(vjust = 0.5),
    size = 0,
    colour = "black",
    fontface = "plain"
  ) +
  
  # Total count at the end of each bar
  geom_text(
    data = total_labels,
    aes(
      x = total,
      y = comparison_label,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.15,
    size = 3, ###
    colour = "black"
  ) +
  
  scale_fill_manual(
    values = category_colours,
    breaks = c(
      "Integrated model only",
      "CpG gain/loss only",
      "Shared between models"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    labels = label_number(
      scale = 1e-3,
      accuracy = 1,
      suffix = "k"
    ),
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.10))
  ) +
  
  labs(
    x = "Number of continent-specific var-CpGs",
    y = NULL
  ) +
  
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        colour = NA
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
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 12,
      colour = "black",
      margin = margin(t = 6)
    ),
    
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.line.x = element_line(
      linewidth = 1,
      colour = "black"
    ),
    axis.ticks.x = element_line(
      linewidth = 1,
      colour = "black"
    ),
    
    panel.grid = element_blank(),
    
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.text = element_text(
      size = 12,
      colour = "black"
    ),
    legend.key.size = unit(0.32, "cm"),
    legend.spacing.x = unit(0.10, "cm"),
    legend.margin = margin(b = 2),
    
    plot.margin = margin(4, 12, 4, 4)
  )

# nature_plot

# ---------------------------
# 5. Save
# ---------------------------
ggsave(
  filename = "CpG_overlap_geno_latent_1vsOthers.pdf",
  width = 5,
  height = 3.6,
  units = "in"
)
