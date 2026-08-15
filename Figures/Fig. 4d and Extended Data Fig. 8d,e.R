library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)
library(grid)

theme_set(theme_classic(base_size = 10))

# ---------------------------
# 1. Read and prepare data
# ---------------------------
data <- fread(
  "CpG_overlap_geno_vs_latent_1vsOthers.tsv",
  sep = "\t",
  header = TRUE
) %>%
  mutate(
    total = overlap + model1_only + model3_only,
    
    comparison_label = recode(
      comparison,
      "AFR_vs_Others" = "AFR vs others",
      "AMR_vs_Others" = "AMR vs others",
      "EAS_vs_Others" = "EAS vs others",
      "EUR_vs_Others" = "EUR vs others",
      "SAS_vs_Others" = "SAS vs others"
    ),
    
    comparison_label = factor(
      comparison_label,
      levels = c(
        "AFR vs others",
        "AMR vs others",
        "EAS vs others",
        "EUR vs others",
        "SAS vs others"
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
    overlap,
    model1_only,
    model3_only
  ) %>%
  pivot_longer(
    cols = c(model3_only, model1_only, overlap),
    names_to = "Category",
    values_to = "Count"
  ) %>%
  mutate(
    Category = factor(
      Category,
      levels = c(
        "model3_only",
        "model1_only",
        "overlap"
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
    size = 2.8,
    colour = "white",
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
    size = 2.8,
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
  
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      size = 8.5,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 9.5,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9.5,
      colour = "black",
      margin = margin(t = 6)
    ),
    
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.line.x = element_line(
      linewidth = 0.35,
      colour = "black"
    ),
    axis.ticks.x = element_line(
      linewidth = 0.30,
      colour = "black"
    ),
    
    panel.grid = element_blank(),
    
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.text = element_text(
      size = 8.5,
      colour = "black"
    ),
    legend.key.size = unit(0.32, "cm"),
    legend.spacing.x = unit(0.10, "cm"),
    legend.margin = margin(b = 2),
    
    plot.margin = margin(4, 12, 4, 4)
  )

# ---------------------------
# 5. Save
# ---------------------------
ggsave(
  filename = "CpG_overlap_geno_latent_1vsOthers.pdf",
  plot = nature_plot,
  width = 6.4,
  height = 3.3,
  units = "in"
)
