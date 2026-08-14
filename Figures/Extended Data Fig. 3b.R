library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

theme_set(theme_classic(base_size = 10))

make_density_plot_obj <- function(prefix, label) {
  
  df_nonref <- read.table(paste0(prefix, ".nonref.log"), header = FALSE)
  df_ref    <- read.table(paste0(prefix, ".ref.log"), header = FALSE)
  
  colnames(df_nonref) <- "Value"
  colnames(df_ref)    <- "Value"
  
  df <- bind_rows(
    df_nonref %>% mutate(Group = "Non-reference"),
    df_ref    %>% mutate(Group = "GRCh38")
  ) %>%
    mutate(
      Value = as.numeric(Value),
      Group = factor(Group, levels = c("GRCh38", "Non-reference"))
    ) %>%
    filter(is.finite(Value))
  
  meds <- df %>%
    group_by(Group) %>%
    summarise(
      med = median(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  group_colors <- c(
    "GRCh38" = "#4D4D4D",
    "Non-reference" = "#7B3294"
  )
  
  ggplot(df, aes(x = Value, fill = Group, color = Group)) +
    geom_density(
      alpha = 0.25,
      linewidth = 0.55,
      adjust = 1
    ) +
    geom_vline(
      data = meds,
      aes(xintercept = med, color = Group),
      linetype = "dashed",
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.06))
    ) +
    labs(
      x = "Methylation level",
      y = "Density",
      title = label,
      fill = NULL,
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.30),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(color = "black", size = 12),
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      legend.position = "top",
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.35, "cm"),
      panel.grid = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}

p_pat <- make_density_plot_obj("HG002.pat", "Paternal")
p_mat <- make_density_plot_obj("HG002.mat", "Maternal")

combined <- (p_pat / p_mat) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(
  "density_hg002_combined_pat_mat.pdf",
  plot = combined,
  width = 4.4,
  height = 3.6
)
