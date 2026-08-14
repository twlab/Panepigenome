library(data.table)
library(ggplot2)
library(scales)

df <- fread("re_con_removeheader_MAF5.txt")
# df$var_residuals<-df$var_residuals/10000
pop_colors <- c(
  AFR = "#0072B2",
  AMR = "#D55E00",
  EAS = "#009E73",
  EUR = "#CC79A7",
  SAS = "#E69F00"
)

df[, group := factor(group, levels = c("AFR", "AMR", "EAS", "EUR", "SAS"))]

summary_stats <- df[
  , .(
    mean_var   = mean(var_residuals, na.rm = TRUE),
    median_var = median(var_residuals, na.rm = TRUE),
    q25        = quantile(var_residuals, 0.25, na.rm = TRUE),
    q75        = quantile(var_residuals, 0.75, na.rm = TRUE),
    N          = .N
  ),
  by = group
]

write.table(
  summary_stats,
  "Variance_con_nonref_geno_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

set.seed(1)
plot_df <- df[
  is.finite(var_residuals),
  .SD[sample(.N, min(.N, 100000))],
  by = group
]

p <- ggplot(plot_df, aes(x = group, y = var_residuals, color = group, fill = group)) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    alpha = 0.25,
    linewidth = 0.35,
    color = NA
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    linewidth = 0.30,
    fill = "white",
    color = "black"
  ) +
  geom_point(
    data = summary_stats,
    aes(x = group, y = median_var),
    inherit.aes = FALSE,
    shape = 95,
    size = 8,
    color = "black"
  ) +
  scale_fill_manual(values = pop_colors) +
  scale_color_manual(values = pop_colors) +
  labs(
    x = NULL,
    y = "Variance"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(5, 6, 5, 5)
  )


ggsave(
  "Variance_con_nonref_geno.pdf",
  plot = p,
  width = 3.8,
  height = 3.1,
  device = "pdf"
)
