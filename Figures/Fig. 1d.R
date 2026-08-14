library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ---------------------------------------------
# Load data
# ---------------------------------------------
df <- read.table("merged.tsv", header = TRUE)

df <- df %>%
  mutate(
    Pop = factor(pop, levels = c("AFR", "AMR", "EAS", "EUR", "SAS")),
    Sex = factor(
      sex,
      levels = c("female", "male"),
      labels = c("Female", "Male")
    )
  )

# ---------------------------------------------
# Variables matching Results section
# ---------------------------------------------
vars <- c("num", "mean_cpg", "mean_len")

var_labels <- c(
  num      = "Number of CGIs",
  mean_cpg = "CpGs per CGI",
  mean_len = "CGI length (bp)"
)

ref_values <- c(
  num      = 51502,
  mean_cpg = 56.869,
  mean_len = 609.194
)

pop_colors <- c(
  AFR = "#0072B2",
  AMR = "#D55E00",
  EAS = "#009E73",
  EUR = "#CC79A7",
  SAS = "#E69F00"
)

wilcox.test(df$mean_len,mu = 609.194)

wilcox.test(df$mean_cpg,mu = 56.869)
mean(df$mean_cpg)
# ---------------------------------------------
# Long format
# ---------------------------------------------
plot_df <- df %>%
  select(Pop, Sex, all_of(vars)) %>%
  pivot_longer(
    cols = all_of(vars),
    names_to = "Metric_raw",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric_raw,
      levels = vars,
      labels = var_labels[vars]
    )
  )

# GRCh38 reference values for dashed lines
ref_df <- tibble(
  Metric_raw = vars,
  Metric = factor(var_labels[vars], levels = var_labels[vars]),
  GRCh38 = as.numeric(ref_values[vars])
)

# ---------------------------------------------
# Plot: compact Fig. 1g
# ---------------------------------------------
p <- ggplot(plot_df, aes(x = Pop, y = Value, fill = Pop)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "black",
    alpha = 0.6,
    linewidth = 0.3
  ) +
  geom_point(
    aes(color = Pop),
    position = position_jitter(width = 0.10, height = 0),
    size = 0.45,
    alpha = 0.3,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 1.5,
    fill = "white",
    color = "black",
    stroke = 0.3
  ) +
  geom_hline(
    data = ref_df,
    aes(yintercept = GRCh38),
    linetype = "dashed",
    color = "#76B7B2",
    linewidth = 1,
    inherit.aes = FALSE,
    alpha = 0.6
  ) +
  facet_grid(
    Metric ~ Sex,
    scales = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(values = pop_colors) +
  scale_color_manual(values = pop_colors) +
  scale_y_continuous(
    labels = function(x) ifelse(abs(x) >= 1000, comma(x), x),
    expand = expansion(mult = c(0.08, 0.10))
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.35, "cm"),
    
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y.left = element_text(size = 12, face = "bold", angle = 90),
    
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.7, "lines")
  )

# p

# ---------------------------------------------
# Save
# ---------------------------------------------
ggsave(
  "Fig1g_CGI_haplotype_resolved_Nature.pdf",
  plot = p,
  width = 6,
  height = 4
)
