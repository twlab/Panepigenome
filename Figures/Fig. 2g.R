library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 14))

df <- read.table(
  "merged_SV_summary.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

df <- df %>%
  filter(!ID %in% c("HG002.hap1", "HG002.hap2")) %>%
  mutate(
    Type = case_when(
      Type == "INS" ~ "CpG gain",
      Type == "DEL" ~ "CpG loss",
      TRUE ~ Type
    ),
    mean_meth = mean / 100,
    median_meth = median / 100,
    Type = factor(Type, levels = c("CpG gain", "CpG loss")),
    TE = factor(TE, levels = c("DNA", "LTR", "LINE", "SINE", "SVA", "Other-RM", "Non-RM"))
  ) %>%
  filter(
    !is.na(Type),
    !is.na(TE),
    !is.na(median_meth),
    !TE %in% c("Other-RM", "Non-RM")
  )

type_colors <- c(
  "CpG gain" = "#B07AA1",
  "CpG loss" = "#76B7B2"
)

summary_df <- df %>%
  group_by(TE, Type) %>%
  summarise(
    mean_meth = mean(median_meth, na.rm = TRUE),
    median_meth = median(median_meth, na.rm = TRUE),
    sd_meth = sd(median_meth, na.rm = TRUE),
    n = sum(!is.na(median_meth)),
    .groups = "drop"
  )

write.table(
  summary_df,
  "SV_TE_CpG_gain_loss_methylation_summary.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pval_df <- df %>%
  select(ID, TE, Type, median_meth) %>%
  group_by(ID, TE, Type) %>%
  summarise(
    median_meth = median(median_meth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Type,
    values_from = median_meth
  ) %>%
  drop_na(`CpG gain`, `CpG loss`) %>%
  group_by(TE) %>%
  summarise(
    n_pair = n(),
    mean_CpG_gain = mean(`CpG gain`, na.rm = TRUE),
    mean_CpG_loss = mean(`CpG loss`, na.rm = TRUE),
    mean_diff = mean_CpG_gain - mean_CpG_loss,
    p_value = wilcox.test(
      `CpG gain`,
      `CpG loss`,
      paired = TRUE,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    p_value_sci = formatC(p_value, format = "e", digits = 2),
    FDR_sci = formatC(FDR, format = "e", digits = 2)
  )

write.table(
  pval_df,
  "SV_TE_CpG_gain_loss_Wilcoxon.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pd <- position_dodge(width = 0.75)

p <- ggplot(df, aes(x = TE, y = median_meth, fill = Type)) +
  geom_violin(
    position = pd,
    width = 0.72,
    trim = TRUE,
    scale = "width",
    color = "black",
    linewidth = 0.25,
    alpha = 0.6
  ) +
  stat_summary(
    aes(group = Type),
    fun = median,
    geom = "point",
    position = pd,
    shape = 23,
    size = 1.8,
    fill = "white",
    color = "black",
    stroke = 0.3
  ) +
  geom_hline(
    yintercept = 0.803,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) +
  scale_fill_manual(values = type_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0.5, 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Methylation level",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"), #angle = 35, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    panel.grid = element_blank()
  )

# p

while (!is.null(dev.list())) dev.off()

ggsave(
  "SV_CpG_gain_loss_TE_methylation_violin.pdf",
  plot = p,
  width = 4.2,
  height = 2.8
)
