library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 14))

# ---------------------------
# Read data
# ---------------------------
df <- read.table(
  "input.log",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(df) <- c("Variant", "Category", "Value")

# ---------------------------
# Clean data
# ---------------------------
df <- df %>%
  mutate(
    Value = Value / 100,
    Category = factor(
      Category,
      levels = c("Benign","Pathogenic", "Risk", "Protective", "Pharmacogenomic")
    ),
    Variant = factor(
      Variant,
      levels = c(
        "variants with CpG gain/loss",
        "variants with CpG gain",
        "variants with CpG loss"
      ),
      labels = c(
        "CpG gain/loss",
        "CpG gain",
        "CpG loss"
      )
    )
  ) %>%
  filter(!is.na(Category), !is.na(Variant), is.finite(Value))

median(df$Value[df$Variant=="CpG gain/loss" & df$Category=='Benign'])
median(df$Value[df$Variant=="CpG gain/loss" & df$Category=='Pathogenic'])
median(df$Value[df$Variant=="CpG gain/loss" & df$Category=='Risk'])
median(df$Value[df$Variant=="CpG gain/loss" & df$Category=='Protective'])
median(df$Value[df$Variant=="CpG gain/loss" & df$Category=='Pharmacogenomic'])

wilcox.test(df$Value[df$Variant=="CpG gain/loss" & df$Category=='Benign'],
            df$Value[df$Variant=="CpG gain/loss" & df$Category==c('Pathogenic','Risk','Protective','Pharmacogenomic')])

# ---------------------------
# Median ± SD summary
# ---------------------------
summary_df <- df %>%
  group_by(Category, Variant) %>%
  summarise(
    median_meth = median(Value, na.rm = TRUE),
    sd_meth     = sd(Value, na.rm = TRUE),
    mean_meth   = mean(Value, na.rm = TRUE),
    n           = sum(!is.na(Value)),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = pmax(median_meth - sd_meth, 0),
    ymax = pmin(median_meth + sd_meth, 1)
  )

write.table(
  summary_df,
  "clinical_variant_methylation_median_SD.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------
# Wilcoxon tests vs Benign within each variant class
# ---------------------------
pval_df <- df %>%
  group_by(Variant) %>%
  group_modify(~ {
    benign <- .x %>% filter(Category == "Benign") %>% pull(Value)
    
    tibble(
      comparison = c("Pathogenic vs Benign", "Risk vs Benign", "Protective vs Benign"),
      Category = factor(
        c("Pathogenic", "Risk", "Protective"),
        levels = c("Pathogenic", "Risk", "Protective", "Benign")
      ),
      p_value = c(
        wilcox.test((.x %>% filter(Category == "Pathogenic") %>% pull(Value)), benign)$p.value,
        wilcox.test((.x %>% filter(Category == "Risk") %>% pull(Value)), benign)$p.value,
        wilcox.test((.x %>% filter(Category == "Protective") %>% pull(Value)), benign)$p.value
      )
    )
  }) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    p_value_sci = formatC(p_value, format = "e", digits = 2),
    FDR_sci = formatC(FDR, format = "e", digits = 2)
  )

write.table(
  pval_df,
  "clinical_variant_methylation_wilcox_vs_benign.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------
# Colors distinct from continental group colors
# ---------------------------
variant_colors <- c(
  "CpG gain/loss" = "#fd7901",
  "CpG gain"      = "#B07AA1",
  "CpG loss"      = "#76B7B2"
)

# ---------------------------
# median ± SD dot plot
# ---------------------------
pd <- position_dodge(width = 0.55)

p <- ggplot(
  summary_df,
  aes(x = Category, y = median_meth, color = Variant)
) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    position = pd,
    width = 0.16,
    linewidth = 0.4,
    alpha = 0.6
  ) +
  geom_point(
    position = pd,
    size = 2.8,
    alpha = 0.95
  ) +
  geom_hline(
    yintercept = 81.405 / 100, ### 
    color = "grey40",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_color_manual(values = variant_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Methylation level",
    color = NULL
  ) +
  theme(
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.4, "cm"),
    panel.grid = element_blank()
  )

# Save PDF
ggsave(
  filename = "clinical_variant_methylation_median_SD_dotplot.pdf",
  plot = p,
  width = 5.4,
  height = 2.6,
  device = "pdf"
)
