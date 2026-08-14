library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 14))

df <- read.table(
  "input",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

df <- df %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Other", "DNA", "LTR", "LINE", "SINE", "SVA")
    ),
    Method = factor(
      Method,
      levels = c("CpG-gain", "CpG-loss"),
      labels = c("CpG gain", "CpG loss")
    ),
    Sample = factor(
      Sample,
      levels = c("SNV", "Indel", "SV")
    )
  ) %>%
  filter(
    !is.na(Group),
    !is.na(Method),
    !is.na(Sample),
    is.finite(Value),
    Sample == "SV"
  )

fold_df <- df %>%
  filter(Group != "Other") %>%
  select(Method, Group, Value) %>%
  group_by(Method, Group) %>%
  summarise(
    Value = sum(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Method,
    values_from = Value
  ) %>%
  mutate(
    fold_change_gain_vs_loss = `CpG gain` / `CpG loss`,
    log2FC_gain_vs_loss = log2(fold_change_gain_vs_loss)
  )

fold_df
# # A tibble: 5 × 5
# Group `CpG gain` `CpG loss` fold_change_gain_vs_loss log2FC_gain_vs_loss
# <fct>      <dbl>      <dbl>                    <dbl>               <dbl>
#   1 DNA       0.0123    0.00813                     1.52               0.602
# 2 LTR       0.0505    0.0376                      1.34               0.427
# 3 LINE      0.0990    0.0564                      1.76               0.812
# 4 SINE      0.248     0.123                       2.01               1.01 
# 5 SVA       0.0666    0.0286                      2.33               1.22 

a<-read.table("TE_ins.txt",header=T)
b<-read.table('TE_del.txt',header=T)

wilcox.test(a$SVA/a$ALL,b$SVA/b$ALL,paired = T)
wilcox.test(a$SINE/a$ALL,b$SINE/b$ALL,paired = T)

feature_colors <- c(
  Other = "grey85",
  DNA   = "#4a72e8",
  LTR   = "#006600",
  LINE  = "#ff6600",
  SINE  = "#cc0000",
  SVA   = "#ea53c4"
)

plot_df <- df %>%
  group_by(Method) %>%
  mutate(
    Fraction = Value / sum(Value, na.rm = TRUE)
  ) %>%
  ungroup()

p <- ggplot(plot_df, aes(x = Method, y = Fraction, fill = Group)) +
  geom_col(
    width = 0.62,
    color = "white",
    linewidth = 0.25,
    alpha = 0.6
  ) +
  scale_fill_manual(values = feature_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "SV-associated var-CpGs (%)",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    panel.grid = element_blank()
  )

# p

ggsave(
  "TE_feature_composition_SV_varCpGs.pdf",
  plot = p,
  width = 3.6,
  height = 2.8
)
