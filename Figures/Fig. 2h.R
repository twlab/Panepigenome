del<-read.table("/scratch/zdong/Projects/PanEpiG/V1-9/PAV_variants/Phased/Merge/INS/Format/CGI/DEL_CpG_island/drivenCGI_line_counts.tsv",header=T)
ins<-read.table("/scratch/zdong/Projects/PanEpiG/V1-9/PAV_variants/Phased/Merge/INS/Format/CGI/INS_CpG_island/drivenCGI_line_counts.tsv",header = T)

re<-cbind(ins,del)
head(re)

library(tidyverse)
library(ggbeeswarm)

theme_set(theme_classic(base_size = 12))

# ---------------------------
# Read input
# ---------------------------

ins$type <- "CGI gain"
del$type <- "CGI loss"

df <- bind_rows(ins, del) %>%
  mutate(
    Line_Count = as.numeric(Line_Count),
    type = factor(type, levels = c("CGI gain", "CGI loss"))
  ) %>%
  filter(is.finite(Line_Count))

# ---------------------------
# Summary
# ---------------------------
summary_df <- df %>%
  group_by(type) %>%
  summarise(
    n = n(),
    mean = mean(Line_Count, na.rm = TRUE),
    median = median(Line_Count, na.rm = TRUE),
    min = min(Line_Count, na.rm = TRUE),
    max = max(Line_Count, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  summary_df,
  "CGI_gain_loss_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---------------------------
# Nature-style colors
# ---------------------------
type_colors <- c(
  "CGI gain" = "#B07AA1",
  "CGI loss" = "#76B7B2"
)

# ---------------------------
# Violin + boxplot + beeswarm
# ---------------------------
p <- ggplot(df, aes(x = type, y = Line_Count)) +
  geom_violin(
    aes(color = type),
    fill = NA,
    linewidth = 0.8,
    trim = TRUE,
    scale = "width",
    show.legend = FALSE
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.35
  ) +
  geom_quasirandom(
    aes(color = type),
    width = 0.18,
    size = 1,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  geom_point(
    data = summary_df,
    aes(x = type, y = median),
    inherit.aes = FALSE,
    shape = 95,
    size = 9,alpha = 0.6,
    color = "black"
  ) +
  scale_color_manual(values = type_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.08)),
    limits = c(600, 1000),
  ) +
  labs(
    x = NULL,
    y = "Number of SV-driven CGIs"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 6, 5, 5)
  )


ggsave(
  "CGI_gain_loss_violin_boxplot_beeswarm.pdf",
  plot = p,
  width = 2.8,
  height = 2.8
)
