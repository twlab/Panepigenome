library(tidyverse)
library(scales)

theme_set(theme_void(base_size = 14))

# ---------------------------
# Read data
# ---------------------------
data <- read.table("cpg.count.pop.tsv", header = TRUE, stringsAsFactors = FALSE)

# ---------------------------
# Summarize CpG counts by haplotype and variant class
# ---------------------------
prop_df <- data %>%
  group_by(id, hap, type) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(
    SNV   = ifelse(is.na(SNV), 0, SNV),
    Indel = ifelse(is.na(Indel), 0, Indel),
    SV    = ifelse(is.na(SV), 0, SV),
    total = SNV + Indel + SV
  ) %>%
  filter(total > 0)

# ---------------------------
# Mean composition across haplotypes
# ---------------------------
composition_df <- prop_df %>%
  summarise(
    SNV   = mean(SNV, na.rm = TRUE),
    Indel = mean(Indel, na.rm = TRUE),
    SV    = mean(SV, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Class",
    values_to = "Count"
  ) %>%
  mutate(
    Fraction = Count / sum(Count),
    Percent = Fraction * 100,
    Count_M = Count / 1e6,
    Class = factor(Class, levels = c("SNV", "SV", "Indel")),
    Label = paste0(Class, "\n", sprintf("%.1f%%", Percent))
  )

# Save composition table
write.table(
  composition_df,
  "varCpG_variant_class_composition.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------
# Donut plot
# ---------------------------
p <- ggplot(composition_df, aes(x = 2, y = Fraction, fill = Class)) +
  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.5
  ) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 3.6,
    color = "black",
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = 0.5,
    y = 0,
    label = paste0(
      "Mean total\n",
      sprintf("%.2fM", sum(composition_df$Count) / 1e6)
    ),
    size = 3.5,
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = c(
      SNV = "#595959",
      SV = "#66C2A5",
      Indel = "#A6761D"
    ),
    labels = c("SNV-CpGs", "SV-CpGs", "Indel-CpGs")
  ) +
  labs(fill = NULL) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

p

ggsave(
  "varCpG_variant_class_donut.pdf",
  plot = p,
  width = 4.2,
  height = 3.5
)
