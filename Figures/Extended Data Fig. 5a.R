library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 14))

# ---------------------------
# Read data
# ---------------------------
data <- read.table(
  "input.log",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# ---------------------------
# Clean variables and combine variant classes
# count is assumed to be a proportion, so convert to percent
# ---------------------------
data <- data %>%
  mutate(
    count_pct = count * 100,
    pop = factor(pop, levels = c("AFR", "AMR", "EAS", "EUR", "SAS")),
    variant_class = case_when(
      type == "SNV" ~ "SNV",
      type %in% c("Indel-insertion", "Indel-deletion") ~ "Indel",
      type %in% c("SV-insertion", "SV-deletion") ~ "SV",
      TRUE ~ NA_character_
    ),
    variant_class = factor(variant_class, levels = c("SNV", "Indel", "SV"))
  ) %>%
  filter(!is.na(variant_class))

# ---------------------------
# Mean ± SD by variant class and population
# ---------------------------
summary_df <- data %>%
  group_by(variant_class, pop) %>%
  summarise(
    mean_pct = mean(count_pct, na.rm = TRUE),
    sd_pct   = sd(count_pct, na.rm = TRUE),
    n        = sum(!is.na(count_pct)),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = pmax(mean_pct - sd_pct, 0),
    ymax = mean_pct + sd_pct
  )

write.table(
  summary_df,
  "mean_SD_percent_variants_with_CpGs_by_variant_class_pop.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------
# Colors
# ---------------------------
population_colors <- c(
  AFR = "#0072B2",
  AMR = "#D55E00",
  EAS = "#009E73",
  EUR = "#CC79A7",
  SAS = "#E69F00"
)

# ---------------------------
#  mean ± SD dot plot
# ---------------------------
pd <- position_dodge(width = 0.55)

p <- ggplot(
  summary_df,
  aes(x = variant_class, y = mean_pct, color = pop)
) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    position = pd,
    width = 0.16,
    linewidth = 0.35,
    alpha = 0.75
  ) +
  geom_point(
    position = pd,
    size = 2.7,
    alpha = 0.95
  ) +
  scale_color_manual(values = population_colors) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Variants altering CpG sites (%)",
    color = NULL
  ) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    panel.grid = element_blank()
  )


ggsave(
  "Prop_variant_CpG_mean_SD_combined_class.pdf",
  plot = p,
  width = 4.4,
  height = 3.3
)
