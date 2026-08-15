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
# Clean variables
# ---------------------------
data <- data %>%
  mutate(
    count_million = count / 1e6,
    pop = factor(pop, levels = c("AFR", "AMR", "EAS", "EUR", "SAS")),
    type = factor(
      type,
      levels = c(
        "SNV-CpGgain",
        "SNV-CpGloss",
        "Indel-insertion",
        "Indel-deletion",
        "SV-insertion",
        "SV-deletion"
      ),
      labels = c(
        "SNV gain",
        "SNV loss",
        "Indel insertion",
        "Indel deletion",
        "SV insertion",
        "SV deletion"
      )
    ),
    variant_class = case_when(
      grepl("^SNV", as.character(type)) ~ "SNV",
      grepl("^Indel", as.character(type)) ~ "Indel",
      grepl("^SV", as.character(type)) ~ "SV",
      TRUE ~ NA_character_
    )
  )

# ---------------------------
# Median ± SD by type and continental group
# ---------------------------
summary_df <- data %>%
  group_by(type, pop) %>%
  summarise(
    median_count = median(count_million, na.rm = TRUE),
    sd_count     = sd(count_million, na.rm = TRUE),
    mean_count   = mean(count_million, na.rm = TRUE),
    n            = sum(!is.na(count_million)),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = pmax(median_count - sd_count, 0),
    ymax = median_count + sd_count
  )

write.table(
  summary_df,
  file = "median_SD_CpG_count_by_type_pop.tsv",
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
# median ± SD dot plot
# ---------------------------
pd <- position_dodge(width = 0.65)

p <- ggplot(
  summary_df,
  aes(x = type, y = median_count, color = pop)
) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    position = pd,
    width = 0.18,
    linewidth = 0.5,
    alpha = 0.6
  ) +
  geom_point(
    position = pd,
    size = 2.5,
    alpha = 0.6
  ) +
  scale_color_manual(values = population_colors) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  labs(
    x = NULL,
    y = "CpG count (millions)",
    color = NULL
  ) +
  theme(
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      vjust = 1,
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.4, "cm"),
    panel.grid = element_blank()
  )

ggsave(
  "CpG_geneticvariants_count_median_SD_by_pop.pdf",
  plot = p,
  width = 5.6,
  height = 3.5
)
