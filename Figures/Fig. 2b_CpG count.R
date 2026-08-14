library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 14))

# ---------------------------
# Read data
# ---------------------------
data <- read.table(
  "input_methylation.log",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# ---------------------------
# Clean variables
# ---------------------------
data <- data %>%
  mutate(
    count = count / 100,
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
    )
  )

wilcox.test(data$count[data$type=='SNV gain'], data$count[data$type=='SNV loss'],paired = T)
shapiro.test(data$count[data$type=='SNV gain'])
median(data$count[data$type=='SNV gain'])
# 0.775
median(data$count[data$type=='SNV loss'])
# 0.775


# ---------------------------
# Median ± SD by type and population
# ---------------------------
summary_df <- data %>%
  group_by(type, pop) %>%
  summarise(
    median_methylation = median(count, na.rm = TRUE),
    sd_methylation     = sd(count, na.rm = TRUE),
    n                  = sum(!is.na(count)),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = pmax(median_methylation - sd_methylation, 0),
    ymax = pmin(median_methylation + sd_methylation, 1)
  )

write.table(
  summary_df,
  file = "median_SD_methylation_by_type_pop.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------
# Colors
# ---------------------------
pop_colors <- c(
  AFR = "#0072B2",
  AMR = "#D55E00",
  EAS = "#009E73",
  EUR = "#CC79A7",
  SAS = "#E69F00"
)

# ---------------------------
# Nature-style median ± SD dot plot
# ---------------------------
pd <- position_dodge(width = 0.65)

p <- ggplot(
  summary_df,
  aes(x = type, y = median_methylation, color = pop)
) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    position = pd,
    width = 0.20,
    linewidth = 0.35,
    alpha = 0.6
  ) +
  geom_point(
    position = pd,
    size = 2.5,
    alpha = 0.6
  ) +
  geom_hline(
    yintercept = 0.803,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) +
  scale_color_manual(values = pop_colors) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0.5, 1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Methylation level",
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
    axis.text = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.4, "cm"),
    panel.grid = element_blank()
  )


# Save PDF
ggsave(
  filename = "CpG_variant_methylation_median_SD_by_pop.pdf",
  plot = p,
  width = 5.6,
  height = 3.5,
  device = "pdf"
)
