library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 12))

df <- tribble(
  ~Category, ~Sample_ID, ~DMR_0.2_5CpG, ~DMR_0.3_5CpG, ~DMR_0.2_10CpG, ~DMR_0.3_10CpG, ~DML_0.2, ~DML_0.3,
  "Instrument", "HG00609", 295, 286, 67, 62, 7144, 7060,
  "Instrument", "HG01192", 12, 1, 7, 1, 6918, 6841,
  "SMRT cell position", "NA18565", 10, 1, 6, 0, 5293, 5210,
  "SMRT cell position", "NA20850", 8, 1, 6, 0, 8123, 8026,
  "Timepoint-Sequel II", "HG04228", 154, 54, 97, 24, 8623, 8569,
  "Timepoint-Revio", "HG03583", 11, 0, 7, 0, 11907, 11714,
  "Timepoint-Revio", "NA19909", 7, 0, 4, 0, 5118, 5038,
  "WGBS duplicate", "HG03516", 1488, 311, 229, 10, 10475, 7898,
  "WGBS duplicate", "HG01978", 2324, 696, 331, 22, 12808, 10856,
  "WGBS duplicate", "HG01952", 2406, 750, 351, 24, 16224, 13694,
  "WGBS duplicate", "HG00741", 2141, 592, 379, 22, 12326, 10115,
  "WGBS duplicate", "HG00621", 1860, 537, 335, 26, 13343, 12080
)

plot_df <- df %>%
  mutate(
    Category = case_when(
      Category %in% c("Timepoint-Sequel II", "Timepoint-Revio") ~ "Timepoint",
      TRUE ~ Category
    )
  ) %>%
  select(
    Category,
    Sample_ID,
    `DMRs, delta beta >= 0.2, >=5 CpGs` = DMR_0.2_5CpG,
    `DMRs, delta beta >= 0.2, >=10 CpGs` = DMR_0.2_10CpG,
    `DMLs, delta beta >= 0.2` = DML_0.2
  ) %>%
  pivot_longer(
    cols = starts_with(c("DMRs", "DMLs")),
    names_to = "Metric",
    values_to = "Count"
  ) %>%
  mutate(
    Category = factor(
      Category,
      levels = c(
        "Instrument",
        "SMRT cell position",
        "Timepoint",
        "WGBS duplicate"
      )
    ),
    Metric = factor(
      Metric,
      levels = c(
        "DMRs, delta beta >= 0.2, >=5 CpGs",
        "DMRs, delta beta >= 0.2, >=10 CpGs",
        "DMLs, delta beta >= 0.2"
      )
    )
  )

p <- ggplot(plot_df, aes(x = Count, y = Category)) +
  geom_point(
    aes(fill = Category),
    shape = 21,
    size = 2.7,
    color = "black",
    stroke = 0.25,
    alpha = 0.85,
    position = position_jitter(height = 0.10, width = 0)
  ) +
  facet_wrap(
    ~Metric,
    scales = "free_x",
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "Instrument" = "#4D4D4D",
      "SMRT cell position" = "#7B3294",
      "Timepoint" = "#0072B2",
      "WGBS duplicate" = "#D55E00"
    )
  ) +
  scale_x_continuous(
    labels = comma,
    expand = expansion(mult = c(0.04, 0.12))
  ) +
  labs(
    x = "Number of differential methylation features",
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth =  0.75, color = "black"),
    axis.ticks = element_line(linewidth = 0.75, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(5, 8, 5, 5)
  )


ggsave(
  "technical_replicate_DMR_DML.pdf",
  plot = p,
  width = 5.2,
  height = 2.4,
  device = "pdf"
)
