library(ggplot2)
library(dplyr)
library(tidyr)

# ---- Load data ----
data <- read.table("/scratch/zdong/Projects/PanEpiG/V1-9/Bed-assembly/Analysis/DNAm/input.log", header = F)

# ---- Ensure proper column names ----
colnames(data)[3:7] <- c("Nonref","Ref","Feature","Pop","continent")

# ---- Reshape to long format for Ref/Nonref ----
data_long <- data %>%
  pivot_longer(cols = c(Nonref, Ref), names_to = "Type", values_to = "Value")

wilcox.test(data_long$Value[data_long$Type=='Nonref' & data_long$Feature=='CGI'],data_long$Value[data_long$Type=='Ref' & data_long$Feature=='CGI'],paired = T)
median(data_long$Value[data_long$Type=='Nonref' & data_long$Feature=='CGI'])
median(data_long$Value[data_long$Type=='Ref' & data_long$Feature=='CGI'])

wilcox.test(data_long$Value[data_long$Type=='Nonref' & data_long$Feature=='promoter'],data_long$Value[data_long$Type=='Ref' & data_long$Feature=='promoter'],paired = T)
median(data_long$Value[data_long$Type=='Nonref' & data_long$Feature=='promoter'])
median(data_long$Value[data_long$Type=='Ref' & data_long$Feature=='promoter'])

# ---- Mean per continent × feature × type ----
summary_df <- data_long %>%
  group_by(continent, Feature, Type) %>%
  summarise(MeanValue = median(Value, na.rm = TRUE), .groups = "drop")

write.table(summary_df,file="summary.txt",sep='\t',quote = F)

# ---- Feature ordering ----
feature_order <- summary_df %>%
  group_by(Feature) %>%
  summarise(OverallMean = mean(MeanValue, na.rm = TRUE)) %>%
  arrange(OverallMean) %>%
  pull(Feature)
summary_df$Feature <- factor(summary_df$Feature, levels = feature_order)


# ---- Feature ordering ----
feature_order <- summary_df %>%
  group_by(Feature) %>%
  summarise(OverallMedian = median(MeanValue, na.rm = TRUE), .groups = "drop") %>%
  arrange(OverallMedian) %>%
  pull(Feature)

summary_df$Feature <- factor(summary_df$Feature, levels = feature_order)

# ---- Nature-style plot ----
continent_cols <- c(
  AFR = "#0072B2",
  AMR = "#D55E00",
  EAS = "#009E73",
  EUR = "#7B3294",
  SAS = "#E69F00"
)

summary_df$Type <- factor(summary_df$Type, levels = c("Ref", "Nonref"))

p <- ggplot(
  summary_df,
  aes(
    x = MeanValue,
    y = Feature,
    color = continent,
    shape = Type
  )
) +
  geom_vline(
    xintercept = c(79.75, 82.55),
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.35
  ) +
  geom_point(
    size = 2.6,
    alpha = 0.9,
    stroke = 0.25,
    position = position_dodge(width = 0.55)
  ) +
  scale_color_manual(values = continent_cols) +
  scale_shape_manual(
    values = c("Ref" = 17, "Nonref" = 16),
    labels = c("GRCh38 CpGs", "Non-reference CpGs")
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = "Median methylation level (%)",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.35, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 8, 5, 5)
  )

ggsave(
  filename = "Fisher_medians_methylation_by_continent_ref_nonref.pdf",
  plot = p,
  width = 8,
  height = 5,
  units = "in"
)
