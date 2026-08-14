library(tidyverse)
library(ggpubr)
library(scales)

theme_set(theme_classic(base_size = 12))

# ---------------------------
# Read and merge data
# ---------------------------
ml <- read.delim(
  "../../ML-coverage/1158/match.log",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  mutate(V1 = sub(".*/", "", V1))

meth <- read.delim(
  "../../Align_methyaltion_files/1158/re_methyaltion_qc_1157.tsv",
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  mutate(File = sub("model\\.pbmm2\\.combined\\.bed", "bam", File))

data <- merge(
  ml[, c("V1", "V4")],
  meth[, c("File", "Count")],
  by.x = "V1",
  by.y = "File"
)

colnames(data) <- c("Sample", "Coverage", "CpG_count")

data <- data %>%
  mutate(
    Coverage = as.numeric(Coverage),
    CpG_count = as.numeric(CpG_count)
  ) %>%
  filter(is.finite(Coverage), is.finite(CpG_count))

# ---------------------------
# Identify CpG-count outliers by IQR
# ---------------------------
q1 <- quantile(data$CpG_count, 0.25, na.rm = TRUE)
q3 <- quantile(data$CpG_count, 0.75, na.rm = TRUE)
iqr <- q3 - q1

lower_bound <- q1 - 1.5 * iqr
upper_bound <- q3 + 1.5 * iqr

data <- data %>%
  mutate(
    OutlierStatus = ifelse(
      CpG_count < lower_bound | CpG_count > upper_bound,
      "Outlier",
      "Non-outlier"
    )
  )

# ---------------------------
# Plot 1: coverage versus CpG count with outliers
# ---------------------------
p1 <- ggplot(data, aes(x = Coverage, y = CpG_count)) +
  geom_point(
    aes(fill = OutlierStatus),
    shape = 21,
    size = 1.8,
    color = "black",
    stroke = 0.15,
    alpha = 0.65
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    color = "black",
    fill = "grey80",
    linewidth = 0.55,
    alpha = 0.45
  ) +
  stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 3.2,
    color = "black",
    parse = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Non-outlier" = "#6baed6",
      "Outlier" = "#888de5"
    )
  ) +
  scale_y_continuous(
    labels = label_number(scale = 1e-6, suffix = "M", accuracy = 1),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  labs(
    x = "Genome coverage",
    y = "CpG call count",
    fill = NULL
  ) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 9),
    panel.grid = element_blank()
  )


# ---------------------------
# Remove low-quality samples
# ---------------------------
out <- read.table("../../lowQ_1158_self.log", stringsAsFactors = FALSE)

data_clean <- data %>%
  filter(!Sample %in% out$V1)

# ---------------------------
# Function to estimate saturation point
# ---------------------------
estimate_saturation <- function(df, threshold_fraction = 100) {
  fit <- loess(CpG_count ~ Coverage, data = df, span = 0.75)
  
  grid <- data.frame(
    Coverage = seq(
      min(df$Coverage, na.rm = TRUE),
      max(df$Coverage, na.rm = TRUE),
      length.out = 5000
    )
  )
  
  grid$pred <- predict(fit, newdata = grid)
  grid$slope <- c(NA, diff(grid$pred) / diff(grid$Coverage))
  
  threshold <- max(df$CpG_count, na.rm = TRUE) / threshold_fraction
  idx <- which(grid$slope < threshold)[1]
  
  if (is.na(idx)) {
    return(data.frame(target_x = NA_real_, target_y = NA_real_))
  }
  
  data.frame(
    target_x = grid$Coverage[idx],
    target_y = grid$pred[idx]
  )
}

target_all <- estimate_saturation(data_clean)

# ---------------------------
# Plot 2: cleaned coverage versus CpG count with saturation point
# ---------------------------
p2 <- ggplot(data_clean, aes(x = Coverage, y = CpG_count)) +
  geom_point(
    size = 1.7,
    alpha = 0.60,
    color = "grey30"
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    color = "black",
    fill = "grey80",
    linewidth = 0.55,
    alpha = 0.45
  ) +
  geom_vline(
    xintercept = target_all$target_x,
    linetype = "dashed",
    color = "#888de5",
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = target_all$target_y,
    linetype = "dashed",
    color = "#888de5",
    linewidth = 0.45
  ) +
  geom_point(
    data = target_all,
    aes(x = target_x, y = target_y),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.7,
    fill = "#888de5",
    color = "black",
    stroke = 0.25
  ) +
  annotate(
    "text",
    x = target_all$target_x,
    y = target_all$target_y,
    label = paste0(
      round(target_all$target_x, 2),
      "×; ",
      label_number(scale = 1e-6, suffix = "M", accuracy = 0.1)(target_all$target_y)
    ),
    hjust = -0.1,
    vjust = -0.8,
    size = 3.0,
    color = "black"
  ) +
  stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 3.2,
    color = "black",
    parse = FALSE
  ) +
  scale_y_continuous(
    labels = label_number(scale = 1e-6, suffix = "M", accuracy = 1),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  labs(
    x = "Genome coverage",
    y = "CpG call count"
  ) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    panel.grid = element_blank()
  )

ggsave(
  "Coverage_CpG_count_saturation_Nature.pdf",
  plot = p2,
  width = 4.2,
  height = 3.5,
  device = "pdf"
)

# ---------------------------
# Add instrument information
# ---------------------------
ins <- read.table(
  "../../ML-coverage/1158/instrument_model.log",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

data_clean <- data_clean %>%
  mutate(
    instrument = case_when(
      Sample %in% ins$filename[ins$instrument_model == "Revio"] ~ "Revio",
      Sample %in% ins$filename[ins$instrument_model == "Sequel II"] ~ "Sequel II",
      TRUE ~ NA_character_
    ),
    instrument = factor(instrument, levels = c("Sequel II", "Revio"))
  ) %>%
  filter(!is.na(instrument))

target_points <- data_clean %>%
  group_by(instrument) %>%
  group_modify(~{
    estimate_saturation(.x)
  }) %>%
  ungroup()

instrument_colors <- c(
  "Sequel II" = "#6baed6",
  "Revio" = "#888de5"
)

# ---------------------------
# Plot 3: instrument-stratified saturation
# ---------------------------
p3 <- ggplot(data_clean, aes(x = Coverage, y = CpG_count, color = instrument)) +
  geom_point(
    size = 1.6,
    alpha = 0.55
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    linewidth = 0.55,
    alpha = 0.25
  ) +
  geom_vline(
    data = target_points,
    aes(xintercept = target_x, color = instrument),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_hline(
    data = target_points,
    aes(yintercept = target_y, color = instrument),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    data = target_points,
    aes(x = target_x, y = target_y, fill = instrument),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.8,
    color = "black",
    stroke = 0.25
  ) +
  scale_color_manual(values = instrument_colors) +
  scale_fill_manual(values = instrument_colors) +
  scale_y_continuous(
    labels = label_number(scale = 1e-6, suffix = "M", accuracy = 1),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  labs(
    x = "Genome coverage",
    y = "CpG count",
    color = NULL,
    fill = NULL
  ) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

ggsave(
  "Coverage_CpG_count_by_instrument.pdf",
  plot = p3,
  width = 4.5,
  height = 3.5,
  device = "pdf"
)
