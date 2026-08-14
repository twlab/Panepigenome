library(dplyr)
library(ggplot2)
library(ggridges)

# -------------------------
# 1. load data
# -------------------------
fileA_path <- "/scratch/zdong/Projects/PanEpiG/V1-9/PAV_variants/Phased/Merge/Reference-basedCpG/Lifover/INS/Freq/final_counts_means_filtered.tsv"  # 替换成你的真实路径
fileB_path <- "/scratch/zdong/Projects/PanEpiG/V1-9/PAV_variants/Phased/Merge/Reference-basedCpG/Lifover/DEL/Freq/final_counts_means_filtered.tsv"

# load data
A <- read.table(fileA_path, stringsAsFactors = FALSE, fill = TRUE,sep="\t")
B <- read.table(fileB_path, stringsAsFactors = FALSE, fill = TRUE,sep="\t")


colnames(A) <- colnames(B) <- c("id",  "x", "y", "value","type")


A$source <- "CpG gain"
B$source <- "CpG loss"

# merge
dat <- bind_rows(A, B)
dat$value<-dat$value/440
dat<-dat[dat$value<=1,]
mean(dat$value)
# [1] 0.04087214
mean(dat$value[dat$type=="SNV"])
# [1] 0.051506
mean(dat$value[dat$type=="IINS" | dat$type=="IDEL"])
# [1] 0.05171943
mean(dat$value[dat$type=="INS" | dat$type=="DEL"])
# [1] 0.02094121

# -------------------------
# 2. class
# -------------------------
dat <- dat %>%
  mutate(
    class = case_when(
      type == "SNV" ~ "SNP",
      type %in% c("IINS", "IDEL") ~ "Indel",
      type %in% c("INS", "DEL") ~ "SV"
    )
  ) %>%
  filter(!is.na(class))

dat$class <- factor(dat$class, levels = c("SNP", "Indel", "SV"))

# -------------------------
# 3. median
# -------------------------
med <- dat %>%
  group_by(class, source) %>%
  summarise(med = median(value), .groups = "drop")

write.table(
  med,
  file = "med_values.txt",  # 输出文件名
  sep = "\t",               # 使用 tab 分隔
  row.names = FALSE,        # 不输出行号
  quote = FALSE             # 不加引号
)


# -------------------------
# 4. plot
# -------------------------
p <- ggplot(dat, aes(x = value, y = class, fill = source)) +
  geom_density_ridges(
    alpha = 0.6,
    scale = 1.3,
    rel_min_height = 0.01,
    color = "grey30",
    size = 0.3
  ) +
  # geom_vline(
  #   data = med,
  #   aes(xintercept = med, color = source),
  #   linetype = "dashed",
  #   size = 0.6,
  #   show.legend = FALSE
  # ) +
  scale_fill_manual(
    values = c("CpG gain" = "#E39C63", "CpG loss" = "#D7D894")
  ) +
  scale_color_manual(
    values = c("CpG gain" = "#E39C63", "CpG loss" = "#D7D894")
  ) +
  theme_classic(base_size = 14) +
  labs(x = "Variant-CpG gain/loss frequency", y = NULL, fill = NULL) +
  theme(
    legend.position = "top",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# -------------------------
# 5. save
# -------------------------
# print(p)

ggsave(
  "ridgeplot_var-cpg_freq.pdf",
  p,
  width = 5.5,
  height = 3.5,
  dpi = 300
)

### count
library(ggridges)
library(dplyr)
library(forcats)
library(tidyverse)

dat <- dat %>%
  mutate(
    class = factor(class, levels = c("SV", "Indel", "SNP")),
    source = factor(source, levels = c("CpG gain", "CpG loss"))
  )

p <- ggplot(
  dat,
  aes(
    x = value,
    fill = source,
    colour = source
  )
) +
  geom_histogram(
    aes(y = after_stat(count / 1e6)),
    bins = 50,
    position = "identity",
    alpha = 0.45,
    linewidth = 0.3
  ) +
  facet_grid(
    rows = vars(class),
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "CpG gain" = "#B07AA1",
      "CpG loss" = "#76B7B2"
    )
  ) +
  scale_colour_manual(
    values = c(
      "CpG gain" = "#B07AA1",
      "CpG loss" = "#76B7B2"
    )
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 2),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    x = "var-CpG gain/loss frequency",
    y = "Count (million)",
    fill = NULL,
    colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black", linewidth = 0.35),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(colour = "black", size = 11),
    strip.background = element_blank(),
    strip.text.y = element_text(
      colour = "black",
      size = 11,
      face = "bold",
      angle = 0
    ),
    panel.spacing.y = unit(0.35, "cm")
  )

ggsave(
  "count_var-cpg_freq.pdf",
  plot = p,
  width = 3.2,
  height = 3.2
)

##### ======================

# pie 
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# dat is your input tibble/data.frame
# Example columns:
# id, x, y, value, type, source, class

# -----------------------------
# 1. Categorize value
# -----------------------------
plot_df <- dat %>%
  mutate(
    class = factor(class, levels = c("SV", "Indel", "SNP")),
    value_group = case_when(
      value < 0.1 ~ "<0.1",
      value > 0.9 ~ ">0.9",
      TRUE ~ "0.1–0.9"
    ),
    value_group = factor(
      value_group,
      levels = c("<0.1", "0.1–0.9", ">0.9")
    )
  ) %>%
  dplyr::count(class, source, value_group, name = "n") %>%
  tidyr::complete(class, source, value_group, fill = list(n = 0)) %>%
  group_by(class, source) %>%
  mutate(
    prop = n / sum(n),
    ymax = cumsum(prop),
    ymin = lag(ymax, default = 0),
    label_pos = (ymax + ymin) / 2,
    label = ifelse(
      prop >= 0.05,
      scales::percent(prop, accuracy = 0.1),
      ""
    )
  ) %>%
  ungroup()

# -----------------------------
# 2. Colors (clean, publication style)
# -----------------------------
pie_cols <- c(
  "<0.1"    = "#FDAE61",   # blue
  "0.1–0.9" = "#67A9CF",   # grey
  ">0.9"    = "#2166AC"    # red
)

# -----------------------------
# 3. Pie plot
# -----------------------------
p <- ggplot(plot_df, aes(x = 1, y = prop, fill = value_group)) +
  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.5
  ) +
  coord_polar(theta = "y") +
  facet_grid(class ~ source) +
  geom_text(
    aes(y = label_pos, label = label),
    color = "black",
    size = 3
  ) +
  scale_fill_manual(values = pie_cols) +
  labs(
    fill = "Value range"
  ) +
  theme_void(base_size = 12) +
  theme(
    strip.text = element_text(size = 11, face = "bold", color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.spacing = unit(1, "lines")
  )


# -----------------------------
# 4. Save
# -----------------------------
ggsave(
  "value_proportion_pie_by_class.pdf",
  plot = p,
  width = 6.5,
  height = 4.8,
  device = cairo_pdf
)


