library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 13))

# ---------------------------
# Read data
# ---------------------------
df <- read.table(
  "pairwise_pop_specific_CpGs.p.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ---------------------------
# Clean and extract pair
# ---------------------------
pop_order <- c("AFR", "AMR", "EAS", "EUR", "SAS")

df2 <- df %>%
  rename(
    file = Content_comparsion,
    Model = Type,
    Count = `Pop-specific CpGs`
  ) %>%
  mutate(
    Model = recode(
      Model,
      "methyaltion difference" = "Shared CpG methylation",
      "methylation difference" = "Shared CpG methylation",
      "CpG gain/loss" = "CpG gain/loss",
      "Intergrated" = "Integrated model",
      "Integrated" = "Integrated model"
    ),
    Model = factor(
      Model,
      levels = c(
        "Shared CpG methylation",
        "CpG gain/loss",
        "Integrated model"
      )
    ),
    pair_raw = str_extract(file, "(AFR|AMR|EAS|EUR|SAS)_vs_(AFR|AMR|EAS|EUR|SAS)"),
    pop1 = str_extract(pair_raw, "^[A-Z]+"),
    pop2 = str_extract(pair_raw, "(?<=_vs_)[A-Z]+"),
    pop1_i = match(pop1, pop_order),
    pop2_i = match(pop2, pop_order),
    PopA = ifelse(pop1_i < pop2_i, pop1, pop2),
    PopB = ifelse(pop1_i < pop2_i, pop2, pop1),
    Count = as.numeric(Count),
    log10_count = log10(Count + 1)
  ) %>%
  filter(!is.na(PopA), !is.na(PopB), !is.na(Model))

# ---------------------------
# Complete all pairwise combinations
# ---------------------------
all_pairs <- expand_grid(PopA = pop_order, PopB = pop_order) %>%
  mutate(
    i = match(PopA, pop_order),
    j = match(PopB, pop_order)
  ) %>%
  filter(i < j) %>%
  select(PopA, PopB)

df_plot <- df2 %>%
  select(Model, PopA, PopB, Count, log10_count) %>%
  right_join(
    expand_grid(
      Model = levels(df2$Model),
      all_pairs
    ),
    by = c("Model", "PopA", "PopB")
  ) %>%
  mutate(
    PopA = factor(PopA, levels = rev(pop_order)),
    PopB = factor(PopB, levels = pop_order),
    Model = factor(
      Model,
      levels = c(
        "Shared CpG methylation",
        "CpG gain/loss",
        "Integrated model"
      )
    ),
    label = ifelse(is.na(Count), "", comma(Count))
  )

# ---------------------------
# Log10 heatmap
# ---------------------------
p <- ggplot(df_plot, aes(x = PopB, y = PopA, fill = log10_count)) +
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = label),
    size = 2.8,
    color = "black"
  ) +
  facet_wrap(
    ~Model,
    nrow = 1
  ) +
  scale_fill_gradient(
    low = "grey95",
    high = "#B2182B",
    na.value = "white",
    name = expression(log[10]~"(CpGs + 1)")
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      size = 11,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank()
  )

# p

ggsave(
  "pairwise_population_specific_CpGs_heatmap_log10.p.pdf",
  plot = p,
  width = 8.4,
  height = 3.4,
  device = cairo_pdf
)
