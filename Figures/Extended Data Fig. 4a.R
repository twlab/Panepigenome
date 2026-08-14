# Load required libraries
library(ggplot2)
library(dplyr)

# Read data
df <- read.table("merged_pop.tsv", header = FALSE)
df <- df[,1:3]
colnames(df) <- c("Sample", "CpG_Count", "Pop")
df$CpG_Count <- df$CpG_Count / 1e6  # convert to millions

# Exclude specific samples (optional)
df <- df[!df$Sample %in% c("HG00272_hap1_hprc_r2_v1.0.1.CpG_sites.bed",
                           "HG00272_hap2_hprc_r2_v1.0.1.CpG_sites.bed"), ]

# === Population and region info ===
pop_info <- data.frame(
  pop = c("CDX", "CHB", "JPT", "KHV", "CHS",
          "BEB", "GIH", "ITU", "PJL", "STU",
          "ASW", "ACB", "ESN", "GWD", "LWK", "MSL", "YRI",
          "GBR", "FIN", "IBS", "TSI", "CEU",
          "CLM", "MXL", "PEL", "PUR", "MKK"),
  region = c("East Asia", "East Asia", "East Asia", "East Asia", "East Asia",
             "South Asia", "South Asia", "South Asia", "South Asia", "South Asia",
             "Africa", "Africa", "Africa", "Africa", "Africa", "Africa", "Africa",
             "Europe", "Europe", "Europe", "Europe", "Europe",
             "America", "America", "America", "America", "Africa")
)

# Merge population info
df <- left_join(df, pop_info, by = c("Pop" = "pop"))

# === Order continents by median CpG ===
region_order <- df %>%
  group_by(region) %>%
  summarise(region_median = median(CpG_Count), .groups = "drop") %>%
  arrange(region_median) %>%
  pull(region)

df$region <- factor(df$region, levels = region_order)

# === Order populations within region by median CpG ===
pop_order <- df %>%
  group_by(region, Pop) %>%
  summarise(median_CpG = median(CpG_Count), .groups = "drop") %>%
  arrange(factor(region, levels = region_order), median_CpG)

df$Pop <- factor(df$Pop, levels = pop_order$Pop)

# === Continent colors ===
region_colors <- c(
  "East Asia" = "#2ca02c",
  "South Asia" = "#ff7f0e",
  "Africa" = "#1f77b4",
  "Europe" = "#9467bd",
  "America" = "#d62728"
)

# === boxplot (no facets) ===
p <- ggplot(df, aes(x = Pop, y = CpG_Count, fill = region)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "black", linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  geom_hline(yintercept = 29.303965, linetype = "dashed", color = "grey50", size = 0.6) +
  scale_fill_manual(values = region_colors) +
  labs(x = NULL, y = "Number of CpGs per haplotype (Million)") +
  coord_cartesian(ylim = c(29, 34)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 13),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3)
  )

# Save figure
ggsave("CpG_Count_boxplot_by_pop_grouped_region.pdf", plot = p, width = 10, height = 4)
