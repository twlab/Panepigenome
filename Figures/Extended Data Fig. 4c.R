# Load required libraries
library(ggplot2)
library(dplyr)

# Read data
df <- read.table("merged_pop.tsv", header = FALSE)
colnames(df) <- c("Sample", "Count","CpG_Count", "Pop", "Sex")

df<-df[-which(df$Pop %in% c("CH","EUR")),]

# Population and region info
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

# Assuming your data is in a data frame called 'df'
medians <- df %>%
  group_by(Sex, Pop,region) %>%         # group by Sex and population/region
  summarise(median_CpG = median(CpG_Count), .groups = "drop")  # calculate median

medians
write.table(medians, "median_CpG_by_sex_population.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

mad_by_sex_region <- df %>%
  group_by(Sex, region) %>%
  summarise(
    MAD_CpG = mad(CpG_Count),
    MAD_SD = mad(CpG_Count) * 1.4826,
    .groups = "drop"
  )

mad_by_sex_region
write.table(mad_by_sex_region, "mad_CpG_by_sex_population.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

mad_re <- medians %>%
  group_by(region) %>%
  summarise(
    MAD_CpG = mad(median_CpG),
    MAD_SD = mad(median_CpG) * 1.4826,
    .groups = "drop"
  )
mad_re
write.table(mad_re , "mad_CpG_re_nosex.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

mad_re <- medians %>%
  group_by(Sex, region) %>%
  summarise(
    MAD_CpG = mad(median_CpG),
    MAD_SD = mad(median_CpG) * 1.4826,
    .groups = "drop"
  )

mad_re 
write.table(mad_re , "mad_CpG_re.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

medians1<-medians[-which(medians$Pop %in% c('ASW','ACB')),]
mad_re <- medians1 %>%
  group_by(Sex, region) %>%
  summarise(
    MAD_CpG = mad(median_CpG),
    MAD_SD = mad(median_CpG) * 1.4826,
    .groups = "drop"
  )

mad_re 
write.table(mad_re , "mad_CpG_re_removeAA.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# --- Order continents by median CpG across all samples ---
region_order <- df %>%
  group_by(region) %>%
  summarise(region_median = median(CpG_Count), .groups = "drop") %>%
  arrange(region_median) %>%
  pull(region)

df$region <- factor(df$region, levels = region_order)

# --- Order populations within each continent by median CpG across all samples ---
pop_order <- df %>%
  group_by(region, Pop) %>%
  summarise(median_CpG = median(CpG_Count), .groups = "drop") %>%
  arrange(factor(region, levels = region_order), median_CpG)

df$Pop <- factor(df$Pop, levels = pop_order$Pop)

# Continent colors
region_colors <- c(
  "East Asia" = "#2ca02c",
  "South Asia" = "#ff7f0e",
  "Africa" = "#1f77b4",
  "Europe" = "#9467bd",
  "America" = "#d62728"
)

# boxplot with facets for Sex
p <- ggplot(df, aes(x = Pop, y = CpG_Count, fill = region)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "black", linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  geom_hline(yintercept = 78, linetype = "dashed", color = "grey50", size = 0.6) +
  scale_fill_manual(values = region_colors) +
  labs(x = NULL, y = "Methylation level") +
  coord_cartesian(ylim = c(50, 100)) +
  facet_wrap(~Sex, scales = "free_x") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 13),
    legend.position = "right",
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing.x = unit(0.5, "lines")
  )

# Save figure
ggsave("CpG_methyaltion_boxplot_by_pop_grouped_region_combined_median.pdf", plot = p, width = 10, height = 4)
