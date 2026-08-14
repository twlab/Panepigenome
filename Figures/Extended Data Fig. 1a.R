# Load required libraries
library(ggplot2)
library(dplyr)
library(ggbeeswarm)  # for stacked/"tree-like" points

# Read data
df <- read.table("merged_density.tsv", header = FALSE)
df$V2<-df$V2/df$V6*1000
df<-df[,1:4]
colnames(df) <- c("Sample", "CpG_Count", "Continent", "Sex")  # assume V4 is Sex
options(scipen = 999)
mean(df$CpG_Count)


wilcox.test(df$CpG_Count,mu=29303965/2934876451*1000)

wilcox.test(df$CpG_Count,mu=33839747/3114496529*1000)


# Optional: exclude specific samples
exclude_samples <- c("HG00272_hap1_hprc_r2_v1.0.1.CpG_sites.bed",
                     "HG00272_hap2_hprc_r2_v1.0.1.CpG_sites.bed")
df <- df[!df$Sample %in% exclude_samples, ]

# Map short codes to full continent names
continent_full <- c("EAS" = "East Asia",
                    "SAS" = "South Asia",
                    "AFR" = "Africa",
                    "EUR" = "Europe",
                    "AMR" = "America")
df$Continent_Full <- recode(df$Continent, !!!continent_full)


# df$CpG_Count <- df$CpG_Count / 1e6  # convert to millions

# Median per continent and sex
medians <- df %>%
  group_by(Sex, Continent_Full) %>%
  summarise(median_CpG = median(CpG_Count), .groups="drop")


diffs <- medians %>%
  pivot_wider(names_from = Sex, values_from = median_CpG) %>%  # make columns female and male
  mutate(female_minus_male = female - male,
         Continent_Abbr = recode(Continent_Full,
                                 "East Asia" = "EAS",
                                 "South Asia" = "SAS",
                                 "Africa" = "AFR",
                                 "Europe" = "EUR",
                                 "America" = "AMR"))  

diffs


# Custom colors
continent_colors1 <- c(
  "EAS" = "#2ca02c", 
  "SAS" = "#ff7f0e", 
  "AFR" = "#1f77b4", 
  "EUR" = "#9467bd", 
  "AMR" = "#d62728"
)

# Continent colors
continent_colors <- c(
  "East Asia" = "#2ca02c", 
  "South Asia" = "#ff7f0e", 
  "Africa" = "#1f77b4", 
  "Europe" = "#9467bd", 
  "America" = "#d62728"
)

# Continent abbreviations
continent_abbr <- c("East Asia" = "EAS",
                    "South Asia" = "SAS",
                    "Africa" = "AFR",
                    "Europe" = "EUR",
                    "America" = "AMR")

# Nature-style bar plot
p <- ggplot(diffs, aes(x = Continent_Abbr, y = female_minus_male, fill = Continent_Abbr)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = female_minus_male), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = continent_colors1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Population",
    y = "CpG Density Difference (Female - Male)"
  )

# Save figure
ggsave("CpG_Density_Difference_Female_Male.pdf", plot = p, width = 7, height = 4)

# Plot: tree-like points, separate panels for Male/Female
p <- ggplot(df, aes(x = Continent_Full, y = CpG_Count, fill = Continent_Full)) +
  geom_quasirandom(method = "smiley", shape = 21, size = 2.2, color = "black",
                   stroke = 0.3, alpha = 0.8) +
  geom_point(data = medians, aes(x = Continent_Full, y = median_CpG),
             color = "black", shape = 95, size = 8) +
  geom_hline(yintercept = 29303965/2934876451*1000, linetype = "dashed", color = "#76B7B2", size = 1.2) +
  geom_hline(yintercept = 33839747/3114496529*1000, linetype = "dashed", color = "#B07AA1", size = 1.2) +
  scale_fill_manual(values = continent_colors) +
  scale_x_discrete(labels = function(x) continent_abbr[x]) +
  labs(x = NULL, y = "CpGs per kilobase") +
  facet_wrap(~Sex) +
  # coord_cartesian(ylim = c(29, 34)) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 11),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.spacing.x = unit(0.4, "lines"),
    strip.background = element_rect(fill = "white", color = NA)
  )

# Save
ggsave("CpG_Count_tree_dots_female_male_facets_per1kb.pdf", plot = p, width = 7, height = 4)
