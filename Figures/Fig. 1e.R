# Load required libraries
library(ggplot2)
library(dplyr)
library(ggbeeswarm)  # for stacked/"tree-like" points

# Read data
df <- read.table("merged.tsv", header = FALSE)
colnames(df) <- c("Sample", "CpG_Count", "Continent", "Sex")  # assume V4 is Sex
options(scipen = 999)
mean(df$CpG_Count)
# 31286564
median(df$CpG_Count)
# 31343583

df$CpG_Count <- df$CpG_Count / 1e6  # convert to millions

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

# Reorder continents by median CpG count (across all sexes)
# df$Continent_Full <- factor(df$Continent_Full,
#                             levels = df %>%
#                               group_by(Continent_Full) %>%
#                               summarise(median_CpG = median(CpG_Count), .groups="drop") %>%
#                               arrange(median_CpG) %>%
#                               pull(Continent_Full))

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

# # A tibble: 5 × 4
# Continent_Full    female      male female_minus_male
# <fct>              <dbl>     <dbl>             <dbl>
#   1 South Asia     32059114. 31771566.           287549 
# 2 East Asia      32106165  31792402.           313764.
# 3 Europe         32081344. 31703084.           378260 
# 4 America        32248782  31722556.           526226.
# 5 Africa         32252188  32001595            250593 

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

# Plot: tree-like points, separate panels for Male/Female
p <- ggplot(df, aes(x = Continent_Full, y = CpG_Count, fill = Continent_Full)) +
  geom_quasirandom(method = "smiley", shape = 21, size = 2.2, color = "black",
                   stroke = 0.3, alpha = 0.6) +
  geom_point(data = medians, aes(x = Continent_Full, y = median_CpG),
             color = "black", shape = 95, size = 8,alpha = 0.6) +
  geom_hline(yintercept = 27.950835, linetype = "dashed", color = "#76B7B2", size = 1.2,alpha = 0.6) +
  geom_hline(yintercept = 29.686360, linetype = "dashed", color = "#B07AA1", size = 1.2,alpha = 0.6) +
  scale_fill_manual(values = continent_colors) +
  scale_x_discrete(labels = function(x) continent_abbr[x]) +
  labs(x = NULL, y = "CpGs per haplotype (Million)") +
  facet_wrap(~Sex) +
  coord_cartesian(ylim = c(27.9, 32.5)) +
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
ggsave("CpG_Count_tree_dots_female_male_facets.pdf", plot = p, width = 6, height = 3)
