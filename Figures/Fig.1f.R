# Load required libraries
library(ggplot2)
library(dplyr)
library(ggbeeswarm)  # for stacked/"tree-like" points

# Read data
df <- read.table("merged.tsv", header = FALSE)
colnames(df) <- c("Sample", "Count","CpG_Count", "Continent", "Sex")  # assume V4 is Sex
options(scipen = 999)
shapiro.test(df$CpG_Count)
# p-value = 0.000000772
median(df$CpG_Count)
# 80.3

summary(df$CpG_Count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 63.90   75.58   80.30   79.28   82.83   89.50 

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

diffs
# # A tibble: 5 × 5
# Continent_Full female  male female_minus_male Continent_Abbr
# <fct>           <dbl> <dbl>             <dbl> <fct>         
#   1 Europe           79.0  78.2             0.800 EUR           
# 2 Africa           79.7  78.4             1.30  AFR           
# 3 East Asia        78.2  79.7            -1.50  EAS           
# 4 South Asia       82.1  80.8             1.30  SAS           
# 5 America          81.2  81.6            -0.350 AMR 

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

library(ggplot2)
library(dplyr)
library(ggbeeswarm)

# --- Nature-style Violin + Boxplot + Beeswarm ---
p <- ggplot(df, aes(x = Continent_Full, y = CpG_Count)) +
  
  # Transparent violin, thin black border (Nature style)
  geom_violin(
    trim = TRUE,
    scale = "width",
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  
  # Internal compact boxplot (thin, subtle)
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    linewidth = 0.35,
    fill = "white",
    color = "black"
  ) +
  
  # Beeswarm (dots) with continent colors → much nicer
  geom_quasirandom(
    aes(color = Continent_Full),
    size = 1.7,
    width = 0.25,
    alpha = 0.6
  ) +
  
  # Median line (highlighted)
  geom_point(
    data = medians,
    aes(x = Continent_Full, y = median_CpG),
    shape = 95,
    size = 9,
    color = "black",alpha = 0.6
  ) +
  
  # Horizontal reference line
  geom_hline(
    yintercept = 78.2,
    linetype = "dashed",
    color = "#76B7B2",
    linewidth = 1.2,alpha = 0.6
  ) +
  
  # Horizontal reference line
  geom_hline(
    yintercept = 78.3,
    linetype = "dashed",
    color = "#B07AA1",
    linewidth = 1.2,alpha = 0.6
  ) +
  
  # Apply fixed continent colors (only for points)
  scale_color_manual(values = continent_colors) +
  
  # X-axis labels: use EAS/SAS/AFR/EUR/AMR
  scale_x_discrete(labels = function(x) continent_abbr[x]) +
  
  # Facet by sex
  facet_wrap(~ Sex) +
  coord_cartesian(ylim = c(50, 100)) +  # y-axis ends at 100
  # Beautiful Nature theme
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    panel.spacing = unit(0.7, "lines"),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Median methylation level"
  )

# Save
ggsave(
  "CpG_methylation_violin_boxplot_female_male_facets.pdf",
  plot = p, width = 6, height = 3
)
