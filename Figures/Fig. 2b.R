library("tidyr")
library("ggplot2")
library(data.table)

df<-fread("re_r2ae_MAF5.txt",header=T)

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = c(r2, r2_pop),
    names_to = "Metric",
    values_to = "R2"
  )

p<-ggplot(df_long, aes(x = R2, fill = Metric)) +
  geom_histogram(
    position = "identity",
    bins = 500,
    alpha = 0.6,
    color = NA
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "Proportion of variance explained",
    y = "Number of var-CpGs"
  ) +
  scale_fill_manual(
    values = c("r2" = "#f4a582", "r2_pop" = "#92c5de"),
    labels = c("Contential group", "Population")
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

ggsave("hist_nonref_ae_beta_maf5.pdf", plot = p, width = 5, height = 4)
