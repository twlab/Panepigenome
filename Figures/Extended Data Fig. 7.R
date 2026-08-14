library(dplyr)
library(ggplot2)

pca_var <- read.table("PCA_variance_explained_variant.txt", header = TRUE)

pc1_var <- round(100 * pca_var[pca_var$PC == "PC1", "variance_explained"], 1)
pc2_var <- round(100 * pca_var[pca_var$PC == "PC3", "variance_explained"], 1)

pop <- read.table(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/Refmeth/continent.log",
  stringsAsFactors = FALSE
)

pc_scores <- read.table("PCA_scores_variant.txt", header = TRUE)

pc_scores$population <- ifelse(
  pc_scores$sample %in% pop$V1,
  pop$V2[match(pc_scores$sample, pop$V1)],
  pc_scores$sample
)

pc_scores$population <- factor(
  pc_scores$population,
  levels = c("AFR", "AMR", "EAS", "EUR", "SAS")
)

continent_colors <- c(
  AFR = "#0072B2",
  AMR = "#D55E00",
  EAS = "#009E73",
  EUR = "#CC79A7",
  SAS = "#E69F00"
)

p <- ggplot(
  pc_scores,
  aes(x = PC1, y = PC3, color = population)
) +
  geom_point(
    size = 1.4,
    alpha = 0.85
  ) +
  stat_ellipse(
    type = "norm",
    linewidth = 0.8,
    level = 0.68
  ) +
  scale_color_manual(values = continent_colors) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC3 (", pc2_var, "%)"),
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

ggsave(
  "PCA_PC1_PC3_con.pdf",
  plot = p,
  width = 6,
  height = 4,
  device = cairo_pdf
)
