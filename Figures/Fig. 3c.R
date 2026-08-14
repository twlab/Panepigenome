library(dplyr)
library(ggplot2)
pca_var<-read.table("PCA_variance_explained_variant.txt",header = T)
# Percentage variance for axis labels
pc1_var <- round(100 * pca_var[pca_var$PC == "PC1", "variance_explained"], 1)
pc2_var <- round(100 * pca_var[pca_var$PC == "PC2", "variance_explained"], 1)

pop <- read.table(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/Refmeth/continent.log",
  stringsAsFactors = FALSE
)
pc_scores<-read.table("PCA_scores_variant.txt",header=T)

# example: pca_scores$sample should be replaced using pop mapping
pc_scores$population <- ifelse(
  pc_scores$sample %in% pop$V1,
  pop$V2[match(pc_scores$sample, pop$V1)],
  pc_scores$sample
)

p <- ggplot(
  pc_scores,
  aes(x = PC1, y = PC2, color = population)
) +
  geom_point(
    size = 1.4,
    alpha = 0.85
  ) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Population"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

p<-p +
  stat_ellipse(
    type = "norm",
    linewidth = 0.8,
    level = 0.68   # ~1 SD, commonly used in genomics
  )

ggsave(
  "PCA_PC1_PC2_con.pdf",
  plot = p,
  width = 6,
  height = 4,
  useDingbats = FALSE
)
