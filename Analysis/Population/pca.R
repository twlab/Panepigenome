library(data.table)
library(ggplot2)

# ===========================
# 1️⃣ Load variant matrix
# ===========================
dong <- fread(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/NN/Saifur/final_best_latent_hidden_4.tsv",
  header = TRUE 
)

# Variant ID
dong[, sample := sapply(strsplit(gene, ":"), function(x) paste(x[1:3], collapse=":"))]

# Keep genotype columns + variant ID
dong <- dong[, 2:ncol(dong)]

# ===========================
# 2️⃣ Missing rate cutoff (≤10%)
# ===========================
dong <- dong[
  rowMeans(is.na(dong[, 1:(ncol(dong)-1), with=FALSE])) <= 0.10
]

colSums(is.na(dong))

# ===========================
# 3️⃣ Metadata
# ===========================
meta <- fread(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/Refmeth/continent.log",
  header = FALSE
)

setnames(meta, c("sample", "population","sex"))
setkey(meta, sample)

# ===========================
# 4️⃣ Build PCA matrix (samples × variants)
# ===========================
latent_mat <- as.matrix(dong[, !("sample"), with=FALSE])
latent_mat <- t(latent_mat)   # samples × variants

# Assign sample IDs
sample_ids <- colnames(dong[, !("sample"), with=FALSE])
rownames(latent_mat) <- sample_ids

# Align metadata
meta <- meta[J(sample_ids)]

# ===========================
# 5️⃣ PCA preprocessing
# ===========================
# Remove zero-variance variants
keep_cols <- apply(latent_mat, 2, var, na.rm = TRUE) > 0
latent_mat <- latent_mat[, keep_cols, drop = FALSE]

# Replace NA / Inf with column means
for (j in seq_len(ncol(latent_mat))) {
  x <- latent_mat[, j]
  bad <- !is.finite(x)
  if (any(bad)) {
    latent_mat[bad, j] <- mean(x[is.finite(x)])
  }
}

# ===========================
# 6️⃣ PCA
# ===========================
pca <- prcomp(
  latent_mat,
  center = TRUE,
  scale. = TRUE   # CORRECT for genotype / variant data
)

# ===========================
# 7️⃣ PCA scores & variance
# ===========================
pc_scores <- as.data.table(pca$x)
pc_scores[, sample := sample_ids]
pc_scores[, population := meta$population]

var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

# Output
fwrite(pc_scores, "PCA_scores_variant.txt", sep = "\t")

pca_var <- data.table(
  PC = paste0("PC", seq_along(var_explained)),
  variance_explained = var_explained,
  cumulative = cumsum(var_explained)
)

fwrite(pca_var, "PCA_variance_explained_variant.txt", sep = "\t")

# ===========================
# 8️⃣ PCA plot (Nature style)
# ===========================
pc1_var <- round(100 * pca_var[PC == "PC1", variance_explained], 1)
pc2_var <- round(100 * pca_var[PC == "PC2", variance_explained], 1)

p <- ggplot(
  pc_scores,
  aes(x = PC1, y = PC2, color = population)
) +
  geom_point(size = 1.4, alpha = 0.85) +
  stat_ellipse(
    type = "norm",
    linewidth = 0.8,
    level = 0.68
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
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

ggsave(
  "PCA_PC1_PC2_variant_NatureStyle.pdf",
  plot = p,
  width = 5.5,
  height = 4.5,
  useDingbats = FALSE
)
