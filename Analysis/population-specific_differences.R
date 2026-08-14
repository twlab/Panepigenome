library(data.table)
library(limma)
library(matrixStats)

# ===============================
# User parameters
# ===============================
MAX_MISSING <- 0.10
DELTA_LATENT_01 <- 0.20   # effect-size threshold (0–1 scale)
EPS <- 1e-6

# ===============================
# Step 1: Load latent data (RAW)
# ===============================
latent_df <- fread(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/NN/Saifur/final_best_latent_hidden_4.tsv"
)

# Assumed format:
# gene | sample | latent
latent_mat <- latent_df

cpg_id <- latent_mat$gene
latent_raw <- as.matrix(latent_mat[, -1])
rownames(latent_raw) <- cpg_id

# ===============================
# Step 2: Load metadata
# ===============================
meta <- fread(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/Refmeth/continent.log",header=F,
  col.names = c("sample", "group", "sex")
)

meta[, sex := fifelse(sex == "female", 0,
                      fifelse(sex == "male",   1, NA_real_))]

# Align samples
common_samples <- intersect(colnames(latent_raw), meta$sample)
latent_raw <- latent_raw[, common_samples, drop = FALSE]
meta<-as.data.frame(meta)
rownames(meta)<-meta$sample
meta <- meta[common_samples,]

stopifnot(identical(colnames(latent_raw), meta$sample))

# ===============================
# Step 3: Missing-value filtering
# ===============================
keep <- rowMeans(is.na(latent_raw)) <= MAX_MISSING
latent_raw <- latent_raw[keep, ]

# ===============================
# Step 4: 0–1 scaling (EFFECT SIZE ONLY)
# ===============================
latent_min <- min(latent_raw, na.rm = TRUE)
latent_max <- max(latent_raw, na.rm = TRUE)

latent_01 <- (latent_raw - latent_min) /
  (latent_max - latent_min)

# ===============================
# Step 5: Pairwise population comparisons
# ===============================
populations <- sort(unique(meta$group))

for (pop1 in populations) {
  
  cat("Processing:", pop1, "vs Others\n")
  
  meta_sub <- meta
  
  latent_sub_raw <- latent_raw[, meta_sub$sample, drop = FALSE]
  
  grp <- ifelse(meta_sub$group == pop1, pop1, "Others")
  grp <- factor(grp, levels = c(pop1, "Others"))
  sex <- meta_sub$sex
  
  design <- model.matrix(~ grp + sex)
  
  fit <- lmFit(latent_sub_raw, design)
  fit <- eBayes(fit)
  
  fitted_raw <- fitted(fit)
  fitted_min <- min(fitted_raw, na.rm = TRUE)
  fitted_max <- max(fitted_raw, na.rm = TRUE)
  
  fitted_01 <- (fitted_raw - fitted_min) /
    (fitted_max - fitted_min)
  
  tab <- topTable(
    fit,
    coef = "grpOthers",
    number = Inf,
    sort.by = "none"
  )
  
  delta_latent_01 <- rowMedians(
    fitted_01[, grp == "Others", drop = FALSE],
    na.rm = TRUE
  ) -
    rowMedians(
      fitted_01[, grp == pop1, drop = FALSE],
      na.rm = TRUE
    )
  
  res <- data.table(
    cpg             = rownames(latent_sub_raw),
    beta_diff_raw   = tab$logFC,
    delta_latent_01 = delta_latent_01,
    p_value         = tab$P.Value,
    fdr             = tab$adj.P.Val
  )
  
  res_sig <- res[
    fdr < 0.05 &
      abs(delta_latent_01) >= DELTA_LATENT_01
  ]
  
  prefix <- sprintf("CpG_latent_limma_%s_vs_Others", pop1)
  
  fwrite(res, paste0(prefix, "_all.txt"), sep = "\t")
  fwrite(res_sig, paste0(prefix, "_sig_deltaLatent01_0.2.txt"), sep = "\t")
}
