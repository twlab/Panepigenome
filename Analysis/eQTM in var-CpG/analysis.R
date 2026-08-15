library(data.table)
library(future)
library(future.apply)
library(mediation)

# ==========================================
# 1. Faster Reading
# ==========================================
meth <- fread("input_meth_maf5%.part_29.csv")
var  <- fread("input_var_maf5%.part_29.csv")
gene <- fread("input_expr_maf5%.part_29.csv")
meta <- fread("/scratch/zdong/Projects/PanEpiG/V1-9/Figures/Clinical/PIEZO1/meta.txt")

# ==========================================
# 2. Transpose Meta and Synchronize Samples
# ==========================================
sample_ids <- colnames(meta)[-1] 
meta_df <- as.data.frame(t(meta[, -1, with=FALSE]))
colnames(meta_df) <- meta$ID    
meta_df$sample_id <- sample_ids

# Identify common samples across all 4 datasets
common_samples <- Reduce(intersect, list(colnames(meth), colnames(var), colnames(gene), meta_df$sample_id))
common_samples <- setdiff(common_samples, c("start", "end", "ID"))

if(length(common_samples) == 0) stop("Zero common samples found! Check sample ID naming.")

rownames(meta_df) <- meta_df$sample_id
meta_df <- meta_df[common_samples, ]
meta_df$sample_id <- NULL 

# Convert data to matrices for fast row-access
meth_mat <- as.matrix(meth[, ..common_samples])
var_mat  <- as.matrix(var[, ..common_samples])
gene_mat <- as.matrix(gene[, ..common_samples])

# Beta -> M-value
beta <- meth_mat / 100
eps  <- 1e-6
meth_mat <- log2((beta + eps) / (1 - beta + eps))
rm(beta); gc()

# ==========================================
# 3. The Optimized Function
# ==========================================
run_one_dt_fast <- function(n, g_vec, m_vec, v_vec, cov_base, id_meth, id_gene) {
  
  # Remove samples with NA in main variables
  keep_idx <- !is.na(g_vec) & !is.na(m_vec) & !is.na(v_vec)
  
  # Requirement: Minimum samples and variance in genotype
  # if(sum(keep_idx) < 10 || var(v_vec[keep_idx], na.rm=TRUE) == 0) return(NULL)
  
  # Build model dataframe
  df <- cbind(cov_base[keep_idx, ], 
              expr = g_vec[keep_idx], 
              meth = m_vec[keep_idx], 
              geno = v_vec[keep_idx])
  
  # --- Models for Filtering ---
  fit.M <- lm(meth ~ geno + ., data = df[, !names(df) %in% c("expr","InferredCov1")])
  fit.Y <- lm(expr ~ geno + meth + ., data = df)
  
  cM <- coef(summary(fit.M))
  cY <- coef(summary(fit.Y))
  
  p_M_geno <- if("geno" %in% rownames(cM)) cM["geno", 4] else 1
  p_Y_meth <- if("meth" %in% rownames(cY)) cY["meth", 4] else 1
  
  # --- Mediation (Only if geno -> meth is promising) ---
  med_res <- list(d0 = NA, d0.p = NA, z0 = NA, z0.p = NA, tau.coef = NA, tau.p = NA)
  
  if (p_M_geno < 0.1) {
    try({
      # First mediation run with 100 bootstrap sims
      med <- mediate(fit.M, fit.Y, treat = "geno", mediator = "meth", boot = TRUE, sims = 100)
      
      # Extract ACME and other mediation results
      med_res <- list(
        d0 = med$d0,
        d0.p = med$d0.p,
        z0 = med$z0,
        z0.p = med$z0.p,
        tau.coef = med$tau.coef,
        tau.p = med$tau.p,
        ACME = med$d0,      # or med$ACME if your version has it
        ACME_p = med$d0.p
      )
      
      # If ACME p-value is 0, redo with more sims
      if (!is.na(med$d0.p) && med$d0.p == 0) {
        med <- mediate(fit.M, fit.Y, treat = "geno", mediator = "meth", boot = TRUE, sims = 10000)
        
        med_res <- list(
          d0 = med$d0,
          d0.p = med$d0.p,
          z0 = med$z0,
          z0.p = med$z0.p,
          tau.coef = med$tau.coef,
          tau.p = med$tau.p,
          ACME = med$d0,
          ACME_p = med$d0.p
        )
      }
      
    }, silent = TRUE)
  }
  # --- Remaining Models ---
  fit.E <- lm(expr ~ meth + ., data = df[, !names(df) %in% c("geno","pc1","pc2","pc3","pc4","pc5")])
  fit.I <- lm(expr ~ geno + meth + geno:meth + ., data = df)
  
  cE <- coef(summary(fit.E))
  cI <- coef(summary(fit.I))
  
  get_val <- function(ct, term, stat) {
    if (!is.null(ct) && term %in% rownames(ct)) ct[term, stat] else NA_real_
  }
  
  return(list(
    id = id_meth, id2 = id_gene,
    E_meth_beta = get_val(cE, "meth", 1), E_meth_p = get_val(cE, "meth", 4),
    M_geno_beta = get_val(cM, "geno", 1), M_geno_p = p_M_geno,
    Y_meth_beta = get_val(cY, "meth", 1), Y_meth_p = p_Y_meth,
    Y_geno_beta = get_val(cY, "geno", 1), Y_geno_p = get_val(cY, "geno", 4),
    I_meth_beta = get_val(cI, "meth", 1), I_meth_p = get_val(cI, "meth", 4),
    I_geno_beta = get_val(cI, "geno", 1), I_geno_p = get_val(cI, "geno", 4),
    I_inter_beta = get_val(cI, "geno:meth", 1), I_inter_p = get_val(cI, "geno:meth", 4),
    ACME_est = med_res$d0, ACME_p = med_res$d0.p,
    ADE_est = med_res$z0, ADE_p = med_res$z0.p,
    Total_est = med_res$tau.coef, Total_p = med_res$tau.p
  ))
}

# ==========================================
# 4. Parallel Setup & Execution
# ==========================================
plan(multisession, workers = 4)
options(future.globals.maxSize = 10 * 1024^3)

idx_chunks <- split(seq_len(nrow(gene_mat)), ceiling(seq_len(nrow(gene_mat)) / 50))

res <- rbindlist(future_lapply(idx_chunks, function(idx_vec) {
  # Inside each worker, run rows and combine to reduce communication overhead
  rbindlist(lapply(idx_vec, function(n) {
    tryCatch({
      as.data.table(run_one_dt_fast(
        n, gene_mat[n,], meth_mat[n,], var_mat[n,], 
        meta_df, paste(meth[n, 1:5], collapse=":"), gene$ID[n]
      ))
    }, error = function(e) NULL)
  }), fill = TRUE)
}))

# ==========================================
# 5. Output
# ==========================================
if(nrow(res) > 0) {
  fwrite(res, "mediation_results_optimized_29.txt", sep = "\t", na = "NA")
  message("Success: Results written to mediation_results_optimized_29.txt")
} else {
  message("Warning: Output is empty. Check sample synchronization or p-value filters.")
}
