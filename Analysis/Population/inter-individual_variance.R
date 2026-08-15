# -------------------------------
# Load packages
# -------------------------------
library(dplyr)

# -------------------------------
# Step 1: Read methylation data
# -------------------------------
df <- read.table(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/NN/Saifur/part_00.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Missing rate cutoff (10%)
df <- df[rowMeans(is.na(df[, 2:ncol(df)])) <= 0.10, ]
he<-read.table("/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/NN/Saifur/header.txt",header=T)
names(df)<-names(he)

# Gene ID
gene_id <- df$gene

# -------------------------------
# Step 2: Metadata
# -------------------------------
# Continental groups
con <- read.table(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/Refmeth/continent.log",
  stringsAsFactors = FALSE
)[, 1:2]
colnames(con) <- c("sample", "group")

# Population groups
pop <- read.table(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/Refmeth/pop.log",
  stringsAsFactors = FALSE
)
pop <- pop[!pop$V2 %in% c("CH", "EUR"), 1:2]
colnames(pop) <- c("sample", "group")

# Expression matrix
expr_mat <- as.matrix(df[, 2:ncol(df)])

# Match sample names
stopifnot(ncol(expr_mat) == nrow(con))
colnames(expr_mat) <- con$sample
expr_mat<-expr_mat[,con$sample]

# -------------------------------
# Step 3: Function
# -------------------------------
compute_group_variance <- function(expr_mat, gene_id, meta, remove_samples = NULL,
                                   eps = 1e-6) {
  
  # Subset samples
  keep <- meta$sample
  if (!is.null(remove_samples)) {
    keep <- setdiff(keep, remove_samples)
  }
  
  expr_mat <- expr_mat[, keep, drop = FALSE]
  meta     <- meta[meta$sample %in% keep, ]
  
  # Precompute sample indices per group
  idx_list <- split(seq_len(nrow(meta)), meta$group)
  
  out <- vector("list", length(idx_list) * nrow(expr_mat))
  k <- 1L
  
  for (g in names(idx_list)) {
    idx <- idx_list[[g]]
    
    for (i in seq_len(nrow(expr_mat))) {
      v <- expr_mat[i, idx]
      v <- v[!is.na(v)]
      
      if (length(v) > 1) {
        var_res <- var(v)
        out[[k]] <- c(
          gene = gene_id[i],
          group = g,
          var_residuals = var_res,
          log_var_residuals = log10(var_res + eps)
        )
        k <- k + 1L
      }
    }
  }
  
  as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
}

# -------------------------------
# Step 4: Run analyses
# -------------------------------
# Continental
con_var <- compute_group_variance(
  expr_mat = expr_mat,
  gene_id  = gene_id,
  meta     = con
)

write.table(
  con_var,
  "con_variance_df_00.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Population
pop_var <- compute_group_variance(
  expr_mat = expr_mat,
  gene_id  = gene_id,
  meta     = pop,
  remove_samples = c("HG002", "HG005")
)

write.table(
  pop_var,
  "pop_variance_df_00.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
