# ===============================
# Load packages
# ===============================
library(dplyr)
library(tidyr)

# ===============================
# Step 1: Read methylation data
# ===============================
df <- read.table(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/NN/Saifur/part_00.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Missing rate cutoff (10%)
df <- df[rowMeans(is.na(df[, 2:ncol(df)])) <= 0.10, ]

# Convert beta → M values
# beta <- df[, 2:ncol(df)] / 100
# eps  <- 1e-6
# df[, 2:ncol(df)] <- log2((beta + eps) / (1 - beta + eps))

# Gene ID
gene_id <- df$V1

# Expression matrix
expr_mat <- df[, 2:ncol(df)]
# CpG / gene ID
he<-read.table("/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/NN/Saifur/header.txt",header=T)
names(expr_mat)<-names(he)[2:ncol(he)]

# ===============================
# Step 2: Metadata
# ===============================
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

# Match samples
long_df <- as.data.frame(expr_mat[,con$sample])
long_df$gene <- gene_id

# ===============================
# Step 3: Long format (CpG × sample)
# ===============================
long_df <- as.data.frame(expr_mat)
long_df$gene <- gene_id

long_df <- pivot_longer(
  long_df,
  cols = -gene,
  names_to = "sample",
  values_to = "y"
)

# ===============================
# Step 4: Compute variance components
# ===============================
compute_Qst_var <- function(df) {
  # σ²_GB = variance of group means
  # σ²_GW = mean within-group variance
  if (length(unique(df$group)) < 2) return(NA_real_)
  
  sigma2_GW <- df %>%
    group_by(group) %>%
    summarize(var_w = var(y, na.rm = TRUE), .groups = "drop") %>%
    summarize(mean(var_w, na.rm = TRUE)) %>%
    pull()
  
  sigma2_GB <- df %>%
    group_by(group) %>%
    summarize(mean_y = mean(y, na.rm = TRUE), .groups = "drop") %>%
    summarize(var(mean_y, na.rm = TRUE)) %>%
    pull()
  
  # Avoid division by zero
  if (is.na(sigma2_GB) || is.na(sigma2_GW) || sigma2_GW == 0) return(NA_real_)
  
  sigma2_GB / (sigma2_GB + 2 * sigma2_GW)
}

# ===============================
# Step 5: Continental Qst (one-vs-rest, separate files)
# ===============================

groups <- sort(unique(con$group))

for (g in groups) {
  
  message("Processing group: ", g)
  
  long_tmp <- long_df %>%
    left_join(con, by = "sample") %>%
    filter(!is.na(y)) %>%
    mutate(
      group_bin = ifelse(group == g, "FOCAL", "REST")
    )
  
  Qst_tmp <- long_tmp %>%
    group_by(gene) %>%
    summarize(
      Qst = compute_Qst_var(
        data.frame(y = y, group = group_bin)
      ),
      .groups = "drop"
    )
  
  out_file <- paste0("Qst_continent_", g, "_vs_rest_00.txt")
  
  write.table(
    Qst_tmp,
    out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
