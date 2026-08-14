# CpG-based FDR
library(data.table)
library(dplyr)

re<-fread("merged_mediation_results.txt")

#### Preprocess ####
### double-check
setDT(re)
re[, col4 := as.numeric(tstrsplit(id, ":", fixed=TRUE)[[4]])]
# re <- re[col4 > 87 & col4 < 353]
# re_filtered <- re

# replace NA with 1 in the mediation test
# re$ACME_p[!is.na(re$M_geno_p) & re$M_geno_p > 0.1] <- 1

#### Two-stage hierarchical FDR ####

# -------------------------------------------------------
# Assume your main table is 're', a data.table
# Columns: p-value columns, transcript ID 'id2', etc.
# Primary p-value columns: primary_p_cols
# Secondary p-value columns: secondary_p_cols
# -------------------------------------------------------

library(data.table)
library(ACAT)
library(qvalue)

# -------------------------------
# Example: 're' is your main data.table with columns:
# - id2 (transcript)
# - p-value columns: primary_p_cols and secondary_p_cols
# -------------------------------
primary_p_cols   <- c("E_meth_p", "M_geno_p", "Y_meth_p", "Y_geno_p")
secondary_p_cols <- c("I_meth_p","I_geno_p","I_inter_p","ACME_p","ADE_p","Total_p")
secondary_fdr_cols <- paste0(secondary_p_cols, "_FDR")

# -------------------------------
# Safe p-value function
# -------------------------------
safe_p <- function(p) {
  p <- as.numeric(p)
  p <- p[is.finite(p)]
  p[p <= 0] <- .Machine$double.eps
  p[p >= 1] <- 1 - .Machine$double.eps
  return(p)
}

# -------------------------------
# Stage 1: Primary hierarchical FDR
# -------------------------------
for (pcol in primary_p_cols) {
  fdr_col <- paste0(pcol, "_FDR")
  message("Processing primary column: ", pcol)
  
  # Split p-values by transcript
  p_list <- split(re[[pcol]], re$id2)
  
  # Transcript-level ACAT
  acat_vector <- sapply(p_list, function(x) {
    x <- safe_p(x)
    if (length(x) == 0) return(NA_real_)
    tryCatch(ACAT::ACAT(x), error=function(e) NA_real_)
  })
  
  gene_level <- data.table(id2 = names(acat_vector), acat_p = acat_vector)
  
  # Gene-level q-values (Storey)
  gene_level[, gene_q := NA_real_]
  valid <- !is.na(gene_level$acat_p)
  if (any(valid)) {
    q_out <- tryCatch(
      qvalue(gene_level$acat_p[valid], pi0.method="bootstrap")$qvalues,
      error=function(e) qvalue(gene_level$acat_p[valid], pi0.method="smoother")$qvalues
    )
    gene_level$gene_q[valid] <- q_out
  }
  
  # Map gene-level q-values to original table
  re[, paste0(pcol, "_gene_q") := gene_level$gene_q[match(id2, gene_level$id2)]]
  
  # CpG-level FDR for significant transcripts
  sig_idx <- which(!is.na(re[[paste0(pcol, "_gene_q")]]) & re[[paste0(pcol, "_gene_q")]] < 0.05)
  re[, (fdr_col) := NA_real_]
  if (length(sig_idx) > 0) {
    valid_sig <- sig_idx[!is.na(re[[pcol]][sig_idx])]
    if (length(valid_sig) > 0) {
      re[valid_sig, (fdr_col) := p.adjust(safe_p(re[[pcol]][valid_sig]), method="BH")]
    }
  }
}

# -------------------------------
# Stage 2: Secondary FDR (only for significant transcripts)
# -------------------------------
# --- IMPROVED STAGE 2 FILTERING (Pair-Level Specificity) ---

# Now, only perform Stage 2 on these specific high-confidence pairs
# -------------------------------
# Stage 2: Secondary FDR with Specific Logic
# -------------------------------

# 1. Identify pairs passing the Mediation "Gatekeeper"
# (SNP affects Methylation AND SNP affects Expression via Methylation)
sig_med_pairs <- re[M_geno_p_FDR < 0.05 & Y_meth_p_FDR < 0.05]

# 2. Identify pairs passing the General "Gatekeeper" 
# (For I_meth, I_geno, I_inter - requiring SNP-Transcript and SNP-Meth links)
sig_gen_pairs <- re[Y_geno_p_FDR < 0.05 & Y_meth_p_FDR < 0.05]

for (fdr_col in secondary_fdr_cols) {
  re[, (fdr_col) := NA_real_]
  pcol_base <- gsub("_FDR$", "", fdr_col)
  
  # Determine which subset of pairs to use based on the variable type
  if (pcol_base %in% c("ACME_p", "ADE_p", "Total_p")) {
    target_dt <- sig_med_pairs
  } else {
    target_dt <- sig_gen_pairs
  }
  
  if (nrow(target_dt) > 0) {
    pvals <- target_dt[[pcol_base]]
    adj <- rep(NA_real_, length(pvals))
    non_na <- !is.na(pvals)
    
    if (any(non_na)) {
      # Global BH correction across the filtered subset
      adj[non_na] <- p.adjust(safe_p(pvals[non_na]), method="BH")
    }
    
    # Map results back to main table
    res_sub <- target_dt[, .(id, id2)]
    res_sub[, temp_adj := adj]
    re[res_sub, (fdr_col) := i.temp_adj, on = .(id, id2)]
  }
}

# Cleanup temporary gene_q columns
re[, grep("_gene_q$", names(re), value=TRUE) := NULL]


re$ADE_p_FDR<-NA
re$Total_p_FDR<-NA

#### Meth-exp and mixed significance: same SV segment, same transcipt, different coordinates ---> different effects  ####
## independent role
ids_E <- with(re, unique(id[!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05]))
ids_Y <- with(re, unique(id[!is.na(Y_meth_p_FDR) & Y_meth_p_FDR < 0.05]))
sum(!is.na(re$E_meth_p_FDR) & re$E_meth_p_FDR < 0.05)

length(ids_E) / length(unique(re$id)) # 0.009824791

length(with(re, unique(id2[!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05]))) # 0.009824791

length(with(re, unique(id2[!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05])))/ length(unique(re$id2))


dong<-with(re, unique(E_meth_beta[!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05]))
length(dong[dong>0])

length(dong[dong<0])


sig_ids <- with(re, id[!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05])
sum(re$E_meth_beta[re$id %in% sig_ids & !is.na(re$E_meth_p_FDR) & re$E_meth_p_FDR < 0.05] > 0, na.rm = TRUE)

sum(re$E_meth_beta[re$id %in% sig_ids & !is.na(re$E_meth_p_FDR) & re$E_meth_p_FDR < 0.05] < 0, na.rm = TRUE)


# multiple sites in the same SVs
# see in sv.R

#### Dependence and independence of meth-exp in genetics ####
# adjust genentic effects
length(ids_Y) / length(unique(re$id)) #  0.009528001

length(with(re, unique(id2[!is.na(Y_meth_p_FDR) & Y_meth_p_FDR < 0.05]))) 

length(with(re, unique(id2[!is.na(Y_meth_p_FDR) & Y_meth_p_FDR < 0.05])))/ length(unique(re$id2))


sum(!is.na(re$Y_meth_p_FDR) & re$Y_meth_p_FDR < 0.05)


subset_re <- re[re$id %in% ids_Y & !is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05, ]
fwrite(subset_re, file = "re_meth_interdependence.txt", sep = "\t", quote = FALSE, row.names = FALSE)

re$id3<-paste(re$id,re$id2,sep = '_')
length(intersect(with(re, id3[!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05]),with(re, id3[!is.na(Y_meth_p_FDR) & Y_meth_p_FDR < 0.05])))


length(intersect(ids_E,ids_Y))

length(intersect(ids_E,ids_Y))/length(ids_E) # 0.8364583 ==> independence/weak effect of epi-var, may be environmental factor or nearby variants 

length(intersect(with(re, unique(id2[!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05])),with(re, unique(id2[!is.na(Y_meth_p_FDR) & Y_meth_p_FDR < 0.05]))))


length(setdiff(ids_E, ids_Y)) # affected by epi-var

length(setdiff(with(re, id3[!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05]),with(re, id3[!is.na(Y_meth_p_FDR) & Y_meth_p_FDR < 0.05])))


length(setdiff(ids_Y, ids_E)) # found only after adjustment for epi-var, one possible reason is GXE interaction buffering methyaltion effect

length(setdiff(with(re, id3[!is.na(Y_meth_p_FDR) & Y_meth_p_FDR < 0.05]),with(re, id3[!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05])))


length(setdiff(with(re, unique(id2[!is.na(Y_meth_p_FDR) & Y_meth_p_FDR < 0.05])),with(re, unique(id2[!is.na(E_meth_p_FDR) & E_meth_p_FDR < 0.05]))))

#### Interaction in methylation part ####
### meth-exp identification affected by Interaction
length(unique(re$id[!is.na(re$I_meth_p_FDR) & re$I_meth_p_FDR < 0.05])) # ===> largely reduced


length(intersect(unique(re$id[!is.na(re$I_meth_p_FDR) & re$I_meth_p_FDR < 0.05]),ids_Y)) 

(length(ids_Y)-length(intersect(unique(re$id[!is.na(re$I_meth_p_FDR) & re$I_meth_p_FDR < 0.05]),ids_Y)))/length(ids_Y)


#  global methylation effect before and after adjust for interactions (test direction)
wilcox.test(
  abs(re$Y_meth_beta[
    !is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05
  ]),
  abs(re$I_meth_beta[
    !is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05
  ]), paired = T
)

# enhancing and buffering effect of pairs
length(unique(
  re$id[
    !is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05 &
      !is.na(re$I_meth_p_FDR)  & re$I_meth_p_FDR  < 0.05 &
      (re$I_meth_beta * re$I_inter_beta > 0) & !is.na(re$I_meth_beta) & !is.na(re$I_inter_beta)
  ]
))
# 0
length(unique(
  re$id[
    !is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05 &
      !is.na(re$I_meth_p_FDR)  & re$I_meth_p_FDR  < 0.05 &
      (re$I_meth_beta * re$I_inter_beta < 0) & !is.na(re$I_meth_beta) & !is.na(re$I_inter_beta)
  ]
))


#### Interaction in genetic part ####

### genetic-exp identification affected by Interaction
ids_Ygeno <- with(re, unique(id[!is.na(Y_geno_p_FDR) & Y_geno_p_FDR < 0.05]))
length(ids_Ygeno)


length(unique(re$id[!is.na(re$I_geno_p_FDR) & re$I_geno_p_FDR < 0.05])) # ===> largely reduced


subset_re <- re[
  re$id %in% ids_Ygeno &
    !is.na(re$Y_geno_p_FDR) &
    re$Y_geno_p_FDR < 0.05 &
    (is.na(re$I_geno_p_FDR) | re$I_geno_p_FDR > 0.05),
]
fwrite(subset_re, file = "re_genetic_interdependence.txt", sep = "\t", quote = FALSE, row.names = FALSE)

length(intersect(unique(re$id[!is.na(re$I_geno_p_FDR) & re$I_geno_p_FDR < 0.05]),ids_Ygeno)) 


(length(ids_Ygeno)-length(intersect(unique(re$id[!is.na(re$I_geno_p_FDR) & re$I_geno_p_FDR < 0.05]),ids_Ygeno)))/length(ids_Ygeno)


## deepseek in interaction
length(unique(re$id[!is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05])) / length(unique(re$id)) # 0.0002558539 => interaction number

length(unique(re$id2[!is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05])) 

write.table(re[!is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05,],file = "re_interaction.csv", sep = "\t", quote = FALSE, row.names = FALSE)


b<-c(length(unique(re$id[!is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05 & !is.na(re$I_geno_p_FDR) & re$I_geno_p_FDR < 0.05])),
     length(unique(re$id[!is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05])), length(with(re, unique(id[!is.na(Y_geno_p_FDR) & Y_geno_p_FDR < 0.05]))), length(unique(re$id)))

dim(b)<-c(2,2)
chisq.test(b)

fisher.test(b)


#  global genetic effect before and after adjust for interactions (test direction)
wilcox.test(
  re$Y_geno_beta[
    !is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05
  ],
  re$I_geno_beta[
    !is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05
  ], paired = T
)
# p-value = 0.01341
# 0.1962508 0.2888000 ===> buffering effect

# enhancing and buffering effect of pairs
length(unique(
  re$id[
    !is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05 &
      !is.na(re$I_geno_p_FDR)  & re$I_geno_p_FDR  < 0.05 &
      (re$I_geno_beta * re$I_inter_beta > 0) & !is.na(re$I_geno_beta) & !is.na(re$I_inter_beta)
  ]
))
# 4
fwrite(re[
  !is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05 &
    !is.na(re$I_geno_p_FDR)  & re$I_geno_p_FDR  < 0.05 &
    (re$I_geno_beta * re$I_inter_beta > 0) & !is.na(re$I_geno_beta) & !is.na(re$I_inter_beta),
], file = "re_hierarchical_fdr_ACAT_4_enhancing.txt", sep = "\t", quote = FALSE, row.names = FALSE)


length(unique(
  re$id[
    !is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05 &
      !is.na(re$I_geno_p_FDR)  & re$I_geno_p_FDR  < 0.05 &
      (re$I_geno_beta * re$I_inter_beta < 0) & !is.na(re$I_geno_beta) & !is.na(re$I_inter_beta)
  ]
))
# 18 ===> more play a role in buffering
fwrite(re[
  !is.na(re$I_inter_p_FDR) & re$I_inter_p_FDR < 0.05 &
    !is.na(re$I_geno_p_FDR)  & re$I_geno_p_FDR  < 0.05 &
    (re$I_geno_beta * re$I_inter_beta < 0) & !is.na(re$I_geno_beta) & !is.na(re$I_inter_beta),
], file = "re_hierarchical_fdr_ACAT_4_buffering.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#### mediation role ####
length(unique(re$id[!is.na(re$ACME_p_FDR) & re$ACME_p_FDR < 0.05])) / length(unique(re$id)) # 9.210742e-05

length(unique(re$id2[!is.na(re$ACME_p_FDR) & re$ACME_p_FDR < 0.05]))

write.table(re[!is.na(re$ACME_p_FDR) & re$ACME_p_FDR < 0.05,],file = "re_mediation.txt", sep = "\t", quote = FALSE, row.names = FALSE)

fwrite(re, file = "re_hierarchical_fdr_ACAT_4.txt", sep = "\t", quote = FALSE, row.names = FALSE)
