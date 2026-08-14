library(data.table)
library(dplyr)
library(ggplot2)

read_bed8 <- function(file, label = NULL) {
  x <- fread(file, header = FALSE, fill = TRUE)
  
  while (ncol(x) < 8) {
    x[, paste0("V", ncol(x) + 1) := NA_character_]
  }
  
  x <- x[, 1:8]
  
  if (!is.null(label)) {
    x[, V8 := label]
  }
  
  x
}

shared <- read_bed8("matched_shared.CREs.bed")
backg  <- read_bed8("matched_background.CREs.bed")

shared_exon   <- read_bed8("matched_shared.exon.bed", "exon")
backg_exon    <- read_bed8("matched_background.exon.bed", "exon")

shared_intron <- read_bed8("matched_shared.intron.bed", "intron")
backg_intron  <- read_bed8("matched_background.intron.bed", "intron")

shared_se     <- read_bed8("matched_shared.SE.bed", "SE")
backg_se      <- read_bed8("matched_background.SE.bed", "SE")

shared <- rbindlist(
  list(shared, shared_exon, shared_intron, shared_se),
  use.names = FALSE,
  fill = TRUE
)

backg <- rbindlist(
  list(backg, backg_exon, backg_intron, backg_se),
  use.names = FALSE,
  fill = TRUE
)

types <- unique(backg$V8)

or <- rep(NA_real_, length(types))
p  <- rep(NA_real_, length(types))
shared_n_vec <- rep(NA_integer_, length(types))
backg_n_vec  <- rep(NA_integer_, length(types))

for (i in seq_along(types)) {
  
  type_i <- types[i]
  
  shared_n <- length(unique(shared$V4[shared$V8 == type_i]))
  backg_n  <- length(unique(backg$V4[backg$V8 == type_i]))
  
  shared_n_vec[i] <- shared_n
  backg_n_vec[i]  <- backg_n
  
  mat <- matrix(
    c(
      shared_n,
      8869 - shared_n,
      backg_n,
      269910 - backg_n
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  ft <- fisher.test(mat)
  
  p[i]  <- ft$p.value
  or[i] <- as.numeric(ft$estimate)
}

res <- data.table(
  type = types,
  shared_n = shared_n_vec,
  backg_n = backg_n_vec,
  odds_ratio = or,
  p_value = p
)

res[, fdr := p.adjust(p_value, method = "fdr")]

res <- res %>%
  mutate(
    feature = case_when(
      type == "1_Active_Promoter" ~ "Active promoter",
      type == "2_Weak_Promoter" ~ "Weak promoter",
      type == "3_Poised_Promoter" ~ "Poised promoter",
      type %in% c("4_Strong_Enhancer", "5_Strong_Enhancer") ~ "Strong enhancer",
      type %in% c("6_Weak_Enhancer", "7_Weak_Enhancer") ~ "Weak enhancer",
      type == "8_Insulator" ~ "Insulator",
      type == "9_Txn_Transition" ~ "Txn transition",
      type == "10_Txn_Elongation" ~ "Txn elongation",
      type == "11_Weak_Txn" ~ "Weak transcription",
      type == "12_Repressed" ~ "Repressed",
      type == "13_Heterochrom/lo" ~ "Heterochromatin",
      type %in% c("14_Repetitive/CNV", "15_Repetitive/CNV") ~ "Repetitive/CNV",
      type == "exon" ~ "Exon",
      type == "intron" ~ "Intron",
      type == "SE" ~ "Super-enhancer",
      TRUE ~ type
    )
  ) %>%
  as.data.table()

res_merged <- res[
  ,
  .(
    shared_n = sum(shared_n),
    backg_n = sum(backg_n)
  ),
  by = feature
]

res_merged[
  ,
  `:=`(
    shared_all = 8869,
    backg_all = 269910
  )
]

res_merged[
  ,
  c("odds_ratio", "p_value") := {
    mat <- matrix(
      c(
        shared_n,
        shared_all,
        backg_n,
        backg_all
      ),
      nrow = 2,
      byrow = TRUE
    )
    
    ft <- fisher.test(mat)
    
    list(
      as.numeric(ft$estimate),
      ft$p.value
    )
  },
  by = feature
]

res_merged[, fdr := p.adjust(p_value, method = "fdr")]
res_merged <- res_merged[order(-odds_ratio)]

fwrite(
  res_merged,
  "CRE_enrichment_odds_ratio_merged_features.tsv",
  sep = "\t"
)

res_merged[, sig := ifelse(fdr < 0.05, "FDR < 0.05", "NS")]
res_merged[, feature := factor(feature, levels = feature[order(odds_ratio)])]

p <- ggplot(res_merged, aes(x = odds_ratio, y = feature)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.45
  ) +
  geom_point(
    aes(fill = sig),
    shape = 21,
    size = 3.2,
    color = "black",
    stroke = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "FDR < 0.05" = "#D55E00",
      "NS" = "grey80"
    )
  ) +
  scale_x_continuous(
    name = "Odds ratio",
    breaks = c(0.5, 1.0, 1.5, 2.0),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  labs(
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10.5, color = "black"),
    axis.text.x = element_text(size = 10.5, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )

# p

ggsave(
  "CRE_enrichment_odds_ratio_merged_features.pdf",
  plot = p,
  width = 5.4,
  height = 4.8,
  device = "pdf"
)
