# ===============================
# Load packages
# ===============================
library(dplyr)
library(tidyr)
library(ggplot2)

# ===============================
# Step 1: Read methylation data
# ===============================
df <- read.table(
  "exp_195_adjusted.bed",
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE,
  comment.char = ""
)

# CpG / gene ID
gene_id <- df$sample

# Expression matrix
expr_mat <- as.matrix(df[, 1:(ncol(df)-1)])

# ===============================
# Step 2: Metadata
# ===============================
# Continental groups
con <- read.table(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/Refmeth/continent.log",
  stringsAsFactors = FALSE
)[, 1:2]
colnames(con) <- c("sample", "group")
rownames(con)<-con$sample

# Match samples
common <- intersect(colnames(expr_mat), rownames(con))
expr_mat <- expr_mat[, common, drop = FALSE]
con <- con[common, , drop = FALSE]
stopifnot(ncol(expr_mat) == nrow(con))

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

a<-long_df[long_df$gene %in% 'ENST00000695812.1_ENSG00000102575.14',]
t.test(a$y[a$sample %in% con$sample[con$group=="AFR"]],a$y[a$sample %in% con$sample[con$group=="EAS"]])
t.test(a$y[a$sample %in% con$sample[con$group=="AFR"]],a$y[a$sample %in% con$sample[con$group=="EUR"]])
t.test(a$y[a$sample %in% con$sample[con$group=="AFR"]],a$y[a$sample %in% con$sample[con$group=="AMR"]])
t.test(a$y[a$sample %in% con$sample[con$group=="AFR"]],a$y[a$sample %in% con$sample[con$group=="SAS"]])

t.test(a$y[a$sample %in% con$sample[con$group=="AFR"]],a$y[a$sample %in% con$sample[con$group!="AFR"]])   

## Prepare data
plot_df <- a %>%
  mutate(group2 = ifelse(sample %in% con$sample[con$group == "AFR"],
                         "AFR", "Non-AFR"))

## t-test
tt <- t.test(y ~ group2, data = plot_df)

pval <- tt$p.value
ci_low <- tt$conf.int[1]
ci_high <- tt$conf.int[2]

p_label <- paste0(
  "P = ", signif(pval, 3),
  "\n95% CI: ", round(ci_low, 3), " to ", round(ci_high, 3)
)

# Make sure group order is fixed
plot_df <- plot_df %>%
  mutate(
    group2 = factor(
      group2,
      levels = unique(group2)
    )
  )

# Safer x-position for the P-value bracket
x_max <- max(plot_df$y, na.rm = TRUE)
x_min <- min(plot_df$y, na.rm = TRUE)
x_range <- x_max - x_min

if (x_range == 0) {
  x_range <- abs(x_max)
}

if (x_range == 0) {
  x_range <- 1
}

x_pos <- x_max + 0.08 * x_range

p <- ggplot(
  plot_df,
  aes(
    x = y,
    y = group2
  )
) +
  
  geom_boxplot(
    width = 0.42,
    outlier.shape = NA,
    linewidth = 0.45,
    fill = "white",
    colour = "black"
  ) +
  
  geom_jitter(
    width = 0,
    height = 0.10,
    size = 1.15,
    alpha = 0.55,
    colour = "grey35"
  ) +
  
  stat_summary(
    fun = median,
    geom = "point",
    shape = 124,
    size = 7,
    colour = "black"
  ) +
  
  # Vertical P-value bracket
  geom_segment(
    x = x_pos,
    xend = x_pos,
    y = 1,
    yend = 2,
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  
  geom_segment(
    x = x_pos,
    xend = x_pos - 0.025 * x_range,
    y = 1,
    yend = 1,
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  
  geom_segment(
    x = x_pos,
    xend = x_pos - 0.025 * x_range,
    y = 2,
    yend = 2,
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  
  annotate(
    "text",
    x = x_pos + 0.04 * x_range,
    y = 1.5,
    label = p_label,
    size = 3.6,
    angle = 90
  ) +
  
  labs(
    x = "ACP5 isoform expression\n(inverse-normalized TPM)",
    y = NULL
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(0.03, 0.20)
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.35,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.35,
      colour = "black"
    ),
    
    legend.position = "none",
    
    plot.margin = margin(
      5,
      18,
      5,
      5
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  )


ggsave(
  "ACP5.pdf",
  plot = p,
  width = 3.6,
  height = 1.9,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)
