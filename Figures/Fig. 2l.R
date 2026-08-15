# ===============================
# Load packages
# ===============================
library(dplyr)
library(tidyr)

# ===============================
# Step 1: Read methylation data
# ===============================
df <- read.table(
  "/scratch/zdong/Projects/PanEpiG/V1-9/Geno-meth-exp/Exp/Exp_adjust2/exp_195_adjusted.bed",
  header = T,
  sep = "\t",
  stringsAsFactors = FALSE,
  comment.char = ""
)

# Missing rate cutoff (10%)
# df <- df[rowMeans(is.na(df[, 5:ncol(df)])) <= 0.10, ]


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

a<-long_df[long_df$gene %in% 'flairiso38094-1_ENSG00000123609.11',] # flairiso151080-1_ENSG00000114650.21
t.test(a$y[a$sample %in% con$sample[con$group=="AFR"]],a$y[a$sample %in% con$sample[con$group=="EAS"]])
t.test(a$y[a$sample %in% con$sample[con$group=="AFR"]],a$y[a$sample %in% con$sample[con$group=="EUR"]])
t.test(a$y[a$sample %in% con$sample[con$group=="AFR"]],a$y[a$sample %in% con$sample[con$group=="AMR"]])
t.test(a$y[a$sample %in% con$sample[con$group=="AFR"]],a$y[a$sample %in% con$sample[con$group=="SAS"]])

t.test(a$y[a$sample %in% con$sample[con$group=="AFR"]],a$y[a$sample %in% con$sample[con$group!="AFR"]])   

library(ggplot2)
library(dplyr)

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

library(ggplot2)
library(dplyr)
library(ggdist)
library(scales)
library(grid)

plot_df <- plot_df %>%
  filter(
    !is.na(group2),
    is.finite(y)
  ) %>%
  mutate(
    group2 = factor(
      group2,
      levels = unique(group2)
    )
  )

if (nlevels(plot_df$group2) != 2) {
  stop("This code expects exactly two groups in group2.")
}

# ------------------------------------------------------------
# P-value bracket position
# ------------------------------------------------------------

y_max <- max(plot_df$y, na.rm = TRUE)
y_min <- min(plot_df$y, na.rm = TRUE)
y_range <- y_max - y_min

if (!is.finite(y_range) || y_range == 0) {
  y_range <- 1
}

bracket_y <- y_max + 0.10 * y_range
arm_height <- 0.025 * y_range
text_y <- bracket_y + 0.045 * y_range

group_levels <- levels(plot_df$group2)

group_colours <- setNames(
  c("#4C78A8", "#E07B54"),
  group_levels
)

p <- ggplot(
  plot_df,
  aes(
    x = group2,
    y = y
  )
) +
  
  # Cloud
  stat_halfeye(
    aes(fill = group2),
    adjust = 0.75,
    width = 0.65,
    justification = -0.25,
    slab_colour = NA,
    slab_alpha = 0.65,
    point_colour = NA,
    interval_colour = NA,
    .width = 0,
    show.legend = FALSE
  ) +
  
  # Rain
  geom_jitter(
    aes(colour = group2),
    width = 0.075,
    height = 0,
    size = 1.05,
    alpha = 0.50,
    show.legend = FALSE
  ) +
  
  # Boxplot
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    linewidth = 0.42,
    fill = "white",
    colour = "black"
  ) +
  
  # Median
  stat_summary(
    fun = median,
    geom = "point",
    shape = 95,
    size = 5.5,
    colour = "black"
  ) +
  
  # P-value bracket
  geom_segment(
    x = 1,
    xend = 2,
    y = bracket_y,
    yend = bracket_y,
    linewidth = 0.40,
    inherit.aes = FALSE
  ) +
  
  geom_segment(
    x = 1,
    xend = 1,
    y = bracket_y,
    yend = bracket_y - arm_height,
    linewidth = 0.40,
    inherit.aes = FALSE
  ) +
  
  geom_segment(
    x = 2,
    xend = 2,
    y = bracket_y,
    yend = bracket_y - arm_height,
    linewidth = 0.40,
    inherit.aes = FALSE
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = text_y,
    label = p_label,
    size = 3.1,
    colour = "black"
  ) +
  
  scale_fill_manual(
    values = group_colours
  ) +
  
  scale_colour_manual(
    values = group_colours
  ) +
  
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(
      mult = c(0.04, 0.22)
    )
  ) +
  
  labs(
    x = NULL,
    y = "NMI5 isoform expression\n(inverse-normalized TPM)"
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.title.y = element_text(
      size = 12,
      colour = "black",
      margin = margin(r = 5)
    ),
    
    axis.line = element_line(
      linewidth = 0.7,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.7,
      colour = "black"
    ),
    
    axis.ticks.length = unit(
      0.10,
      "cm"
    ),
    
    legend.position = "none",
    
    panel.grid = element_blank(),
    
    plot.margin = margin(
      t = 14,
      r = 4,
      b = 4,
      l = 4
    )
  )


ggsave(
  filename = "NMI_flairiso38094-1_raincloud.pdf",
  plot = p,
  width = 2.1,
  height = 2.6
)
