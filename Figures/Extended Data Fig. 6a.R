library(tidyverse)
library(scales)
library(grid)

df <- read.delim("input.txt", header = TRUE, check.names = FALSE)

df2 <- df %>%
  rowwise() %>%
  mutate(
    var_overlap    = round(`var-cpg_overlapped`),
    var_nonoverlap = round(`var-cpg_all` - `var-cpg_overlapped`),
    cpg_overlap    = round(`cpg_overlapped`),
    cpg_nonoverlap = round(`cpg_all` - `cpg_overlapped`),
    
    fold_enrichment = (`var-cpg_overlapped` / `var-cpg_all`) /
      (`cpg_overlapped` / `cpg_all`),
    
    fisher_p = fisher.test(
      matrix(
        c(var_overlap, var_nonoverlap,
          cpg_overlap, cpg_nonoverlap),
        nrow = 2,
        byrow = TRUE
      )
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(fisher_p, method = "BH"),
    neglog10_FDR = -log10(pmax(FDR, 1e-300)),
    
    Method = factor(Method, levels = c("CpG gain", "CpG loss")),
    Sample = factor(Sample, levels = c("SNV", "Indel", "SV")),
    Group  = factor(Group, levels = c("Gene", "Promoter", "CGI", "Shore", "Shelf"))
  )

write.table(df2,file='enrichment.txt',sep="\t",quote = F)

p <- ggplot(df2, aes(x = Sample, y = Group)) +
  geom_point(
    aes(size = fold_enrichment, fill = neglog10_FDR),
    shape = 21,
    colour = "black",
    stroke = 0.35,
    alpha = 0.95
  ) +
  facet_wrap(~ Method, nrow = 1) +
  scale_size_area(
    max_size = 11,
    limits = c(0, 5.5),
    breaks = c(1, 2, 3, 4, 5),
    name = "Fold\nenrichment"
  ) +
  scale_fill_gradientn(
    colours = c("#F2F4F7", "#F3B6B1", "#D95F5F", "#A50F15"),
    name = expression(-log[10]~FDR),
    limits = c(0, max(df2$neglog10_FDR, na.rm = TRUE)),
    oob = squish
  ) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 0.9) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 12,
      face = "bold",
      colour = "black",
      margin = margin(b = 6)
    ),
    axis.text.x = element_text(
      size = 11,
      colour = "black",
      margin = margin(t = 3)
    ),
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.35
    ),
    panel.spacing.x = unit(0.35, "in"),
    legend.position = "right",
    legend.title = element_text(size = 12, colour = "black"),
    legend.text = element_text(size = 10, colour = "black"),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.35, "cm"),
    plot.margin = margin(5, 8, 5, 5)
  )

# p

ggsave(
  "varCpG_enrichment_FDR_bubble.pdf",
  p,
  width = 6.2,
  height = 3.4,
  units = "in"
)
