ben<-read.table("matched_meth_Benign_Likely_benign.bed.filtered")
median(ben$V4)


df<-read.table("PBR.bed.filtered")
median(df$V4)

df2<-read.table("../../../../Drug/ClinPGx/overlap_Both/DNAm/matched_meth.bed.filtered")
median(df2$V4)

df3<-rbind(df,df2)
median(df3$V4)

shapiro.test(df3$V4)
wilcox.test(df3$V4,ben$V4)


library(tidyverse)
library(scales)

theme_set(theme_classic(base_size = 12))

# ============================================================
# 1. Read data
# ============================================================

ben <- read.table(
  "matched_meth_Benign_Likely_benign.bed.filtered",
  header = FALSE,
  stringsAsFactors = FALSE
)

df <- read.table(
  "PBR.bed.filtered",
  header = FALSE,
  stringsAsFactors = FALSE
)

df2 <- read.table(
  "../../../../Drug/ClinPGx/overlap_Both/DNAm/matched_meth.bed.filtered",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Combine non-benign clinical and pharmacogenomic var-CpGs
df3 <- rbind(df, df2)

# ============================================================
# 2. Prepare methylation data
# ============================================================

meth_df <- bind_rows(
  tibble(
    type = "Non-benign clinical or pharmacogenomic",
    methylation = as.numeric(df3$V4) / 100
  ),
  tibble(
    type = "Benign clinical",
    methylation = as.numeric(ben$V4) / 100
  )
) %>%
  filter(
    !is.na(methylation),
    is.finite(methylation),
    methylation >= 0,
    methylation <= 1
  ) %>%
  mutate(
    type = factor(
      type,
      levels = c(
        "Non-benign clinical or pharmacogenomic",
        "Benign clinical"
      )
    )
  )

# ============================================================
# 3. Summary statistics
# ============================================================

summary_df <- meth_df %>%
  group_by(type) %>%
  summarise(
    n = n(),
    mean = mean(methylation),
    median = median(methylation),
    sd = sd(methylation),
    .groups = "drop"
  )

print(summary_df)

# ============================================================
# 4. Wilcoxon rank-sum test
# ============================================================

wilcox_result <- wilcox.test(
  methylation ~ type,
  data = meth_df,
  exact = FALSE
)

p_value <- wilcox_result$p.value

format_p <- function(p) {
  
  if (is.na(p)) {
    return("P = NA")
  }
  
  if (p < 2.2e-16) {
    return("P < 2.2 × 10\u207b\u00b9\u2076")
  }
  
  if (p < 0.001) {
    exponent <- floor(log10(p))
    coefficient <- p / 10^exponent
    
    superscript_digits <- c(
      "0" = "\u2070",
      "1" = "\u00b9",
      "2" = "\u00b2",
      "3" = "\u00b3",
      "4" = "\u2074",
      "5" = "\u2075",
      "6" = "\u2076",
      "7" = "\u2077",
      "8" = "\u2078",
      "9" = "\u2079",
      "-" = "\u207b"
    )
    
    exponent_text <- paste0(
      superscript_digits[
        strsplit(as.character(exponent), "")[[1]]
      ],
      collapse = ""
    )
    
    return(
      paste0(
        "P = ",
        sprintf("%.2f", coefficient),
        " × 10",
        exponent_text
      )
    )
  }
  
  paste0("P = ", sprintf("%.3f", p))
}

p_label <- format_p(p_value)

# ============================================================
# 5. Annotation text
# ============================================================

annotation_df <- summary_df %>%
  mutate(
    annotation = paste0(
      type,
      ": n = ",
      comma(n),
      ", median = ",
      sprintf("%.3f", median)
    )
  )

annotation_text <- paste(
  annotation_df$annotation,
  collapse = "\n"
)

# ============================================================
# 6. colours
# ============================================================

type_colours <- c(
  "Non-benign clinical or pharmacogenomic" = "#9E1B32",
  "Benign clinical" = "#3B6F8F"
)

# ============================================================
# 7. ECDF plot
# ============================================================

p <- ggplot(
  meth_df,
  aes(
    x = methylation,
    colour = type
  )
) +
  
  stat_ecdf(
    geom = "step",
    linewidth = 1.05,
    pad = TRUE
  ) +
  
  # Median lines
  geom_vline(
    data = summary_df,
    aes(
      xintercept = median,
      colour = type
    ),
    linetype = "dashed",
    linewidth = 0.65,
    show.legend = FALSE
  ) +
  
  annotate(
    "text",
    x = 0.025,
    y = 0.97,
    label = annotation_text,
    hjust = 0,
    vjust = 1,
    size = 3.15,
    lineheight = 1.15,
    colour = "black"
  ) +
  
  annotate(
    "text",
    x = 0.975,
    y = 0.08,
    label = p_label,
    hjust = 1,
    vjust = 0,
    size = 3.4,
    colour = "black"
  ) +
  
  scale_colour_manual(
    values = type_colours,
    breaks = c(
      "Non-benign clinical or pharmacogenomic",
      "Benign clinical"
    ),
    labels = c(
      "Non-benign clinical or\npharmacogenomic",
      "Benign clinical"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = "Methylation level",
    y = "Cumulative proportion"
  ) +
  
  guides(
    colour = guide_legend(
      override.aes = list(
        linewidth = 1.2,
        linetype = "solid"
      )
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      margin = margin(t = 6)
    ),
    
    axis.title.y = element_text(
      margin = margin(r = 6)
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
      0.13,
      "cm"
    ),
    
    panel.grid = element_blank(),
    
    legend.position = c(0.97, 0.52),
    legend.justification = c(1, 0.5),
    
    legend.text = element_text(
      size = 10.5,
      colour = "black",
      lineheight = 0.95
    ),
    
    legend.key.width = unit(
      0.65,
      "cm"
    ),
    
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      colour = NA
    ),
    
    plot.margin = margin(
      t = 6,
      r = 7,
      b = 5,
      l = 5
    )
  )


# ============================================================
# 8. Save
# ============================================================

ggsave(
  filename = "clinical_benign_methylation_ECDF.pdf",
  plot = p,
  width = 3.2,
  height = 3.2
)
