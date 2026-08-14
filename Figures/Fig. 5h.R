library(data.table)
library(ggplot2)

dt <- fread(
  "input.log",
  header = TRUE,
  sep = "\t"
)

dt[, tissue_label := gsub("_", " ", tissue)]

dt[, type := factor(
  type,
  levels = c("mQTL", "sQTL", "eQTL")
)]

# Order tissues by mQTL fold enrichment
ord_m <- dt[
  type == "mQTL"
][
  order(FE),
  tissue_label
]

ord_all <- c(
  ord_m,
  setdiff(unique(dt$tissue_label), ord_m)
)

dt[, tissue_label := factor(
  tissue_label,
  levels = ord_all
)]

# Range of QTL enrichments within each tissue
seg <- dt[, .(
  FE_min = min(FE, na.rm = TRUE),
  FE_max = max(FE, na.rm = TRUE)
), by = tissue_label]

# ---------------------------
# Rotated plot
# ---------------------------
p2 <- ggplot() +
  
  # Vertical range connecting QTL types within each tissue
  geom_segment(
    data = seg,
    aes(
      x = tissue_label,
      xend = tissue_label,
      y = FE_min,
      yend = FE_max
    ),
    colour = "grey70",
    linewidth = 0.9,
    lineend = "round"
  ) +
  
  # QTL-specific enrichment points
  geom_point(
    data = dt,
    aes(
      x = tissue_label,
      y = FE,
      fill = type
    ),
    shape = 21,
    size = 3.6,
    colour = "black",
    stroke = 0.28
  ) +
  
  scale_fill_manual(
    values = c(
      mQTL = "#D55E00",
      sQTL = "#0072B2",
      eQTL = "#009E73"
    )
  ) +
  
  scale_y_continuous(
    limits = c(1, 2),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = NULL,
    y = "Fold enrichment",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x = element_text(
      size = 11,
      colour = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.title.y = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    
    legend.position = "top",
    
    legend.text = element_text(
      size = 10
    ),
    
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    plot.margin = margin(
      5.5,
      12,
      18,
      5.5
    )
  )


ggsave(
  "lollipop_QTL_rotated.pdf",
  plot = p2,
  width = 10.2,
  height = 5,
  units = "in",
  useDingbats = FALSE
)
