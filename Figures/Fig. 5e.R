library(tidyverse)
library(scales)
library(ggrepel)
library(patchwork)
library(grid)

theme_set(
  theme_classic(
    base_size = 9,
    base_family = "Arial"
  )
)

# ============================================================
# 1. Association data
# ============================================================

df <- tribble(
  ~position, ~FDR,
  ">27868100:77:79:82:SV_ins", 0.341417822,
  ">27868100:108:110:82:SV_ins", 0.330728739,
  ">27868100:114:116:82:SV_ins", 0.902123618,
  ">27868100:133:135:82:SV_ins", 0.734890506,
  ">27868100:145:147:82:SV_ins", 0.513217221,
  ">27868100:158:160:82:SV_ins", 0.866262694,
  ">27868100:168:170:82:SV_ins", 0.986412393,
  ">27868100:173:175:82:SV_ins", 0.329056187,
  ">27868100:180:182:82:SV_ins", 0.120468373,
  ">27868100:182:184:82:SV_ins", 0.121579309,
  ">27868100:189:191:82:SV_ins", 0.339709917,
  ">27868100:197:199:82:SV_ins", 0.258690548,
  ">27868100:209:211:82:SV_ins", 0.061750087,
  ">27868100:214:216:82:SV_ins", 0.094353857,
  ">27868100:226:228:82:SV_ins", 0.084544883,
  ">27868100:251:253:82:SV_ins", 0.457594787,
  ">27868100:254:256:82:SV_ins", 0.935170319,
  ">27868100:263:265:82:SV_ins", 0.329584272,
  ">27868100:284:286:82:SV_ins", 0.008901145,
  ">27868100:304:306:82:SV_ins", 0.178877883,
  ">27868100:316:318:82:SV_ins", 0.187122293,
  ">27868100:352:354:82:SV_ins", 0.548228001,
  ">27868100:356:358:82:SV_ins", 0.190372532,
  ">27868100:364:366:82:SV_ins", 0.160185367,
  ">27868100:388:390:82:SV_ins", 0.013395645,
  ">27868100:412:414:82:SV_ins", 0.035312718,
  ">27868100:427:429:82:SV_ins", 0.006646012
)

plot_df <- df %>%
  separate(
    position,
    into = c(
      "anchor",
      "start",
      "end",
      "sv_id",
      "variant_type"
    ),
    sep = ":",
    remove = FALSE
  ) %>%
  mutate(
    anchor = as.numeric(str_remove(anchor, "^>")),
    start = as.numeric(start),
    end = as.numeric(end),
    midpoint = (start + end) / 2,
    logFDR = -log10(FDR),
    
    status = factor(
      if_else(FDR < 0.05, "FDR < 0.05", "Not significant"),
      levels = c("Not significant", "FDR < 0.05")
    ),
    
    CpG_label = as.character(round(midpoint))
  )

fdr_threshold <- -log10(0.05)

point_colours <- c(
  "Not significant" = "#BDBDBD",
  "FDR < 0.05"      = "#B2182B"
)

# ============================================================
# 2. Regional Manhattan-style association panel
# ============================================================

p_association <- ggplot(
  plot_df,
  aes(
    x = midpoint,
    y = logFDR
  )
) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 0,
    ymax = fdr_threshold,
    fill = "#F7F7F7",
    colour = NA
  ) +
  
  geom_hline(
    yintercept = fdr_threshold,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "#666666"
  ) +
  
  geom_point(
    aes(fill = status),
    shape = 21,
    size = 2.8,
    colour = "black",
    stroke = 0.25
  ) +
  
  geom_text_repel(
    data = plot_df %>%
      filter(FDR < 0.05),
    aes(label = CpG_label),
    size = 2.5,
    colour = "black",
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.18,
    segment.colour = "#777777",
    segment.linewidth = 0.25,
    direction = "both",
    max.overlaps = Inf,
    seed = 12
  ) +
  
  annotate(
    "text",
    x = min(plot_df$midpoint),
    y = fdr_threshold,
    label = "FDR = 0.05",
    hjust = 0,
    vjust = -0.55,
    size = 2.5,
    colour = "#555555"
  ) +
  
  scale_fill_manual(
    values = point_colours,
    breaks = c("FDR < 0.05", "Not significant"),
    name = NULL
  ) +
  
  scale_x_continuous(
    breaks = pretty_breaks(n = 7),
    expand = expansion(mult = c(0.025, 0.055))
  ) +
  
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.18))
  ) +
  
  labs(
    x = "Position within the inserted sequence (bp)",
    y = expression(-log[10]("FDR"))
  ) +
  
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 2.5)
    )
  ) +
  
  theme_classic(
    base_size = 9,
    base_family = "Arial"
  ) +
  
  theme(
    axis.text = element_text(
      size = 8,
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 9,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      margin = margin(t = 5)
    ),
    
    axis.title.y = element_text(
      margin = margin(r = 5)
    ),
    
    axis.line = element_line(
      linewidth = 0.35,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.3,
      colour = "black"
    ),
    
    legend.position = "none",
    legend.justification = "left",
    legend.box.just = "left",
    
    legend.text = element_text(
      size = 8,
      colour = "black"
    ),
    
    legend.key.width = unit(0.35, "cm"),
    legend.key.height = unit(0.30, "cm"),
    
    panel.grid = element_blank(),
    
    plot.margin = margin(
      t = 3,
      r = 6,
      b = 2,
      l = 3
    )
  )

# ============================================================
# 3. GRCh38 locus information
# ============================================================

# ============================================================
# 3. GRCh38 locus information using the real transcript model
# ============================================================

chromosome <- "chr13"
insertion_position <- 22913826

library(data.table)
library(dplyr)
library(stringr)

gene_exons <- fread(
  "LINC00621_exons_GRCh38.bed",
  header = FALSE,
  sep = "\t",
  col.names = c(
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand"
  )
) %>%
  as_tibble() %>%
  filter(
    chrom == "chr13",
    str_detect(name, "^ENST00000668623\\.1_exon_")
  ) %>%
  mutate(
    start = as.numeric(start) + 1,
    end = as.numeric(end)
  )

# Check that exon records were loaded
print(gene_exons)

if (nrow(gene_exons) == 0) {
  stop(
    "No ENST00000668623.1 exons were found. ",
    "Check the BED filename and transcript identifier."
  )
}

# Real transcript boundaries and strand
gene_start  <- min(gene_exons$start, na.rm = TRUE)
gene_end    <- max(gene_exons$end, na.rm = TRUE)
gene_strand <- unique(gene_exons$strand)

if (length(gene_strand) != 1) {
  stop("More than one strand value was found in the exon BED file.")
}

# Plot the full transcript plus flanking sequence
locus_padding <- 3000

window_start <- min(
  gene_start,
  insertion_position
) - locus_padding

window_end <- max(
  gene_end,
  insertion_position
) + locus_padding

# Transcript-direction arrows
arrow_spacing <- 3500

arrow_positions <- seq(
  gene_start + 1200,
  gene_end - 1200,
  by = arrow_spacing
)

gene_arrows <- tibble(
  x = arrow_positions
)

if (gene_strand == "-") {
  gene_arrows <- gene_arrows %>%
    mutate(xend = x - 800)
} else {
  gene_arrows <- gene_arrows %>%
    mutate(xend = x + 800)
}

# Genomic-coordinate ticks
coordinate_breaks <- pretty(
  c(window_start, window_end),
  n = 5
)

# Check whether insertion is within the transcript
insertion_context <- ifelse(
  insertion_position >= gene_start &
    insertion_position <= gene_end,
  "within transcript",
  "outside transcript"
)

message(
  "Transcript: chr13:",
  comma(gene_start),
  "-",
  comma(gene_end),
  "; strand ",
  gene_strand,
  "; insertion is ",
  insertion_context
)
# ============================================================
# 4. UCSC-style locus track
# ============================================================

# ============================================================
# 4. UCSC-style locus track
# ============================================================

p_locus <- ggplot() +
  
  # Track backgrounds
  annotate(
    "rect",
    xmin = window_start,
    xmax = window_end,
    ymin = 2.55,
    ymax = 3.45,
    fill = "#FAFAFA",
    colour = NA
  ) +
  annotate(
    "rect",
    xmin = window_start,
    xmax = window_end,
    ymin = 1.55,
    ymax = 2.45,
    fill = "white",
    colour = NA
  ) +
  annotate(
    "rect",
    xmin = window_start,
    xmax = window_end,
    ymin = 0.55,
    ymax = 1.45,
    fill = "#FAFAFA",
    colour = NA
  ) +
  
  # Track separators
  annotate(
    "segment",
    x = window_start,
    xend = window_end,
    y = 2.5,
    yend = 2.5,
    linewidth = 0.25,
    colour = "#D0D0D0"
  ) +
  annotate(
    "segment",
    x = window_start,
    xend = window_end,
    y = 1.5,
    yend = 1.5,
    linewidth = 0.25,
    colour = "#D0D0D0"
  ) +
  
  # Genomic ruler
  annotate(
    "segment",
    x = window_start,
    xend = window_end,
    y = 3,
    yend = 3,
    linewidth = 0.4,
    colour = "black"
  ) +
  
  geom_segment(
    data = tibble(position = coordinate_breaks),
    aes(
      x = position,
      xend = position,
      y = 2.92,
      yend = 3.08
    ),
    linewidth = 0.3,
    colour = "black"
  ) +
  
  # Transcript intron line
  annotate(
    "segment",
    x = gene_start,
    xend = gene_end,
    y = 2,
    yend = 2,
    linewidth = 0.55,
    colour = "#2C3E8C"
  ) +
  
  # Direction arrows along the transcript
  geom_segment(
    data = gene_arrows,
    aes(
      x = x,
      xend = xend,
      y = 2,
      yend = 2
    ),
    linewidth = 0.45,
    colour = "#2C3E8C",
    arrow = arrow(
      length = unit(0.07, "cm"),
      type = "closed"
    )
  ) +
  
  # Real exons from ENST00000668623.1
  geom_rect(
    data = gene_exons,
    aes(
      xmin = start,
      xmax = end,
      ymin = 1.75,
      ymax = 2.25
    ),
    inherit.aes = FALSE,
    fill = "#2C3E8C",
    colour = "#2C3E8C",
    linewidth = 0.2
  ) +
  
  # Transcript label
  annotate(
    "text",
    x = gene_start,
    y = 2.34,
    label = "LINC00621",
    parse = TRUE,
    hjust = 0,
    size = 3.2,
    colour = "#2C3E8C"
  ) +
  
  # Insertion breakpoint
  annotate(
    "segment",
    x = insertion_position,
    xend = insertion_position,
    y = 0.78,
    yend = 1.30,
    linewidth = 1,
    colour = "#B2182B"
  ) +
  
  annotate(
    "point",
    x = insertion_position,
    y = 1.30,
    shape = 25,
    size = 3.2,
    fill = "#B2182B",
    colour = "#B2182B"
  ) +
  
  annotate(
    "text",
    x = insertion_position,
    y = 0.72,
    label = "Insertion\nchr13:22,913,826",
    size = 2.5,
    colour = "#B2182B",
    vjust = 1
  ) +
  
  # Dashed guide through tracks
  annotate(
    "segment",
    x = insertion_position,
    xend = insertion_position,
    y = 0.55,
    yend = 3.45,
    linewidth = 0.3,
    linetype = "dashed",
    colour = "#B2182B"
  ) +
  
  # Track labels
  annotate(
    "text",
    x = window_start,
    y = 2,
    label = "Gene",
    hjust = 1.08,
    size = 2.6,
    colour = "black"
  ) +
  
  annotate(
    "text",
    x = window_start,
    y = 1,
    label = "SV",
    hjust = 1.08,
    size = 2.6,
    colour = "black"
  ) +
  
  scale_x_continuous(
    limits = c(window_start, window_end),
    breaks = coordinate_breaks,
    labels = label_number(
      scale = 1e-6,
      accuracy = 0.001,
      suffix = " Mb"
    ),
    expand = c(0, 0),
    position = "top"
  ) +
  
  scale_y_continuous(
    limits = c(0.40, 3.55),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  
  labs(
    x = "chr13 position (GRCh38)",
    y = NULL
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_void(
    base_size = 9,
    base_family = "Arial"
  ) +
  
  theme(
    axis.text.x.top = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.title.x.top = element_text(
      size = 11,
      colour = "black",
      margin = margin(b = 4)
    ),
    
    axis.ticks.x.top = element_line(
      linewidth = 0.3,
      colour = "black"
    ),
    
    plot.margin = margin(
      t = 6,
      r = 6,
      b = 3,
      l = 55
    )
  )

# ============================================================
# 5. Combine association and locus panels
# ============================================================

combined_plot <- p_locus / p_association +
  plot_layout(
    heights = c(1.35, 2.4)
  )


# ============================================================
# 6. Save
# ============================================================

ggsave(
  filename =
    "LINC00621_insertion_varCpG_regional_association_UCSC.pdf",
  plot = combined_plot,
  width = 6,
  height = 4.6,
  units = "in"
)
