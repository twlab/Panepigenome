library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

theme_set(theme_classic(base_size = 9))

# ============================================================
# 1. Read data
# ============================================================

self <- fread(
  "../../G-M_var-cpg/Matrixeqtls/MAF5/Permutation_analysis_itself/results_with_qvalue_significant.txt"
)

other <- fread(
  "../../G-M_var-cpg/Matrixeqtls/MAF5/Permutation_analysis_pair/results_with_qvalue_significant.txt"
)

gg <- fread(
  "../../G-M_var-cpg-g/MAF5/Permutation_analysis/results_with_qvalue_significant.txt"
)

# ============================================================
# 2. Define sets
# ============================================================

A <- unique(na.omit(self$id))
B <- unique(na.omit(other$id))
C <- unique(na.omit(gg$id))

A <- A[A != ""]
B <- B[B != ""]
C <- C[C != ""]

all_ids <- unique(c(A, B, C))

fwrite(
  data.table(ID = all_ids),
  "all_var-cpg_methbyvariants.txt",
  sep = "\t",
  col.names = FALSE
)

# ============================================================
# 3. Membership and intersection counts
# ============================================================

membership_df <- tibble(
  id = all_ids,
  A = all_ids %in% A,
  B = all_ids %in% B,
  C = all_ids %in% C
)

summary_table <- membership_df %>%
  mutate(
    category = case_when(
      A & !B & !C ~ "A only",
      !A & B & !C ~ "B only",
      !A & !B & C ~ "C only",
      A & B & !C ~ "A & B",
      A & !B & C ~ "A & C",
      !A & B & C ~ "B & C",
      A & B & C ~ "A & B & C",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(category)) %>%
  group_by(category) %>%
  summarise(
    n = dplyr::n(),
    .groups = "drop"
  )

fwrite(
  as.data.table(summary_table),
  "varCpG_methylation_association_intersections.tsv",
  sep = "\t"
)

# ============================================================
# 4. Intersection order
# ============================================================

category_order <- summary_table %>%
  arrange(desc(n)) %>%
  pull(category)

summary_table <- summary_table %>%
  mutate(
    category = factor(category, levels = category_order),
    category_number = as.numeric(category),
    count_label = case_when(
      n >= 1e6 ~ paste0(number(n / 1e6, accuracy = 0.1), "m"),
      n >= 1e3 ~ paste0(number(n / 1e3, accuracy = 0.1), "k"),
      TRUE ~ comma(n)
    )
  ) %>%
  arrange(category_number)

n_intersections <- length(category_order)

# ============================================================
# 5. Matrix data
# ============================================================

matrix_df <- expand_grid(
  category = category_order,
  set_code = c("A", "B", "C")
) %>%
  mutate(
    present = case_when(
      category == "A only" & set_code == "A" ~ TRUE,
      category == "B only" & set_code == "B" ~ TRUE,
      category == "C only" & set_code == "C" ~ TRUE,
      category == "A & B" & set_code %in% c("A", "B") ~ TRUE,
      category == "A & C" & set_code %in% c("A", "C") ~ TRUE,
      category == "B & C" & set_code %in% c("B", "C") ~ TRUE,
      category == "A & B & C" ~ TRUE,
      TRUE ~ FALSE
    ),
    category = factor(category, levels = category_order),
    category_number = as.numeric(category),
    set_code = factor(
      set_code,
      levels = c("C", "B", "A")
    ),
    set_number = as.numeric(set_code)
  )

segment_df <- matrix_df %>%
  filter(present) %>%
  group_by(category_number) %>%
  summarise(
    ymin = min(set_number),
    ymax = max(set_number),
    .groups = "drop"
  ) %>%
  filter(ymin < ymax)

# ============================================================
# 6. Set labels and total set sizes
# ============================================================

set_df <- tibble(
  set_code = c("C", "B", "A"),
  set_number = c(1, 2, 3),
  set_label = c(
    "Nearby non-CpG-altering\nvariants",
    "Nearby var-CpG\ncopy number",
    "Self var-CpG\ncopy number"
  ),
  set_size = c(
    length(C),
    length(B),
    length(A)
  )
) %>%
  mutate(
    count_label = case_when(
      set_size >= 1e6 ~ paste0(
        number(set_size / 1e6, accuracy = 0.1),
        "m"
      ),
      set_size >= 1e3 ~ paste0(
        number(set_size / 1e3, accuracy = 0.1),
        "k"
      ),
      TRUE ~ comma(set_size)
    )
  )

# ============================================================
# 7. Intersection-size bar plot
# Numeric x-axis
# ============================================================

p_bar <- ggplot(
  summary_table,
  aes(
    x = category_number,
    y = n
  )
) +
  geom_col(
    width = 0.7,
    fill = "#4D4D4D",
    colour = NA
  ) +
  geom_text(
    aes(label = count_label),
    vjust = -0.35,
    size = 3.2,
    colour = "black"
  ) +
  scale_x_continuous(
    limits = c(0.5, n_intersections + 0.5),
    breaks = seq_len(n_intersections),
    labels = NULL,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    labels = label_number(
      scale_cut = cut_short_scale()
    ),
    breaks = pretty_breaks(n = 4),
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    x = NULL,
    y = "Associated var-CpGs"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 12,
      colour = "black",
      margin = margin(r = 6)
    ),
    axis.line.y = element_line(
      linewidth = 1,
      colour = "black"
    ),
    axis.ticks.y = element_line(
      linewidth = 1,
      colour = "black"
    ),
    panel.grid = element_blank(),
    plot.margin = margin(3, 3, 0, 3)
  )

# ============================================================
# 8. Set-label panel
# ============================================================

p_labels <- ggplot(
  set_df,
  aes(
    x = 1,
    y = set_number,
    label = set_label
  )
) +
  geom_text(
    hjust = 1,
    size = 3.2,
    colour = "black",
    lineheight = 0.95
  ) +
  scale_x_continuous(
    limits = c(0, 1.02),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, 3.5),
    breaks = 1:3,
    expand = c(0, 0)
  ) +
  theme_void() +
  theme(
    plot.margin = margin(0, 2, 3, 2)
  )

# ============================================================
# 9. Dot matrix
# Numeric x-axis matching p_bar exactly
# ============================================================

p_matrix <- ggplot() +
  geom_segment(
    data = segment_df,
    aes(
      x = category_number,
      xend = category_number,
      y = ymin,
      yend = ymax
    ),
    linewidth = 1,
    colour = "#333333",
    lineend = "round"
  ) +
  geom_point(
    data = matrix_df %>% filter(!present),
    aes(
      x = category_number,
      y = set_number
    ),
    shape = 16,
    size = 3.2,
    colour = "#D9D9D9"
  ) +
  geom_point(
    data = matrix_df %>% filter(present),
    aes(
      x = category_number,
      y = set_number
    ),
    shape = 16,
    size = 3.2,
    colour = "#333333"
  ) +
  scale_x_continuous(
    limits = c(0.5, n_intersections + 0.5),
    breaks = seq_len(n_intersections),
    labels = NULL,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, 3.5),
    breaks = 1:3,
    labels = NULL,
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(0, 3, 3, 3)
  )

# ============================================================
# 10. Total set-size horizontal bars
# ============================================================

set_df <- set_df %>%
  mutate(
    set_row = factor(
      set_number,
      levels = c(1, 2, 3)
    )
  )

p_set_size <- ggplot(
  set_df,
  aes(
    x = set_size,
    y = set_row
  )
) +
  geom_col(
    width = 0.7,
    fill = "#737373",
    colour = NA,
    orientation = "y"
  ) +
  
  geom_text(
    aes(label = count_label),
    hjust = -0.15,
    size = 3.2,
    colour = "black"
  ) +
  
  scale_x_continuous(
    limits = c(
      0,
      max(set_df$set_size, na.rm = TRUE) * 1.28
    ),
    labels = label_number(
      scale_cut = cut_short_scale()
    ),
    breaks = pretty_breaks(n = 3),
    expand = c(0, 0)
  ) +
  
  scale_y_discrete(
    limits = levels(set_df$set_row),
    labels = NULL,
    expand = expansion(add = 0.5)
  ) +
  
  labs(
    x = "Total associated var-CpGs",
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      size = 12,
      colour = "black",
      margin = margin(t = 4)
    ),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    
    axis.line.x = element_line(
      linewidth = 1,
      colour = "black"
    ),
    
    axis.ticks.x = element_line(
      linewidth = 1,
      colour = "black"
    ),
    
    panel.grid = element_blank(),
    
    plot.margin = margin(
      t = 0,
      r = 8,
      b = 3,
      l = 3
    )
  )

# ============================================================
# 11. Assemble
#
# The same three column widths must be used in both rows.
# ============================================================

column_widths <- c(2.15, 4.8, 1.85)

top_row <- plot_spacer() + p_bar + plot_spacer() +
  plot_layout(widths = column_widths)

bottom_row <- p_labels + p_matrix + p_set_size +
  plot_layout(widths = column_widths)

p <- top_row / bottom_row +
  plot_layout(
    heights = c(3.0, 1.25)
  )

# ============================================================
# 12. Save
# ============================================================

ggsave(
  "varCpG_methylation_association_upset_with_set_sizes.pdf",
  plot = p,
  width = 6.6,
  height = 3.8
)
