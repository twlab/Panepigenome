library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)
library(grid)

theme_set(theme_classic(base_size = 10))

make_density_plot_obj <- function(prefix, label) {
  
  # --------------------------------------------------
  # 1. Find gain and loss files
  # Expected filenames contain:
  # SNV / Indel / SV
  # gain / loss
  # --------------------------------------------------
  variant_files <- list.files(
    path = ".",
    pattern = paste0(
      "^",
      gsub("\\.", "\\\\.", prefix),
      ".*(gain|loss).*\\.log$"
    ),
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  ref_file <- paste0(prefix, ".ref.log")
  
  if (length(variant_files) == 0) {
    stop("No gain or loss files found for: ", prefix)
  }
  
  if (!file.exists(ref_file)) {
    stop("Reference file not found: ", ref_file)
  }
  
  # --------------------------------------------------
  # 2. Read gain/loss files and infer variant class
  # --------------------------------------------------
  df_variant <- rbindlist(
    lapply(variant_files, function(f) {
      
      file_name <- basename(f)
      
      variant_class <- case_when(
        grepl("SNV", file_name, ignore.case = TRUE) ~ "SNV",
        grepl("indel", file_name, ignore.case = TRUE) ~ "Indel",
        grepl("SV", file_name, ignore.case = TRUE) ~ "SV",
        TRUE ~ NA_character_
      )
      
      event_type <- case_when(
        grepl("gain", file_name, ignore.case = TRUE) ~ "CpG gain",
        grepl("loss", file_name, ignore.case = TRUE) ~ "CpG loss",
        TRUE ~ NA_character_
      )
      
      x <- fread(
        f,
        header = FALSE,
        select = 1,
        col.names = "Value",
        showProgress = FALSE
      )
      
      x[, `:=`(
        Variant_class = variant_class,
        Event = event_type,
        Source_file = file_name
      )]
      
      x
    }),
    use.names = TRUE,
    fill = TRUE
  )
  
  # Remove files whose class/event could not be inferred
  df_variant <- df_variant[
    !is.na(Variant_class) & !is.na(Event)
  ]
  
  # --------------------------------------------------
  # 3. Read GRCh38 reference values
  # --------------------------------------------------
  df_ref <- fread(
    ref_file,
    header = FALSE,
    select = 1,
    col.names = "Value",
    showProgress = FALSE
  )
  
  df_ref[, `:=`(
    Variant_class = "GRCh38",
    Event = "Reference",
    Source_file = basename(ref_file)
  )]
  
  # --------------------------------------------------
  # 4. Combine and format
  # --------------------------------------------------
  df <- rbindlist(
    list(df_ref, df_variant),
    use.names = TRUE,
    fill = TRUE
  ) %>%
    as.data.frame() %>%
    mutate(
      Value = as.numeric(Value),
      Variant_class = factor(
        Variant_class,
        levels = c("GRCh38", "SNV", "Indel", "SV")
      ),
      Event = factor(
        Event,
        levels = c("Reference", "CpG gain", "CpG loss")
      )
    ) %>%
    filter(
      is.finite(Value),
      Value >= 0,
      Value <= 100
    )
  
  # --------------------------------------------------
  # 5. Median values
  # --------------------------------------------------
  medians <- df %>%
    group_by(Variant_class, Event) %>%
    summarise(
      median_value = median(Value, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  print(medians)
  
  # --------------------------------------------------
  # 6. Colours and line types
  # --------------------------------------------------
  class_colours <- c(
    "GRCh38" = "#000000",
    "SNV"    = "#595959",
    "SV"     = "#66C2A5",
    "Indel"  = "#A6761D"
  )
  
  event_linetypes <- c(
    "Reference" = "solid",
    "CpG gain"  = "solid",
    "CpG loss"  = "dashed"
  )
  
  # --------------------------------------------------
  # 7. Density plot
  # --------------------------------------------------
  ggplot(
    df,
    aes(
      x = Value,
      colour = Variant_class,
      linetype = Event,
      group = interaction(Variant_class, Event)
    )
  ) +
    geom_density(
      linewidth = 0.75,
      adjust = 1,
      na.rm = TRUE
    ) +
    geom_vline(
      data = medians,
      aes(
        xintercept = median_value,
        colour = Variant_class,
        linetype = Event
      ),
      linewidth = 0.35,
      alpha = 0.65,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = class_colours,
      breaks = c("GRCh38", "SNV", "Indel", "SV"),
      name = NULL
    ) +
    scale_linetype_manual(
      values = event_linetypes,
      breaks = c("CpG gain", "CpG loss"),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.01))
    ) +
    coord_cartesian(xlim = c(0, 100)) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.06))
    ) +
    labs(
      x = "Methylation level",
      y = "Density",
      title = label
    ) +
    guides(
      colour = guide_legend(
        order = 1,
        override.aes = list(linetype = "solid", linewidth = 0.9)
      ),
      linetype = guide_legend(
        order = 2,
        override.aes = list(colour = "black", linewidth = 0.9)
      )
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      axis.text = element_text(colour = "black", size = 9),
      axis.title = element_text(colour = "black", size = 10),
      plot.title = element_text(
        size = 10.5,
        hjust = 0.5,
        face = "bold"
      ),
      legend.position = "top",
      legend.box = "vertical",
      legend.text = element_text(size = 8.5),
      legend.key.width = unit(0.9, "cm"),
      legend.key.height = unit(0.30, "cm"),
      panel.grid = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}

p_pat <- make_density_plot_obj(
  prefix = "HG002.pat",
  label = "Paternal"
)

p_mat <- make_density_plot_obj(
  prefix = "HG002.mat",
  label = "Maternal"
)

combined <- (p_pat / p_mat) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(
  "density_HG002_pat_mat_by_variant.pdf",
  plot = combined,
  width = 4.2,
  height = 3.5
)
