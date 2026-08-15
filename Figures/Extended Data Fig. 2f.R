# Load and prep data
ins <- read.table("../../ML-coverage/1158/instrument_model.log", sep = "\t", header = TRUE)
data <- read.table("re_methyaltion_qc_1157.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# data$File <- sub("\\.5mc\\.jasmine\\.hifi_reads\\.model\\.pbmm2\\.combined\\.bed$", "", data$File)
data$File <- sub("\\.model\\.pbmm2\\.combined\\.bed$", ".bam", data$File)
# data$File <- gsub("\\.pm_v1\\.4\\.0", "", data$File)


lowq<-read.table("../../lowQ_1158_self.log")
# lowq$V1<-gsub("\\.pm_v1\\.4\\.0", "", lowq$V1)
data<-data[-which(data$File %in% lowq$V1),]

data$ins <- NA
data$ins[data$File %in% ins$filename[ins$instrument_model == "Revio"]] <- "Revio"
data$ins[data$File %in% ins$filename[ins$instrument_model == "Sequel II"]] <- "Sequel II"

#### === Boxplot for lab === ###
library(tidyverse)
library(ggpubr)

# === Load Data ===
ins <- read.table("lab_model.log", sep = "\t", header = TRUE)
# data <- read.table("re_methyaltion_qc_1151.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data$File <- sub("\\.5mc\\.jasmine\\.hifi_reads\\.bam$", "", data$File)
data$File <- sub("\\.model\\.pbmm2\\.combined\\.bed$", ".bam", data$File)
data$File <- gsub("\\.pm_v1\\.4\\.0", "", data$File)
data$File <- gsub("\\.ccs_v6\\.0\\.0", "", data$File)

data$File <- sub("\\.ccs_reprocessed.jasmine.5mc", "", data$File)
data$File <- sub("fixed_kinetics", "5mc", data$File)
data$File <- sub("m54329U_200610_234222", "m54329U_200610_234222.hifi_reads", data$File)
data$File <- sub("m54329U_200612_200443", "m54329U_200612_200443.hifi_reads", data$File)
data$File <- sub("m54329U_200614_021746", "m54329U_200614_021746.hifi_reads", data$File)
data$File <- sub("m54329U_200615_084313", "m54329U_200615_084313.hifi_reads", data$File)

data <- data %>%
  mutate(
    File = sub("\\.5mc\\.jasmine\\.hifi_reads\\.bam$", "", File),
    File = sub("\\.model\\.pbmm2\\.combined\\.bed$", ".bam", File),
    File = gsub("\\.pm_v1\\.4\\.0", "", File),
    File = gsub("\\.ccs_v6\\.0\\.0", "", File),
    File = sub("\\.ccs_reprocessed.jasmine.5mc", "", File),
    File = sub("fixed_kinetics", "5mc", File),
    File = sub("m54329U_200610_234222", "m54329U_200610_234222.hifi_reads", File),
    File = sub("m54329U_200612_200443", "m54329U_200612_200443.hifi_reads", File),
    File = sub("m54329U_200614_021746", "m54329U_200614_021746.hifi_reads", File),
    File = sub("m54329U_200615_084313", "m54329U_200615_084313.hifi_reads", File),
    File = sub("HG02984.m64055e_220517_182131.5mc.hifi_reads.bam", "HG02984.m64055e_220517_182131.fixed_kinetics.hifi_reads.bam", File),
    File = sub("m54329U_210507_224951-bc1009_BAK8A_OA.5mc.hifi_reads.bam", "m54329U_210507_224951-bc1009_BAK8A_OA.fixed_kinetics.hifi_reads.bam", File),
    File = sub("m54329U_210507_224951-bc1010_BAK8A_OA.5mc.hifi_reads.bam", "m54329U_210507_224951-bc1010_BAK8A_OA.fixed_kinetics.hifi_reads.bam", File),           
    File = sub("m54329U_220203_054640-bc1018.5mc.hifi_reads.bam", "m54329U_220203_054640-bc1018.fixed_kinetics.hifi_reads.bam", File),           
    # File = sub("m54329U_220903_190900-bc2053.5mc.hifi_reads.bam", "m54329U_220903_190900-bc2055.5mc.hifi_reads.bam", File),
    # File = sub("m54329U_220905_144015-bc2053.5mc.hifi_reads.bam", "m54329U_220905_144015-bc2055.5mc.hifi_reads.bam", File),
    # File = sub("m54329U_220907_100858-bc2053.5mc.hifi_reads.bam", "m54329U_220907_100858-bc2055.5mc.hifi_reads.bam", File),
    File = sub("m64043_210901_164337-bc1022.5mc.hifi_reads.bam", "m64043_210901_164337-bc1022.fixed_kinetics.hifi_reads.bam", File),
    File = sub("m64043_210903_163032-bc1022.5mc.hifi_reads.bam", "m64043_210903_163032-bc1022.fixed_kinetics.hifi_reads.bam", File),
    File = sub("m64043_220726_203720-bc1010.5mc.hifi_reads.bam", "m64043_220726_203720-bc1010.fixed_kinetics.hifi_reads.bam", File),
    File = sub("m64043_220728_173215-bc1016.5mc.hifi_reads.bam", "m64043_220728_173215-bc1016.fixed_kinetics.hifi_reads.bam", File),
    File = sub("m64076_200211_192227.hifi_reads.bam", "m64076_200211_192227.hifi_reads.ccs_v6.0.0.pm_v1.4.0.bam", File),
    File = sub("m64076_210502_044702.hifi_reads.bam", "m64076_210502_044702.hifi_reads.pm_v1.4.0.bam", File),
    File = sub("m64076_211215_225159.hifi_reads.bam", "m64076_211215_225159.hifi_reads.pm_v1.4.0.bam", File),
    File = sub("m64076_220216_013707-bc1016.5mc.hifi_reads.bam", "m64076_220216_013707-bc1016.fixed_kinetics.hifi_reads.bam", File),
    File = sub("m64076_220216_013707-bc1016.5mc.hifi_reads.bam", "m64076_220216_013707-bc1016.fixed_kinetics.hifi_reads.bam", File),
    # File = sub("m84036_230317_175945_s2.hifi_reads.bc2012.bam", "m84036_230317_175945_s2.hifi_reads.bc2013.bam", File),
    File = sub("m84081_231222_081401_s1.hifi_reads.bc2022.bam", "m84081_231222_081401_s1.hifi_reads.bc2022.pm_v1.4.0.bam", File),
    File = sub("m84081_231222_081401_s1.hifi_reads.bc2027.bam", "m84081_231222_081401_s1.hifi_reads.bc2027.pm_v1.4.0.bam", File),
  )

ins$filename <- sub("fixed_kinetics", "5mc", ins$filename)

data<-data[data$ins=="Revio",]
# === Merge Instrument Info ===
data <- data %>%
  left_join(ins %>% dplyr::select(filename, production), by = c("File" = "filename"))
names(data)[ncol(data)]<-"Lab"

write.table(data$File[is.na(data$Lab)],
            file = "missing_lab_files.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# === Filter and Clean ===
data <- data %>%
  # filter(Outlier == 0) %>%
  mutate(
    ML = as.numeric(Median),
    # Coverage = as.numeric(Coverage),
    Lab = as.character(Lab)
  )

# === Optional: Group lab names ===
data <- data %>%
  mutate(
    LabGroup = case_when(
      str_detect(Lab, "WUSTL") ~ "WUSTL",
      str_detect(Lab, "UW") ~ "UW",
      str_detect(Lab, "RU") ~ "RU",
      str_detect(Lab, "AMED") ~ "AMED",
      str_detect(Lab, "UCSC") ~ "UCSC",
      str_detect(Lab, "NISC") ~ "NISC",
      str_detect(Lab, "NISC") ~ "NISC",
      str_detect(Lab, "HTGM") ~ "HTGM",
      str_detect(Lab, "PacBio") ~ "PacBio",
      TRUE ~ "Other"
    )
  )

# === Filter for plotting (at least 2 samples per group) ===
df_filtered <- data %>%
  filter(!is.na(LabGroup)) %>%
  group_by(LabGroup) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  mutate(LabGroup = factor(LabGroup))

# === Generate Comparisons ===
lab_levels <- levels(df_filtered$LabGroup)
comparisons <- if (length(lab_levels) >= 2) combn(lab_levels, 2, simplify = FALSE) else list()

# === Define Colors ===
palette <- RColorBrewer::brewer.pal(max(3, length(lab_levels)), "Set2")
names(palette) <- lab_levels
palette <- c(
  Other = "#66C2A5",
  RU    = "#FC8D62",
  UCSC  = "#8DA0CB",
  UW    = "#E78AC3",
  WUSTL = "#A6D854"
)

# === Plot ===
p <- ggplot(df_filtered, aes(x = LabGroup, y = ML, fill = LabGroup)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = comparisons,
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  scale_fill_manual(values = palette) +
  labs(
    x = NULL,
    y = "Methyaltion median value"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave("methylation_lab_Revio.pdf", plot = p, width = 5, height = 4)

mean(data$Median[data$LabGroup=="AMED"],na.rm=T)
mean(data$Median[data$LabGroup=="HTGM"],na.rm=T)
mean(data$Median[data$LabGroup=="NISC"],na.rm=T)
mean(data$Median[data$LabGroup=="Other"],na.rm=T)
mean(data$Median[data$LabGroup=="PacBio"],na.rm=T)
mean(data$Median[data$LabGroup=="RU"],na.rm=T)
mean(data$Median[data$LabGroup=="UCSC"],na.rm=T)
mean(data$Median[data$LabGroup=="UW"],na.rm=T)
mean(data$Median[data$LabGroup=="WUSTL"],na.rm=T)

median(data$Median[data$LabGroup=="NISC"],na.rm=T)
median(data$Median[data$LabGroup=="Other"],na.rm=T)
median(data$Median[data$LabGroup=="PacBio"],na.rm=T)
median(data$Median[data$LabGroup=="RU"],na.rm=T)
median(data$Median[data$LabGroup=="UCSC"],na.rm=T)
median(data$Median[data$LabGroup=="UW"],na.rm=T)
median(data$Median[data$LabGroup=="WUSTL"],na.rm=T)
median(data$Median[data$LabGroup!=c("NISC","PacBio")],na.rm=T)

wilcox.test(data$Median[data$LabGroup!=c("NISC","PacBio")],data$Median[data$LabGroup=="NISC"])

wilcox.test(data$Median[data$LabGroup!=c("PacBio")],data$Median[data$LabGroup=="PacBio"])


median(data$Median[data$LabGroup!=c("PacBio")])

median(data$Median[data$LabGroup==c("PacBio")])


df_test <- data %>%
  mutate(
    Median = as.numeric(Median),
    CompareGroup = case_when(
      LabGroup == "PacBio" ~ "PacBio",
      !LabGroup %in% c("NISC", "PacBio") ~ "Other labs",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(CompareGroup),
    is.finite(Median)
  ) %>%
  mutate(
    CompareGroup = factor(CompareGroup, levels = c("Other labs", "PacBio"))
  )

p <- ggplot(df_test, aes(x = CompareGroup, y = Median)) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    fill = "grey88",
    color = "black",
    linewidth = 0.35
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    shape = 21,
    size = 1.6,
    fill = "grey35",
    color = "black",
    stroke = 0.15,
    alpha = 0.65
  ) +
  stat_compare_means(
    comparisons = list(c("Other labs", "PacBio")),
    method = "wilcox.test",
    label = "p.format",
    size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Median methylation level (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      size = 10,
      color = "black",
      angle = 30,
      hjust = 1
    ),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    panel.grid = element_blank(),
    legend.position = "none"
  )


ggsave(
  "PacBio_vs_other_labs_methylation_boxplot.pdf",
  plot = p,
  width = 1.6,
  height = 3.0,
  device = "pdf"
)
