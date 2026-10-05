## Cayman Descriptive Analysis & QC Plots
library(tidyverse)
library(ggsci)
library(ggpubr)
library(rstatix)
library(lme4)
library(lmerTest)

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family) +
     theme(plot.title = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA, fill = NA),
           plot.background = element_rect(colour = NA, fill = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold", size = rel(0.8)),
           axis.title.y = element_text(angle=90, vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(size = rel(0.7)),
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing  = unit(0, "cm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")))
}

df_raw <- rio::import("data/shotgun/cayman_results/oct2026_results/families_rpkm_table.tsv") |> select(-HELIBA_103370, -HELIFU_103370)
stats <- rio::import("data/shotgun/cayman_results/oct2026_results/sample_statistics.tsv") |>
  # S103370's baseline CAZy profile essentially failed (43 families, 0.13% CAZy
  # reads vs a cohort median of 350 / 3.7%, despite normal overall read depth);
  # excluded here to match the family-table exclusion used in scripts 2-4.
  filter(!sample %in% c("HELIBA_103370", "HELIFU_103370")) |>
  mutate(timepoint = case_when(str_detect(sample, "HELIBA") ~ "Baseline",
                                str_detect(sample, "HELIFU") ~ "Follow-up"),
          timepoint = as.factor(timepoint)) |>
  rename(sampleID = sample) |> select(-timepoint)
names(stats)
str(stats$sampleID)
clinical <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
names(clinical)
str(clinical$sampleID)
stats <- left_join(stats, clinical, by = "sampleID") |> filter(!is.na(timepoint))
head(stats)[1:5,1:10]

# Create QC output directory
dir.create("results/4_functional_change/cayman", showWarnings = FALSE, recursive = TRUE)
qc_dir <- "results/4_functional_change/cayman/qc"
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Histogram of pct_cazy_reads
gg_pct_cazy_hist <- gghistogram(stats, x = "pct_cazy_reads", bins = 30, fill = "#69b3a2", color = "black") +
  theme_Publication() +
  facet_wrap(~timepoint) +
  labs(title = "% CAZy Reads per Sample", x = "% CAZy Reads", y = "Count")
ggsave(file.path(qc_dir, "pct_cazy_reads_histogram.pdf"), gg_pct_cazy_hist, width = 6, height = 4)

# 2. Boxplot of pct_cazy_reads by timepoint (if timepoint column exists)
gg_pct_cazy_box <- ggboxplot(stats, x = "timepoint", y = "pct_cazy_reads", fill = "timepoint", width = 0.4) +
    theme_Publication() +
    stat_compare_means() +
    scale_fill_jco() +
    labs(title = "% CAZy Reads by Timepoint", x = "Timepoint", y = "% CAZy Reads")
  ggsave(file.path(qc_dir, "pct_cazy_reads_by_timepoint_boxplot.pdf"), gg_pct_cazy_box, width = 4, height = 6)

pwc_pct <- stats %>%
  group_by(timepoint) %>%
  wilcox_test(pct_cazy_reads ~ EthnicityTot) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "timepoint", dodge = 0.8)

(gg_pct_cazy_box <- ggboxplot(stats, x = "timepoint", y = "pct_cazy_reads", fill = "EthnicityTot", width = 0.4) +
    theme_Publication() +
    stat_pvalue_manual(pwc_pct, label = "p.adj.signif", tip.length = 0, size = 5) +
    scale_fill_jco() +
    labs(title = "% CAZy Reads by Timepoint", x = "Timepoint", y = "% CAZy Reads", fill = ""))
ggsave(file.path(qc_dir, "pct_cazy_reads_by_timepoint_ethnicity_boxplot.pdf"), gg_pct_cazy_box, width = 5, height = 6)

# 3. Scatterplot of pct_cazy_reads vs filtered_reads with correlation

  cor_pctcazy_total <- cor.test(stats$filtered_reads, stats$pct_cazy_reads, method = "spearman")
  subtitle_pctcazy_total <- paste0("Spearman's rho = ", round(cor_pctcazy_total$estimate, 3),
                                   ", p = ", formatC(cor_pctcazy_total$p.value, format = "e", digits = 2))
  gg_pct_cazy_scatter <- ggplot(stats, aes(x = filtered_reads, y = pct_cazy_reads)) +
    geom_point(color = "#0072B2", alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    theme_Publication() +
    labs(title = "% CAZy Reads vs Filtered Reads", x = "Filtered Reads", y = "% CAZy Reads", subtitle = subtitle_pctcazy_total)
  ggsave(file.path(qc_dir, "pct_cazy_reads_vs_filtered_reads.pdf"), gg_pct_cazy_scatter, width = 6, height = 4)

# 4. Histogram of family richness (n_families; catalog size = 466 families)
gg_richness_hist <- gghistogram(stats, x = "n_families", bins = 30, fill = "#E69F00", color = "black") +
  theme_Publication() +
  facet_wrap(~timepoint) +
  labs(title = "Family richness", x = "Family richness (n_families)", y = "Count")
ggsave(file.path(qc_dir, "richness_histogram.pdf"), gg_richness_hist, width = 6, height = 4)

# 5. Boxplot of family richness by timepoint and ethnicity
pwc_richness <- stats %>%
  group_by(timepoint) %>%
  wilcox_test(n_families ~ EthnicityTot) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "timepoint", dodge = 0.8)

# LMM: ethnicity x timepoint interaction on family richness (repeated measures, 1|ID)
lmm_richness   <- lmerTest::lmer(n_families ~ EthnicityTot * timepoint + (1 | ID), data = stats)
lmm_richness_p <- summary(lmm_richness)$coefficients["EthnicityTotSouth-Asian Surinamese:timepointfollow-up", "Pr(>|t|)"]
lmm_richness_label <- paste0("LMM: Ethnicity × Timepoint p = ",
                              ifelse(lmm_richness_p < 0.001,
                                     formatC(lmm_richness_p, format = "e", digits = 2),
                                     formatC(lmm_richness_p, format = "f", digits = 3)))

(gg_richness_box <- ggboxplot(stats, x = "timepoint", y = "n_families", fill = "EthnicityTot", width = 0.4) +
    theme_Publication() +
    stat_pvalue_manual(pwc_richness, label = "p.adj.signif", tip.length = 0, size = 5) +
    scale_fill_jco() +
    labs(title = "Family richness", x = "Timepoint", y = "Number of families detected", fill = "",
         subtitle = lmm_richness_label) +
    theme(plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic")))
ggsave(file.path(qc_dir, "richness_by_timepoint_ethnicity_boxplot.pdf"), gg_richness_box, width = 5, height = 6)
