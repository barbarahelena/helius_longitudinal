## VFDB virulence factor comparison of Odoribacter splanchnicus bins
## Dutch vs South-Asian Surinamese, baseline and follow-up
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

#### Libraries ####
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggsci)
library(ggthemes)
library(grid)
library(ggpubr)

#### Theme ####
theme_Publication <- function(base_size = 14, base_family = "sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title        = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
      text              = element_text(),
      panel.background  = element_rect(colour = NA, fill = NA),
      plot.background   = element_rect(colour = NA, fill = NA),
      panel.border      = element_rect(colour = NA),
      axis.title        = element_text(face = "bold", size = rel(0.8)),
      axis.title.y      = element_text(angle = 90, vjust = 2),
      axis.title.x      = element_text(vjust = -0.2),
      axis.text         = element_text(size = rel(0.7)),
      axis.line         = element_line(colour = "black"),
      axis.ticks        = element_line(),
      panel.grid.major  = element_line(colour = "#f0f0f0"),
      panel.grid.minor  = element_blank(),
      legend.key        = element_rect(colour = NA),
      legend.position   = "bottom",
      legend.key.size   = unit(0.2, "cm"),
      legend.spacing    = unit(0, "cm"),
      strip.background  = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text        = element_text(face = "bold"),
      plot.caption      = element_text(size = rel(0.5), face = "italic"),
      plot.subtitle     = element_text(size = 8, hjust = 0.5, face = "italic")
    ))
}

#### Paths ####
vfdb_file   <- "data/shotgun/odoribacter_annotation/all_vfdb_results.txt"
fasta_file  <- "data/shotgun/odoribacter_annotation/VFDB_setB_pro.fas"
trans_file  <- "data/shotgun/odoribacter_annotation/bin_translation_table.tsv"
clin_file   <- "data/clinicaldata/clinicaldata_long.RDS"
batch_files <- c(
  "data/shotgun/odoribacter_annotation/bins_odoribacter_batch1.csv",
  "data/shotgun/odoribacter_annotation/bins_odoribacter_batch2.csv",
  "data/shotgun/odoribacter_annotation/bins_odoribacter_batch3.csv"
)
results_dir <- "results/3_species_change/4_odoribacter_anno"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

#### 1. Load bin translation table ####
trans <- read.delim(trans_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE) %>%
  rename(locus_prefix = locus_tag_prefix)

cat("Bins in translation table:", nrow(trans), "\n")

#### 2. Parse VFDB FASTA headers → annotation lookup ####
# Header format:
# >VFG037170(gb|WP_...) (gene_symbol) protein name [VF_name (VFxxxxx) - VF_category (VFCxxxxx)] [Organism]
#
# We extract: vfg_id, gene_symbol, vf_name, vf_id, vf_category, vfc_id, organism

raw_headers <- readLines(fasta_file) %>%
  keep(startsWith, ">") %>%
  sub("^>", "", .)

vfdb_anno <- tibble(header = raw_headers) %>%
  mutate(
    vfg_id      = str_extract(header, "^VFG\\d+"),
    gene_symbol = str_match(header, "\\)\\s+\\((\\w[^)]+)\\)")[, 2],
    # VF name and VF ID: [VF_name (VFxxxxx) - ...]
    vf_name     = str_match(header, "\\[([^\\[\\]]+)\\s+\\(VF\\d+\\)")[, 2],
    vf_id       = str_extract(header, "VF\\d+"),
    # VF category and VFC ID: [... - VF_category (VFCxxxxx)]
    vf_category = str_match(header, "-\\s+([^\\[\\]]+?)\\s+\\(VFC\\d+\\)")[, 2],
    vfc_id      = str_extract(header, "VFC\\d+"),
    organism    = str_match(header, "\\[([^\\[\\]]+)\\]\\s*$")[, 2]
  ) %>%
  dplyr::select(-header) %>%
  filter(!is.na(vfg_id)) %>%
  distinct(vfg_id, .keep_all = TRUE)  # one row per VFG ID

cat("VFG IDs in FASTA:", nrow(vfdb_anno), "\n")
cat("VF categories:", n_distinct(vfdb_anno$vf_category), "\n")
cat("VF names:", n_distinct(vfdb_anno$vf_name), "\n")

#### 3. Parse merged DIAMOND output ####
# Merged file: repeated # header blocks (one per bin).
# Bin identity recovered from query locus prefix (PREFIX_NNNNN → PREFIX).

col_names <- c("query", "subject", "pident", "length", "mismatch",
               "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

vfdb_hits <- read.table(
  vfdb_file,
  header       = FALSE,
  sep          = "\t",
  comment.char = "#",
  quote        = "",
  col.names    = col_names,
  stringsAsFactors = FALSE,
  fill         = TRUE
)

cat("Total VFDB hits loaded:", nrow(vfdb_hits), "\n")

#### 4. Filter hits and extract VFG IDs ####
# Thresholds: ≥30% identity and bitscore ≥50 (conservative for remote homologs).
MIN_PIDENT   <- 30
MIN_BITSCORE <- 50

vfdb_hits <- vfdb_hits %>%
  filter(pident >= MIN_PIDENT, bitscore >= MIN_BITSCORE) %>%
  mutate(
    locus_prefix = sub("_.*", "", query),
    vfg_id       = str_extract(subject, "^VFG\\d+")
  )

cat("Hits after filtering (pident >=", MIN_PIDENT, ", bitscore >=", MIN_BITSCORE, "):",
    nrow(vfdb_hits), "\n")

#### 5. Join DIAMOND hits with VFDB annotations ####
vfdb_hits <- vfdb_hits %>%
  left_join(vfdb_anno, by = "vfg_id")

unmatched_vfg <- sum(is.na(vfdb_hits$vf_category))
if (unmatched_vfg > 0)
  warning(unmatched_vfg, " hits have no VF category annotation (check FASTA coverage)")

cat("Hits with VF category annotation:", sum(!is.na(vfdb_hits$vf_category)),
    "/", nrow(vfdb_hits), "\n")

#### 6. Per-bin VF category proportions ####
# Denominator = total VFDB hits per bin (after filtering).
# This reflects the VF repertoire composition of each bin.

total_hits_per_bin <- vfdb_hits %>%
  filter(!is.na(vf_category)) %>%
  count(locus_prefix, name = "total_vf_hits")

vf_long <- vfdb_hits %>%
  filter(!is.na(vf_category)) %>%
  count(locus_prefix, vf_category, name = "n_hits") %>%
  left_join(total_hits_per_bin, by = "locus_prefix") %>%
  mutate(proportion = n_hits / total_vf_hits)

cat("Unique VF categories in data:", n_distinct(vf_long$vf_category), "\n")
cat("Unique bins with VF hits:", n_distinct(vf_long$locus_prefix), "\n")

#### 7. Load bin depths (presence/absence per sample) ####
depth_raw <- map_dfr(batch_files, function(f) {
  read.csv(f, check.names = FALSE) %>%
    dplyr::select(bin, starts_with("Depth "))
}) %>%
  mutate(bin_name = sub("\\.fa$", "", bin)) %>%
  dplyr::select(bin_name, starts_with("Depth "))

depth_long <- depth_raw %>%
  pivot_longer(cols = starts_with("Depth "),
               names_to  = "sampleID",
               values_to = "depth") %>%
  mutate(sampleID = sub("^Depth ", "", sampleID),
         depth    = replace_na(depth, 0))

depth_own <- depth_long %>%
  inner_join(trans %>% dplyr::select(bin_name, locus_prefix, subject_id),
             by = "bin_name") %>%
  filter(sampleID == paste0("HELIBA_", subject_id) |
         sampleID == paste0("HELIFU_",  subject_id))

#### 8. Join with clinical data ####
clin <- readRDS(clin_file) %>%
  filter(EthnicityTot %in% c("Dutch", "South-Asian Surinamese")) %>%
  droplevels() %>%
  dplyr::select(sampleID, ID, EthnicityTot, timepoint)

sample_map <- trans %>%
  dplyr::select(locus_prefix, subject_id) %>%
  mutate(
    sampleID_BA = paste0("HELIBA_", subject_id),
    sampleID_FU = paste0("HELIFU_", subject_id)
  ) %>%
  pivot_longer(cols = c(sampleID_BA, sampleID_FU),
               names_to  = NULL,
               values_to = "sampleID") %>%
  mutate(subject_id = as.character(subject_id))

bin_clin <- sample_map %>%
  inner_join(clin, by = "sampleID") %>%
  left_join(depth_own %>% dplyr::select(locus_prefix, sampleID, depth),
            by = c("locus_prefix", "sampleID")) %>%
  mutate(
    depth        = replace_na(depth, 0),
    present      = depth > 0,
    timepoint    = factor(timepoint, levels = c("baseline", "follow-up")),
    EthnicityTot = factor(EthnicityTot,
                          levels = c("Dutch", "South-Asian Surinamese"))
  )

cat("\nAll samples with clinical data (present AND absent):\n")
print(table(bin_clin$EthnicityTot, bin_clin$timepoint))

#### 9. Build analysis dataset ####
# Only bins where depth > 0 AND that have at least one VFDB hit.
# Complete grid: every present bin × every VF category, missing filled with 0.

present_bins <- bin_clin %>%
  filter(present, locus_prefix %in% unique(vf_long$locus_prefix)) %>%
  dplyr::select(locus_prefix, sampleID, EthnicityTot, timepoint, depth)

cat("\nPresent bins with VF hits (used for analysis):\n")
print(
  present_bins %>%
    distinct(locus_prefix, EthnicityTot, timepoint) %>%
    count(EthnicityTot, timepoint)
)

vf_per_bin <- expand.grid(
  locus_prefix = unique(present_bins$locus_prefix),
  vf_category  = unique(vf_long$vf_category),
  stringsAsFactors = FALSE
) %>%
  left_join(vf_long %>% dplyr::select(locus_prefix, vf_category, proportion),
            by = c("locus_prefix", "vf_category")) %>%
  mutate(proportion = replace_na(proportion, 0))

vf_df <- present_bins %>%
  left_join(vf_per_bin, by = "locus_prefix", relationship = "many-to-many")

#### 10. Statistical analysis ####
run_stats_vf <- function(df) {
  features <- unique(df$vf_category)

  results <- map_dfr(features, function(feat) {
    sub_df <- df %>% filter(vf_category == feat)

    ba_data  <- sub_df %>% filter(timepoint == "baseline")
    dutch_ba <- ba_data$proportion[ba_data$EthnicityTot == "Dutch"]
    sas_ba   <- ba_data$proportion[ba_data$EthnicityTot == "South-Asian Surinamese"]
    wx_ba <- tryCatch(
      wilcox.test(dutch_ba, sas_ba, exact = FALSE),
      error = function(e) list(statistic = NA, p.value = NA)
    )

    fu_data  <- sub_df %>% filter(timepoint == "follow-up")
    dutch_fu <- fu_data$proportion[fu_data$EthnicityTot == "Dutch"]
    sas_fu   <- fu_data$proportion[fu_data$EthnicityTot == "South-Asian Surinamese"]
    wx_fu <- tryCatch(
      wilcox.test(dutch_fu, sas_fu, exact = FALSE),
      error = function(e) list(statistic = NA, p.value = NA)
    )

    lmm_res <- tryCatch({
      mod <- lmer(
        proportion ~ EthnicityTot * timepoint + (1 | locus_prefix),
        data    = sub_df,
        REML    = FALSE,
        control = lmerControl(optimizer = "bobyqa")
      )
      coef_tab <- summary(mod)$coefficients

      eth_row <- grep("^EthnicityTot", rownames(coef_tab))
      eth_row <- eth_row[!grepl("timepoint", rownames(coef_tab)[eth_row])]
      if (length(eth_row) == 0) eth_row <- NA_integer_

      int_row <- grep("EthnicityTot.*timepoint|timepoint.*EthnicityTot",
                      rownames(coef_tab))
      if (length(int_row) == 0) int_row <- NA_integer_

      extract <- function(r, col) if (!is.na(r)) coef_tab[r, col] else NA_real_

      list(
        lmm_eth_estimate = extract(eth_row, "Estimate"),
        lmm_eth_se       = extract(eth_row, "Std. Error"),
        lmm_eth_pval     = extract(eth_row, "Pr(>|t|)"),
        lmm_int_estimate = extract(int_row, "Estimate"),
        lmm_int_se       = extract(int_row, "Std. Error"),
        lmm_int_pval     = extract(int_row, "Pr(>|t|)")
      )
    }, error = function(e) {
      list(lmm_eth_estimate = NA_real_, lmm_eth_se = NA_real_, lmm_eth_pval = NA_real_,
           lmm_int_estimate = NA_real_, lmm_int_se = NA_real_, lmm_int_pval = NA_real_)
    })

    tibble(
      vf_category          = feat,
      wilcox_baseline_stat = wx_ba$statistic,
      wilcox_baseline_p    = wx_ba$p.value,
      wilcox_fu_stat       = wx_fu$statistic,
      wilcox_fu_p          = wx_fu$p.value,
      lmm_eth_estimate     = lmm_res$lmm_eth_estimate,
      lmm_eth_se           = lmm_res$lmm_eth_se,
      lmm_eth_pval         = lmm_res$lmm_eth_pval,
      lmm_int_estimate     = lmm_res$lmm_int_estimate,
      lmm_int_se           = lmm_res$lmm_int_se,
      lmm_int_pval         = lmm_res$lmm_int_pval
    )
  }) %>%
    mutate(
      wilcox_baseline_fdr = p.adjust(wilcox_baseline_p, method = "BH"),
      wilcox_fu_fdr       = p.adjust(wilcox_fu_p,       method = "BH"),
      lmm_eth_fdr         = p.adjust(lmm_eth_pval,      method = "BH"),
      lmm_int_fdr         = p.adjust(lmm_int_pval,      method = "BH")
    )

  return(results)
}

cat("\nRunning VF category statistics...\n")
vf_results <- run_stats_vf(vf_df)

#### 11. Save results ####
write.csv(vf_results, file.path(results_dir, "vfdb_category_stats.csv"), row.names = FALSE)

cat("\nSignificant VF categories (any comparison, FDR < 0.05):\n")
sig_vf <- vf_results %>%
  filter(wilcox_baseline_fdr < 0.05 | wilcox_fu_fdr < 0.05 |
         lmm_eth_fdr < 0.05 | lmm_int_fdr < 0.05)
print(sig_vf)

#### 12. Visualisation ####
jco_cols <- pal_jco()(2)
names(jco_cols) <- c("Dutch", "South-Asian Surinamese")

## --- 12a. Boxplots for significant VF categories ---
if (nrow(sig_vf) > 0) {
  plot_vf_box <- vf_df %>%
    filter(vf_category %in% sig_vf$vf_category)

  p_vf_box <- ggplot(plot_vf_box,
                     aes(x = EthnicityTot, y = proportion, fill = EthnicityTot)) +
    geom_boxplot(outlier.size = 0.8, width = 0.5) +
    scale_fill_manual(values = jco_cols, name = "Ethnicity") +
    facet_grid(vf_category ~ timepoint, scales = "free_y") +
    labs(
      title    = "VF category proportions",
      subtitle = "Significant differences (FDR < 0.05), Odoribacter splanchnicus",
      x = "", y = "Proportion of VFDB hits",
      caption  = "Proportion = hits in VF category / total VFDB hits per bin"
    ) +
    theme_Publication() +
    theme(axis.text.x  = element_text(angle = 30, hjust = 1),
          strip.text.y = element_text(angle = 0, size = rel(0.7)))

  ggsave(
    file.path(results_dir, "vfdb_category_boxplot.pdf"),
    plot   = p_vf_box,
    width  = 8, height = max(4, nrow(sig_vf) * 1.5 + 2)
  )

  ## Figure panel: baseline only, no legend
  pl_fig3_F <- ggplot(plot_vf_box %>% filter(timepoint == "baseline"),
                      aes(x = EthnicityTot, y = proportion, fill = EthnicityTot)) +
    geom_boxplot(outlier.size = 0.8, width = 0.5) +
    stat_compare_means(comparisons = list(c("Dutch", "South-Asian Surinamese")),
                       method = "wilcox.test", label = "p.signif", tip.length = 0) +
    scale_fill_manual(values = jco_cols, guide = "none") +
    facet_wrap(~ vf_category, scales = "free_y") +
    scale_y_continuous(expand = expansion(add = c(0, 0.010))) +
    labs(title = "VFDB: Odoribacter splanchnicus", x = "", y = "Proportion of VFDB hits") +
    theme_Publication() +
    theme(strip.text  = element_text(size = rel(0.8)))

  cat("VF category boxplot saved.\n")
} else {
  cat("No significantly different VF categories (FDR < 0.05).\n")
}

## --- 12b. LMM ethnicity heatmap (top 20 by absolute effect) ---
vf_lmm_heat <- vf_results %>%
  filter(!is.na(lmm_eth_estimate)) %>%
  slice_max(abs(lmm_eth_estimate), n = 20) %>%
  mutate(
    sig_label   = case_when(
      lmm_eth_fdr < 0.001 ~ "***",
      lmm_eth_fdr < 0.01  ~ "**",
      lmm_eth_fdr < 0.05  ~ "*",
      TRUE                 ~ ""
    ),
    vf_category = factor(vf_category,
                         levels = vf_results %>%
                           filter(!is.na(lmm_eth_estimate)) %>%
                           slice_max(abs(lmm_eth_estimate), n = 20) %>%
                           arrange(lmm_eth_estimate) %>%
                           pull(vf_category))
  )

lim <- max(abs(vf_lmm_heat$lmm_eth_estimate), na.rm = TRUE)

p_vf_heat <- ggplot(vf_lmm_heat,
                    aes(x = 1, y = vf_category, fill = lmm_eth_estimate)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3, vjust = 0.5) +
  scale_fill_gradient2(
    low      = jco_cols[["Dutch"]],
    mid      = "white",
    high     = jco_cols[["South-Asian Surinamese"]],
    midpoint = 0,
    limits   = c(-lim, lim),
    name     = "LMM estimate\n(SAS vs Dutch)"
  ) +
  scale_x_continuous(breaks = NULL) +
  labs(
    title    = "Top 20 VF categories — LMM ethnicity effect",
    subtitle = paste0("Odoribacter splanchnicus bins  |  blue = higher in Dutch, ",
                      "yellow = higher in SAS\n* FDR<0.05  ** FDR<0.01  *** FDR<0.001"),
    x = "", y = "",
    caption  = "LMM: proportion ~ EthnicityTot * timepoint + (1|bin); estimate = SAS vs Dutch main effect"
  ) +
  theme_Publication() +
  theme(axis.text.y = element_text(size = rel(0.75)))

ggsave(
  file.path(results_dir, "vfdb_category_heatmap_lmm.pdf"),
  plot   = p_vf_heat,
  width  = 7, height = 9
)

## --- 12c. Baseline heatmap — Wilcoxon + mean proportion difference ---
vf_baseline_means <- vf_df %>%
  filter(timepoint == "baseline") %>%
  distinct(locus_prefix, vf_category, EthnicityTot, proportion) %>%
  group_by(vf_category, EthnicityTot) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = EthnicityTot, values_from = mean_prop, values_fill = 0) %>%
  mutate(diff_Dutch_SAS = Dutch - `South-Asian Surinamese`)

vf_ba_heat <- vf_results %>%
  filter(!is.na(wilcox_baseline_p)) %>%
  slice_max(abs(wilcox_baseline_stat), n = 20) %>%
  left_join(vf_baseline_means %>% dplyr::select(vf_category, diff_Dutch_SAS),
            by = "vf_category") %>%
  mutate(
    sig_label   = case_when(
      wilcox_baseline_fdr < 0.001 ~ "***",
      wilcox_baseline_fdr < 0.01  ~ "**",
      wilcox_baseline_fdr < 0.05  ~ "*",
      TRUE                        ~ ""
    ),
    vf_category = factor(vf_category,
                         levels = vf_results %>%
                           filter(!is.na(wilcox_baseline_p)) %>%
                           slice_max(abs(wilcox_baseline_stat), n = 20) %>%
                           left_join(vf_baseline_means %>%
                                       dplyr::select(vf_category, diff_Dutch_SAS),
                                     by = "vf_category") %>%
                           arrange(diff_Dutch_SAS) %>%
                           pull(vf_category))
  )

lim_ba <- max(abs(vf_ba_heat$diff_Dutch_SAS), na.rm = TRUE)

p_vf_heat_ba <- ggplot(vf_ba_heat,
                       aes(x = 1, y = vf_category, fill = diff_Dutch_SAS)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3, vjust = 0.5) +
  scale_fill_gradient2(
    low      = jco_cols[["South-Asian Surinamese"]],
    mid      = "white",
    high     = jco_cols[["Dutch"]],
    midpoint = 0,
    limits   = c(-lim_ba, lim_ba),
    name     = "Mean proportion\nDutch − SAS (baseline)"
  ) +
  scale_x_continuous(breaks = NULL) +
  labs(
    title    = "Top 20 VF categories — ethnicity difference at baseline",
    subtitle = paste0("Odoribacter splanchnicus bins  |  blue = higher in Dutch, ",
                      "yellow = higher in SAS\n* Wilcoxon FDR<0.05  ** FDR<0.01  *** FDR<0.001"),
    x = "", y = "",
    caption  = "Ranked by Wilcoxon W statistic; colour = mean proportion Dutch − SAS at baseline"
  ) +
  theme_Publication() +
  theme(axis.text.y = element_text(size = rel(0.75)))

ggsave(
  file.path(results_dir, "vfdb_category_heatmap_baseline.pdf"),
  plot   = p_vf_heat_ba,
  width  = 7, height = 9
)

## --- 12d. Top VF names per category (interpretation aid) ---
vf_name_summary <- vfdb_hits %>%
  filter(!is.na(vf_category), !is.na(vf_name)) %>%
  count(vf_category, vf_name, name = "n_hits") %>%
  group_by(vf_category) %>%
  slice_max(n_hits, n = 5) %>%
  arrange(vf_category, desc(n_hits))

write.csv(vf_name_summary,
          file.path(results_dir, "vfdb_top_vf_names_per_category.csv"),
          row.names = FALSE)

cat("\nAll plots saved. Results written to:", results_dir, "\n")

#### 13. Per-VF-name analysis ####
# VF names are more granular than categories; many will be rare across bins.
# We apply a minimum prevalence filter: only analyse VF names detected (≥1 hit)
# in at least MIN_PREV bins.

MIN_PREV <- 5   # minimum number of bins a VF name must be present in

## 13a. Count presence per VF name (across all bins with any VF hits)
vfname_prevalence <- vfdb_hits %>%
  filter(!is.na(vf_name), locus_prefix %in% unique(present_bins$locus_prefix)) %>%
  distinct(locus_prefix, vf_name) %>%
  count(vf_name, name = "n_bins")

prevalent_vf_names <- vfname_prevalence %>%
  filter(n_bins >= MIN_PREV) %>%
  pull(vf_name)

cat("\nVF names with prevalence >=", MIN_PREV, "bins:", length(prevalent_vf_names), "\n")

## 13b. Per-bin VF name proportions (same denominator: total VFDB hits per bin)
vfname_long <- vfdb_hits %>%
  filter(!is.na(vf_name), locus_prefix %in% unique(present_bins$locus_prefix)) %>%
  filter(vf_name %in% prevalent_vf_names) %>%
  count(locus_prefix, vf_name, name = "n_hits") %>%
  left_join(total_hits_per_bin, by = "locus_prefix") %>%
  mutate(proportion = n_hits / total_vf_hits)

## 13c. Complete grid: every present bin × every prevalent VF name
vfname_per_bin <- expand.grid(
  locus_prefix = unique(present_bins$locus_prefix),
  vf_name      = prevalent_vf_names,
  stringsAsFactors = FALSE
) %>%
  left_join(vfname_long %>% dplyr::select(locus_prefix, vf_name, proportion),
            by = c("locus_prefix", "vf_name")) %>%
  mutate(proportion = replace_na(proportion, 0))

vfname_df <- present_bins %>%
  left_join(vfname_per_bin, by = "locus_prefix", relationship = "many-to-many")

## 13d. Statistics (same structure as run_stats_vf, iterating over vf_name)
run_stats_vfname <- function(df) {
  features <- unique(df$vf_name)

  results <- map_dfr(features, function(feat) {
    sub_df <- df %>% filter(vf_name == feat)

    ba_data  <- sub_df %>% filter(timepoint == "baseline")
    dutch_ba <- ba_data$proportion[ba_data$EthnicityTot == "Dutch"]
    sas_ba   <- ba_data$proportion[ba_data$EthnicityTot == "South-Asian Surinamese"]
    wx_ba <- tryCatch(
      wilcox.test(dutch_ba, sas_ba, exact = FALSE),
      error = function(e) list(statistic = NA, p.value = NA)
    )

    fu_data  <- sub_df %>% filter(timepoint == "follow-up")
    dutch_fu <- fu_data$proportion[fu_data$EthnicityTot == "Dutch"]
    sas_fu   <- fu_data$proportion[fu_data$EthnicityTot == "South-Asian Surinamese"]
    wx_fu <- tryCatch(
      wilcox.test(dutch_fu, sas_fu, exact = FALSE),
      error = function(e) list(statistic = NA, p.value = NA)
    )

    lmm_res <- tryCatch({
      mod <- lmer(
        proportion ~ EthnicityTot * timepoint + (1 | locus_prefix),
        data    = sub_df,
        REML    = FALSE,
        control = lmerControl(optimizer = "bobyqa")
      )
      coef_tab <- summary(mod)$coefficients

      eth_row <- grep("^EthnicityTot", rownames(coef_tab))
      eth_row <- eth_row[!grepl("timepoint", rownames(coef_tab)[eth_row])]
      if (length(eth_row) == 0) eth_row <- NA_integer_

      int_row <- grep("EthnicityTot.*timepoint|timepoint.*EthnicityTot",
                      rownames(coef_tab))
      if (length(int_row) == 0) int_row <- NA_integer_

      extract <- function(r, col) if (!is.na(r)) coef_tab[r, col] else NA_real_

      list(
        lmm_eth_estimate = extract(eth_row, "Estimate"),
        lmm_eth_se       = extract(eth_row, "Std. Error"),
        lmm_eth_pval     = extract(eth_row, "Pr(>|t|)"),
        lmm_int_estimate = extract(int_row, "Estimate"),
        lmm_int_se       = extract(int_row, "Std. Error"),
        lmm_int_pval     = extract(int_row, "Pr(>|t|)")
      )
    }, error = function(e) {
      list(lmm_eth_estimate = NA_real_, lmm_eth_se = NA_real_, lmm_eth_pval = NA_real_,
           lmm_int_estimate = NA_real_, lmm_int_se = NA_real_, lmm_int_pval = NA_real_)
    })

    tibble(
      vf_name              = feat,
      wilcox_baseline_stat = wx_ba$statistic,
      wilcox_baseline_p    = wx_ba$p.value,
      wilcox_fu_stat       = wx_fu$statistic,
      wilcox_fu_p          = wx_fu$p.value,
      lmm_eth_estimate     = lmm_res$lmm_eth_estimate,
      lmm_eth_se           = lmm_res$lmm_eth_se,
      lmm_eth_pval         = lmm_res$lmm_eth_pval,
      lmm_int_estimate     = lmm_res$lmm_int_estimate,
      lmm_int_se           = lmm_res$lmm_int_se,
      lmm_int_pval         = lmm_res$lmm_int_pval
    )
  }) %>%
    mutate(
      wilcox_baseline_fdr = p.adjust(wilcox_baseline_p, method = "BH"),
      wilcox_fu_fdr       = p.adjust(wilcox_fu_p,       method = "BH"),
      lmm_eth_fdr         = p.adjust(lmm_eth_pval,      method = "BH"),
      lmm_int_fdr         = p.adjust(lmm_int_pval,      method = "BH")
    )

  return(results)
}

cat("Running VF name statistics (", length(prevalent_vf_names), " names)...\n")
vfname_results <- run_stats_vfname(vfname_df)

## 13e. Add VF category annotation to results
# Use the most common category per VF name (some VF names map to >1 category
# in the FASTA, which would create duplicate rows via left_join).
vfname_results <- vfname_results %>%
  left_join(
    vfdb_anno %>%
      filter(!is.na(vf_name), !is.na(vf_category)) %>%
      count(vf_name, vf_category, name = "n") %>%
      group_by(vf_name) %>%
      slice_max(n, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      dplyr::select(vf_name, vf_category),
    by = "vf_name"
  )

## 13f. Save results
write.csv(vfname_results,
          file.path(results_dir, "vfdb_name_stats.csv"),
          row.names = FALSE)

cat("\nSignificant VF names (any comparison, FDR < 0.05):\n")
sig_vfname <- vfname_results %>%
  filter(wilcox_baseline_fdr < 0.05 | wilcox_fu_fdr < 0.05 |
         lmm_eth_fdr < 0.05 | lmm_int_fdr < 0.05)
print(sig_vfname %>% dplyr::select(vf_name, vf_category, lmm_eth_estimate,
                                    lmm_eth_fdr, wilcox_baseline_fdr, wilcox_fu_fdr))

## --- 13g. LMM ethnicity heatmap (top 25 VF names by absolute effect) ---
top_n_vfname <- 25

vfname_lmm_heat <- vfname_results %>%
  filter(!is.na(lmm_eth_estimate)) %>%
  slice_max(abs(lmm_eth_estimate), n = top_n_vfname, with_ties = FALSE) %>%
  mutate(
    sig_label = case_when(
      lmm_eth_fdr < 0.001 ~ "***",
      lmm_eth_fdr < 0.01  ~ "**",
      lmm_eth_fdr < 0.05  ~ "*",
      TRUE                 ~ ""
    ),
    vf_name = factor(vf_name,
                     levels = vfname_results %>%
                       filter(!is.na(lmm_eth_estimate)) %>%
                       slice_max(abs(lmm_eth_estimate), n = top_n_vfname, with_ties = FALSE) %>%
                       arrange(lmm_eth_estimate) %>%
                       pull(vf_name))
  )

lim_vfname <- max(abs(vfname_lmm_heat$lmm_eth_estimate), na.rm = TRUE)

if (nrow(vfname_lmm_heat) > 0 && !is.infinite(lim_vfname)) {
  p_vfname_heat <- ggplot(vfname_lmm_heat,
                          aes(x = 1, y = vf_name, fill = lmm_eth_estimate)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_label), size = 3, vjust = 0.5) +
    scale_fill_gradient2(
      low      = jco_cols[["Dutch"]],
      mid      = "white",
      high     = jco_cols[["South-Asian Surinamese"]],
      midpoint = 0,
      limits   = c(-lim_vfname, lim_vfname),
      name     = "LMM estimate\n(SAS vs Dutch)"
    ) +
    scale_x_continuous(breaks = NULL) +
    labs(
      title    = paste0("Top ", top_n_vfname, " VF names - LMM ethnicity effect"),
      subtitle = paste0("Odoribacter splanchnicus bins  |  blue = higher in Dutch, ",
                        "yellow = higher in SAS\n* FDR<0.05  ** FDR<0.01  *** FDR<0.001\n",
                        "VF names in >=", MIN_PREV, " bins"),
      x = "", y = "",
      caption  = "LMM: proportion ~ EthnicityTot * timepoint + (1|bin); estimate = SAS vs Dutch main effect"
    ) +
    theme_Publication() +
    theme(axis.text.y = element_text(size = rel(0.7)))

  ggsave(
    file.path(results_dir, "vfdb_name_heatmap_lmm.pdf"),
    plot   = p_vfname_heat,
    width  = 8, height = 10
  )
  cat("VF name LMM heatmap saved.\n")
} else {
  cat("No VF names with estimable LMM — LMM heatmap not produced.\n")
}

## --- 13h. Baseline heatmap (top 25 VF names by Wilcoxon W, colour = mean diff) ---
vfname_baseline_means <- vfname_df %>%
  filter(timepoint == "baseline") %>%
  distinct(locus_prefix, vf_name, EthnicityTot, proportion) %>%
  group_by(vf_name, EthnicityTot) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = EthnicityTot, values_from = mean_prop, values_fill = 0) %>%
  mutate(diff_Dutch_SAS = Dutch - `South-Asian Surinamese`)

vfname_ba_heat <- vfname_results %>%
  filter(!is.na(wilcox_baseline_p)) %>%
  slice_max(abs(wilcox_baseline_stat), n = top_n_vfname, with_ties = FALSE) %>%
  left_join(vfname_baseline_means %>% dplyr::select(vf_name, diff_Dutch_SAS),
            by = "vf_name") %>%
  mutate(
    sig_label = case_when(
      wilcox_baseline_fdr < 0.001 ~ "***",
      wilcox_baseline_fdr < 0.01  ~ "**",
      wilcox_baseline_fdr < 0.05  ~ "*",
      TRUE                        ~ ""
    ),
    vf_name = factor(vf_name,
                     levels = vfname_results %>%
                       filter(!is.na(wilcox_baseline_p)) %>%
                       slice_max(abs(wilcox_baseline_stat), n = top_n_vfname, with_ties = FALSE) %>%
                       left_join(vfname_baseline_means %>%
                                   dplyr::select(vf_name, diff_Dutch_SAS),
                                 by = "vf_name") %>%
                       arrange(diff_Dutch_SAS) %>%
                       pull(vf_name))
  )

lim_ba_vfname <- max(abs(vfname_ba_heat$diff_Dutch_SAS), na.rm = TRUE)

if (nrow(vfname_ba_heat) > 0 && !is.infinite(lim_ba_vfname)) {
  p_vfname_heat_ba <- ggplot(vfname_ba_heat,
                             aes(x = 1, y = vf_name, fill = diff_Dutch_SAS)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_label), size = 3, vjust = 0.5) +
    scale_fill_gradient2(
      low      = jco_cols[["South-Asian Surinamese"]],
      mid      = "white",
      high     = jco_cols[["Dutch"]],
      midpoint = 0,
      limits   = c(-lim_ba_vfname, lim_ba_vfname),
      name     = "Mean proportion\nDutch − SAS (baseline)"
    ) +
    scale_x_continuous(breaks = NULL) +
    labs(
      title    = paste0("Top ", top_n_vfname, " VF names - ethnicity difference at baseline"),
      subtitle = paste0("Odoribacter splanchnicus bins  |  blue = higher in Dutch, ",
                        "yellow = higher in SAS\n* Wilcoxon FDR<0.05  ** FDR<0.01  *** FDR<0.001\n",
                        "VF names in >=", MIN_PREV, " bins"),
      x = "", y = "",
      caption  = "Ranked by Wilcoxon W statistic; colour = mean proportion Dutch − SAS at baseline"
    ) +
    theme_Publication() +
    theme(axis.text.y = element_text(size = rel(0.7)))

  ggsave(
    file.path(results_dir, "vfdb_name_heatmap_baseline.pdf"),
    plot   = p_vfname_heat_ba,
    width  = 8, height = 10
  )
  cat("VF name baseline heatmap saved.\n")
} else {
  cat("No VF names with estimable baseline stats — baseline heatmap not produced.\n")
}

cat("\nVF name analysis complete. Results written to:", results_dir, "\n")
