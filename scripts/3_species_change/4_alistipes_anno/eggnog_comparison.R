## Functional annotation comparison of Alistipes putredinis bins
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
      plot.title   = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
      text         = element_text(),
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
anno_file      <- "data/shotgun/alistipes_annotation/eggnog/all_eggnog_results.annotations"
trans_file     <- "data/shotgun/alistipes_annotation/bin_translation_table.tsv"
clin_file      <- "data/clinicaldata_long.RDS"
results_dir    <- "results/3_species_change/4_alistipes_anno"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

#### 1. Load bin translation table ####
trans <- read.delim(trans_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE) %>%
  rename(bin_name       = bin_name,
         locus_prefix   = locus_tag_prefix,
         subject_id     = subject_id,
         best_sample    = best_sample)

cat("Bins in translation table:", nrow(trans), "\n")

#### 2. Parse eggNOG annotation file ####
# Merged file: many blocks of ## comments + repeated #query header lines (one per bin).
# Strategy: read column names from the FIRST #query line, then read the whole file
# with comment.char = "#" (skips both ## and all #query repetitions).
raw_lines <- readLines(anno_file, n = 500)   # enough to find first header
header_idx <- which(startsWith(raw_lines, "#query"))
if (length(header_idx) == 0)
  stop("No header line (#query) found in the first 500 lines of annotations file")

col_names <- strsplit(sub("^#", "", raw_lines[header_idx[1]]), "\t")[[1]]

anno <- read.table(
  anno_file,
  header       = FALSE,
  sep          = "\t",
  comment.char = "#",    # drops ALL lines starting with # (## comments + #query headers)
  quote        = "",
  col.names    = col_names,
  stringsAsFactors = FALSE,
  fill         = TRUE
)

cat("Total annotated genes loaded:", nrow(anno), "\n")

#### 3. Extract locus_tag prefix from query ####
# #query format: PREFIX_NNNNN  (e.g. FNBCDG_00001 → prefix = FNBCDG)
anno <- anno %>%
  mutate(locus_prefix = sub("_.*", "", query))

cat("Unique locus prefixes in annotations:", n_distinct(anno$locus_prefix), "\n")
cat("Unique locus prefixes in translation table:", n_distinct(trans$locus_prefix), "\n")
unmatched <- setdiff(anno$locus_prefix, trans$locus_prefix)
if (length(unmatched) > 0)
  warning(length(unmatched), " locus prefix(es) in annotations not found in translation table: ",
          paste(head(unmatched, 5), collapse = ", "))

#### 4. Calculate per-bin functional profiles ####

## Helper: clean a delimited field — return NA for "-" or blank entries
clean_field <- function(x) {
  x <- trimws(x)
  ifelse(x == "" | x == "-", NA_character_, x)
}

## Total annotated genes per bin (denominator):
## genes where at least KEGG_Module or COG_category is not "-"/blank
anno_filt <- anno %>%
  mutate(
    kegg_clean = clean_field(KEGG_Module),
    cog_clean  = clean_field(COG_category)
  ) %>%
  filter(!is.na(kegg_clean) | !is.na(cog_clean))

total_per_bin <- anno_filt %>%
  count(locus_prefix, name = "total_annotated")

## 4a. KEGG Module proportions per bin
kegg_long <- anno_filt %>%
  filter(!is.na(kegg_clean)) %>%
  mutate(module = strsplit(kegg_clean, ",")) %>%
  unnest(module) %>%
  mutate(module = trimws(module)) %>%
  filter(module != "" & module != "-") %>%
  count(locus_prefix, module, name = "n_genes") %>%
  left_join(total_per_bin, by = "locus_prefix") %>%
  mutate(proportion = n_genes / total_annotated)

cat("Unique KEGG modules:", n_distinct(kegg_long$module), "\n")

#### 4c. Fetch KEGG module descriptions ####
# Download all module descriptions once; cache locally to avoid repeat API calls.
kegg_lookup_file <- file.path(results_dir, "kegg_module_lookup.rds")

if (file.exists(kegg_lookup_file)) {
  module_lookup <- readRDS(kegg_lookup_file)
  cat("Loaded KEGG module lookup from cache.\n")
} else {
  cat("Fetching KEGG module descriptions from REST API (requires internet)...\n")
  # KEGG REST: https://rest.kegg.jp/list/module  returns "md:M00001\tDescription\n"
  raw <- tryCatch(
    readLines("https://rest.kegg.jp/list/module", warn = FALSE),
    error = function(e) { warning("KEGG API unavailable: ", e$message); character(0) }
  )
  if (length(raw) > 0) {
    module_lookup <- tibble(raw = raw) %>%
      separate(raw, into = c("module", "description"), sep = "\t", extra = "merge") %>%
      mutate(module = sub("^md:", "", module))
    saveRDS(module_lookup, kegg_lookup_file)
    cat("KEGG module lookup cached to", kegg_lookup_file, "\n")
  } else {
    module_lookup <- tibble(module = character(), description = character())
    cat("Warning: KEGG lookup empty — module IDs will be used as labels.\n")
  }
}

# Attach descriptions to kegg_long; keep bare module ID where lookup misses
kegg_long <- kegg_long %>%
  left_join(module_lookup, by = "module") %>%
  mutate(module_label = if_else(is.na(description), module,
                                paste0(module, ": ", description)))

cat("Modules with description:", sum(!is.na(kegg_long$description)), "/",
    nrow(kegg_long), "\n")

#### 5. Load bin depths (presence / absence per sample) ####
# The batch CSVs contain coverage depth of every bin across all samples.
# We only need each bin's depth in its OWN subject's two samples (HELIBA_ / HELIFU_).

batch_files <- c(
  "data/shotgun/alistipes_annotation/bins_alistipes_batch1.csv",
  "data/shotgun/alistipes_annotation/bins_alistipes_batch2.csv",
  "data/shotgun/alistipes_annotation/bins_alistipes_batch3.csv"
)

depth_raw <- map_dfr(batch_files, function(f) {
  read.csv(f, check.names = FALSE) %>%
    dplyr::select(bin, starts_with("Depth "))
}) %>%
  mutate(bin_name = sub("\\.fa$", "", bin)) %>%   # strip .fa suffix
  dplyr::select(bin_name, starts_with("Depth "))

depth_long <- depth_raw %>%
  pivot_longer(
    cols      = starts_with("Depth "),
    names_to  = "sampleID",
    values_to = "depth"
  ) %>%
  mutate(sampleID = sub("^Depth ", "", sampleID),
         depth    = replace_na(depth, 0))

## Keep only each bin's own subject's HELIBA / HELIFU samples
depth_own <- depth_long %>%
  inner_join(trans %>% dplyr::select(bin_name, locus_prefix, subject_id),
             by = "bin_name") %>%
  filter(
    sampleID == paste0("HELIBA_", subject_id) |
    sampleID == paste0("HELIFU_",  subject_id)
  )

cat("\nDepth summary per own-subject sample (depth > 0 = present):\n")
cat("  Present at baseline: ",
    sum(depth_own$depth > 0 & grepl("HELIBA", depth_own$sampleID)), "\n")
cat("  Present at follow-up:",
    sum(depth_own$depth > 0 & grepl("HELIFU",  depth_own$sampleID)), "\n")

#### 6. Join with clinical data ####

## Clinical long: columns include sampleID, ID, EthnicityTot, timepoint
clin <- readRDS(clin_file) %>%
  filter(EthnicityTot %in% c("Dutch", "South-Asian Surinamese")) %>%
  droplevels() %>%
  dplyr::select(sampleID, ID, EthnicityTot, timepoint)

## Build sample look-up: each subject → HELIBA_ and HELIFU_ sample IDs
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

## Join with clinical and attach depth (depth = 0 → absent)
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

cat("\nAll samples with annotation + clinical data (present AND absent):\n")
print(table(bin_clin$EthnicityTot, bin_clin$timepoint))

#### 7. Build analysis datasets (present samples only, zeros filled per feature) ####
# Only samples where depth > 0 (Alistipes present) are included.
# Within those, bins that LACK a given module/category get proportion = 0
# so that per-feature presence-absence is captured in the statistics.

present_bins <- bin_clin %>%
  filter(present) %>%
  dplyr::select(locus_prefix, sampleID, EthnicityTot, timepoint, depth)

## module_label lookup (one row per module, independent of locus_prefix)
module_labels <- kegg_long %>%
  distinct(module, module_label)

## Complete KEGG: every present bin × every module, fill missing with 0
kegg_per_bin <- expand.grid(
  locus_prefix = unique(present_bins$locus_prefix),
  module       = unique(kegg_long$module),
  stringsAsFactors = FALSE
) %>%
  left_join(kegg_long %>% dplyr::select(locus_prefix, module, proportion),
            by = c("locus_prefix", "module")) %>%
  mutate(proportion = replace_na(proportion, 0)) %>%
  left_join(module_labels, by = "module")

kegg_df <- present_bins %>%
  left_join(kegg_per_bin, by = "locus_prefix", relationship = "many-to-many")

## Sanity check: n present subjects per ethnicity × timepoint
cat("\nSanity check — present subjects per ethnicity × timepoint:\n")
print(
  present_bins %>%
    distinct(locus_prefix, EthnicityTot, timepoint) %>%
    count(EthnicityTot, timepoint)
)

#### 8. Statistical analysis ####

run_stats <- function(df, feature_col) {
  features <- unique(df[[feature_col]])

  results <- map_dfr(features, function(feat) {
    sub_df <- df %>% filter(.data[[feature_col]] == feat)

    ## --- Wilcoxon (abundance including zeros) at baseline ---
    ba_data  <- sub_df %>% filter(timepoint == "baseline")
    dutch_ba <- ba_data$proportion[ba_data$EthnicityTot == "Dutch"]
    sas_ba   <- ba_data$proportion[ba_data$EthnicityTot == "South-Asian Surinamese"]

    wx_ba <- tryCatch(
      wilcox.test(dutch_ba, sas_ba, exact = FALSE),
      error = function(e) list(statistic = NA, p.value = NA)
    )

    ## --- Wilcoxon at follow-up ---
    fu_data  <- sub_df %>% filter(timepoint == "follow-up")
    dutch_fu <- fu_data$proportion[fu_data$EthnicityTot == "Dutch"]
    sas_fu   <- fu_data$proportion[fu_data$EthnicityTot == "South-Asian Surinamese"]

    wx_fu <- tryCatch(
      wilcox.test(dutch_fu, sas_fu, exact = FALSE),
      error = function(e) list(statistic = NA, p.value = NA)
    )

    ## --- LMM: ethnicity * timepoint, (1|subject) ---
    lmm_res <- tryCatch({
      mod <- lmer(
        proportion ~ EthnicityTot * timepoint + (1 | locus_prefix),
        data    = sub_df,
        REML    = FALSE,
        control = lmerControl(optimizer = "bobyqa")
      )
      coef_tab <- summary(mod)$coefficients

      ## Main effect of ethnicity (SAS vs Dutch, positive = higher in SAS)
      eth_row <- grep("^EthnicityTot", rownames(coef_tab))
      eth_row <- eth_row[!grepl("timepoint", rownames(coef_tab)[eth_row])]
      if (length(eth_row) == 0) eth_row <- NA_integer_

      ## Interaction term (ethnicity × timepoint)
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
      feature                  = feat,
      wilcox_baseline_stat     = wx_ba$statistic,
      wilcox_baseline_p        = wx_ba$p.value,
      wilcox_fu_stat           = wx_fu$statistic,
      wilcox_fu_p              = wx_fu$p.value,
      lmm_eth_estimate         = lmm_res$lmm_eth_estimate,
      lmm_eth_se               = lmm_res$lmm_eth_se,
      lmm_eth_pval             = lmm_res$lmm_eth_pval,
      lmm_int_estimate         = lmm_res$lmm_int_estimate,
      lmm_int_se               = lmm_res$lmm_int_se,
      lmm_int_pval             = lmm_res$lmm_int_pval
    )
  })

  ## BH correction within each comparison type
  results <- results %>%
    mutate(
      wilcox_baseline_fdr  = p.adjust(wilcox_baseline_p,  method = "BH"),
      wilcox_fu_fdr        = p.adjust(wilcox_fu_p,        method = "BH"),
      lmm_eth_fdr          = p.adjust(lmm_eth_pval,       method = "BH"),
      lmm_int_fdr          = p.adjust(lmm_int_pval,       method = "BH")
    )

  return(results)
}

cat("\nRunning KEGG module statistics...\n")
kegg_results <- run_stats(kegg_df, "module")

#### 9. Save results tables ####
kegg_results <- kegg_results %>%
  left_join(module_labels %>% rename(feature = module), by = "feature") %>%
  dplyr::relocate(module_label, .after = feature)

write.csv(kegg_results, file.path(results_dir, "kegg_module_stats.csv"), row.names = FALSE)

cat("\nSignificant KEGG modules (any comparison, FDR < 0.05):\n")
sig_kegg <- kegg_results %>%
  filter(wilcox_baseline_fdr < 0.05 | wilcox_fu_fdr < 0.05 |
         lmm_eth_fdr < 0.05 | lmm_int_fdr < 0.05)
print(sig_kegg)

#### 10. Visualisation ####

## JCO colour palette: Dutch = first (blue), SAS = second (yellow)
jco_cols <- pal_jco()(2)
names(jco_cols) <- c("Dutch", "South-Asian Surinamese")

## --- 9a. KEGG module heatmap / dot plot ---
if (nrow(sig_kegg) > 0) {

  plot_kegg <- kegg_df %>%
    filter(module %in% sig_kegg$feature) %>%
    group_by(module_label, EthnicityTot, timepoint) %>%
    summarise(
      mean_prop = mean(proportion, na.rm = TRUE),
      se_prop   = sd(proportion, na.rm = TRUE) / sqrt(n()),
      .groups   = "drop"
    )

  ## Dot plot: x = ethnicity, y = module_label, size = mean proportion, colour = ethnicity
  p_kegg_dot <- ggplot(plot_kegg,
                       aes(x = EthnicityTot, y = module_label,
                           size   = mean_prop,
                           colour = EthnicityTot)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(name = "Mean proportion", range = c(2, 10)) +
    scale_colour_manual(values = jco_cols, name = "Ethnicity") +
    facet_wrap(~timepoint) +
    labs(
      title    = "Significantly different KEGG modules",
      subtitle = "Alistipes putredinis bins — Dutch vs South-Asian Surinamese",
      x = "", y = "",
      caption  = "FDR < 0.05 in at least one comparison"
    ) +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(
    file.path(results_dir, "kegg_modules_dotplot.pdf"),
    plot   = p_kegg_dot,
    width  = 10, height = max(4, nrow(sig_kegg) * 0.5 + 3)
  )

  ## Boxplot per significant module (faceted by timepoint)
  plot_kegg_box <- kegg_df %>%
    filter(module %in% sig_kegg$feature)

  p_kegg_box <- ggplot(plot_kegg_box,
                       aes(x = EthnicityTot, y = proportion,
                           fill = EthnicityTot)) +
    geom_boxplot(outlier.size = 0.8, width = 0.5) +
    scale_fill_manual(values = jco_cols, name = "Ethnicity") +
    scale_colour_manual(values = jco_cols, guide = "none") +
    facet_grid(module_label ~ timepoint, scales = "free_y") +
    labs(
      title    = "KEGG module proportions",
      subtitle = "Significant differences (FDR < 0.05), Alistipes putredinis",
      x = "", y = "Proportion of annotated genes",
      caption  = "Proportion = genes in module / total annotated genes per bin"
    ) +
    theme_Publication() +
    theme(axis.text.x  = element_text(angle = 30, hjust = 1),
          strip.text.y = element_text(angle = 0, size = rel(0.6)))

  ggsave(
    file.path(results_dir, "kegg_modules_boxplot.pdf"),
    plot   = p_kegg_box,
    width  = 7, height = max(4, nrow(sig_kegg) * 1.5 + 2)
  )

  plot_kegg_box <- plot_kegg_box |> mutate(module_label = fct_recode(module_label,
   "Ascorbate biosynthesis" = "M00129: Ascorbate biosynthesis, animals, glucose-1P => ascorbate",
   "Glucuronate pathway" = "M00014: Glucuronate pathway (uronate pathway)"
  ))

  ## Figure panel: baseline only, no legend
  pl_fig3_E <- ggplot(plot_kegg_box %>% filter(timepoint == "baseline"),
                      aes(x = EthnicityTot, y = proportion, fill = EthnicityTot)) +
    geom_boxplot(outlier.size = 0.8, width = 0.5) +
    stat_compare_means(comparisons = list(c("Dutch", "South-Asian Surinamese")),
                       method = "wilcox.test", label = "p.signif", tip.length = 0) +
    scale_fill_manual(values = jco_cols, guide = "none") +
    facet_wrap(~ module_label, scales = "free_y", labeller = label_wrap_gen(width = 45)) +
    scale_y_continuous(expand = expansion(add = c(0, 0.00025))) +
    labs(title = "KEGG pathways: Alistipes putredinis", x = "", y = "Proportion of annotated genes") +
    theme_Publication() +
    theme(strip.text  = element_text(size = rel(0.8)))

  cat("KEGG plots saved.\n")
} else {
  cat("No significantly different KEGG modules (FDR < 0.05) — no plot produced.\n")
}

## --- 9c. Ethnicity heatmap — LMM ethnicity main effect (SAS vs Dutch) ---
# Colour = LMM estimate for EthnicityTotSouth-Asian Surinamese (positive = higher in SAS).
# Ranked by absolute LMM ethnicity estimate. Stars from lmm_eth_fdr.

kegg_lmm_heat <- kegg_results %>%
  filter(!is.na(lmm_eth_estimate)) %>%
  slice_max(abs(lmm_eth_estimate), n = 25) %>%
  mutate(
    sig_label    = case_when(
      lmm_eth_fdr < 0.001 ~ "***",
      lmm_eth_fdr < 0.01  ~ "**",
      lmm_eth_fdr < 0.05  ~ "*",
      TRUE                 ~ ""
    ),
    module_label = factor(module_label,
                          levels = kegg_results %>%
                            filter(!is.na(lmm_eth_estimate)) %>%
                            slice_max(abs(lmm_eth_estimate), n = 25) %>%
                            arrange(lmm_eth_estimate) %>%
                            pull(module_label))
  )

lim <- max(abs(kegg_lmm_heat$lmm_eth_estimate), na.rm = TRUE)

p_heatmap <- ggplot(kegg_lmm_heat,
                    aes(x = 1, y = module_label, fill = lmm_eth_estimate)) +
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
    title    = "Top 25 KEGG modules — LMM ethnicity effect",
    subtitle = paste0("Alistipes putredinis bins  |  blue = higher in Dutch, ",
                      "yellow = higher in SAS\n* FDR<0.05  ** FDR<0.01  *** FDR<0.001"),
    x = "", y = "",
    caption  = "LMM: proportion ~ EthnicityTot * timepoint + (1|bin); estimate = SAS vs Dutch main effect"
  ) +
  theme_Publication() +
  theme(axis.text.y = element_text(size = rel(0.65)))

ggsave(
  file.path(results_dir, "kegg_modules_heatmap_lmm.pdf"),
  plot   = p_heatmap,
  width  = 7, height = 10
)

write.csv(kegg_results, file.path(results_dir, "kegg_module_summary.csv"),
          row.names = FALSE)

cat("LMM ethnicity heatmap saved.\n")

## --- 9d. Baseline heatmap — Wilcoxon ethnicity effect at baseline ---
# At baseline each bin is observed once, so no random effect is appropriate.
# Rank by Wilcoxon baseline statistic; colour by mean proportion difference
# (Dutch − SAS) at baseline to show direction.

kegg_baseline_means <- kegg_df %>%
  filter(timepoint == "baseline") %>%
  distinct(locus_prefix, module, module_label, EthnicityTot, proportion) %>%
  group_by(module, module_label, EthnicityTot) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = EthnicityTot, values_from = mean_prop, values_fill = 0) %>%
  mutate(diff_Dutch_SAS = Dutch - `South-Asian Surinamese`)

kegg_ba_heat <- kegg_results %>%
  filter(!is.na(wilcox_baseline_p)) %>%
  slice_max(abs(wilcox_baseline_stat), n = 25) %>%
  left_join(kegg_baseline_means %>% dplyr::select(module, diff_Dutch_SAS),
            by = c("feature" = "module")) %>%
  mutate(
    sig_label    = case_when(
      wilcox_baseline_fdr < 0.001 ~ "***",
      wilcox_baseline_fdr < 0.01  ~ "**",
      wilcox_baseline_fdr < 0.05  ~ "*",
      TRUE                        ~ ""
    ),
    module_label = factor(module_label,
                          levels = kegg_results %>%
                            filter(!is.na(wilcox_baseline_p)) %>%
                            slice_max(abs(wilcox_baseline_stat), n = 25) %>%
                            left_join(kegg_baseline_means %>%
                                        dplyr::select(module, diff_Dutch_SAS),
                                      by = c("feature" = "module")) %>%
                            arrange(diff_Dutch_SAS) %>%
                            pull(module_label))
  )

lim_ba <- max(abs(kegg_ba_heat$diff_Dutch_SAS), na.rm = TRUE)

p_heatmap_ba <- ggplot(kegg_ba_heat,
                       aes(x = 1, y = module_label, fill = diff_Dutch_SAS)) +
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
    title    = "Top 25 KEGG modules — ethnicity difference at baseline",
    subtitle = paste0("Alistipes putredinis bins  |  blue = higher in Dutch, ",
                      "yellow = higher in SAS\n* Wilcoxon FDR<0.05  ** FDR<0.01  *** FDR<0.001"),
    x = "", y = "",
    caption  = "Ranked by Wilcoxon W statistic; colour = mean proportion Dutch − SAS at baseline"
  ) +
  theme_Publication() +
  theme(axis.text.y = element_text(size = rel(0.65)))

ggsave(
  file.path(results_dir, "kegg_modules_heatmap_baseline.pdf"),
  plot   = p_heatmap_ba,
  width  = 7, height = 10
)

write.csv(
  kegg_results %>%
    dplyr::select(feature, module_label, wilcox_baseline_stat,
                  wilcox_baseline_p, wilcox_baseline_fdr) %>%
    left_join(kegg_baseline_means, by = c("feature" = "module")),
  file.path(results_dir, "kegg_module_wilcox_baseline.csv"),
  row.names = FALSE
)

cat("Baseline heatmap saved.\n")

cat("\nDone. Results written to:", results_dir, "\n")
