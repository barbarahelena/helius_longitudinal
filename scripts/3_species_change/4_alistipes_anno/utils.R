## Shared utilities — Alistipes/Odoribacter annotation analyses
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggsci)
library(ggthemes)
library(grid)
library(ggpubr)

#### Presence-calling threshold ####
# Minimum independently re-mapped depth at a timepoint required to call a MAG
# "present" there. Depths in (0, 1x) are trace-level signal (single/few reads,
# plausibly cross-mapping from a closely related co-occurring strain) that is
# not reliably distinguishable from background — see
# results/3_species_change/4_alistipes_anno/alistipes_depth_wide.csv, which
# shows a clean bimodal split (<1x vs >8x) with no bins in between once one
# timepoint is genuinely colonised.
PRESENCE_THRESHOLD <- 1

#### MAG quality as a covariate ####
# Gene-content measures depend on assembly quality: an incomplete MAG makes
# genes look absent, and contamination from a co-occurring organism adds genes
# that are not the genome's own. Quality is not evenly distributed — Dutch MAGs
# are less complete (median 92.7% vs 95.5%, p = 6.6e-4) and more contaminated
# (7.20% vs 4.34%, p = 1.8e-7) than South-Asian Surinamese ones, so for the
# ethnicity comparison it is a genuine confounder. It does not differ by clade
# (p = 0.67 and 0.74).
#
# Rather than discard the ~72% of MAGs that fall below a high-quality cutoff,
# completeness and contamination are carried as covariates in the models. Use
# this to attach them to a per-bin table keyed on locus_prefix.
add_bin_quality <- function(df, trans, batch_files) {
  q <- map_dfr(batch_files, function(f) {
    read.csv(f, check.names = FALSE) %>%
      dplyr::select(bin, Completeness, Contamination)
  }) %>%
    mutate(bin_name = sub("\\.fa$", "", bin)) %>%
    inner_join(trans %>% dplyr::select(bin_name, locus_prefix), by = "bin_name") %>%
    dplyr::select(locus_prefix, Completeness, Contamination)
  if (any(duplicated(q$locus_prefix)))
    stop("duplicate locus prefixes in the batch quality tables")
  df %>% left_join(q, by = "locus_prefix")
}

#### Gene-content denominator ####
# Proportions are expressed per predicted protein-coding gene, rather than per
# gene that eggNOG could assign a KEGG module or COG category (only ~76% of
# them, excluding a quarter of each genome for reasons unrelated to virulence
# factors).
#
# Two sources, in order of preference:
#  1. Bakta's own CDS count, if bakta_cds_counts.tsv is present. This is exactly
#     the set of proteins DIAMOND searched against VFDB, so numerator and
#     denominator come from the same gene prediction. Produced on Snellius by
#     alistipes_bins_annotation/extract_bakta_cds_counts.sh.
#  2. Otherwise CheckM2's Total_Coding_Sequences from the batch tables, which is
#     CheckM2's internal prodigal call. Close to Bakta's but not the same caller.
# Whichever is used is printed, so it is never ambiguous which denominator a run
# was built on.
BAKTA_CDS_FILE <- "data/shotgun/alistipes_annotation/bakta_cds_counts.tsv"

load_gene_counts <- function(trans, batch_files, file = BAKTA_CDS_FILE) {
  if (file.exists(file)) {
    cds <- read.delim(file, stringsAsFactors = FALSE)
    if (!all(c("locus_prefix", "cds_count") %in% names(cds)))
      stop(file, " must have columns 'locus_prefix' and 'cds_count'")
    if (any(duplicated(cds$locus_prefix)))
      stop("duplicate locus prefixes in ", file)
    cat("Gene-content denominator: Bakta CDS counts from", file, "\n")
    return(cds %>%
             dplyr::select(locus_prefix, total_cds = cds_count) %>%
             dplyr::filter(!is.na(total_cds), total_cds > 0))
  }
  cat("Gene-content denominator: CheckM2 (prodigal) Total_Coding_Sequences.\n",
      "  For Bakta's own counts, run extract_bakta_cds_counts.sh on Snellius\n",
      "  and copy the result to ", file, "\n", sep = "")
  map_dfr(batch_files, function(f) {
    read.csv(f, check.names = FALSE) %>%
      dplyr::select(bin, Total_Coding_Sequences)
  }) %>%
    mutate(bin_name = sub("\\.fa$", "", bin)) %>%
    inner_join(trans %>% dplyr::select(bin_name, locus_prefix), by = "bin_name") %>%
    dplyr::select(locus_prefix, total_cds = Total_Coding_Sequences) %>%
    dplyr::filter(!is.na(total_cds), total_cds > 0)
}

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

#### Colours ####
jco_palette <- function() {
  cols <- pal_jco()(2)
  names(cols) <- c("Dutch", "South-Asian Surinamese")
  cols
}

#### Data loading ####
# Loads bin depths from batch CSVs, joins with clinical data, returns bin_clin.
# bin_clin has one row per bin × sample (both timepoints), with columns:
#   locus_prefix, sampleID, EthnicityTot, timepoint, depth, present
load_bin_clin <- function(trans, batch_files, clin_file) {
  depth_raw <- map_dfr(batch_files, function(f) {
    read.csv(f, check.names = FALSE) %>%
      dplyr::select(bin, starts_with("Depth "))
  }) %>%
    mutate(bin_name = sub("\\.fa$", "", bin)) %>%
    dplyr::select(bin_name, starts_with("Depth "))

  depth_long <- depth_raw %>%
    pivot_longer(cols = starts_with("Depth "),
                 names_to = "sampleID", values_to = "depth") %>%
    mutate(sampleID = sub("^Depth ", "", sampleID),
           depth    = replace_na(depth, 0))

  depth_own <- depth_long %>%
    inner_join(trans %>% dplyr::select(bin_name, locus_prefix, subject_id),
               by = "bin_name") %>%
    filter(sampleID == paste0("HELIBA_", subject_id) |
           sampleID == paste0("HELIFU_",  subject_id))

  clin <- readRDS(clin_file) %>%
    filter(EthnicityTot %in% c("Dutch", "South-Asian Surinamese")) %>%
    droplevels() %>%
    dplyr::select(sampleID, ID, EthnicityTot, timepoint)

  sample_map <- trans %>%
    dplyr::select(locus_prefix, subject_id) %>%
    mutate(sampleID_BA = paste0("HELIBA_", subject_id),
           sampleID_FU = paste0("HELIFU_", subject_id)) %>%
    pivot_longer(cols = c(sampleID_BA, sampleID_FU),
                 names_to = NULL, values_to = "sampleID") %>%
    mutate(subject_id = as.character(subject_id))

  sample_map %>%
    inner_join(clin, by = "sampleID") %>%
    left_join(depth_own %>% dplyr::select(locus_prefix, sampleID, depth),
              by = c("locus_prefix", "sampleID")) %>%
    mutate(
      depth        = replace_na(depth, 0),
      present      = depth >= PRESENCE_THRESHOLD,
      timepoint    = factor(timepoint, levels = c("baseline", "follow-up")),
      EthnicityTot = factor(EthnicityTot,
                            levels = c("Dutch", "South-Asian Surinamese"))
    )
}

#### Statistics ####
# Bin-level comparison of a gene-content measure between ethnic groups.
#
# Gene content is a property of the assembled genome, and each participant has
# exactly one MAG, so it has no baseline and follow-up value: expanding to
# bin x timepoint rows duplicates each bin's value and leaves zero within-bin
# variance. A mixed model with timepoint and (1 | locus_prefix) is then
# degenerate — the timepoint terms are exactly zero and the random intercept
# absorbs the rest, leaving no residual degrees of freedom. So this takes one
# row per bin and fits a plain linear model.
#
# Completeness and contamination are adjusted for because both bias gene-content
# measures and both differ between ethnic groups (see add_bin_quality).
# The unadjusted Wilcoxon is kept alongside as a descriptive comparison.
# Expects one row per bin per feature, with columns: locus_prefix, EthnicityTot,
# proportion, Completeness, Contamination. Returns one row per feature.
run_bin_stats <- function(df, feature_col) {
  stopifnot(!any(duplicated(df[c("locus_prefix", feature_col)])))
  features <- unique(df[[feature_col]])

  results <- map_dfr(features, function(feat) {
    sub_df <- df %>% filter(.data[[feature_col]] == feat)
    dutch <- sub_df$proportion[sub_df$EthnicityTot == "Dutch"]
    sas   <- sub_df$proportion[sub_df$EthnicityTot == "South-Asian Surinamese"]

    wx <- tryCatch(wilcox.test(dutch, sas, exact = FALSE),
                   error = function(e) list(statistic = NA_real_, p.value = NA_real_))

    fit <- tryCatch({
      mod <- lm(proportion ~ EthnicityTot + Completeness + Contamination, data = sub_df)
      ct  <- summary(mod)$coefficients
      row <- grep("^EthnicityTot", rownames(ct))
      if (length(row) != 1) stop("expected one ethnicity term")
      list(estimate = ct[row, "Estimate"],
           se       = ct[row, "Std. Error"],
           pval     = ct[row, "Pr(>|t|)"])
    }, error = function(e) list(estimate = NA_real_, se = NA_real_, pval = NA_real_))

    tibble(
      feature        = feat,
      n_dutch        = length(dutch),
      n_sas          = length(sas),
      median_dutch   = median(dutch),
      median_sas     = median(sas),
      wilcox_stat    = unname(wx$statistic),
      wilcox_p       = wx$p.value,
      adj_estimate   = fit$estimate,   # SAS vs Dutch, adjusted for MAG quality
      adj_se         = fit$se,
      adj_p          = fit$pval
    )
  })

  results %>%
    mutate(wilcox_fdr = p.adjust(wilcox_p, method = "BH"),
           adj_fdr    = p.adjust(adj_p,    method = "BH"))
}

# Runs Wilcoxon (baseline + follow-up) and optionally LMM (ethnicity main effect +
# ethnicity × timepoint interaction) for each level of feature_col.
# Set lmm = FALSE to skip the LMM (faster, Wilcoxon only).
# Returns one row per feature with raw p-values and BH-adjusted FDR.
# The result column is always named "feature"; rename downstream if needed.
run_stats <- function(df, feature_col, lmm = TRUE) {
  features <- unique(df[[feature_col]])

  results <- map_dfr(features, function(feat) {
    sub_df <- df %>% filter(.data[[feature_col]] == feat)

    ba_data  <- sub_df %>% filter(timepoint == "baseline")
    dutch_ba <- ba_data$proportion[ba_data$EthnicityTot == "Dutch"]
    sas_ba   <- ba_data$proportion[ba_data$EthnicityTot == "South-Asian Surinamese"]
    wx_ba <- tryCatch(wilcox.test(dutch_ba, sas_ba, exact = FALSE),
                      error = function(e) list(statistic = NA, p.value = NA))

    fu_data  <- sub_df %>% filter(timepoint == "follow-up")
    dutch_fu <- fu_data$proportion[fu_data$EthnicityTot == "Dutch"]
    sas_fu   <- fu_data$proportion[fu_data$EthnicityTot == "South-Asian Surinamese"]
    wx_fu <- tryCatch(wilcox.test(dutch_fu, sas_fu, exact = FALSE),
                      error = function(e) list(statistic = NA, p.value = NA))

    if (lmm) {
      lmm_res <- tryCatch({
        # Completeness and contamination are adjusted for rather than filtered
        # on: they differ systematically between ethnic groups and both bias
        # gene-content measures. Added only when present, so callers that do
        # not supply them still get the unadjusted model.
        fml <- if (all(c("Completeness", "Contamination") %in% names(sub_df))) {
          proportion ~ EthnicityTot * timepoint + Completeness + Contamination +
            (1 | locus_prefix)
        } else {
          proportion ~ EthnicityTot * timepoint + (1 | locus_prefix)
        }
        mod <- lmer(
          fml,
          data = sub_df, REML = FALSE,
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
    }

    row <- tibble(
      feature              = feat,
      n_dutch_baseline     = length(dutch_ba),
      n_sas_baseline       = length(sas_ba),
      wilcox_baseline_stat = wx_ba$statistic,
      wilcox_baseline_p    = wx_ba$p.value,
      wilcox_fu_stat       = wx_fu$statistic,
      wilcox_fu_p          = wx_fu$p.value
    )

    if (lmm) {
      row <- row %>% mutate(
        lmm_eth_estimate = lmm_res$lmm_eth_estimate,
        lmm_eth_se       = lmm_res$lmm_eth_se,
        lmm_eth_pval     = lmm_res$lmm_eth_pval,
        lmm_int_estimate = lmm_res$lmm_int_estimate,
        lmm_int_se       = lmm_res$lmm_int_se,
        lmm_int_pval     = lmm_res$lmm_int_pval
      )
    }
    row
  })

  fdr_cols <- results %>%
    mutate(
      wilcox_baseline_fdr = p.adjust(wilcox_baseline_p, method = "BH"),
      wilcox_fu_fdr       = p.adjust(wilcox_fu_p,       method = "BH")
    )

  if (lmm) {
    fdr_cols <- fdr_cols %>%
      mutate(
        lmm_eth_fdr = p.adjust(lmm_eth_pval, method = "BH"),
        lmm_int_fdr = p.adjust(lmm_int_pval, method = "BH")
      )
  }
  fdr_cols
}
