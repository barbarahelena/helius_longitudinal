## Shared utilities — Alistipes/Odoribacter annotation analyses
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

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
      present      = depth > 0,
      timepoint    = factor(timepoint, levels = c("baseline", "follow-up")),
      EthnicityTot = factor(EthnicityTot,
                            levels = c("Dutch", "South-Asian Surinamese"))
    )
}

#### Statistics ####
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
        mod <- lmer(
          proportion ~ EthnicityTot * timepoint + (1 | locus_prefix),
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
