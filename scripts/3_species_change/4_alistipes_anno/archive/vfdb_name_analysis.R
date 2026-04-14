## VFDB per-VF-name granular analysis — Alistipes putredinis bins
## Secondary/exploratory analysis below VF category level.
## Run vfdb_comparison.R first (requires: vfdb_hits, present_bins,
##   total_hits_per_bin, vf_df, jco_cols, results_dir from that script).
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

source("scripts/3_species_change/4_alistipes_anno/utils.R")

MIN_PREV <- 5   # minimum bins a VF name must appear in

#### Prevalence filter ####
vfname_prevalence <- vfdb_hits %>%
  filter(!is.na(vf_name), locus_prefix %in% unique(present_bins$locus_prefix)) %>%
  distinct(locus_prefix, vf_name) %>%
  count(vf_name, name = "n_bins")

prevalent_vf_names <- vfname_prevalence %>%
  filter(n_bins >= MIN_PREV) %>%
  pull(vf_name)

cat("VF names with prevalence >=", MIN_PREV, "bins:", length(prevalent_vf_names), "\n")

#### Per-bin VF name proportions ####
vfname_long <- vfdb_hits %>%
  filter(!is.na(vf_name), locus_prefix %in% unique(present_bins$locus_prefix),
         vf_name %in% prevalent_vf_names) %>%
  count(locus_prefix, vf_name, name = "n_hits") %>%
  left_join(total_hits_per_bin, by = "locus_prefix") %>%
  mutate(proportion = n_hits / total_vf_hits)

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

#### Statistics ####
cat("Running VF name statistics (", length(prevalent_vf_names), " names)...\n")
vfname_results <- run_stats(vfname_df, "vf_name") %>%
  rename(vf_name = feature) %>%
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

write.csv(vfname_results,
          file.path(results_dir, "vfdb_name_stats.csv"),
          row.names = FALSE)

sig_vfname <- vfname_results %>%
  filter(wilcox_baseline_fdr < 0.05 | wilcox_fu_fdr < 0.05 |
           lmm_eth_fdr < 0.05 | lmm_int_fdr < 0.05)

cat("\nSignificant VF names (FDR < 0.05):\n")
print(sig_vfname %>% dplyr::select(vf_name, vf_category, lmm_eth_estimate,
                                    lmm_eth_fdr, wilcox_baseline_fdr, wilcox_fu_fdr))

#### LMM ethnicity heatmap (top 25 VF names) ####
top_n <- 25
jco_cols <- jco_palette()

vfname_lmm_heat <- vfname_results %>%
  filter(!is.na(lmm_eth_estimate)) %>%
  slice_max(abs(lmm_eth_estimate), n = top_n, with_ties = FALSE) %>%
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
                       slice_max(abs(lmm_eth_estimate), n = top_n, with_ties = FALSE) %>%
                       arrange(lmm_eth_estimate) %>%
                       pull(vf_name))
  )

lim <- max(abs(vfname_lmm_heat$lmm_eth_estimate), na.rm = TRUE)

if (nrow(vfname_lmm_heat) > 0 && !is.infinite(lim)) {
  p_vfname_heat_lmm <- ggplot(vfname_lmm_heat,
                               aes(x = 1, y = vf_name, fill = lmm_eth_estimate)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_label), size = 3, vjust = 0.5) +
    scale_fill_gradient2(low = jco_cols[["Dutch"]], mid = "white",
                         high = jco_cols[["South-Asian Surinamese"]],
                         midpoint = 0, limits = c(-lim, lim),
                         name = "LMM estimate\n(SAS vs Dutch)") +
    scale_x_continuous(breaks = NULL) +
    labs(title    = paste0("Top ", top_n, " VF names — LMM ethnicity effect"),
         subtitle = paste0("blue = higher in Dutch, yellow = higher in SAS\n",
                           "VF names in >= ", MIN_PREV, " bins"),
         x = "", y = "") +
    theme_Publication() +
    theme(axis.text.y = element_text(size = rel(0.7)))

  ggsave(file.path(results_dir, "vfdb_name_heatmap_lmm.pdf"),
         plot = p_vfname_heat_lmm, width = 8, height = 10)
}

#### Baseline heatmap (top 25 VF names by Wilcoxon W) ####
vfname_baseline_means <- vfname_df %>%
  filter(timepoint == "baseline") %>%
  distinct(locus_prefix, vf_name, EthnicityTot, proportion) %>%
  group_by(vf_name, EthnicityTot) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = EthnicityTot, values_from = mean_prop, values_fill = 0) %>%
  mutate(diff_Dutch_SAS = Dutch - `South-Asian Surinamese`)

vfname_ba_heat <- vfname_results %>%
  filter(!is.na(wilcox_baseline_p)) %>%
  slice_max(abs(wilcox_baseline_stat), n = top_n, with_ties = FALSE) %>%
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
                       slice_max(abs(wilcox_baseline_stat), n = top_n, with_ties = FALSE) %>%
                       left_join(vfname_baseline_means %>%
                                   dplyr::select(vf_name, diff_Dutch_SAS),
                                 by = "vf_name") %>%
                       arrange(diff_Dutch_SAS) %>%
                       pull(vf_name))
  )

lim_ba <- max(abs(vfname_ba_heat$diff_Dutch_SAS), na.rm = TRUE)

if (nrow(vfname_ba_heat) > 0 && !is.infinite(lim_ba)) {
  p_vfname_heat_ba <- ggplot(vfname_ba_heat,
                              aes(x = 1, y = vf_name, fill = diff_Dutch_SAS)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_label), size = 3, vjust = 0.5) +
    scale_fill_gradient2(low = jco_cols[["South-Asian Surinamese"]], mid = "white",
                         high = jco_cols[["Dutch"]], midpoint = 0,
                         limits = c(-lim_ba, lim_ba),
                         name = "Mean proportion\nDutch − SAS (baseline)") +
    scale_x_continuous(breaks = NULL) +
    labs(title    = paste0("Top ", top_n, " VF names — ethnicity difference at baseline"),
         subtitle = paste0("blue = higher in Dutch, yellow = higher in SAS\n",
                           "VF names in >= ", MIN_PREV, " bins"),
         x = "", y = "") +
    theme_Publication() +
    theme(axis.text.y = element_text(size = rel(0.7)))

  ggsave(file.path(results_dir, "vfdb_name_heatmap_baseline.pdf"),
         plot = p_vfname_heat_ba, width = 8, height = 10)
}

cat("\nVF name analysis done. Results written to:", results_dir, "\n")
