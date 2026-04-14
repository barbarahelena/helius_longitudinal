## Exploratory / QC plots — Alistipes putredinis annotation
## Run eggnog_comparison.R and vfdb_comparison.R first to populate the environment.
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

source("scripts/3_species_change/4_alistipes_anno/utils.R")

results_dir <- "results/3_species_change/4_alistipes_anno"
jco_cols    <- jco_palette()

#### KEGG — full boxplot (baseline + follow-up) ####
p_kegg_box <- ggplot(
    kegg_df %>% filter(module %in% sig_kegg$module),
    aes(x = EthnicityTot, y = proportion, fill = EthnicityTot)
  ) +
  geom_boxplot(outlier.size = 0.8, width = 0.5) +
  scale_fill_manual(values = jco_cols, name = "Ethnicity") +
  facet_grid(module_label ~ timepoint, scales = "free_y") +
  labs(title    = "KEGG module proportions",
       subtitle = "Significant differences (FDR < 0.05), Alistipes putredinis",
       x = "", y = "Proportion of annotated genes",
       caption  = "Proportion = genes in module / total annotated genes per bin") +
  theme_Publication() +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1),
        strip.text.y = element_text(angle = 0, size = rel(0.6)))

ggsave(file.path(results_dir, "qc_kegg_boxplot_all_timepoints.pdf"),
       plot = p_kegg_box,
       width = 7, height = max(4, nrow(sig_kegg) * 1.5 + 2))

#### KEGG — dot plot (significant modules, both timepoints) ####
plot_kegg_dot <- kegg_df %>%
  filter(module %in% sig_kegg$module) %>%
  group_by(module_label, EthnicityTot, timepoint) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop")

p_kegg_dot <- ggplot(
    plot_kegg_dot,
    aes(x = EthnicityTot, y = module_label, size = mean_prop, colour = EthnicityTot)
  ) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(name = "Mean proportion", range = c(2, 10)) +
  scale_colour_manual(values = jco_cols, name = "Ethnicity") +
  facet_wrap(~ timepoint) +
  labs(title    = "Significantly different KEGG modules",
       subtitle = "Alistipes putredinis — Dutch vs South-Asian Surinamese",
       x = "", y = "", caption = "FDR < 0.05 in at least one comparison") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(results_dir, "qc_kegg_dotplot.pdf"),
       plot = p_kegg_dot,
       width = 10, height = max(4, nrow(sig_kegg) * 0.5 + 3))

#### KEGG — LMM ethnicity heatmap (top 25 by absolute effect) ####
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

p_kegg_heat_lmm <- ggplot(kegg_lmm_heat,
                           aes(x = 1, y = module_label, fill = lmm_eth_estimate)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3, vjust = 0.5) +
  scale_fill_gradient2(low = jco_cols[["Dutch"]], mid = "white",
                       high = jco_cols[["South-Asian Surinamese"]],
                       midpoint = 0, limits = c(-lim, lim),
                       name = "LMM estimate\n(SAS vs Dutch)") +
  scale_x_continuous(breaks = NULL) +
  labs(title    = "Top 25 KEGG modules — LMM ethnicity effect",
       subtitle = "blue = higher in Dutch, yellow = higher in SAS",
       x = "", y = "") +
  theme_Publication() +
  theme(axis.text.y = element_text(size = rel(0.65)))

ggsave(file.path(results_dir, "qc_kegg_heatmap_lmm.pdf"),
       plot = p_kegg_heat_lmm, width = 7, height = 10)

#### KEGG — baseline heatmap (top 25 by Wilcoxon W) ####
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
            by = "module") %>%
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
                                      by = "module") %>%
                            arrange(diff_Dutch_SAS) %>%
                            pull(module_label))
  )

lim_ba <- max(abs(kegg_ba_heat$diff_Dutch_SAS), na.rm = TRUE)

p_kegg_heat_ba <- ggplot(kegg_ba_heat,
                         aes(x = 1, y = module_label, fill = diff_Dutch_SAS)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3, vjust = 0.5) +
  scale_fill_gradient2(low = jco_cols[["South-Asian Surinamese"]], mid = "white",
                       high = jco_cols[["Dutch"]], midpoint = 0,
                       limits = c(-lim_ba, lim_ba),
                       name = "Mean proportion\nDutch − SAS (baseline)") +
  scale_x_continuous(breaks = NULL) +
  labs(title    = "Top 25 KEGG modules — ethnicity difference at baseline",
       subtitle = "blue = higher in Dutch, yellow = higher in SAS",
       x = "", y = "") +
  theme_Publication() +
  theme(axis.text.y = element_text(size = rel(0.65)))

ggsave(file.path(results_dir, "qc_kegg_heatmap_baseline.pdf"),
       plot = p_kegg_heat_ba, width = 7, height = 10)

#### VFDB — full boxplot (baseline + follow-up) ####
p_vf_box <- ggplot(
    vf_df %>% filter(vf_category %in% sig_vf$vf_category),
    aes(x = EthnicityTot, y = proportion, fill = EthnicityTot)
  ) +
  geom_boxplot(outlier.size = 0.8, width = 0.5) +
  scale_fill_manual(values = jco_cols, name = "Ethnicity") +
  facet_grid(vf_category ~ timepoint, scales = "free_y") +
  labs(title    = "VF category proportions",
       subtitle = "Significant differences (FDR < 0.05), Alistipes putredinis",
       x = "", y = "Proportion of VFDB hits",
       caption  = "Proportion = hits in VF category / total VFDB hits per bin") +
  theme_Publication() +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1),
        strip.text.y = element_text(angle = 0, size = rel(0.7)))

ggsave(file.path(results_dir, "qc_vfdb_boxplot_all_timepoints.pdf"),
       plot = p_vf_box,
       width = 8, height = max(4, nrow(sig_vf) * 1.5 + 2))

#### VFDB — LMM ethnicity heatmap (top 20) ####
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

lim_vf <- max(abs(vf_lmm_heat$lmm_eth_estimate), na.rm = TRUE)

p_vf_heat_lmm <- ggplot(vf_lmm_heat,
                         aes(x = 1, y = vf_category, fill = lmm_eth_estimate)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3, vjust = 0.5) +
  scale_fill_gradient2(low = jco_cols[["Dutch"]], mid = "white",
                       high = jco_cols[["South-Asian Surinamese"]],
                       midpoint = 0, limits = c(-lim_vf, lim_vf),
                       name = "LMM estimate\n(SAS vs Dutch)") +
  scale_x_continuous(breaks = NULL) +
  labs(title    = "Top 20 VF categories — LMM ethnicity effect",
       subtitle = "blue = higher in Dutch, yellow = higher in SAS",
       x = "", y = "") +
  theme_Publication() +
  theme(axis.text.y = element_text(size = rel(0.75)))

ggsave(file.path(results_dir, "qc_vfdb_heatmap_lmm.pdf"),
       plot = p_vf_heat_lmm, width = 7, height = 9)

#### VFDB — baseline heatmap (top 20 by Wilcoxon W) ####
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

lim_ba_vf <- max(abs(vf_ba_heat$diff_Dutch_SAS), na.rm = TRUE)

p_vf_heat_ba <- ggplot(vf_ba_heat,
                       aes(x = 1, y = vf_category, fill = diff_Dutch_SAS)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3, vjust = 0.5) +
  scale_fill_gradient2(low = jco_cols[["South-Asian Surinamese"]], mid = "white",
                       high = jco_cols[["Dutch"]], midpoint = 0,
                       limits = c(-lim_ba_vf, lim_ba_vf),
                       name = "Mean proportion\nDutch − SAS (baseline)") +
  scale_x_continuous(breaks = NULL) +
  labs(title    = "Top 20 VF categories — ethnicity difference at baseline",
       subtitle = "blue = higher in Dutch, yellow = higher in SAS",
       x = "", y = "") +
  theme_Publication() +
  theme(axis.text.y = element_text(size = rel(0.75)))

ggsave(file.path(results_dir, "qc_vfdb_heatmap_baseline.pdf"),
       plot = p_vf_heat_ba, width = 7, height = 9)

cat("\nQC plots saved to:", results_dir, "\n")
