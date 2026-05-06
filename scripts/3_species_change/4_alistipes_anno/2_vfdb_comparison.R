## VFDB virulence factor annotation — Alistipes putredinis bins
## Dutch vs South-Asian Surinamese, baseline and follow-up
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

source("scripts/3_species_change/4_alistipes_anno/utils.R")

#### Paths ####
vfdb_file   <- "data/shotgun/alistipes_annotation/all_vfdb_results.txt"
fasta_file  <- "data/shotgun/alistipes_annotation/VFDB_setB_pro.fas"
trans_file  <- "data/shotgun/alistipes_annotation/bin_translation_table.tsv"
clin_file   <- "data/clinicaldata/clinicaldata_long.RDS"
batch_files <- c(
  "data/shotgun/alistipes_annotation/bins_alistipes_batch1.csv",
  "data/shotgun/alistipes_annotation/bins_alistipes_batch2.csv",
  "data/shotgun/alistipes_annotation/bins_alistipes_batch3.csv"
)
results_dir <- "results/3_species_change/4_alistipes_anno"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

#### 1. Load bin translation table ####
trans <- read.delim(trans_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE) %>%
  rename(locus_prefix = locus_tag_prefix)

cat("Bins in translation table:", nrow(trans), "\n")

#### 2. Parse VFDB FASTA headers → annotation lookup ####
# Header format:
# >VFG037170(gb|WP_...) (gene_symbol) protein [VF_name (VFxxxxx) - VF_category (VFCxxxxx)] [Organism]
raw_headers <- readLines(fasta_file) %>%
  keep(startsWith, ">") %>%
  sub("^>", "", .)

vfdb_anno <- tibble(header = raw_headers) %>%
  mutate(
    vfg_id      = str_extract(header, "^VFG\\d+"),
    gene_symbol = str_match(header, "\\)\\s+\\((\\w[^)]+)\\)")[, 2],
    vf_name     = str_match(header, "\\[([^\\[\\]]+)\\s+\\(VF\\d+\\)")[, 2],
    vf_id       = str_extract(header, "VF\\d+"),
    vf_category = str_match(header, "-\\s+([^\\[\\]]+?)\\s+\\(VFC\\d+\\)")[, 2],
    vfc_id      = str_extract(header, "VFC\\d+"),
    organism    = str_match(header, "\\[([^\\[\\]]+)\\]\\s*$")[, 2]
  ) %>%
  dplyr::select(-header) %>%
  filter(!is.na(vfg_id)) %>%
  distinct(vfg_id, .keep_all = TRUE)

cat("VFG IDs in FASTA:", nrow(vfdb_anno), "\n")

#### 3. Parse and filter DIAMOND output ####
col_names <- c("query", "subject", "pident", "length", "mismatch",
               "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

MIN_PIDENT   <- 30
MIN_BITSCORE <- 50

vfdb_hits <- read.table(
  vfdb_file,
  header = FALSE, sep = "\t", comment.char = "#",
  quote = "", col.names = col_names, stringsAsFactors = FALSE, fill = TRUE
) %>%
  filter(pident >= MIN_PIDENT, bitscore >= MIN_BITSCORE) %>%
  mutate(locus_prefix = sub("_.*", "", query),
         vfg_id       = str_extract(subject, "^VFG\\d+")) %>%
  left_join(vfdb_anno, by = "vfg_id")

cat("Hits after filtering (pident >=", MIN_PIDENT, ", bitscore >=", MIN_BITSCORE, "):",
    nrow(vfdb_hits), "\n")

#### 4. Per-bin VF category proportions ####
total_hits_per_bin <- vfdb_hits %>%
  filter(!is.na(vf_category)) %>%
  count(locus_prefix, name = "total_vf_hits")

vf_long <- vfdb_hits %>%
  filter(!is.na(vf_category)) %>%
  count(locus_prefix, vf_category, name = "n_hits") %>%
  left_join(total_hits_per_bin, by = "locus_prefix") %>%
  mutate(proportion = n_hits / total_vf_hits)

cat("Unique VF categories:", n_distinct(vf_long$vf_category), "\n")

#### 5. Join with clinical data ####
bin_clin <- load_bin_clin(trans, batch_files, clin_file)

baseline_bins <- bin_clin %>%
  filter(present, timepoint == "baseline",
         locus_prefix %in% unique(vf_long$locus_prefix)) %>%
  distinct(locus_prefix)

present_bins <- bin_clin %>%
  filter(present, locus_prefix %in% baseline_bins$locus_prefix) %>%
  dplyr::select(locus_prefix, sampleID, EthnicityTot, timepoint, depth)

cat("\nPresent bins with VF hits per ethnicity × timepoint:\n")
print(present_bins %>%
        distinct(locus_prefix, EthnicityTot, timepoint) %>%
        count(EthnicityTot, timepoint))

#### 6. Build analysis dataset ####
# Complete grid: every present bin × every VF category; missing proportion = 0
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

#### 7. Prevalence filter ####
MIN_PREV <- 0.10

prevalent_vf <- vf_df %>%
  filter(timepoint == "baseline") %>%
  group_by(vf_category) %>%
  summarise(prev = mean(proportion > 0), .groups = "drop") %>%
  filter(prev >= MIN_PREV) %>%
  pull(vf_category)

cat("\nVF categories passing prevalence filter (>=", MIN_PREV * 100, "% of baseline bins):",
    length(prevalent_vf), "of", n_distinct(vf_df$vf_category), "\n")

vf_df_filt <- vf_df %>% filter(vf_category %in% prevalent_vf)

#### 8. Statistics ####
cat("\nRunning VF category statistics...\n")
vf_results <- run_stats(vf_df_filt, "vf_category") %>%
  rename(vf_category = feature)

write.csv(vf_results, file.path(results_dir, "vfdb_category_stats.csv"),
          row.names = FALSE)

## Top VF names per category (interpretation aid)
vf_name_summary <- vfdb_hits %>%
  filter(!is.na(vf_category), !is.na(vf_name)) %>%
  count(vf_category, vf_name, name = "n_hits") %>%
  group_by(vf_category) %>%
  slice_max(n_hits, n = 5) %>%
  arrange(vf_category, desc(n_hits))

write.csv(vf_name_summary,
          file.path(results_dir, "vfdb_top_vf_names_per_category.csv"),
          row.names = FALSE)

sig_vf <- vf_results %>%
  filter(lmm_eth_fdr < 0.05 | lmm_int_fdr < 0.05 |
           wilcox_baseline_fdr < 0.05 | wilcox_fu_fdr < 0.05)

cat("\nTesting", n_distinct(vf_df_filt$vf_category), "VF categories (N Dutch:",
    unique(vf_results$n_dutch_baseline), "| N SAS:",
    unique(vf_results$n_sas_baseline), ")\n")
cat("Significant VF categories (FDR < 0.05 in any test):", nrow(sig_vf), "\n")
print(sig_vf %>% dplyr::select(vf_category, n_dutch_baseline, n_sas_baseline,
                                lmm_eth_estimate, lmm_eth_fdr,
                                lmm_int_estimate, lmm_int_fdr,
                                wilcox_baseline_fdr, wilcox_fu_fdr))

#### 9. Figure panel pl_fig3_F ####
jco_cols <- jco_palette()

plot_vf_box <- vf_df_filt %>% filter(vf_category %in% sig_vf$vf_category)

pl_fig3_F <- ggplot(
    plot_vf_box %>% filter(timepoint == "baseline"),
    aes(x = EthnicityTot, y = proportion, fill = EthnicityTot)
  ) +
  geom_boxplot(outlier.size = 0.8, width = 0.5) +
  stat_compare_means(comparisons = list(c("Dutch", "South-Asian Surinamese")),
                     method = "wilcox.test", label = "p.format", tip.length = 0) +
  scale_fill_manual(values = jco_cols, guide = "none") +
  facet_wrap(~ vf_category, scales = "free_y") +
  scale_y_continuous(expand = expansion(add = c(0, 0.010))) +
  labs(title = "VFDB: Alistipes putredinis",
       x = "", y = "Proportion of VFDB hits") +
  theme_Publication() +
  theme(strip.text = element_text(size = rel(0.8)))

ggsave(file.path(results_dir, "vfdb_sig_baseline.pdf"),
       plot = pl_fig3_F,
       width = 8, height = max(5, ceiling(nrow(sig_vf) / 3) * 2.8 + 1))

## Full VFDB boxplot (both timepoints)
p_vf_box <- ggplot(
    vf_df_filt %>% filter(vf_category %in% sig_vf$vf_category),
    aes(x = EthnicityTot, y = proportion, fill = EthnicityTot)
  ) +
  geom_boxplot(outlier.size = 0.8, width = 0.5) +
  stat_compare_means(comparisons = list(c("Dutch", "South-Asian Surinamese")),
                     method = "wilcox.test", label = "p.format", tip.length = 0) +
  scale_fill_manual(values = jco_cols, name = "Ethnicity") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
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
       width = 8, height = max(8, nrow(sig_vf) * 4 + 3))

## VFDB baseline heatmap (top 20 by Wilcoxon W)
vf_baseline_means <- vf_df_filt %>%
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

cat("\nDone. Results written to:", results_dir, "\n")
