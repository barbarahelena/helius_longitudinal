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
# Same denominator as the clade comparisons in 3_draw_tree.R: CDS predicted by
# Bakta per bin. Dividing by the bin's total VF hits instead would make the
# categories compositional (they would sum to 1), so a rise in one category
# would force a fall in the others and the categories could not be tested
# independently.
total_per_bin <- load_gene_counts(trans, batch_files)

vf_long <- vfdb_hits %>%
  filter(!is.na(vf_category)) %>%
  count(locus_prefix, vf_category, name = "n_hits") %>%
  inner_join(total_per_bin, by = "locus_prefix") %>%
  mutate(proportion = n_hits / total_cds)

cat("Unique VF categories:", n_distinct(vf_long$vf_category), "\n")

#### 5. Join with clinical data ####
bin_clin <- load_bin_clin(trans, batch_files, clin_file)

# One row per bin. Gene content is a property of the assembled genome and each
# participant has one MAG, so it has no per-timepoint value; a bin is included
# if it was detected at either timepoint.
bins_present <- bin_clin %>%
  filter(present, locus_prefix %in% unique(vf_long$locus_prefix)) %>%
  distinct(locus_prefix, EthnicityTot)
stopifnot(!any(duplicated(bins_present$locus_prefix)))

cat("\nBins with VF hits, per ethnicity:\n")
print(count(bins_present, EthnicityTot))

#### 6. Build analysis dataset ####
# Complete grid: every bin × every VF category; missing proportion = 0.
# Completeness and contamination travel with each bin so the model can adjust
# for them — both bias gene-content measures and both differ by ethnicity.
vf_df <- expand.grid(
  locus_prefix = bins_present$locus_prefix,
  vf_category  = unique(vf_long$vf_category),
  stringsAsFactors = FALSE
) %>%
  left_join(vf_long %>% dplyr::select(locus_prefix, vf_category, proportion),
            by = c("locus_prefix", "vf_category")) %>%
  mutate(proportion = replace_na(proportion, 0)) %>%
  left_join(bins_present, by = "locus_prefix") %>%
  add_bin_quality(trans, batch_files)

cat("\nMAG quality by ethnicity (adjusted for, not filtered on):\n")
print(vf_df %>%
        distinct(locus_prefix, EthnicityTot, Completeness, Contamination) %>%
        group_by(EthnicityTot) %>%
        summarise(n = n(),
                  median_completeness  = round(median(Completeness), 1),
                  median_contamination = round(median(Contamination), 2),
                  .groups = "drop"))

#### 7. Prevalence filter ####
MIN_PREV <- 0.10

prevalent_vf <- vf_df %>%
  group_by(vf_category) %>%
  summarise(prev = mean(proportion > 0), .groups = "drop") %>%
  filter(prev >= MIN_PREV) %>%
  pull(vf_category)

cat("\nVF categories passing prevalence filter (>=", MIN_PREV * 100, "% of bins):",
    length(prevalent_vf), "of", n_distinct(vf_df$vf_category), "\n")

vf_df_filt <- vf_df %>% filter(vf_category %in% prevalent_vf)

#### 8. Statistics ####
cat("\nRunning VF category statistics (bin-level, adjusted for MAG quality)...\n")
vf_results <- run_bin_stats(vf_df_filt, "vf_category") %>%
  rename(vf_category = feature)
print(vf_results %>%
        dplyr::select(vf_category, n_dutch, n_sas, wilcox_p, wilcox_fdr, adj_p, adj_fdr) %>%
        arrange(adj_p))

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

# Significance is taken from the quality-adjusted model. The unadjusted Wilcoxon
# is reported alongside it, but MAG completeness and contamination differ
# between the two groups, so it is descriptive rather than a test of ethnicity.
sig_vf <- vf_results %>% filter(adj_fdr < 0.05)

cat("\nTesting", n_distinct(vf_df_filt$vf_category), "VF categories (N Dutch:",
    unique(vf_results$n_dutch), "| N SAS:", unique(vf_results$n_sas), ")\n")
cat("Significant VF categories (adjusted FDR < 0.05):", nrow(sig_vf), "\n")
print(sig_vf %>% dplyr::select(vf_category, n_dutch, n_sas,
                               median_dutch, median_sas,
                               adj_estimate, adj_p, adj_fdr,
                               wilcox_p, wilcox_fdr))
cat("Categories significant before adjustment but not after:",
    sum(vf_results$wilcox_fdr < 0.05 & vf_results$adj_fdr >= 0.05, na.rm = TRUE), "\n")

#### 9. Figure panel pl_fig3_F ####
jco_cols <- jco_palette()

# The two boxplots below show only the categories reaching FDR < 0.05. When none
# do, facet_wrap has nothing to lay out, so skip them rather than fail: a null
# result is a valid outcome, not a broken run. The heatmap below covers all
# categories and is drawn either way.
plot_vf_box <- vf_df_filt %>% filter(vf_category %in% sig_vf$vf_category)

if (nrow(sig_vf) == 0) {
  cat("\nNo VF category reached FDR < 0.05; skipping the significant-category",
      "boxplots (vfdb_sig_baseline.pdf, qc_vfdb_boxplot_all_timepoints.pdf).\n")
  pl_fig3_F <- NULL
} else {

pl_fig3_F <- ggplot(
    plot_vf_box,
    aes(x = EthnicityTot, y = proportion, fill = EthnicityTot)
  ) +
  geom_boxplot(outlier.size = 0.8, width = 0.5) +
  stat_compare_means(comparisons = list(c("Dutch", "South-Asian Surinamese")),
                     method = "wilcox.test", label = "p.format", tip.length = 0) +
  scale_fill_manual(values = jco_cols, guide = "none") +
  facet_wrap(~ vf_category, scales = "free_y") +
  scale_y_continuous(expand = expansion(add = c(0, 0.010))) +
  labs(title    = "VFDB: Alistipes putredinis",
       subtitle = "Significant after adjustment for MAG completeness and contamination",
       x = "", y = "Proportion of predicted genes",
       caption  = "Bracket shows the unadjusted Wilcoxon p; significance is from the adjusted model") +
  theme_Publication() +
  theme(strip.text = element_text(size = rel(0.8)))

ggsave(file.path(results_dir, "vfdb_sig_categories.pdf"),
       plot = pl_fig3_F,
       width = 8, height = max(5, ceiling(nrow(sig_vf) / 3) * 2.8 + 1))

}  # end: significant categories present

## VFDB heatmap (top 20 by unadjusted Wilcoxon W)
vf_baseline_means <- vf_df_filt %>%
  distinct(locus_prefix, vf_category, EthnicityTot, proportion) %>%
  group_by(vf_category, EthnicityTot) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = EthnicityTot, values_from = mean_prop, values_fill = 0) %>%
  mutate(diff_Dutch_SAS = Dutch - `South-Asian Surinamese`)

# Stars mark the quality-adjusted FDR, so the heatmap and the boxplots above
# agree on what counts as significant.
vf_heat_base <- vf_results %>%
  filter(!is.na(wilcox_p)) %>%
  slice_max(abs(wilcox_stat), n = 20) %>%
  left_join(vf_baseline_means %>% dplyr::select(vf_category, diff_Dutch_SAS),
            by = "vf_category")

vf_ba_heat <- vf_heat_base %>%
  mutate(
    sig_label   = case_when(
      adj_fdr < 0.001 ~ "***",
      adj_fdr < 0.01  ~ "**",
      adj_fdr < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    vf_category = factor(vf_category,
                         levels = vf_heat_base %>%
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
                       name = "Mean proportion\nDutch − SAS") +
  scale_x_continuous(breaks = NULL) +
  labs(title    = "Top 20 VF categories — ethnicity difference",
       subtitle = "blue = higher in Dutch, yellow = higher in SAS",
       x = "", y = "") +
  theme_Publication() +
  theme(axis.text.y = element_text(size = rel(0.75)))

ggsave(file.path(results_dir, "qc_vfdb_heatmap_baseline.pdf"),
       plot = p_vf_heat_ba, width = 7, height = 9)

cat("\nDone. Results written to:", results_dir, "\n")
