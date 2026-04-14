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

present_bins <- bin_clin %>%
  filter(present, locus_prefix %in% unique(vf_long$locus_prefix)) %>%
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

#### 7. Statistics ####
cat("\nRunning VF category statistics...\n")
vf_results <- run_stats(vf_df, "vf_category") %>%
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
  filter(wilcox_baseline_fdr < 0.05 | wilcox_fu_fdr < 0.05 |
           lmm_eth_fdr < 0.05 | lmm_int_fdr < 0.05)

cat("\nSignificant VF categories (FDR < 0.05 in any comparison):", nrow(sig_vf), "\n")
print(sig_vf %>% dplyr::select(vf_category, lmm_eth_estimate,
                                lmm_eth_fdr, wilcox_baseline_fdr))

#### 8. Figure panel pl_fig3_F ####
jco_cols <- jco_palette()

plot_vf_box <- vf_df %>% filter(vf_category %in% sig_vf$vf_category)

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

cat("\nDone. Results written to:", results_dir, "\n")
