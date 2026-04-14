## eggNOG functional annotation — Alistipes putredinis bins
## Dutch vs South-Asian Surinamese, baseline and follow-up
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

source("scripts/3_species_change/4_alistipes_anno/utils.R")

#### Paths ####
anno_file   <- "data/shotgun/alistipes_annotation/eggnog/all_eggnog_results.annotations"
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

#### 2. Parse eggNOG annotation file ####
# Merged file: blocks of ## comments + repeated #query header lines per bin.
# Read column names from the first #query line, then skip all # lines.
raw_lines  <- readLines(anno_file, n = 500)
header_idx <- which(startsWith(raw_lines, "#query"))
if (length(header_idx) == 0)
  stop("No header line (#query) found in the first 500 lines of annotations file")

col_names <- strsplit(sub("^#", "", raw_lines[header_idx[1]]), "\t")[[1]]

anno <- read.table(
  anno_file,
  header = FALSE, sep = "\t", comment.char = "#",
  quote = "", col.names = col_names, stringsAsFactors = FALSE, fill = TRUE
) %>%
  mutate(locus_prefix = sub("_.*", "", query))

cat("Total annotated genes:", nrow(anno), "\n")

#### 3. Per-bin KEGG module proportions ####
clean_field <- function(x) {
  x <- trimws(x)
  ifelse(x == "" | x == "-", NA_character_, x)
}

anno_filt <- anno %>%
  mutate(kegg_clean = clean_field(KEGG_Module),
         cog_clean  = clean_field(COG_category)) %>%
  filter(!is.na(kegg_clean) | !is.na(cog_clean))

total_per_bin <- anno_filt %>%
  count(locus_prefix, name = "total_annotated")

kegg_long <- anno_filt %>%
  filter(!is.na(kegg_clean)) %>%
  mutate(module = strsplit(kegg_clean, ",")) %>%
  unnest(module) %>%
  mutate(module = trimws(module)) %>%
  filter(module != "" & module != "-") %>%
  count(locus_prefix, module, name = "n_genes") %>%
  left_join(total_per_bin, by = "locus_prefix") %>%
  mutate(proportion = n_genes / total_annotated)

#### 4. Fetch KEGG module descriptions (cached) ####
kegg_lookup_file <- file.path(results_dir, "kegg_module_lookup.rds")

if (file.exists(kegg_lookup_file)) {
  module_lookup <- readRDS(kegg_lookup_file)
} else {
  cat("Fetching KEGG module descriptions from REST API...\n")
  raw <- tryCatch(
    readLines("https://rest.kegg.jp/list/module", warn = FALSE),
    error = function(e) { warning("KEGG API unavailable: ", e$message); character(0) }
  )
  if (length(raw) > 0) {
    module_lookup <- tibble(raw = raw) %>%
      separate(raw, into = c("module", "description"), sep = "\t", extra = "merge") %>%
      mutate(module = sub("^md:", "", module))
    saveRDS(module_lookup, kegg_lookup_file)
  } else {
    module_lookup <- tibble(module = character(), description = character())
    warning("KEGG lookup empty — module IDs will be used as labels.")
  }
}

kegg_long <- kegg_long %>%
  left_join(module_lookup, by = "module") %>%
  mutate(module_label = if_else(is.na(description), module,
                                paste0(module, ": ", description)))

module_labels <- kegg_long %>% distinct(module, module_label)

#### 5. Join with clinical data ####
bin_clin <- load_bin_clin(trans, batch_files, clin_file)

present_bins <- bin_clin %>%
  filter(present) %>%
  dplyr::select(locus_prefix, sampleID, EthnicityTot, timepoint, depth)

cat("\nPresent bins per ethnicity × timepoint:\n")
print(present_bins %>%
        distinct(locus_prefix, EthnicityTot, timepoint) %>%
        count(EthnicityTot, timepoint))

#### 6. Build analysis dataset ####
# Complete grid: every present bin × every module; missing proportion = 0
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

#### 7. Statistics ####
cat("\nRunning KEGG module statistics...\n")
kegg_results <- run_stats(kegg_df, "module") %>%
  rename(module = feature) %>%
  left_join(module_labels, by = "module") %>%
  dplyr::relocate(module_label, .after = module)

write.csv(kegg_results, file.path(results_dir, "kegg_module_stats.csv"),
          row.names = FALSE)

sig_kegg <- kegg_results %>%
  filter(wilcox_baseline_fdr < 0.05 | wilcox_fu_fdr < 0.05 |
           lmm_eth_fdr < 0.05 | lmm_int_fdr < 0.05)

cat("\nSignificant KEGG modules (FDR < 0.05 in any comparison):", nrow(sig_kegg), "\n")
print(sig_kegg %>% dplyr::select(module, module_label, lmm_eth_estimate,
                                  lmm_eth_fdr, wilcox_baseline_fdr))

#### 8. Figure panel pl_fig3_E ####
jco_cols <- jco_palette()

plot_kegg_box <- kegg_df %>%
  filter(module %in% sig_kegg$module) %>%
  mutate(module_label = fct_recode(module_label,
    "Ascorbate biosynthesis"  = "M00129: Ascorbate biosynthesis, animals, glucose-1P => ascorbate",
    "Glucuronate pathway"     = "M00014: Glucuronate pathway (uronate pathway)"
  ))

pl_fig3_E <- ggplot(
    plot_kegg_box %>% filter(timepoint == "baseline"),
    aes(x = EthnicityTot, y = proportion, fill = EthnicityTot)
  ) +
  geom_boxplot(outlier.size = 0.8, width = 0.5) +
  stat_compare_means(comparisons = list(c("Dutch", "South-Asian Surinamese")),
                     method = "wilcox.test", label = "p.format", tip.length = 0) +
  scale_fill_manual(values = jco_cols, guide = "none") +
  facet_wrap(~ module_label, scales = "free_y",
             labeller = label_wrap_gen(width = 45)) +
  scale_y_continuous(expand = expansion(add = c(0, 0.00025))) +
  labs(title = "KEGG pathways: Alistipes putredinis",
       x = "", y = "Proportion of annotated genes") +
  theme_Publication() +
  theme(strip.text = element_text(size = rel(0.8)))

cat("\nDone. Results written to:", results_dir, "\n")
