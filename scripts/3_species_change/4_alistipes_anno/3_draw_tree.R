## Phylogenetic tree — Alistipes putredinis bins
## Annotated with ethnicity, bootstrap, VFDB motility, KEGG modules, bin quality
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl
##
## Required packages (install once if needed):
##   BiocManager::install(c("ggtree", "treeio"))
##   install.packages("ggnewscale")

source("scripts/3_species_change/4_alistipes_anno/utils.R")

library(ggtree)
library(treeio)
library(ggnewscale)
library(ggsci)

#### Paths ####
tree_file   <- "data/shotgun/alistipes_annotation/alistipes_new.treefile"
vfdb_file   <- "data/shotgun/alistipes_annotation/all_vfdb_results.txt"
fasta_file  <- "data/shotgun/alistipes_annotation/VFDB_setB_pro.fas"
anno_file   <- "data/shotgun/alistipes_annotation/eggnog/all_eggnog_results.annotations"
trans_file  <- "data/shotgun/alistipes_annotation/bin_translation_table.tsv"
clin_file   <- "data/clinicaldata/clinicaldata_long.RDS"
diet_file   <- "data/clinicaldata_long_pcdiet.RDS"
batch_files <- c(
  "data/shotgun/alistipes_annotation/bins_alistipes_batch1.csv",
  "data/shotgun/alistipes_annotation/bins_alistipes_batch2.csv",
  "data/shotgun/alistipes_annotation/bins_alistipes_batch3.csv"
)
results_dir <- "results/3_species_change/4_alistipes_anno"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

#### Constants ####
MIN_PIDENT      <- 30
MIN_BITSCORE    <- 50
TARGET_MODULES  <- c("M00129", "M00014")
MIN_CLADE_SIZE  <- 5   # clades smaller than this are excluded from PLOT_CLADES (too few bins for group-wise stats)
N_TOP_VF        <- 2   # number of top VFDB categories to display
# PLOT_CLADES is set dynamically in section 8a, after clade calling, from all
# clades with >= MIN_CLADE_SIZE tips — see that section for why.

#### 1. Load tree ####
tree <- read.iqtree(tree_file)

#### 2. Load translation sheet ####
trans <- read.delim(trans_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  rename(locus_prefix = locus_tag_prefix)

cat("Bins in translation table:", nrow(trans), "\n")

#### 3. Subject ID, ethnicity, timepoint per bin ####
# Determine dominant sample (max depth) per bin → subject_id + timepoint
depth_raw <- map_dfr(batch_files, function(f) {
  df <- read.csv(f, check.names = FALSE)
  if (nrow(df) == 0) return(NULL)
  df %>%
    dplyr::select(bin, starts_with("Depth ")) %>%
    mutate(bin_name = sub("\\.fa$", "", bin),
           across(starts_with("Depth "), as.numeric)) %>%
    dplyr::select(bin_name, starts_with("Depth "))
})

bin_dominant <- depth_raw %>%
  pivot_longer(cols = starts_with("Depth "),
               names_to = "sampleID", values_to = "depth") %>%
  mutate(sampleID = sub("^Depth ", "", sampleID),
         depth    = replace_na(depth, 0)) %>%
  filter(depth > 0) %>%
  group_by(bin_name) %>%
  slice_max(depth, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    timepoint = if_else(str_starts(sampleID, "HELIBA"), "baseline", "follow-up")
  ) %>%
  dplyr::select(bin_name, sampleID, timepoint)

# Join via sampleID directly — clinical data ID is "S100001" format, not numeric
clin_all <- readRDS(clin_file) %>%
  filter(EthnicityTot %in% c("Dutch", "South-Asian Surinamese")) %>%
  dplyr::select(sampleID, ID, EthnicityTot, Age, Sex, BMI, FUtime)

clin_eth <- clin_all %>% distinct(sampleID, .keep_all = TRUE)

bin_clin_meta <- bin_dominant %>%
  left_join(clin_eth, by = "sampleID")

cat("Bins with ethnicity assigned:",
    sum(!is.na(bin_clin_meta$EthnicityTot)), "of", nrow(bin_clin_meta), "\n")

#### 4. Depth at own subject's baseline and follow-up ####
subject_ids <- bin_dominant %>%
  mutate(subject_id = str_extract(sampleID, "\\d+$")) %>%
  dplyr::select(bin_name, subject_id)

depth_tp <- depth_raw %>%
  pivot_longer(cols = starts_with("Depth "),
               names_to = "sampleID", values_to = "depth") %>%
  mutate(sampleID = sub("^Depth ", "", sampleID),
         depth    = replace_na(depth, 0)) %>%
  inner_join(subject_ids, by = "bin_name") %>%
  filter(sampleID == paste0("HELIBA_", subject_id) |
         sampleID == paste0("HELIFU_",  subject_id)) %>%
  mutate(tp = if_else(str_starts(sampleID, "HELIBA"), "depth_baseline", "depth_followup")) %>%
  dplyr::select(bin_name, tp, depth) %>%
  pivot_wider(names_from = tp, values_from = depth, values_fill = 0)

cat("Bins with depth at both timepoints:", nrow(depth_tp), "\n")

#### 5. KEGG module proportions per bin ("M00129", "M00014") ####
raw_lines  <- readLines(anno_file, n = 500)
header_idx <- which(startsWith(raw_lines, "#query"))
if (length(header_idx) == 0)
  stop("No #query header found in eggnog annotations file")
eg_col_names <- strsplit(sub("^#", "", raw_lines[header_idx[1]]), "\t")[[1]]

anno <- read.table(
  anno_file,
  header = FALSE, sep = "\t", comment.char = "#",
  quote = "", col.names = eg_col_names, stringsAsFactors = FALSE, fill = TRUE
) %>%
  mutate(locus_prefix = sub("_.*", "", query))

clean_field <- function(x) {
  x <- trimws(x)
  ifelse(x == "" | x == "-", NA_character_, x)
}

total_per_bin <- anno %>%
  mutate(kegg_clean = clean_field(KEGG_Module),
         cog_clean  = clean_field(COG_category)) %>%
  filter(!is.na(kegg_clean) | !is.na(cog_clean)) %>%
  count(locus_prefix, name = "total_annotated")

kegg_per_bin <- anno %>%
  mutate(KEGG_Module = trimws(KEGG_Module)) %>%
  filter(KEGG_Module != "", KEGG_Module != "-") %>%
  mutate(module = strsplit(KEGG_Module, ",")) %>%
  unnest(module) %>%
  mutate(module = trimws(module)) %>%
  filter(module %in% TARGET_MODULES) %>%
  count(locus_prefix, module, name = "n_genes") %>%
  left_join(total_per_bin, by = "locus_prefix") %>%
  mutate(proportion = n_genes / total_annotated) %>%
  dplyr::select(locus_prefix, module, proportion) %>%
  pivot_wider(names_from = module, values_from = proportion, values_fill = 0)

cat("Bins with KEGG annotation:", nrow(kegg_per_bin), "\n")

#### 6. VFDB category proportions per bin (all categories) ####
raw_headers <- readLines(fasta_file) %>%
  keep(startsWith, ">") %>%
  sub("^>", "", .)

vfdb_anno <- tibble(header = raw_headers) %>%
  mutate(
    vfg_id      = str_extract(header, "^VFG\\d+"),
    vf_category = str_match(header, "-\\s+([^\\[\\]]+?)\\s+\\(VFC\\d+\\)")[, 2]
  ) %>%
  dplyr::select(vfg_id, vf_category) %>%
  filter(!is.na(vfg_id)) %>%
  distinct(vfg_id, .keep_all = TRUE)

vf_col_names <- c("query", "subject", "pident", "length", "mismatch",
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

vfdb_hits <- read.table(
  vfdb_file,
  header = FALSE, sep = "\t", comment.char = "#",
  quote = "", col.names = vf_col_names, stringsAsFactors = FALSE, fill = TRUE
) %>%
  filter(pident >= MIN_PIDENT, bitscore >= MIN_BITSCORE) %>%
  mutate(
    locus_prefix = sub("_.*", "", query),
    vfg_id       = str_extract(subject, "^VFG\\d+")
  ) %>%
  left_join(vfdb_anno, by = "vfg_id")

vf_per_bin <- vfdb_hits %>%
  filter(!is.na(vf_category)) %>%
  count(locus_prefix, vf_category, name = "n_hits") %>%
  left_join(total_per_bin, by = "locus_prefix") %>%
  mutate(proportion = n_hits / total_annotated) %>%
  dplyr::select(locus_prefix, vf_category, proportion) %>%
  pivot_wider(names_from = vf_category, values_from = proportion, values_fill = 0)

# All VFDB category column names (used throughout)
vf_cat_cols <- setdiff(names(vf_per_bin), "locus_prefix")
cat("VF categories detected:", length(vf_cat_cols), ":",
    paste(vf_cat_cols, collapse = ", "), "\n")

#### 7. Assemble tip metadata ####
# Key: bin_name (= tree tip label)
tip_meta <- trans %>%
  left_join(bin_clin_meta, by = "bin_name") %>%
  left_join(vf_per_bin,    by = "locus_prefix") %>%
  left_join(kegg_per_bin,  by = "locus_prefix") %>%
  left_join(depth_tp,      by = "bin_name") %>%
  mutate(
    EthnicityTot = factor(EthnicityTot,
                          levels = c("Dutch", "South-Asian Surinamese")),
    timepoint    = factor(timepoint, levels = c("baseline", "follow-up")),
    across(all_of(c(vf_cat_cols, TARGET_MODULES)), ~ replace_na(.x, 0)),
    depth_baseline = replace_na(depth_baseline, 0),
    depth_followup = replace_na(depth_followup, 0)
  )

cat("Tips in metadata before filter:", nrow(tip_meta), "\n")
cat("Tips in tree:", length(tree@phylo$tip.label), "\n")

tip_meta <- tip_meta %>%
  filter(bin_name %in% tree@phylo$tip.label)

cat("Tips in metadata after filter:", nrow(tip_meta), "\n")

#### 8. Draw tree ####
eth_colors <- jco_palette()  # Dutch = blue, SAS = orange/yellow

# ---- 8a. Call clades from cophenetic (patristic) distances ----
# Clades are called with a threshold h on the maximum pairwise cophenetic
# distance within a group, applied TOP-DOWN on the tree topology itself:
# starting at the root, descend into a node's children only if its subtree's
# internal max distance exceeds h; otherwise the whole subtree is one clade.
# Every group returned is therefore a genuine clade (or a singleton) BY
# CONSTRUCTION — unlike hclust(..., method = "average") on the same distance
# matrix, which can merge tips that are close in patristic distance without
# being each other's closest relatives on the tree. That matters here: at
# h = 0.02, average-linkage/cutree put 83 of 180 tips in one "clade" whose
# tips' MRCA was the root of the whole tree (i.e. not a clade at all), and a
# second 13-tip group was likewise non-monophyletic. Both problems are fixed
# by construction below (see the diagnostic comparison printed further down).
phylo_tree <- tree@phylo
D          <- ape::cophenetic.phylo(phylo_tree)

call_monophyletic_clades <- function(phylo, D, h) {
  root         <- setdiff(phylo$edge[, 1], phylo$edge[, 2])
  n_tip        <- length(phylo$tip.label)
  get_children <- function(node) phylo$edge[phylo$edge[, 1] == node, 2]
  get_tips     <- function(node) {
    if (node <= n_tip) return(phylo$tip.label[node])
    ape::extract.clade(phylo, node)$tip.label
  }
  clusters <- list()
  stack    <- list(root)
  while (length(stack) > 0) {
    node  <- stack[[1]]
    stack <- stack[-1]
    tips  <- get_tips(node)
    max_d <- if (length(tips) == 1) 0 else max(D[tips, tips])
    if (max_d <= h) {
      clusters[[length(clusters) + 1]] <- tips
    } else {
      stack <- c(stack, as.list(get_children(node)))
    }
  }
  setNames(
    unlist(lapply(seq_along(clusters), function(i) rep(i, length(clusters[[i]])))),
    unlist(clusters)
  )
}

# H_CUT: a fine sweep of h (0.001 steps, 0.015-0.05) found no natural break
# near 0.02 for the old average-linkage method. Under the monophyletic caller
# the same sweep shows a clean, wide plateau at 9 clades for h in
# [0.029, 0.040] (12 consecutive grid points with no change) — by far the
# most stable partition in the whole range — so that plateau is used as the
# primary, data-driven cutoff. h = 0.02 (finer-grained, 26 clades) is kept
# below purely as a sensitivity comparison, alongside the rest of the sweep.
H_CUT             <- 0.03    # primary: start of the [0.029, 0.040] plateau
H_CUT_SENSITIVITY <- 0.02    # for comparison only, see clade_h_sensitivity.csv

tip_clusters <- call_monophyletic_clades(phylo_tree, D, H_CUT)
cat("Monophyletic clades detected at h =", H_CUT, ":", n_distinct(tip_clusters), "\n")

# Diagnostic: how different is this from the old average-linkage/cutree call
# at the same h? (quantifies the effect of the monophyly fix)
hc_old           <- hclust(as.dist(D), method = "average")
tip_clusters_old <- cutree(hc_old, h = H_CUT)
cat("For comparison, average-linkage cutree at the same h gives:",
    n_distinct(tip_clusters_old), "clusters (not all guaranteed monophyletic)\n")

# For each cluster find its MRCA node (used for geom_hilight). With the
# monophyletic caller, this MRCA's descendant tips equal the cluster exactly.
clade_nodes <- map_dfr(sort(unique(tip_clusters)), function(cl) {
  tips <- names(tip_clusters)[tip_clusters == cl]
  nd <- if (length(tips) == 1) {
    which(tree@phylo$tip.label == tips)   # tip node index directly
  } else {
    ape::getMRCA(tree@phylo, tips)
  }
  tibble(node = nd, n_tips = length(tips), cluster = cl)
}) %>%
  filter(!is.na(node))

# Readable "Clade I, II, ..." labels are given only to clusters big enough
# for group-wise stats (n_tips >= MIN_CLADE_SIZE), numbered by descending
# size so Clade I is the largest — this keeps the numbering contiguous
# (I, II, III, ...) instead of skipping labels for excluded clusters.
# Smaller clusters get a plain descriptive label instead of a roman numeral.
clade_nodes <- clade_nodes %>% arrange(desc(n_tips))
named_idx   <- which(clade_nodes$n_tips >= MIN_CLADE_SIZE)
clade_nodes$clade <- NA_character_
clade_nodes$clade[named_idx]  <- paste0("Clade ", as.roman(seq_along(named_idx)))
clade_nodes$clade[-named_idx] <- paste0("Unclassified (cluster ", clade_nodes$cluster[-named_idx],
                                         ", n=", clade_nodes$n_tips[-named_idx], ")")
clade_nodes <- clade_nodes %>% arrange(node)

# Sanity check: every multi-tip cluster must be monophyletic by construction.
mono_check <- clade_nodes %>%
  filter(n_tips > 1) %>%
  rowwise() %>%
  mutate(monophyletic = ape::is.monophyletic(
    phylo_tree, names(tip_clusters)[tip_clusters == cluster])) %>%
  ungroup()
if (!all(mono_check$monophyletic)) {
  stop("Non-monophyletic cluster produced by call_monophyletic_clades() — this should not happen.")
}
cat("Monophyly check: all", nrow(mono_check), "multi-tip clades confirmed monophyletic.\n")

cat("Detected clades (all, before any size filtering):\n")
print(clade_nodes %>% dplyr::select(clade, n_tips, node))

# PLOT_CLADES: every clade with enough tips for meaningful group-wise stats,
# not a hand-picked subset — so no clade is silently dropped from the figures
# or tables. Smaller clades are reported here (and in the CSV below) but
# excluded from group comparisons that need reasonable group sizes.
PLOT_CLADES <- clade_nodes %>% filter(n_tips >= MIN_CLADE_SIZE) %>% pull(clade)
cat("\nPLOT_CLADES (n_tips >=", MIN_CLADE_SIZE, "):", paste(PLOT_CLADES, collapse = ", "), "\n")
cat("Excluded from PLOT_CLADES (too few tips for group-wise stats):",
    paste(setdiff(clade_nodes$clade, PLOT_CLADES), collapse = ", "), "\n")

# Full ethnicity composition of ALL detected clades, not just PLOT_CLADES —
# exported so nothing is left unaccounted for in the text/methods.
clade_eth_all <- tibble(bin_name = names(tip_clusters), cluster = tip_clusters) %>%
  left_join(clade_nodes %>% dplyr::select(cluster, clade, n_tips), by = "cluster") %>%
  left_join(tip_meta %>% dplyr::select(bin_name, EthnicityTot), by = "bin_name") %>%
  count(clade, n_tips, EthnicityTot) %>%
  arrange(desc(n_tips))
cat("\nEthnicity composition of ALL detected clades:\n")
print(clade_eth_all, n = Inf)
write.csv(clade_eth_all, file.path(results_dir, "clade_ethnicity_all_clades.csv"), row.names = FALSE)

# Hand-picked high-contrast palette for clade highlights; extended via
# interpolation if there are more clades than base colors (was fixed at 10,
# but the monophyletic caller can produce more clades than that).
clade_pal_base <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
                     "#76B7B2", "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC")
clade_pal <- if (nrow(clade_nodes) <= length(clade_pal_base)) {
  clade_pal_base[seq_len(nrow(clade_nodes))]
} else {
  grDevices::colorRampPalette(clade_pal_base)(nrow(clade_nodes))
}
clade_colors <- setNames(clade_pal, clade_nodes$clade)

# ---- 8a-ii. Assign clades and rank VFDB categories by clade difference ----
clade_map <- tibble(
  bin_name = names(tip_clusters),
  cluster  = tip_clusters
) %>%
  left_join(clade_nodes %>% dplyr::select(cluster, clade), by = "cluster")

tip_meta_clades <- tip_meta %>%
  left_join(clade_map %>% dplyr::select(bin_name, clade), by = "bin_name")

# Run KW for all VFDB categories across PLOT_CLADES only (baseline bins only)
# — restricted to clades with enough tips (MIN_CLADE_SIZE) so the comparison
# isn't diluted by singleton/near-singleton groups.
tip_meta_focal <- tip_meta_clades %>%
  filter(timepoint == "baseline", clade %in% PLOT_CLADES)

vf_kw <- map_dfr(vf_cat_cols, function(cat) {
  vals  <- tip_meta_focal[[cat]]
  grp   <- factor(tip_meta_focal$clade)
  if (all(is.na(vals)) || var(vals, na.rm = TRUE) == 0)
    return(tibble(vf_category = cat, kw_p = 1, kw_stat = NA_real_))
  kw <- kruskal.test(vals ~ grp)
  tibble(vf_category = cat, kw_p = kw$p.value, kw_stat = kw$statistic)
}) %>%
  arrange(kw_p)

cat("\nVF categories ranked by clade difference (KW p-value):\n")
print(vf_kw, n = Inf)

TOP_VF_CATS <- vf_kw %>%
  filter(kw_p < 0.05) %>%
  slice_head(n = N_TOP_VF) %>%
  pull(vf_category)

# Fallback: take top N_TOP_VF even if none reach p < 0.05
if (length(TOP_VF_CATS) < N_TOP_VF) {
  TOP_VF_CATS <- vf_kw %>%
    slice_head(n = N_TOP_VF) %>%
    pull(vf_category)
  cat("Note: fewer than", N_TOP_VF, "categories reached KW p < 0.05; using top",
      N_TOP_VF, "by test statistic.\n")
}

cat("Top", N_TOP_VF, "VF categories for plots:", paste(TOP_VF_CATS, collapse = ", "), "\n")

# ---- 8a-iii. Sensitivity of clade calling and downstream findings to h ----
# Re-derive clades at each h with the same monophyletic caller, then check
# whether (i) the number of clades that are majority-Dutch vs majority-SAS
# among PLOT_CLADES-sized clusters is stable, and (ii) the top VF categories
# above stay significantly different across clades (Kruskal-Wallis).
H_SWEEP <- c(0.02, 0.025, 0.03, 0.035, 0.04)

sensitivity_results <- map_dfr(H_SWEEP, function(h) {
  cl_h  <- call_monophyletic_clades(phylo_tree, D, h)
  map_h <- tibble(bin_name = names(cl_h), cluster = cl_h)
  meta_h <- tip_meta %>%
    inner_join(map_h, by = "bin_name") %>%
    filter(!is.na(EthnicityTot))

  big_clusters <- meta_h %>%
    count(cluster) %>%
    filter(n >= MIN_CLADE_SIZE) %>%
    pull(cluster)

  eth_by_cluster <- meta_h %>%
    filter(cluster %in% big_clusters) %>%
    count(cluster, EthnicityTot) %>%
    group_by(cluster) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  n_dutch_enriched <- eth_by_cluster %>%
    filter(EthnicityTot == "Dutch", prop > 0.5) %>% nrow()
  n_sas_enriched <- eth_by_cluster %>%
    filter(EthnicityTot == "South-Asian Surinamese", prop > 0.5) %>% nrow()

  vf_h <- map_dfr(TOP_VF_CATS, function(cat) {
    sub <- meta_h %>% filter(cluster %in% big_clusters, timepoint == "baseline")
    if (n_distinct(sub$cluster) < 2 || var(sub[[cat]], na.rm = TRUE) == 0) {
      return(tibble(vf_category = cat, kw_p = NA_real_))
    }
    tibble(vf_category = cat,
           kw_p = kruskal.test(sub[[cat]] ~ factor(sub$cluster))$p.value)
  }) %>%
    pivot_wider(names_from = vf_category, values_from = kw_p, names_prefix = "kw_p_")

  tibble(h = h, n_clusters_total = n_distinct(cl_h),
         n_clusters_ge_min = length(big_clusters),
         n_dutch_enriched  = n_dutch_enriched,
         n_sas_enriched    = n_sas_enriched) %>%
    bind_cols(vf_h)
})

cat("\nSensitivity of clade calling across h:\n")
print(sensitivity_results, n = Inf)
write.csv(sensitivity_results, file.path(results_dir, "clade_h_sensitivity.csv"), row.names = FALSE)

# ---- 8b. Base tree with clade highlights ----
# Pre-compute tip layout here so we can restrict each highlight to the exact
# angular (y) range of its own tips — prevents the largest clade from wrapping
# around the entire circle.
p_base <- ggtree(tree, layout = "circular", size = 0.35)

tree_layout <- fortify(tree)
tip_layout  <- tree_layout %>% filter(isTip)
max_tip_x   <- max(tip_layout$x, na.rm = TRUE)

clade_y_ranges <- tibble(
  label   = names(tip_clusters),
  cluster = tip_clusters
) %>%
  left_join(clade_nodes %>% dplyr::select(cluster, clade), by = "cluster") %>%
  left_join(tip_layout %>% dplyr::select(label, y), by = "label") %>%
  filter(!is.na(clade)) %>%
  group_by(clade) %>%
  summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop")

for (i in seq_len(nrow(clade_y_ranges))) {
  p_base <- p_base +
    annotate("rect",
             xmin  = 0,
             xmax  = max_tip_x * 1.006,
             ymin  = clade_y_ranges$ymin[i],
             ymax  = clade_y_ranges$ymax[i],
             fill  = clade_colors[clade_y_ranges$clade[i]],
             alpha = 0.2)
}

p_base <- p_base %<+%
  (tip_meta %>% rename(label = bin_name)) +
  geom_tippoint(
    aes(color = EthnicityTot, shape = timepoint),
    size = 2.5, na.rm = TRUE
  ) +
  scale_color_manual(
    values   = eth_colors,
    name     = "Ethnicity",
    na.value = "grey70"
  ) +
  scale_shape_manual(
    values   = c("baseline" = 16, "follow-up" = 17),
    name     = "Timepoint",
    na.value = 1
  ) +
  theme_tree() +
  theme(
    legend.position   = "right",
    legend.box        = "vertical",
    legend.text       = element_text(size = 14),
    legend.title      = element_text(size = 15, face = "bold"),
    legend.key.size   = unit(0.7, "cm"),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA)
  )

# Clade labels — positioned at the mean y (angle) of each clade's tips,
# pushed outward past the heatmap rings.
clade_nodes_focal <- clade_nodes %>%
  filter(clade %in% PLOT_CLADES) %>%
  arrange(clade) %>%
  mutate(short_label = sub("^Clade ", "", clade))

# Compute mean y (= angular position) of tips per clade
clade_tip_angles <- tibble(
  label   = names(tip_clusters),
  cluster = tip_clusters
) %>%
  left_join(clade_nodes %>% dplyr::select(cluster, clade), by = "cluster") %>%
  left_join(tip_layout %>% dplyr::select(label, y), by = "label") %>%
  filter(clade %in% PLOT_CLADES) %>%
  group_by(clade) %>%
  summarise(mean_y = mean(y, na.rm = TRUE), .groups = "drop")

# ---- 8d. Annotation strips via gheatmap ----
all_tips <- tree@phylo$tip.label

# Helper: build a matrix aligned to tree tips (rownames = tip labels)
align_to_tips <- function(df) {
  df %>%
    tibble::column_to_rownames("bin_name") %>%
    .[all_tips[all_tips %in% rownames(.)], , drop = FALSE]
}

short_names <- function(cats) abbreviate(cats, minlength = 1)

# All ring widths and offsets are expressed in tree branch-length units.
# gheatmap offset = distance from tips (branch-length units); width = fraction of tree x-range.
gap    <- max_tip_x * 0.02   # small gap between rings (branch-length units)
w_vfdb <- 0.15               # VFDB ring: 15% of tree x-range

o_vfdb <- max_tip_x * 0.02

# Clade labels sit just outside the outermost ring
# o_vfdb is in branch-length units; convert to fraction of max_tip_x before adding
label_x_final <- max_tip_x * (1 + o_vfdb / max_tip_x + w_vfdb + 0.12)

clade_label_df <- clade_nodes_focal %>%
  left_join(clade_tip_angles, by = "clade") %>%
  mutate(x_pos = label_x_final,
         color = clade_colors[clade])

p_base <- p_base +
  new_scale_color() +
  geom_text(
    data        = clade_label_df,
    aes(x = x_pos, y = mean_y, label = short_label, color = color),
    inherit.aes = FALSE,
    fontface    = "bold",
    size        = 6
  ) +
  scale_color_identity()

# --- Ring 1: VFDB top categories (z-scored per bin) ---
heat_vfdb <- tip_meta %>%
  dplyr::select(bin_name, all_of(TOP_VF_CATS)) %>%
  align_to_tips() %>%
  scale() %>%
  as.data.frame()
colnames(heat_vfdb) <- short_names(TOP_VF_CATS)

p1 <- gheatmap(
  p_base, heat_vfdb,
  offset = o_vfdb, width = w_vfdb,
  colnames_angle = 95, colnames_offset_y = 0.25, font.size = 6
) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "VFDB\n(z-score)", na.value = "grey92"
  ) +
  ggtitle("Tree of Alistipes putredinis bins") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

ggsave(
  file.path(results_dir, "alistipes_tree.pdf"),
  p1,
  width  = 10,
  height = 9
)
cat("\nPer-sample tree saved to:", file.path(results_dir, "alistipes_tree.pdf"), "\n")

#### 9. Per-clade summary ####
## 9a. Chi-square: clade × ethnicity association — restricted to PLOT_CLADES
clade_eth_tab <- table(
  clade    = tip_meta_clades$clade[tip_meta_clades$clade %in% PLOT_CLADES],
  ethnicity = tip_meta_clades$EthnicityTot[tip_meta_clades$clade %in% PLOT_CLADES]
)

cat("\nClade × ethnicity table:\n")
print(clade_eth_tab)

chisq_res <- chisq.test(clade_eth_tab)
cat("\nChi-square test (clade × ethnicity):\n")
print(chisq_res)

## 9b. Per-clade functional summary (counts = bins)
clade_summary <- tip_meta_clades %>%
  group_by(clade) %>%
  summarise(
    n_total   = n(),
    n_dutch   = sum(EthnicityTot == "Dutch",                  na.rm = TRUE),
    n_sas     = sum(EthnicityTot == "South-Asian Surinamese", na.rm = TRUE),
    pct_dutch = round(n_dutch / n_total * 100, 1),
    pct_sas   = round(n_sas   / n_total * 100, 1),
    across(all_of(TOP_VF_CATS),
           ~ round(mean(.x, na.rm = TRUE), 5),
           .names = "mean_{.col}"),
    .groups = "drop"
  ) %>%
  arrange(clade)

cat("\nPer-clade summary:\n")
print(clade_summary, n = Inf)

write.csv(clade_summary,
          file.path(results_dir, "clade_functional_summary.csv"),
          row.names = FALSE)

cat("\nClade summary written to:", results_dir, "\n")

# Save tip metadata with clade assignments for downstream scripts
saveRDS(tip_meta_clades, file.path(results_dir, "tip_meta_clades.RDS"))
cat("tip_meta_clades saved to:", file.path(results_dir, "tip_meta_clades.RDS"), "\n")

# Read and save bin quality data with clade + ethnicity for downstream scripts
bin_quality_clade <- map_dfr(batch_files, function(f) {
  read.csv(f, check.names = FALSE) %>%
    dplyr::select(bin, Completeness, Contamination) %>%
    mutate(bin_name = sub("\\.fa$", "", bin)) %>%
    dplyr::select(bin_name, Completeness, Contamination)
}) %>%
  left_join(clade_map %>% dplyr::select(bin_name, clade), by = "bin_name") %>%
  left_join(bin_clin_meta %>% dplyr::select(bin_name, EthnicityTot), by = "bin_name")
saveRDS(bin_quality_clade, file.path(results_dir, "bin_quality_clade.RDS"))
cat("bin_quality_clade saved to:", file.path(results_dir, "bin_quality_clade.RDS"), "\n")

#### 10. Per-clade plots ####
jco_cols <- jco_palette()

## 10a. Clade distribution per ethnicity (panel B — flipped from original)
# Shows: for each ethnicity, what proportion of its bins fall in each clade.
# Chi-square based on present-sample (depth > 0 at baseline) counts.
sample_counts_clade <- tip_meta_clades %>%
  filter(!is.na(EthnicityTot), depth_baseline > 0, clade %in% PLOT_CLADES) %>%
  count(clade, EthnicityTot) %>%
  pivot_wider(names_from = EthnicityTot, values_from = n, values_fill = 0L)

eth_chisq_mat <- sample_counts_clade %>%
  tibble::column_to_rownames("clade") %>%
  as.matrix()

eth_chisq <- chisq.test(eth_chisq_mat)
cat("\nChi-square test (ethnicity × clade, present-sample counts, PLOT_CLADES only):\n")
print(eth_chisq)

# Stacked bar: proportion of each ethnicity per clade, ordered by clade name
eth_clade_long <- tip_meta_clades %>%
  filter(!is.na(EthnicityTot), clade %in% PLOT_CLADES) %>%
  count(clade, EthnicityTot) %>%
  group_by(clade) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    EthnicityTot = factor(EthnicityTot, levels = c("Dutch", "South-Asian Surinamese")),
    clade        = factor(clade, levels = rev(PLOT_CLADES))   # rev so Clade I is at top after coord_flip
  )

eth_pal <- rev(pal_simpsons()(2))

p_eth <- ggplot(eth_clade_long, aes(x = clade, y = prop, fill = EthnicityTot)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = eth_pal) +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() +
  labs(title    = "Ethnicity composition per clade",
       # subtitle = paste0("Chi-square p = ", signif(eth_chisq$p.value, 3)),
       x = "", y = "Proportion of bins", fill = "") +
  theme_Publication()

## 10b. Functional proportions per clade (boxplot, one panel per top VF category)
# Colour/pair helpers shared across functional plots
clade_levels    <- PLOT_CLADES
clade_pairs_all <- combn(clade_levels, 2, simplify = FALSE)
clade_fill_cols <- clade_colors[PLOT_CLADES]

# Restricted to baseline bins only (avoids pseudoreplication from paired timepoints)
func_long <- tip_meta_clades %>%
  filter(timepoint == "baseline") %>%
  dplyr::select(clade, all_of(TOP_VF_CATS)) %>%
  pivot_longer(cols = -clade,
               names_to = "feature", values_to = "proportion")

func_long_filt <- func_long %>%
  filter(clade %in% PLOT_CLADES) %>%
  mutate(clade = factor(clade, levels = PLOT_CLADES))

# Kruskal-Wallis per feature; pairwise Wilcoxon (BH) only where KW p < 0.05
sig_pairs <- map_dfr(unique(func_long_filt$feature), function(feat) {
  sub  <- func_long_filt %>% filter(feature == feat)
  kw_p <- kruskal.test(proportion ~ clade, data = sub)$p.value
  if (kw_p >= 0.05) return(NULL)
  map_dfr(clade_pairs_all, function(pair) {
    d <- sub %>% filter(clade %in% pair)
    p <- wilcox.test(proportion ~ clade, data = d, exact = FALSE)$p.value
    tibble(feature = feat, group1 = pair[1], group2 = pair[2], p_raw = p)
  }) %>%
    mutate(p_adj = p.adjust(p_raw, method = "BH")) %>%
    filter(p_adj < 0.05)
})

cat("\nSignificant pairwise comparisons (Wilcoxon BH, after KW p<0.05):\n")
print(sig_pairs)

# Build one plot per feature so comparisons stay within the right facet,
# then assemble with ggarrange
feat_plots <- purrr::map(unique(func_long_filt$feature), function(feat) {
  df    <- func_long_filt %>% filter(feature == feat)
  pairs <- sig_pairs %>%
    filter(feature == feat) %>%
    { purrr::map2(.$group1, .$group2, c) }

  p <- ggplot(df, aes(x = clade, y = proportion, fill = clade)) +
    geom_boxplot(outlier.size = 0.6, width = 0.6, alpha = 0.7) +
    scale_fill_manual(values = clade_fill_cols, guide = "none") +
    labs(title = feat, x = "", y = "Proportion of annotated genes") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (length(pairs) > 0)
    p <- p + stat_compare_means(comparisons = pairs, method = "wilcox.test",
                                label = "{p.format}", tip.length = 0)
  p
})

## 10c. Horizontal boxplot: baseline depth per clade, Dutch vs SAS
abund_clade <- tip_meta_clades %>%
  filter(!is.na(EthnicityTot), clade %in% PLOT_CLADES, depth_baseline > 0) %>%
  mutate(
    EthnicityTot = factor(EthnicityTot, levels = c("South-Asian Surinamese", "Dutch")),
    clade        = factor(clade, levels = rev(PLOT_CLADES))
  )

wx_abund_clade <- map_dfr(PLOT_CLADES, function(cl) {
  d <- abund_clade %>% filter(clade == cl)
  if (n_distinct(d$EthnicityTot) < 2) return(NULL)
  p <- wilcox.test(depth_baseline ~ EthnicityTot, data = d, exact = FALSE)$p.value
  tibble(clade = cl, p = p)
}) %>%
  mutate(p_label = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  ))

cat("\nWilcoxon baseline depth per clade (Dutch vs SAS):\n")
print(wx_abund_clade)

x_max_abund <- max(abund_clade$depth_baseline, na.rm = TRUE)

p_abund_clade <- ggplot(abund_clade,
                        aes(x = depth_baseline, y = clade, fill = EthnicityTot)) +
  geom_boxplot(outlier.size = 0.5, width = 0.55, alpha = 0.8,
               position = position_dodge(0.65)) +
  geom_text(
    data = wx_abund_clade %>%
      mutate(clade = factor(clade, levels = rev(PLOT_CLADES))),
    aes(x = 100, y = clade, label = p_label),
    inherit.aes = FALSE,
    size = 6, fontface = "bold", hjust = 0
  ) +
  scale_fill_manual(values = rev(jco_cols[1:2]), name = "") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(
    title = "Baseline abundance per clade",
    x     = "Sequencing depth (baseline)",
    y     = ""
  ) +
  theme_Publication() +
  theme(legend.position = "bottom")

# Standalone version: labels A–D (makes sense without the tree)
p_clade_standalone <- ggarrange(
  p_eth,
  p_abund_clade,
  ggarrange(plotlist = feat_plots, ncol = 2, nrow = 1, labels = c("C", "D")),
  nrow    = 3,
  labels  = c("A", "B", ""),
  heights = c(1, 1, 1)
)

ggsave(
  file.path(results_dir, "clade_summary_plots.pdf"),
  p_clade_standalone,
  width  = 10,
  height = 15
)

# Combined version: labels B–E (panel A = tree in tree_withplots.pdf)
p_clade_summary <- ggarrange(
  p_eth,
  p_abund_clade,
  ggarrange(plotlist = feat_plots, ncol = 2, nrow = 1, labels = c("D", "E")),
  nrow    = 3,
  labels  = c("B", "C", ""),
  heights = c(1, 1, 1)
)

cat("\nClade summary plots saved to:", results_dir, "\n")

# Combined: tree (top, full width) + summary panels (bottom row)
# Portrait layout — tree takes top half, summary panels take bottom half
p_combined <- ggarrange(
  p1,
  p_clade_summary,
  nrow    = 2,
  labels  = c("A", ""),
  heights = c(1.4, 1)
)
ggsave(file.path(results_dir, "tree_withplots.pdf"),
  p_combined,
  width  = 12,
  height = 24)

cat("\nCombined figure saved to:", results_dir, "\n")
