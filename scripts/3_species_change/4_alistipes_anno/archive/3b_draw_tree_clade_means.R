## ARCHIVE — Clade-mean heatmap tree for Alistipes putredinis bins
## Each tip is coloured by its clade's mean VFDB z-score rather than its own value.
## Moved here from 3_draw_tree.R because the per-tip heatmap (p1) is preferred.
##
## Requires: run 3_draw_tree.R first so that tip_meta_clades.RDS exists in results_dir.
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

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
results_dir <- "results/3_species_change/4_alistipes_anno"

#### Constants ####
MIN_PIDENT  <- 30
MIN_BITSCORE <- 50
PLOT_CLADES  <- c("Clade I", "Clade II", "Clade III", "Clade IV")
N_TOP_VF     <- 2

clade_pal <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
               "#76B7B2", "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC")

#### Load data saved by 3_draw_tree.R ####
tip_meta_clades <- readRDS(file.path(results_dir, "tip_meta_clades.RDS"))

# Reconstruct clade_colors
all_clades   <- sort(unique(tip_meta_clades$clade))
clade_colors <- setNames(clade_pal[seq_along(all_clades)], all_clades)

# TOP_VF_CATS — read from saved functional summary to stay in sync
clade_summary <- read.csv(file.path(results_dir, "clade_functional_summary.csv"),
                          check.names = FALSE)
TOP_VF_CATS <- grep("^mean_", names(clade_summary), value = TRUE) %>%
  sub("^mean_", "", .)

cat("Using TOP_VF_CATS:", paste(TOP_VF_CATS, collapse = ", "), "\n")

#### Rebuild tip metadata (needed for tree layout) ####
tree <- read.iqtree(tree_file)

# Re-derive cluster assignments
H_CUT    <- 0.02
hc           <- hclust(as.dist(ape::cophenetic.phylo(tree@phylo)), method = "average")
tip_clusters <- cutree(hc, h = H_CUT)

clade_nodes <- map_dfr(sort(unique(tip_clusters)), function(cl) {
  tips <- names(tip_clusters)[tip_clusters == cl]
  nd <- if (length(tips) == 1) which(tree@phylo$tip.label == tips) else
    ape::getMRCA(tree@phylo, tips)
  tibble(node = nd, n_tips = length(tips), cluster = cl,
         clade = paste0("Clade ", as.roman(cl)))
}) %>%
  filter(!is.na(node)) %>%
  arrange(node)

clade_colors <- setNames(clade_pal[seq_len(nrow(clade_nodes))], clade_nodes$clade)

clade_map <- tibble(bin_name = names(tip_clusters), cluster = tip_clusters) %>%
  left_join(clade_nodes %>% dplyr::select(cluster, clade), by = "cluster")

#### Re-build p_base (minimal, just for gheatmap) ####
eth_colors <- jco_palette()

p_base <- ggtree(tree, layout = "circular", size = 0.35)
tree_layout <- fortify(tree)
tip_layout  <- tree_layout %>% filter(isTip)
max_tip_x   <- max(tip_layout$x, na.rm = TRUE)

# Clade highlights (all clades)
clade_y_ranges <- tibble(label = names(tip_clusters), cluster = tip_clusters) %>%
  left_join(clade_nodes %>% dplyr::select(cluster, clade), by = "cluster") %>%
  left_join(tip_layout %>% dplyr::select(label, y), by = "label") %>%
  filter(!is.na(clade)) %>%
  group_by(clade) %>%
  summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop")

for (i in seq_len(nrow(clade_y_ranges))) {
  p_base <- p_base +
    annotate("rect",
             xmin = 0, xmax = max_tip_x * 1.006,
             ymin = clade_y_ranges$ymin[i], ymax = clade_y_ranges$ymax[i],
             fill = clade_colors[clade_y_ranges$clade[i]], alpha = 0.7)
}

tip_meta <- tip_meta_clades   # already has all columns needed for %<+%

p_base <- p_base %<+%
  (tip_meta %>% rename(label = bin_name)) +
  geom_tippoint(aes(color = EthnicityTot, shape = timepoint), size = 2.5, na.rm = TRUE) +
  scale_color_manual(values = eth_colors, name = "Ethnicity", na.value = "grey70") +
  scale_shape_manual(values = c("baseline" = 16, "follow-up" = 17),
                     name = "Timepoint", na.value = 1) +
  theme_tree() +
  theme(legend.position = "right", legend.box = "vertical",
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

# Clade labels
clade_nodes_focal <- clade_nodes %>%
  filter(clade %in% PLOT_CLADES) %>%
  arrange(clade) %>%
  mutate(short_label = as.character(as.roman(cluster)))

clade_tip_angles <- tibble(label = names(tip_clusters), cluster = tip_clusters) %>%
  left_join(clade_nodes %>% dplyr::select(cluster, clade), by = "cluster") %>%
  left_join(tip_layout %>% dplyr::select(label, y), by = "label") %>%
  filter(clade %in% PLOT_CLADES) %>%
  group_by(clade) %>%
  summarise(mean_y = mean(y, na.rm = TRUE), .groups = "drop")

w_vfdb        <- 0.15
o_vfdb        <- max_tip_x * 0.02
label_x_final <- max_tip_x * (1 + o_vfdb / max_tip_x + w_vfdb + 0.12)

clade_label_df <- clade_nodes_focal %>%
  left_join(clade_tip_angles, by = "clade") %>%
  mutate(x_pos = label_x_final, color = clade_colors[clade])

p_base <- p_base +
  new_scale_color() +
  geom_text(data = clade_label_df,
            aes(x = x_pos, y = mean_y, label = short_label, color = color),
            inherit.aes = FALSE, fontface = "bold", size = 4) +
  scale_color_identity()

#### Clade-mean heatmap ####
all_tips     <- tree@phylo$tip.label
align_to_tips <- function(df) {
  df %>% tibble::column_to_rownames("bin_name") %>%
    .[all_tips[all_tips %in% rownames(.)], , drop = FALSE]
}
short_names <- function(cats) abbreviate(cats, minlength = 1)

tip_meta_focal <- tip_meta_clades %>% filter(timepoint == "baseline")

clade_vf_means <- tip_meta_focal %>%
  group_by(clade) %>%
  summarise(across(all_of(TOP_VF_CATS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

heat_vfdb_clade <- tip_meta_clades %>%
  left_join(clade_vf_means, by = "clade") %>%
  dplyr::select(bin_name, all_of(paste0(TOP_VF_CATS, ".y"))) %>%
  rename_with(~ str_remove(.x, "\\.y$"), ends_with(".y")) %>%
  align_to_tips() %>%
  scale() %>%
  as.data.frame()

colnames(heat_vfdb_clade) <- short_names(TOP_VF_CATS)

p2 <- gheatmap(
  p_base, heat_vfdb_clade,
  offset = o_vfdb, width = w_vfdb,
  colnames_angle = 95, colnames_offset_y = 0.25, font.size = 3,
  color = NA
) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "VFDB clade mean\n(z-score)", na.value = "grey92"
  )

ggsave(
  file.path(results_dir, "alistipes_tree_clade_means.pdf"),
  p2,
  width  = 10,
  height = 9
)
cat("Clade-mean tree saved to:", file.path(results_dir, "alistipes_tree_clade_means.pdf"), "\n")
