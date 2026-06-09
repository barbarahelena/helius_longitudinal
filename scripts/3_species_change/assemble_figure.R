## Figure 3 — Species-level dynamics (shotgun, Dutch vs SAS)
##
## Panels:
##   A — LMM forest plot: species with significant ethnicity x timepoint interaction
##   B — Ethnicity composition per clade (barplot)
##   C — Baseline abundance per clade, Dutch vs SAS (boxplot)
##   D — VFDB boxplot 1
##   E — VFDB boxplot 2
##   F — Alistipes phylogenetic tree (full width)
##   G — Strain stability vs follow-up time
##   H — Strain stability vs Shannon index (baseline)
##   I — Strain stability vs Bray-Curtis dissimilarity
##   J — Strain sharing percentage comparison between Dutch and SAS
##
## ML AUROC panels (formerly A & B) moved to suppl_figure_ml_auroc.pdf

library(ggpubr)
library(ggplotify)

## Source subscripts (run analyses and define panel objects) ----
# source("scripts/3_species_change/2_mlmodels/3_ml_processing/process_eth_time_models_shotgun.R")  # produces suppl_figure_ml_auroc.pdf
source("scripts/3_species_change/3_species/1_lmm_sg.R")                                          # defines pl_fig3_C                                        # defines pl_fig3_D (heatmap)
source("scripts/3_species_change/4_alistipes_anno/3_draw_tree.R")                                 # defines p1 (tree), p_eth (barplot), feat_plots (boxplots)
source("scripts/3_species_change/5_strain_stability/strainsharing_plot.R")                        # defines pl_fig3_G, pl_fig3_H, pl_fig3_I, pl_fig3_J, pl_fig3_K

## Assemble Figure 3 ----
# Stack ethnicity bar + abundance boxplot into one column
eth_abund_col <- ggarrange(
  p_eth, p_abund_clade,
  nrow   = 2,
  labels = c("B", "C"),
  common.legend = TRUE,
  legend = "bottom"
)

# Panel A: aplot combines forest + heatmap with proper y-axis alignment;
# convert to grob so ggarrange can place it
panel_A <- as_ggplot(as.grob(pl_forest_heatmap))

# Row 1: [forest + heatmap] | [eth bar / abundance] | VF boxplot 1 | VF boxplot 2
top_row <- ggarrange(
  panel_A, eth_abund_col, feat_plots[[1]], feat_plots[[2]],
  ncol   = 4,
  widths = c(2.5, 1.2, 1, 1),
  labels = c("A", "", "D", "E")
)

# Row 2: tree alone, full width
tree_row <- ggarrange(p1, labels = "F")

fig3 <- ggarrange(
  top_row, tree_row,
  nrow    = 2,
  heights = c(1.0, 1.6)
)

dir.create("results/3_species_change", showWarnings = FALSE, recursive = TRUE)
ggsave(fig3, filename = "results/3_species_change/figure3.pdf",
       width = 18, height = 16, device = cairo_pdf)

## Supplementary figure: strain stability panels (G–J) ----
suppl_strain <- ggarrange(
  pl_fig3_H, pl_fig3_G, pl_fig3_I, pl_fig3_J,
  nrow   = 2,
  ncol   = 2,
  labels = c("A", "B", "C", "D")
)

ggsave(suppl_strain, filename = "results/3_species_change/suppl_figure_strain_stability.pdf",
       width = 12, height = 14, device = cairo_pdf)
