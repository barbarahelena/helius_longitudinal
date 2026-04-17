## Figure 3 — Species-level dynamics (shotgun, Dutch vs SAS)
##
## Panels:
##   A — ROC ethnicity prediction at baseline
##   B — ROC ethnicity prediction at follow-up
##   C — LMM forest plot: species with significant ethnicity x timepoint interaction
##   D — Heatmap: species-clinical variable associations
##   E — KEGG pathways boxplots, Alistipes putredinis, baseline only
##   F — VFDB category boxplots, Alistipes putredinis, baseline only
##   G — Strain stability vs Shannon index (baseline)
##   H — Strain stability vs follow-up time
##   I — Strain stability vs Bray-Curtis dissimilarity
##   J — Strain sharing percentage comparison between Dutch and SAS
##   K — Strain sharing connected dot plot + abundance density + n subjects

library(ggpubr)

## Source subscripts (run analyses and define panel objects) ----
source("scripts/3_species_change/2_mlmodels/3_ml_processing/process_eth_time_models_shotgun.R")  # defines pl_fig3_A, pl_fig3_B
source("scripts/3_species_change/3_species/1_lmm_sg.R")                                          # defines pl_fig3_C
source("scripts/3_species_change/3_species/2_heatmap.R")                                         # defines pl_fig3_D (heatmap)
source("scripts/3_species_change/4_alistipes_anno/2_eggnog_comparison.R")                        # defines pl_fig3_E
source("scripts/3_species_change/4_alistipes_anno/3_vfdb_comparison.R")                          # defines pl_fig3_F
source("scripts/3_species_change/5_strain_stability/strainsharing_plot.R")                       # defines pl_fig3_G, pl_fig3_H, pl_fig3_I, pl_fig3_J, pl_fig3_K

## Assemble Figure 3 ----
top_row    <- ggarrange(pl_fig3_A, pl_fig3_B, pl_fig3_C,
                        ncol = 3, widths = c(1.0, 1.0, 1.6),
                        labels = c("A", "B", "C"))
mid_row    <- ggarrange(pl_fig3_D, 
                        ggarrange(pl_fig3_E, ggarrange(pl_fig3_F, NULL, labels = c("F", "")),
                            nrow = 2),
                        ncol = 2, widths = c(2.0, 2.0),
                        labels = c("D", "E"))

scatter_row <- ggarrange(pl_fig3_H, pl_fig3_G, pl_fig3_I, pl_fig3_J,
                         nrow = 1, labels = c("G", "H", "I", "J"))

fig3 <- ggarrange(top_row, mid_row, scatter_row,
                  nrow = 3, heights = c(0.6, 0.7, 0.5))

dir.create("results/3_species_change", showWarnings = FALSE, recursive = TRUE)
ggsave(fig3, filename = "results/3_species_change/figure3.pdf",
       width = 15, height = 18, device = cairo_pdf)
