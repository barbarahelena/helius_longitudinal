## Figure 1 — Longitudinal microbiome change across ethnicities (16S)
##
## Panels:
##   A — Violin: follow-up time by ethnicity
##   B — PCoA: community composition baseline vs follow-up
##   C — Scatter: Bray-Curtis dissimilarity vs follow-up time
##   D — Violin: Bray-Curtis dissimilarity by ethnicity (confounder-adjusted)
##   E — Violin: Shannon by ethnicity × timepoint
##   F — Scatter: baseline Shannon vs Bray-Curtis dissimilarity (all ethnicities)
##   G — Heatmap: core OTUs across ethnic groups (baseline)

library(ggpubr)
library(ggplotify)

## Source subscripts (generates individual PDFs and defines panel objects) ----
source("scripts/1_longitudinal_change/3_alphadiversity.R")
source("scripts/1_longitudinal_change/2c_ordination.R")
source("scripts/1_longitudinal_change/4_heatmap.R")

## Recreate all ComplexHeatmap objects fresh ----
## Legend uses reference semantics and is mutated when drawn, so reuse causes
ht_base_fig <- make_core_heatmap(base_mats_union$median, base_mats_union$prevalence,
                                  col_fun_abund, row_labels_vec, title = "Baseline",
                                  show_legend = FALSE, name = "Baseline_abund")
ht_fu_fig   <- make_core_heatmap(fu_mats_union$median, fu_mats_union$prevalence,
                                  col_fun_abund, row_labels_vec, title = "Follow-up",
                                  show_legend = TRUE, name = "Median rel.ab. log10")

heatmap_gg <- as.ggplot(function() {
    draw(ht_base_fig + ht_fu_fig,
         annotation_legend_side = "right",
         padding = unit(c(5, 2, 5, 10), "mm"))
})

## Assemble Figure 1
left_panels <- ggarrange(
    ggarrange(pl_fig1_A, pl_fig1_B, widths = c(1, 1.1), labels = c("A", "B")),
    ggarrange(pl_fig1_C, pl_fig1_D, nrow = 1, labels = c("C", "D"), widths = c(1, 1.25)),
    ggarrange(pl_fig1_E, pl_fig1_F, widths = c(1.25, 1), labels = c("E", "F")),
    nrow = 3, heights = c(1, 1, 1.4))

fig1 <- ggarrange(left_panels, heatmap_gg,
                  ncol = 2, widths = c(1, 1),
                  labels = c("", "G"))

ggsave(fig1, filename = "results/1_longitudinal_change/figure1.pdf",
          width = 20, height = 13, device = "pdf")
