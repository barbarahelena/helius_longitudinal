## Supplementary Figure 6 — Shotgun longitudinal microbiome change (mirrors Fig 1 A–F)
##
## Panels:
##   A — Violin: follow-up time by ethnicity (shotgun)
##   B — PCoA: community composition baseline vs follow-up, per ethnicity (shotgun)
##   C — Scatter: Bray-Curtis dissimilarity vs follow-up time (shotgun)
##   D — Violin: Bray-Curtis dissimilarity by ethnicity (confounder-adjusted, shotgun)
##   E — Violin: Shannon by ethnicity × timepoint (shotgun)
##   F — Scatter: baseline Shannon vs Bray-Curtis dissimilarity (shotgun)
##
## Panels A–D are defined at the end of 1_ordination_sg.R
## Panels E–F are defined at the end of 2_alphadiversity_sg.R

library(ggpubr)

## Source panel scripts ----
source("scripts/3_species_change/1_comparison_16s/1_ordination_sg.R")
source("scripts/3_species_change/1_comparison_16s/2_alphadiversity_sg.R")

## Assemble ----
suppl_fig6 <- ggarrange(
    ggarrange(pl_sfig6_A, pl_sfig6_B, labels = c("A", "B"), widths = c(1, 1.5)),
    ggarrange(pl_sfig6_C, pl_sfig6_D, nrow = 1, labels = c("C", "D"), widths = c(1, 1.25)),
    ggarrange(pl_sfig6_E, pl_sfig6_F, widths = c(1.25, 1), labels = c("E", "F")),
    nrow = 3, heights = c(1.4, 1.4, 1.4)
)
resultsfolder <- "results/3_species_change/1_comparison_16s"
dir.create(resultsfolder, showWarnings = FALSE, recursive = TRUE)
ggsave(suppl_fig6,
       filename = file.path(resultsfolder, "suppl_fig6.pdf"),
       width = 13, height = 15, device = "pdf")
