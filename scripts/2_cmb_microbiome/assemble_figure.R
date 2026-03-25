## Figure 2 — Cardiometabolic disease and microbiome instability (16S)
##
## Panels:
##   A — Combined prevalence: DM, HT, LLD at baseline + follow-up by ethnicity (faceted bar)
##   B — Overall Bray-Curtis effect sizes: DM, HTN, Dyslipidemia (dot-whisker)
##   C — Per-ethnicity DM effect on Bray-Curtis (interaction forest plot)
##   D — Per-ethnicity HT effect on Bray-Curtis (interaction forest plot)
##   E — Per-ethnicity LLD effect on Bray-Curtis (interaction forest plot)
##   F — Overall Shannon effect sizes: DM, HTN, Dyslipidemia (dot-whisker)
##   G — Per-ethnicity DM effect on Shannon (interaction forest plot)
##   H — Per-ethnicity HT effect on Shannon (interaction forest plot)
##   I — Per-ethnicity LLD effect on Shannon (interaction forest plot)

library(ggpubr)

## Source subscripts (generates supplementary PDFs and defines panel objects) ----
source("scripts/2_cmb_microbiome/betadiversity_cmb_16s.R")    # defines pl_fig2_prev, pl_fig2_A, pl_fig2_Aint, pl_fig2_Aint_ht, pl_fig2_Aint_lld
source("scripts/2_cmb_microbiome/alphadiversity_cmb_16s.R")   # defines pl_shan_A, pl_shan_Aint_dm, pl_shan_Aint_ht, pl_shan_Aint_lld

## Assemble Figure 2 ----
bray_row <- ggarrange(pl_fig2_A, pl_fig2_Aint, pl_fig2_Aint_ht, pl_fig2_Aint_lld,
                      ncol = 4, labels = c("B", "C", "D", "E"))
shan_row <- ggarrange(pl_shan_A, pl_shan_Aint_dm, pl_shan_Aint_ht, pl_shan_Aint_lld,
                      ncol = 4, labels = c("F", "G", "H", "I"))

fig2 <- ggarrange(ggarrange(pl_fig2_prev, NULL, widths = c(1.0, 0.3)), bray_row, shan_row,
                  nrow = 3, heights = c(1.0, 1.0, 1.0), labels = c("A", "", ""))

ggsave(fig2, filename = "results/2_cmb_microbiome/figure2.pdf",
       width = 18, height = 14, device = "pdf")
