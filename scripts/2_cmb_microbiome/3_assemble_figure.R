## Figure 2 — Cardiometabolic disease and microbiome instability (16S)
##
## Panels:
##   A — Combined prevalence: DM, HT, LLD at baseline + follow-up by ethnicity (faceted bar)
##   B — Overall Bray-Curtis effect sizes: DM, HTN, Dyslipidemia (dot-whisker)
##   C — Overall Shannon effect sizes: DM, HTN, Dyslipidemia (dot-whisker)
##   D — Bray-Curtis by diabetes status, all 6 ethnicities (violin, ordered by effect size)
##   E — Bray-Curtis by dyslipidemia status, all 6 ethnicities (violin, same ethnicity order)

library(ggpubr)

## Source subscripts (generates supplementary PDFs and defines panel objects) ----
source("scripts/2_cmb_microbiome/1_betadiversity_cmb_16s.R")    # defines pl_fig2_prev, pl_fig2_A, pl_fig2_D_violin, pl_fig2_E_violin
source("scripts/2_cmb_microbiome/2_alphadiversity_cmb_16s.R")   # defines pl_shan_A

## Assemble Figure 2 ----
bc_col  <- ggarrange(pl_fig2_A, pl_shan_A, nrow = 2, labels = c("B", "C"), common.legend = TRUE, legend = "bottom")
row1    <- ggarrange(pl_fig2_prev, bc_col, ncol = 2, widths = c(2, 1),
                     labels = c("A", ""))

fig2 <- ggarrange(
    row1,
    pl_fig2_D_violin,
    pl_fig2_E_violin,
    nrow = 3,
    heights = c(1.0, 1.0, 1.0),
    labels = c("", "D", "E")
)

ggsave(fig2, filename = "results/2_cmb_microbiome/figure2.pdf",
       width = 12, height = 13, device = "pdf")
