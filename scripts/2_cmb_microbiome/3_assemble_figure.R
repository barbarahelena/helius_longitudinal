## Figure 2 — Cardiometabolic disease and microbiome instability (16S)
##
## Panels:
##   A — Combined prevalence: DM, HT, Dyslipidemia, MetSyn at baseline + follow-up by ethnicity
##   B — Extended forest plot: multi-domain predictors of Bray-Curtis instability
##
## Note: alpha-diversity panels are assembled as a supplementary figure
##       in scripts/2_cmb_microbiome/2_alphadiversity_cmb_16s.R

library(ggpubr)

## Source beta-diversity subscript ----
source("scripts/2_cmb_microbiome/1_betadiversity_cmb_16s.R")   # defines pl_fig2_prev, pl_bc_extended

## Assemble Figure 2 ----
fig2 <- ggarrange(
    pl_fig2_prev,
    pl_bc_combined,
    nrow   = 2,
    labels = c("A", "B"),
    heights = c(1, 2.5)
)

ggsave(fig2, filename = "results/2_cmb_microbiome/figure2.pdf",
       width = 11, height = 12, device = "pdf")
