## Figure 1 — Longitudinal microbiome change across ethnicities (16S)
##
## Panels:
##   A — Violin: follow-up time by ethnicity
##   B — PCoA: community composition baseline vs follow-up
##   C — Scatter: Bray-Curtis dissimilarity vs follow-up time
##   D — Violin: Bray-Curtis dissimilarity by ethnicity (unadjusted)
##   E — Violin: Bray-Curtis dissimilarity by ethnicity (confounder-adjusted)
##   F — Violin: Shannon by ethnicity × timepoint
##   G — Scatter: baseline Shannon vs Bray-Curtis dissimilarity (all ethnicities)

library(ggpubr)

## Source subscripts (generates individual PDFs and defines panel objects) ----
source("scripts/1_longitudinal_change/3_alphadiversity.R")
source("scripts/1_longitudinal_change/2c_ordination.R")

## Assemble Figure 1
fig1 <- ggarrange(ggarrange(pl_fig1_A, pl_fig1_B, widths = c(1.3, 1), labels = c("A", "B")),
                  ggarrange(pl_fig1_C, pl_fig1_D, pl_fig1_E, nrow = 1, labels = c("C", "D", "E")),
                  ggarrange(pl_fig1_F, pl_fig1_G, widths = c(1.25, 1), labels = c("F", "G")),
                  nrow = 3, heights = c(1, 1, 1.25))

ggsave(fig1, filename = "results/1_longitudinal_change/figure1.pdf", 
          width = 13, height = 15, device = "pdf")
