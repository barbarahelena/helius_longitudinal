## Figure 4 — Functional shifts: CAZymes & metabolic pathways (shotgun)
##
## Panels:
##   A — CAZyme forest plot: 54 families with ethnicity × timepoint interactions
##         (adjusted LMM, all p < 0.05; coloured by FDR significance)
##   B — GH13, GH5, GH32 emmeans trajectories
##         (predicted marginal means ± 95 % CI by ethnicity over time)
##   C — GutSMASH metabolic pathway interactions
##         (forest plot of FDR-significant ethnicity × timepoint shifts;
##          falls back to top-3 trajectory panels if < 2 FDR-sig pathways)

library(ggpubr)
library(tidyverse)

## ── 1. Source analysis scripts ───────────────────────────────────────────────

source("scripts/4_functional_change/cayman/3_cayman_longitudinal_lmm.R")

# Stash CAZyme objects before gutsmash sourcing can overwrite shared names
pl_4A_base      <- p_forest_adj       # forest: adjusted LMM, all p < 0.05
cayman_emm      <- plist_emm          # list of emmeans panels per family
cayman_sig      <- statres_adj_sig    # data.frame: $family column gives order

source("scripts/4_functional_change/gutsmash/3_gutsmash_longitudinal_lmm.R")

gutsmash_plist   <- plist             # trajectory plots for FDR-sig pathways
gutsmash_statres <- statres           # full LMM results (family, estimate, ci, padj)

## ── 2. Panel A: CAZyme forest plot ───────────────────────────────────────────
pl_4A <- pl_4A_base +
  labs(
    title    = "Differential CAZyme dynamics by ethnicity",
    subtitle = "Adjusted LMM: ethnicity × timepoint interaction (p < 0.05)\nAdjusted for baseline age, sex, BMI, smoking, PPI"
  )

## ── 3. Panel B: GH13, GH5, GH32 emmeans trajectories ────────────────────────
target_fam <- c("GH13", "GH5", "GH32")
emm_idx    <- match(target_fam, cayman_sig$family)

# For any family not found in top-20 emmeans list, fall back to next available
emm_idx_clean <- mapply(function(idx, fallback) {
  if (!is.na(idx) && idx <= length(cayman_emm)) idx
  else min(fallback, length(cayman_emm))
}, emm_idx, seq_along(target_fam))

pl_B_panels <- lapply(seq_along(target_fam), function(k) {
  idx  <- emm_idx_clean[k]
  fam  <- target_fam[k]
  base <- cayman_emm[[idx]]
  # Override title to show the requested family name prominently
  base +
    labs(title = fam) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = rel(0.9), hjust = 0.5))
})

pl_4B <- ggarrange(
  plotlist      = pl_B_panels,
  ncol          = 3,
  common.legend = TRUE,
  legend        = "bottom"
)

## ── 4. Panel C: GutSMASH pathway interactions ────────────────────────────────
gs_fdr <- gutsmash_statres %>%
  filter(padj < 0.05) %>%
  mutate(
    label = str_remove_all(family, "[\\[\\]']"),
    label = str_trunc(label, 55),
    label = fct_reorder(label, estimate)
  )

if (nrow(gs_fdr) >= 2) {
  # Forest plot (mirrors 4A style)
  pl_4C <- ggplot(gs_fdr, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(
      aes(xmin = conflow, xmax = confhigh),
      color = "#E64B35", size = 0.45, linewidth = 0.55
    ) +
    theme_Publication() +
    labs(
      x        = "Estimate (Ethnicity × Timepoint interaction)",
      y        = "",
      title    = "GutSMASH metabolic pathway shifts",
      subtitle = "LMM ethnicity × timepoint, FDR < 0.05"
    )
} else {
  # Fallback: up to 3 top trajectory panels
  n_gs   <- min(3, length(gutsmash_plist))
  pl_4C  <- ggarrange(
    plotlist      = gutsmash_plist[seq_len(n_gs)],
    ncol          = n_gs,
    common.legend = TRUE,
    legend        = "bottom"
  )
}

## ── 5. Assemble Figure 4 ─────────────────────────────────────────────────────
#
#  Layout (landscape, 18 × 14 in):
#
#   ┌──────────────────┬─────────────────────────────────────┐
#   │                  │  B: GH13 | GH5 | GH32              │
#   │  A: CAZyme       │  emmeans trajectories               │
#   │  forest plot     ├─────────────────────────────────────┤
#   │  (p < 0.05)      │  C: GutSMASH pathway interactions   │
#   └──────────────────┴─────────────────────────────────────┘

right_col <- ggarrange(
  pl_4B, pl_4C,
  nrow    = 2,
  labels  = c("B", "C"),
  heights = c(1.0, 1.0)
)

fig4 <- ggarrange(
  pl_4A,
  right_col,
  ncol   = 2,
  labels = c("A", ""),
  widths = c(1.0, 1.8)
)

dir.create("results/4_functional_change", showWarnings = FALSE, recursive = TRUE)
ggsave(
  fig4,
  filename = "results/4_functional_change/figure4.pdf",
  width    = 18,
  height   = 14,
  device   = cairo_pdf
)
