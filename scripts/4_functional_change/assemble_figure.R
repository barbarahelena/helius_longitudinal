## Figure 4 — Functional shifts: HUMAnN pathways & CAZyme families
##
## Panels:
##   A  HUMAnN forest + cross-sectional heatmap (FDR < 0.05)
##   B  HUMAnN violin: Gluconeogenesis III
##   C  HUMAnN violin: L-lysine biosynthesis II
##   D  CAZyme Canberra PCoA — baseline by ethnicity
##   (–) CAZyme Canberra PCoA — follow-up by ethnicity
##   E  CAZyme forest + cross-sectional heatmap (FDR < 0.05)
##   F  CAZyme violin: GH13_16
##   G  CAZyme violin: GH5_26
##   H  GAG-to-DF ratio

library(ggpubr)
library(tidyverse)
library(ggsci)

## ── Helper: render an aplot composite to a ggplot panel ───────────────────────
## patchwork combines forest + heatmap with guides="keep" so each panel's legend
## stays in its configured position; grid.grabExpr captures the rendered grob.

aplot_to_gg <- function(ap) {
  if (is.null(ap)) return(ggplot() + theme_void())
  tryCatch(
    ggpubr::as_ggplot(grid::grid.grabExpr({
      grid::grid.newpage()
      print(ap)
    }, warn = 0)),
    error = function(e) {
      message("aplot conversion failed: ", conditionMessage(e))
      ggplot() + theme_void()
    }
  )
}

## ── 1. Source HUMAnN — snapshot before CAZyme overwrites shared names ─────────

if (exists("pl_combined")) rm(pl_combined)
source("scripts/4_functional_change/humann/2_lmm_humann.R")

pl_humann_combined <- if (exists("pl_combined")) pl_combined else NULL
df_clin_humann     <- df_clin
statres_humann     <- statres
pseudocount_humann <- pseudocount

## ── 2. Source CAZyme LMM ──────────────────────────────────────────────────────

if (exists("pl_combined")) rm(pl_combined)
source("scripts/4_functional_change/cayman/3_cayman_longitudinal_lmm.R")

pl_cazyme_combined <- if (exists("pl_combined")) pl_combined else NULL
# make_boxviolin (CAZyme version), dftot_adj, statres_adj remain in environment

## ── 3. Source CAZyme ordination — snapshot Canberra PCoA panels ───────────────

source("scripts/4_functional_change/cayman/4_cayman_ordination.R")

pl_can_bl <- if (exists("ethcan_bl")) ethcan_bl else ggplot() + theme_void()
pl_can_fu <- if (exists("ethcan_fu")) ethcan_fu else ggplot() + theme_void()

## ── 4. Source CAZyme ratios — build GAG violin immediately ────────────────────

source("scripts/4_functional_change/cayman/5_cayman_ratios.R")
# dftot (ratios), wilcox_res now in environment

p_gag_vln <- make_ratio_vln(
  dftot, "log10_GAG_DF", "GAG / DF (log₁₀ ratio)", "GAG-to-DF ratio",
  lmm_df$pval_interact[lmm_df$ratio == "log10_GAG_DF"],
  show_pval = FALSE
) + theme(plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"))

## ── 5. Convert aplot composites → ggplot panels ───────────────────────────────

pl_A <- aplot_to_gg(pl_humann_combined)
pl_E <- aplot_to_gg(pl_cazyme_combined)

## ── 6. Build named HUMAnN violin panels (B, C) ────────────────────────────────

find_col <- function(partial, cols) {
  m <- grep(partial, cols, value = TRUE, ignore.case = TRUE)
  if (length(m) == 0) stop("Column not found: ", partial)
  m[1]
}

make_humann_vln <- function(partial_name) {
  nm         <- find_col(partial_name, names(df_clin_humann))
  pval       <- statres_humann$pval[statres_humann$pathway == nm]
  if (length(pval) == 0) pval <- NA_real_
  nm_clean   <- str_trunc(str_remove(nm, "^[A-Z0-9_-]+: "), 40)
  pval_label <- formatC(pval, format = "e", digits = 2)
  df_tmp     <- df_clin_humann
  df_tmp$mb  <- log10(df_tmp[[nm]] * 100 + pseudocount_humann)
  ggplot(df_tmp, aes(x = timepoint, y = mb, fill = EthnicityTot)) +
    geom_violin(colour = NA, aes(alpha = timepoint)) +
    geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
    facet_wrap(~EthnicityTot) +
    scale_fill_jco(guide = "none") +
    scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
    theme_Publication() +
    labs(x        = "",
         y        = "log10(abundance % + pseudocount)",
         title    = nm_clean,
         subtitle = paste0("Ethnicity × Timepoint p=", pval_label)) +
    theme(plot.title    = element_text(face = "bold", size = rel(0.75), hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"))
}

pl_B <- make_humann_vln("Gluconeogenesis III")
pl_C <- make_humann_vln("L-lysine biosynthesis II")

## ── 7. Build named CAZyme violin panels (F, G) ────────────────────────────────

get_cazyme_pval <- function(fam) {
  row <- statres_adj[statres_adj$family == fam, ]
  if (nrow(row) == 0) NA_real_ else row$pval[1]
}

make_cazyme_vln <- function(fam_name) {
  make_boxviolin(dftot_adj, fam_name, get_cazyme_pval(fam_name), "log10(CPM + 1)",
                 show_pval = FALSE) +
    labs(title = fam_name) +
    theme(plot.title    = element_text(face = "bold", size = rel(0.9), hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"))
}

pl_F <- make_cazyme_vln("GH13_16")
pl_G <- make_cazyme_vln("GH5_26")

## ── 8. Assemble rows ──────────────────────────────────────────────────────────

# Row 1 — HUMAnN forest+heatmap (A) + two CAZyme PCoAs (B)
top_row <- ggarrange(
  pl_A, pl_can_bl, pl_can_fu,
  ncol   = 3,
  labels = c("A", "B", ""),
  widths = c(2, 1, 1)
)

# Row 2 — CAZyme forest+heatmap (C) + GAG ratio violin (D)
bottom_row <- ggarrange(
  pl_E, ggarrange(p_gag_vln, NULL, nrow = 2, heights = c(1.5, 0.5)),
  ncol   = 2,
  labels = c("C", "D"),
  widths = c(2, 1)
)

## ── 9. Final assembly ─────────────────────────────────────────────────────────

fig4 <- ggarrange(
  top_row,
  bottom_row,
  nrow    = 2,
  heights = c(1, 1.5)
)

## ── 10. Save ──────────────────────────────────────────────────────────────────

dir.create("results/4_functional_change", showWarnings = FALSE, recursive = TRUE)
ggsave(
  fig4,
  filename = "results/4_functional_change/figure4.pdf",
  width    = 18,
  height   = 12,
  device   = cairo_pdf
)

## ── 11. Supplementary figure: violin panels B, C, F, G ───────────────────────

suppl_fig4 <- ggarrange(
  pl_B, pl_C, pl_F, pl_G,
  nrow   = 2,
  ncol   = 2,
  labels = c("A", "B", "C", "D")
)

ggsave(
  suppl_fig4,
  filename = "results/4_functional_change/suppl_figure_violins.pdf",
  width    = 12,
  height   = 12,
  device   = cairo_pdf
)
