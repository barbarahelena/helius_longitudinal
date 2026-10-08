## Figure 4 — Functional shifts: HUMAnN pathways & CAZyme families
##
## Panels (as rendered in figure4.pdf):
##   A  HUMAnN forest + cross-sectional heatmap (FDR < 0.05)
##   B  CAZyme Canberra PCoA — baseline by ethnicity
##   (–) CAZyme Canberra PCoA — follow-up by ethnicity (unlabelled, below B)
##   C  CAZyme forest + cross-sectional heatmap (FDR < 0.05)
##   D  GAG-to-DF ratio
##   E  Mucin-to-DF ratio
##   F  CAZyme family richness
##
## Supplementary violins (suppl_figure_violins.pdf):
##   A  HUMAnN violin: Gluconeogenesis III
##   B  HUMAnN violin: L-lysine biosynthesis II
##   C+ CAZyme violins: top FDR-significant families by padj (target_fam,
##      computed in 3_cayman_longitudinal_lmm.R) — not hardcoded, so this
##      tracks whichever families are currently most significant

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

## ── 4. Source CAZyme ratios — reuse its Mucin/GAG violin panels ───────────────

# Drop objects from earlier runs so a failure inside the ratios script errors
# below instead of silently reusing a stale qc_plots (e.g. an old richness panel)
if (exists("qc_plots")) rm(qc_plots)
source("scripts/4_functional_change/cayman/5_cayman_ratios.R")
# dftot (ratios), wilcox_res now in environment; p_mucin, p_gag built with
# show_pval = TRUE by default (Dutch-vs-SAS Wilcoxon p-value shown per timepoint facet)

p_mucin_vln <- p_mucin + theme(plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"))
p_gag_vln   <- p_gag   + theme(plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"))

# Family richness violin, built in 5_cayman_ratios.R alongside the technical checks
p_rich_vln  <- qc_plots[[match("richness", qc_vars)]] +
  theme(plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"))

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
    scale_alpha_manual(values = c("baseline" = 0.60, "follow-up" = 0.90), guide = "none") +
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

## ── 7. Build CAZyme violin panels for the top FDR-significant families ───────
## Uses target_fam (top 3 by padj among the top-20-by-|estimate| set), already
## computed in 3_cayman_longitudinal_lmm.R — avoids hardcoding family names
## that drift out of date as the underlying data/results change.

get_cazyme_pval <- function(fam) {
  row <- statres_adj[statres_adj$family == fam, ]
  if (nrow(row) == 0) NA_real_ else row$pval[1]
}

make_cazyme_vln <- function(fam_name) {
  make_boxviolin(dftot_adj, fam_name, get_cazyme_pval(fam_name), "log10(RPKM + 1)",
                 show_pval = FALSE) +
    labs(title = fam_name) +
    theme(plot.title    = element_text(face = "bold", size = rel(0.9), hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"))
}

pl_cazyme_top <- lapply(target_fam, make_cazyme_vln)

## ── 8. Assemble rows ──────────────────────────────────────────────────────────

# Row 1 — HUMAnN forest+heatmap (A) + two CAZyme PCoAs stacked (B)
pcoa_col <- ggarrange(pl_can_bl, pl_can_fu, nrow = 2, labels = c("B", ""))

top_row <- ggarrange(
  pl_A, pcoa_col,
  ncol   = 2,
  labels = c("A", ""),
  widths = c(1.8, 1)
)

# Row 2 — CAZyme forest+heatmap (C) + GAG/DF, Mucin/DF ratio and richness violins (D, E, F)
ratio_col <- ggarrange(p_gag_vln, p_mucin_vln, p_rich_vln, nrow = 3, heights = c(1, 1, 1),
                       labels = c("D", "E", "F"))

bottom_row <- ggarrange(
  pl_E, ratio_col,
  ncol   = 2,
  labels = c("C", ""),
  widths = c(2, 1)
)

## ── 9. Final assembly ─────────────────────────────────────────────────────────

fig4 <- ggarrange(
  top_row,
  bottom_row,
  nrow    = 2,
  heights = c(1.5, 1.8)
)

## ── 10. Save ──────────────────────────────────────────────────────────────────

dir.create("results/4_functional_change", showWarnings = FALSE, recursive = TRUE)
ggsave(
  fig4,
  filename = "results/4_functional_change/figure4.pdf",
  width    = 14,
  height   = 17,
  device   = cairo_pdf
)

## ── 11. Supplementary figure: HUMAnN violins + top CAZyme violins ────────────

suppl_panels <- c(list(pl_B, pl_C), pl_cazyme_top)
suppl_ncol   <- 3
suppl_nrow   <- ceiling(length(suppl_panels) / suppl_ncol)

suppl_fig4 <- ggarrange(
  plotlist = suppl_panels,
  nrow     = suppl_nrow,
  ncol     = suppl_ncol,
  labels   = LETTERS[seq_along(suppl_panels)]
)

ggsave(
  suppl_fig4,
  filename = "results/4_functional_change/suppl_figure_violins.pdf",
  width    = 5 * suppl_ncol,
  height   = 5 * suppl_nrow,
  device   = cairo_pdf
)
