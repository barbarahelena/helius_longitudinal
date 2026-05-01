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
## aplot::insert_right() aligns axes via patchwork/grid; grid.grabExpr captures
## the rendered output as a vector grob; as_ggplot wraps it for ggarrange.

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

plab_gag_bl <- formatC(wilcox_res$p_GAG_DF[wilcox_res$timepoint == "baseline"],
                        format = "e", digits = 2)
plab_gag_fu <- formatC(wilcox_res$p_GAG_DF[wilcox_res$timepoint == "follow-up"],
                        format = "e", digits = 2)

p_gag_vln <- ggplot(dftot, aes(x = timepoint, y = log10_GAG_DF, fill = EthnicityTot)) +
  geom_violin(colour = NA, aes(alpha = timepoint)) +
  geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
  facet_wrap(~EthnicityTot) +
  scale_fill_jco(guide = "none") +
  scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
  theme_Publication() +
  labs(x        = "",
       y        = "GAG / DF (log₁₀ ratio)",
       title    = "GAG-to-DF ratio",
       subtitle = paste0("Baseline p=", plab_gag_bl, "  Follow-up p=", plab_gag_fu))

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
    theme(plot.title = element_text(face = "bold", size = rel(0.75), hjust = 0.5))
}

pl_B <- make_humann_vln("Gluconeogenesis III")
pl_C <- make_humann_vln("L-lysine biosynthesis II")

## ── 7. Build named CAZyme violin panels (F, G) ────────────────────────────────

get_cazyme_pval <- function(fam) {
  row <- statres_adj[statres_adj$family == fam, ]
  if (nrow(row) == 0) NA_real_ else row$pval[1]
}

make_cazyme_vln <- function(fam_name) {
  make_boxviolin(dftot_adj, fam_name, get_cazyme_pval(fam_name), "log10(CPM + 1)") +
    labs(title = fam_name) +
    theme(plot.title = element_text(face = "bold", size = rel(0.9), hjust = 0.5))
}

pl_F <- make_cazyme_vln("GH13_16")
pl_G <- make_cazyme_vln("GH5_26")

## ── 8. Assemble rows ──────────────────────────────────────────────────────────

# Row 1 — HUMAnN: forest+heatmap (wide) + two pathway violins
humann_row <- ggarrange(
  pl_A, pl_B, pl_C,
  ncol   = 3,
  labels = c("A", "B", "C"),
  widths = c(2, 1, 1)
)

# Row 2 — CAZyme ordination: two Canberra PCoAs (D) + forest+heatmap (E)
cazyme_ordin_row <- ggarrange(
  pl_can_bl, pl_can_fu, pl_E,
  ncol   = 3,
  labels = c("D", "", "E"),
  widths = c(1, 1, 2)
)

# Row 3 — CAZyme details: two family violins (F, G) + GAG ratio (H)
cazyme_vln_row <- ggarrange(
  pl_F, pl_G, p_gag_vln,
  ncol   = 3,
  labels = c("F", "G", "H"),
  widths = c(1, 1, 1)
)

## ── 9. Final assembly ─────────────────────────────────────────────────────────

fig4 <- ggarrange(
  humann_row,
  cazyme_ordin_row,
  cazyme_vln_row,
  nrow    = 3,
  heights = c(1.2, 1.2, 1.0)
)

## ── 10. Save ──────────────────────────────────────────────────────────────────

dir.create("results/4_functional_change", showWarnings = FALSE, recursive = TRUE)
ggsave(
  fig4,
  filename = "results/4_functional_change/figure4.pdf",
  width    = 18,
  height   = 15,
  device   = cairo_pdf
)
