## Figure 4 — Functional shifts: CAZymes & metabolic pathways (shotgun)
##
## Panels:
##   A   — CAZyme forest plot: families with FDR < 0.05 ethnicity × timepoint
##           interactions (adjusted LMM; coloured by FDR significance)
##   B–D — GH13, GH5, GH32 box-violin plots
##           (raw data by ethnicity × timepoint; annotated with LMM p)
##   E   — GutSMASH metabolic pathway forest (FDR < 0.05)
##   F–H — Top-3 FDR-significant GutSMASH pathway box-violin plots

library(ggpubr)
library(tidyverse)

## ── 1. Source analysis scripts ───────────────────────────────────────────────

source("scripts/4_functional_change/cayman/3_cayman_longitudinal_lmm.R")

# Stash CAZyme objects before gutsmash sourcing can overwrite shared names
cayman_statres_adj <- statres_adj      # full adjusted LMM results (has padj)
cayman_data_adj    <- dftot_adj        # individual-level data for box-violins
cayman_sig         <- statres_adj_sig  # pval < 0.05 significant families

source("scripts/4_functional_change/gutsmash/3_gutsmash_longitudinal_lmm.R")

gutsmash_statres <- statres   # full LMM results (has padj)
gutsmash_data    <- df_clin   # individual-level data for box-violins

## ── Helper: box-violin panel ─────────────────────────────────────────────────
make_boxviolin <- function(df, family_name, pval, y_label) {
  df$mb <- log10(df[[family_name]] + 1)
  pval_label <- formatC(pval, format = "e", digits = 2)

  df$eth_label <- gsub("South-Asian Surinamese", "South-Asian\nSurinamese", as.character(df$EthnicityTot))
  df$eth_label <- factor(df$eth_label, levels = gsub("South-Asian Surinamese", "South-Asian\nSurinamese",
                                                      levels(df$EthnicityTot)))

  ggplot(df, aes(x = eth_label, y = mb, fill = EthnicityTot)) +
    geom_violin(colour = NA, alpha = 0.75) +
    geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
    facet_wrap(~timepoint) +
    scale_fill_jco(guide = "none") +
    theme_Publication() +
    labs(
      x        = "",
      y        = y_label,
      subtitle = paste0("Ethnicity \u00d7 Timepoint p=", pval_label)
    )
}

## ── 2. Panel A: CAZyme forest plot (top 20 FDR < 0.05 by |estimate|) ─────────
statres_adj_q <- cayman_statres_adj %>%
  filter(padj < 0.05) %>%
  slice_max(order_by = abs(estimate), n = 20) %>%
  mutate(family    = fct_reorder(family, estimate),
         direction = ifelse(estimate > 0, "SAS", "Dutch"))

pl_4A <- ggplot(statres_adj_q, aes(x = estimate, y = family, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(xmin = conflow, xmax = confhigh), size = 0.4, linewidth = 0.5) +
  scale_colour_manual(values = c("Dutch" = "#2166AC", "SAS" = "#E6B800"),
                      name = "", labels = c("Dutch" = "More increase in Dutch", "SAS" = "More increase in SAS")) +
  theme_Publication() +
  labs(
    x        = "Interaction effect (\u00b1 95% CI)",
    y        = "",
    title    = "Differential CAZyme dynamics by ethnicity",
    subtitle = "Adjusted LMM: ethnicity \u00d7 timepoint interaction (FDR < 0.05)\nAdjusted for baseline age, sex, BMI, smoking, PPI"
  )

## ── 3. Panels B–D: top-3 FDR-significant CAZyme box-violin plots ─────────────
target_fam <- statres_adj_q %>%
  arrange(padj) %>%
  slice_head(n = 3) %>%
  pull(family) %>%
  as.character()

pl_B_panels <- lapply(target_fam, function(fam) {
  pval_row <- cayman_statres_adj %>% filter(family == fam)
  pval <- if (nrow(pval_row) > 0) pval_row$pval[1] else NA_real_

  make_boxviolin(cayman_data_adj, fam, pval, y_label = "log10(RPKM + 1)") +
    labs(title = fam) +
    theme(plot.title = element_text(face = "bold", size = rel(0.9), hjust = 0.5))
})

pl_4BCD <- ggarrange(
  plotlist = pl_B_panels,
  ncol     = 3,
  labels   = c("B", "C", "D")
)

## ── 4. Panel E: GutSMASH forest plot (FDR < 0.05) ────────────────────────────
gs_fdr <- gutsmash_statres %>%
  filter(padj < 0.05) %>%
  mutate(
    label = str_remove_all(family, "[\\[\\]']"),
    label = str_trunc(label, 55),
    label = fct_reorder(label, estimate)
  )

if (nrow(gs_fdr) >= 2) {
  gs_fdr <- gs_fdr %>% mutate(direction = ifelse(estimate > 0, "SAS", "Dutch"))

  pl_4E <- ggplot(gs_fdr, aes(x = estimate, y = label, colour = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(
      aes(xmin = conflow, xmax = confhigh),
      size = 0.45, linewidth = 0.55
    ) +
    scale_colour_manual(values = c("Dutch" = "#2166AC", "SAS" = "#E6B800"),
                        name = "", labels = c("Dutch" = "More increase in Dutch", "SAS" = "More increase in SAS")) +
    theme_Publication() +
    labs(
      x        = "Interaction effect (\u00b1 95% CI)",
      y        = "",
      title    = "GutSMASH metabolic pathway shifts",
      subtitle = "LMM ethnicity \u00d7 timepoint, FDR < 0.05"
    )
} else {
  # Fallback: single-panel placeholder if fewer than 2 FDR-sig pathways
  pl_4E <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No FDR-significant\nGutSMASH pathways",
             hjust = 0.5, vjust = 0.5, size = 4)
}

## ── 5. Panels F–H: top-3 FDR-significant GutSMASH box-violin plots ───────────
gs_sig_top <- gutsmash_statres %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = 3)

if (nrow(gs_sig_top) > 0) {
  pl_D_panels <- lapply(seq_len(nrow(gs_sig_top)), function(i) {
    nm       <- gs_sig_top$family[i]
    pval     <- gs_sig_top$pval[i]
    nm_clean <- str_trunc(str_remove_all(nm, "[\\[\\]']"), 40)

    make_boxviolin(gutsmash_data, nm, pval,
                   y_label = "log10(abundance % + 0.01)") +
      labs(title = nm_clean) +
      theme(plot.title = element_text(face = "bold", size = rel(0.75),
                                      hjust = 0.5))
  })

  pl_4FGH <- ggarrange(
    plotlist = pl_D_panels,
    ncol     = length(pl_D_panels),
    labels   = c("F", "G", "H")[seq_len(nrow(gs_sig_top))]
  )
} else {
  pl_4FGH <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No FDR-significant\nGutSMASH pathways",
             hjust = 0.5, vjust = 0.5, size = 4)
}

## ── 6. Assemble Figure 4 ─────────────────────────────────────────────────────
#
#  Layout (portrait, 14 × 18 in):
#
#   ┌──────────────────┬──────────────────┐
#   │  A: CAZyme       │  E: GutSMASH     │  <- forest plots side by side
#   │  forest (top 20) │  forest (FDR)    │
#   ├──────────────────┴──────────────────┤
#   │  B: GH13   │  C: GH5   │  D: GH32  │  <- CAZyme violins, full width
#   ├────────────────────────────────────┤
#   │  F: pw1    │  G: pw2   │  H: pw3   │  <- GutSMASH violins, full width
#   └────────────────────────────────────┘

forest_row <- ggarrange(
  pl_4A, pl_4E,
  ncol   = 2,
  labels = c("A", "E")
)

fig4 <- ggarrange(
  forest_row,
  pl_4BCD,
  pl_4FGH,
  nrow    = 3,
  heights = c(1.4, 1.0, 1.0),
  labels  = c("", "", "")
)

dir.create("results/4_functional_change", showWarnings = FALSE, recursive = TRUE)
ggsave(
  fig4,
  filename = "results/4_functional_change/figure4.pdf",
  width    = 12,
  height   = 12,
  device   = cairo_pdf
)

