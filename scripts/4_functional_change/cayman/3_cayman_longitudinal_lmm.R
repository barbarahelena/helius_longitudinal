## Cayman Longitudinal Analysis with Linear Mixed Models
library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(lmerTest)
library(patchwork)

# Theme ------------------------------------------------------------------------
theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family) +
    theme(plot.title        = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
          text              = element_text(),
          panel.background  = element_rect(colour = NA, fill = NA),
          plot.background   = element_rect(colour = NA, fill = NA),
          panel.border      = element_rect(colour = NA),
          axis.title        = element_text(face = "bold", size = rel(0.8)),
          axis.title.y      = element_text(angle = 90, vjust = 2),
          axis.title.x      = element_text(vjust = -0.2),
          axis.text         = element_text(size = rel(0.7)),
          axis.line         = element_line(colour = "black"),
          axis.ticks        = element_line(),
          panel.grid.major  = element_line(colour = "#f0f0f0"),
          panel.grid.minor  = element_blank(),
          legend.key        = element_rect(colour = NA),
          legend.position   = "bottom",
          legend.key.size   = unit(0.2, "cm"),
          legend.spacing    = unit(0, "cm"),
          strip.background  = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
          strip.text        = element_text(face = "bold"),
          plot.caption      = element_text(size = rel(0.5), face = "italic"),
          plot.subtitle     = element_text(size = 8, hjust = 0.5, face = "italic")))
}

# Data import ------------------------------------------------------------------
df_raw   <- rio::import("data/shotgun/cayman_results/families_cpm_table.tsv") |>
  dplyr::select(-HELIBA_103370, -HELIFU_103370)
clinical <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
stats    <- rio::import("data/shotgun/cayman_results/sample_statistics.tsv") |>
  dplyr::rename(sampleID = sample)
dir.create("results/4_functional_change/cayman/longitudinal", showWarnings = FALSE, recursive = TRUE)

# Prevalence filtering (>5%) and data preparation ------------------------------
rownames(df_raw)    <- df_raw$family
family_prevalence   <- rowSums(df_raw > 0) / ncol(df_raw) * 100
prevalent_families  <- names(family_prevalence)[family_prevalence > 5]

df_raw        <- df_raw %>% filter(family %in% prevalent_families)
df_raw$family <- NULL
df            <- as.data.frame(t(as.matrix(df_raw)))
df$sampleID   <- rownames(df)
gene_families <- setdiff(colnames(df), "sampleID")

dftot <- left_join(df, clinical) |> droplevels()
write.csv2(data.frame(family = prevalent_families),
           "results/4_functional_change/cayman/longitudinal/prevalent_families_list.csv",
           row.names = FALSE)

# Gene family-level LMMs ------------------------------------------------------
# Model: CAZy ~ EthnicityTot * timepoint + FUtime + (1|ID)
# Covariates (Age, BMI, Sex, Smoking, PPI) used only to restrict to complete-case
# sample; not included in model because participants were matched on these variables.
dftot_adj <- dftot

statres_adj <- data.frame()
for (gf in gene_families) {
  print(gf)
  dftot_adj$mb <- log10(dftot_adj[[gf]] + 1)
  tryCatch({
    model_adj <- lmer(mb ~ EthnicityTot * timepoint + FUtime + (1|ID), data = dftot_adj)
    res <- summary(model_adj)
    ci  <- confint(model_adj, method = "Wald")
    interaction_row <- grep("EthnicityTotSouth-Asian Surinamese:timepointfollow-up",
                            rownames(res$coefficients))
    ci_row <- grep("EthnicityTotSouth-Asian Surinamese:timepointfollow-up", rownames(ci))
    statres_adj <- rbind(statres_adj, data.frame(
      family   = gf,
      estimate = res$coefficients[interaction_row, 1],
      conflow  = ifelse(length(ci_row) > 0, ci[ci_row, 1], NA),
      confhigh = ifelse(length(ci_row) > 0, ci[ci_row, 2], NA),
      pval     = res$coefficients[interaction_row, 5]
    ))
  }, error = function(e) NULL)
}
statres_adj <- as.data.frame(statres_adj) %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv2(statres_adj,
           "results/4_functional_change/cayman/longitudinal/lmm_ethnicity_timepoint_adjusted.csv",
           row.names = FALSE)

# Figure 4 panels (A, B–D) -----------------------------------------------------
make_boxviolin <- function(df, family_name, pval, y_label, show_pval = TRUE) {
  df$mb        <- log10(df[[family_name]] + 1)
  df$timepoint <- factor(df$timepoint, levels = c("baseline", "follow-up"))
  df           <- df %>%
    group_by(EthnicityTot, ID) %>%
    filter(n() == 2) %>%
    ungroup() %>%
    arrange(EthnicityTot, ID, timepoint)

  # Paired Wilcoxon computed manually per facet to avoid ggpubr facet+paired bug
  pval_annot <- df %>%
    group_by(EthnicityTot) %>%
    group_modify(~ {
      bl <- .x$mb[.x$timepoint == "baseline"]
      fu <- .x$mb[.x$timepoint == "follow-up"]
      p  <- tryCatch(wilcox.test(bl, fu, paired = TRUE)$p.value, error = function(e) NA_real_)
      data.frame(p = p, y_pos = max(.x$mb, na.rm = TRUE) + diff(range(.x$mb, na.rm = TRUE)) * 0.08)
    }) %>%
    ungroup() %>%
    mutate(label = ifelse(p < 0.05,
                          paste0("p=", ifelse(p < 0.001,
                                              formatC(p, format = "e", digits = 2),
                                              formatC(p, format = "f", digits = 3))),
                          ""))

  pval_label <- formatC(pval, format = "e", digits = 2)

  ggplot(df, aes(x = timepoint, y = mb, fill = EthnicityTot)) +
    geom_violin(colour = NA, aes(alpha = timepoint)) +
    geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
    { if (show_pval) geom_text(data = pval_annot,
                               aes(x = 1.5, y = y_pos, label = label),
                               inherit.aes = FALSE, size = 3) } +
    facet_wrap(~EthnicityTot) +
    scale_fill_jco(guide = "none") +
    scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
    theme_Publication() +
    labs(x = "", y = y_label,
         subtitle = paste0("Ethnicity \u00d7 Timepoint p=", pval_label))
}

# Panel A: forest plot (FDR < 0.05, top 20 by |estimate|)
statres_adj_q <- statres_adj %>%
  filter(padj < 0.05) %>%
  slice_max(order_by = abs(estimate), n = 20) %>%
  mutate(family    = fct_reorder(family, estimate),
         direction = ifelse(estimate > 0, "SAS", "Dutch"))

pl_4A <- ggplot(statres_adj_q, aes(x = estimate, y = family, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(xmin = conflow, xmax = confhigh), size = 0.4, linewidth = 0.5) +
  scale_colour_manual(values = c("Dutch" = "#2166AC", "SAS" = "#E6B800"),
                      name = "", labels = c("Dutch" = "Greater positive change in Dutch",
                                            "SAS"   = "Greater positive change in SAS")) +
  theme_Publication() +
  labs(x        = "Interaction effect (\u00b1 95% CI)",
       y        = "",
       title    = "Differential CAZyme changes by ethnicity",
       subtitle = "Adjusted LMM: ethnicity \u00d7 timepoint interaction (FDR < 0.05)")

# Panels B–D: top-3 FDR-significant CAZyme box-violin plots
target_fam <- statres_adj_q %>%
  arrange(padj) %>%
  slice_head(n = 3) %>%
  pull(family) %>%
  as.character()

pl_B_panels <- lapply(target_fam, function(fam) {
  pval_row <- statres_adj %>% filter(family == fam)
  pval     <- if (nrow(pval_row) > 0) pval_row$pval[1] else NA_real_
  make_boxviolin(dftot_adj, fam, pval, y_label = "log10(CPM + 1)") +
    labs(title = fam) +
    theme(plot.title = element_text(face = "bold", size = rel(0.9), hjust = 0.5))
})

pl_4BCD <- ggarrange(plotlist = pl_B_panels, ncol = 3, labels = c("B", "C", "D"))

# ---------------------------------------------------------------------------
# Overlap: cross-sectional significant vs LMM interaction significant
# ---------------------------------------------------------------------------
cs_results <- read.csv2("results/4_functional_change/cayman/crossectional/crosssectional_results_combined.csv",
                        stringsAsFactors = FALSE) |>
    mutate(across(starts_with("padj_"), as.numeric))

lmm_annotated <- statres_adj |>
    mutate(
        sig_lmm_interaction = padj < 0.05,
        lmm_direction = ifelse(estimate > 0, "SAS increases more", "Dutch increases more")
    ) |>
    dplyr::rename(gene_family = family)

overlap <- cs_results |>
    left_join(lmm_annotated |> dplyr::select(gene_family, sig_lmm_interaction, lmm_direction, estimate, padj),
              by = "gene_family") |>
    mutate(
        pattern = case_when(
            sig_lmm_interaction & sig_baseline & sig_followup &
                cs_baseline_direction == "SAS higher"  & lmm_direction == "SAS increases more" ~ "divergence (SAS)",
            sig_lmm_interaction & sig_baseline & sig_followup &
                cs_baseline_direction == "Dutch higher" & lmm_direction == "Dutch increases more" ~ "divergence (Dutch)",
            sig_lmm_interaction & sig_baseline & sig_followup &
                cs_baseline_direction == "SAS higher"  & lmm_direction == "Dutch increases more" ~ "convergence (SAS shrinks)",
            sig_lmm_interaction & sig_baseline & sig_followup &
                cs_baseline_direction == "Dutch higher" & lmm_direction == "SAS increases more" ~ "convergence (Dutch shrinks)",
            sig_lmm_interaction & sig_baseline & !sig_followup &
                cs_baseline_direction == "SAS higher"  & lmm_direction == "Dutch increases more" ~ "convergence (SAS shrinks)",
            sig_lmm_interaction & sig_baseline & !sig_followup &
                cs_baseline_direction == "Dutch higher" & lmm_direction == "SAS increases more" ~ "convergence (Dutch shrinks)",
            sig_lmm_interaction & sig_baseline & !sig_followup ~ "divergence (lost significance)",
            sig_lmm_interaction & !sig_baseline & sig_followup &
                cs_followup_direction == "SAS higher"  & lmm_direction == "SAS increases more" ~ "emergence (SAS)",
            sig_lmm_interaction & !sig_baseline & sig_followup &
                cs_followup_direction == "Dutch higher" & lmm_direction == "Dutch increases more" ~ "emergence (Dutch)",
            sig_lmm_interaction & !sig_baseline & sig_followup ~ "emergence (other)",
            sig_lmm_interaction & !sig_baseline & !sig_followup ~ "LMM only",
            TRUE ~ NA_character_
        )
    )

cat("\n--- Cross-sectional vs LMM interaction overlap (Cayman) ---\n")
cat("LMM interaction-significant gene families (FDR < 0.05):", sum(overlap$sig_lmm_interaction, na.rm = TRUE), "\n")
cat("Of those, also sig at baseline in cross-sectional: ",
    sum(overlap$sig_lmm_interaction & overlap$sig_baseline,  na.rm = TRUE), "\n")
cat("Of those, also sig at follow-up in cross-sectional: ",
    sum(overlap$sig_lmm_interaction & overlap$sig_followup,  na.rm = TRUE), "\n")
cat("Of those, sig at both timepoints in cross-sectional:",
    sum(overlap$sig_lmm_interaction & overlap$sig_baseline & overlap$sig_followup, na.rm = TRUE), "\n")
cat("\nPattern classification (LMM-significant gene families):\n")
print(table(overlap$pattern[overlap$sig_lmm_interaction], useNA = "ifany"))

write.csv2(overlap,
           "results/4_functional_change/cayman/longitudinal/crosssectional_vs_lmm_overlap.csv",
           row.names = FALSE)

# ---------------------------------------------------------------------------
# Forest plot + cross-sectional heatmap (aplot)
# ---------------------------------------------------------------------------
if (nrow(statres_adj_q) >= 2) {
    fam_vec <- as.character(statres_adj_q$family)

    # log10 fold change: log10(SAS median) − log10(Dutch median) = log10(SAS/Dutch)
    heatmap_diff <- dftot_adj |>
        dplyr::select(sampleID, EthnicityTot, timepoint, all_of(fam_vec)) |>
        pivot_longer(cols = all_of(fam_vec), names_to = "gene_family", values_to = "abund") |>
        mutate(log_abund = log10(abund + 1)) |>
        group_by(gene_family, timepoint, EthnicityTot) |>
        summarise(median_log = median(log_abund, na.rm = TRUE), .groups = "drop") |>
        pivot_wider(names_from = EthnicityTot, values_from = median_log) |>
        mutate(log2fc = `South-Asian Surinamese` - Dutch)

    cs_sig_mask <- overlap |>
        filter(gene_family %in% fam_vec) |>
        dplyr::select(gene_family, sig_baseline, sig_followup) |>
        mutate(across(c(sig_baseline, sig_followup), as.logical)) |>
        left_join(
            read.csv2("results/4_functional_change/cayman/crossectional/crosssectional_results_combined.csv",
                      stringsAsFactors = FALSE) |>
                dplyr::select(gene_family, padj_baseline, padj_followup),
            by = "gene_family"
        )

    heatmap_data <- heatmap_diff |>
        left_join(cs_sig_mask, by = "gene_family") |>
        mutate(
            diff_display = case_when(
                timepoint == "baseline"  & sig_baseline  ~ log2fc,
                timepoint == "follow-up" & sig_followup  ~ log2fc,
                TRUE ~ NA_real_
            ),
            padj_cs = case_when(
                timepoint == "baseline"  ~ padj_baseline,
                timepoint == "follow-up" ~ padj_followup
            ),
            star = case_when(
                !is.na(diff_display) & padj_cs < 0.001 ~ "***",
                !is.na(diff_display) & padj_cs < 0.01  ~ "**",
                !is.na(diff_display) & padj_cs < 0.05  ~ "*",
                TRUE ~ ""
            )
        ) |>
        left_join(
            statres_adj_q |>
                mutate(gene_family = as.character(family)) |>
                dplyr::select(gene_family, family),
            by = "gene_family"
        ) |>
        mutate(timepoint = factor(timepoint, levels = c("baseline", "follow-up")))

    abs_lim <- max(abs(heatmap_data$diff_display), na.rm = TRUE)

    pl_heatmap <- ggplot(heatmap_data,
                         aes(x = timepoint, y = family, fill = diff_display)) +
        geom_tile(color = "white", linewidth = 0.4) +
        geom_text(aes(label = star), color = "black", size = 2.5, vjust = 0.75) +
        scale_fill_gradient2(
            low      = "#2166AC",
            mid      = "white",
            high     = "#E6B800",
            na.value = "grey93",
            limits   = c(-abs_lim, abs_lim),
            name     = "log₁₀ FC\n(SAS / Dutch)",
            guide    = guide_colorbar(barheight = unit(6, "cm"),
                                      barwidth  = unit(0.5, "cm"))
        ) +
        theme_Publication() +
        theme(
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y  = element_blank(),
            axis.text.x  = element_text(angle = 45, hjust = 1, size = rel(1.0)),
            legend.position = "right"
        ) +
        labs(x = "", y = "")

    pl_combined <- (pl_4A | pl_heatmap) +
        plot_layout(widths = c(4, 1), guides = "keep")

    combined_height <- max(4, nrow(statres_adj_q) * 0.3 + 2)
    cairo_pdf("results/4_functional_change/cayman/longitudinal/forest_heatmap_cayman.pdf",
              width = 9, height = combined_height)
    print(pl_combined)
    dev.off()
}

# ===========================================================================
# SUPPLEMENT: all FDR-significant families — forest plot + violin plots
# ===========================================================================
dir.create("results/4_functional_change/cayman/longitudinal/supplement",
           showWarnings = FALSE, recursive = TRUE)

statres_supp <- statres_adj %>%
  filter(padj < 0.05) %>%
  mutate(family    = fct_reorder(family, estimate),
         direction = ifelse(estimate > 0, "SAS", "Dutch"))

# Supplement forest plot: all FDR-significant families (no n=20 cap)
pl_supp_forest <- ggplot(statres_supp,
                         aes(x = estimate, y = family, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(xmin = conflow, xmax = confhigh), size = 0.4, linewidth = 0.5) +
  scale_colour_manual(values = c("Dutch" = "#2166AC", "SAS" = "#E6B800"),
                      name   = "",
                      labels = c("Dutch" = "Greater positive change in Dutch",
                                 "SAS"   = "Greater positive change in SAS")) +
  theme_Publication() +
  labs(x        = "Interaction effect (± 95% CI)",
       y        = "",
       title    = "Differential CAZyme dynamics by ethnicity",
       subtitle = "All FDR-significant gene families (padj < 0.05)")

supp_forest_height <- max(5, nrow(statres_supp) * 0.28 + 2)
ggsave(
  "results/4_functional_change/cayman/longitudinal/supplement/supplement_forest_all.pdf",
  pl_supp_forest, width = 8, height = supp_forest_height, device = cairo_pdf
)

# Supplement violin plots: all FDR-significant families, 9 per page
supp_families <- as.character(statres_supp$family)

pdf("results/4_functional_change/cayman/longitudinal/supplement/supplement_violins_all.pdf",
    width = 18, height = 15)
letter_idx <- 1
for (i in seq(1, length(supp_families), by = 9)) {
  batch <- supp_families[i:min(i + 8, length(supp_families))]
  plots <- lapply(batch, function(fam) {
    pval_row <- statres_adj %>% filter(family == fam)
    pval     <- if (nrow(pval_row) > 0) pval_row$pval[1] else NA_real_
    make_boxviolin(dftot_adj, fam, pval, "log10(CPM + 1)") +
      labs(title = fam) +
      theme(plot.title = element_text(face = "bold", size = rel(0.9), hjust = 0.5))
  })
  batch_labels <- LETTERS[letter_idx:(letter_idx + length(batch) - 1)]
  letter_idx   <- letter_idx + length(batch)
  print(ggarrange(plotlist = plots, ncol = 3, nrow = 3, labels = batch_labels))
}
dev.off()
