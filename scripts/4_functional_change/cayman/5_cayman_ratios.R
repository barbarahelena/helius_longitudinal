## Cayman Mucin/DF and GAG/DF ratios by ethnicity
library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ComplexHeatmap)
library(circlize)
library(Cairo)
library(ggsci)

theme_Publication <- function(base_size = 14, base_family = "sans") {
  library(grid)
  library(ggthemes)
  suppressWarnings(theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title       = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
      text             = element_text(),
      panel.background = element_rect(colour = NA, fill = NA),
      plot.background  = element_rect(colour = NA, fill = NA),
      panel.border     = element_rect(colour = NA),
      axis.title       = element_text(face = "bold", size = rel(0.8)),
      axis.title.y     = element_text(angle = 90, vjust = 2),
      axis.title.x     = element_text(vjust = -0.2),
      axis.text        = element_text(size = rel(0.7)),
      axis.line        = element_line(colour = "black"),
      axis.ticks       = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key       = element_rect(colour = NA),
      legend.position  = "bottom",
      legend.key.size  = unit(0.2, "cm"),
      legend.spacing   = unit(0, "cm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 8, hjust = 0.5, face = "italic")
    ))
}

ETH_DUTCH <- "Dutch"
ETH_SAS   <- "South-Asian Surinamese"
eth_colors <- c("#2166AC", "#E6B800")
names(eth_colors) <- c(ETH_DUTCH, ETH_SAS)

dir.create("results/4_functional_change/cayman/ratios", showWarnings = FALSE, recursive = TRUE)

# --- Load data ------------------------------------------------------------
# Substrate group sums and ratios (Mucin/DF, GAG/DF) are now precomputed
# upstream in the pipeline, not calculated here.
df_sums <- rio::import("data/shotgun/cayman_results/oct2026_results/substrates_rpkm_table.tsv") |>
  # S103370's baseline CAZy profile essentially failed (43 families detected,
  # 0.13% CAZy reads), which distorts its substrate ratios (e.g. log10_GAG_DF
  # swings from -2.77 at baseline to -0.99 at follow-up); excluded here to
  # match the family-table exclusion used in scripts 1-4.
  dplyr::filter(!sample %in% c("HELIBA_103370", "HELIFU_103370")) |>
  dplyr::rename(
    sampleID       = sample,
    ratio_Mucin_DF = mucin_df_ratio,
    ratio_GAG_DF   = gag_df_ratio,
    log10_Mucin_DF = log10_mucin_df_ratio,
    log10_GAG_DF   = log10_gag_df_ratio
  )

clinical <- readRDS("data/clinicaldata/clinicaldata_long.RDS")

# --- Merge with clinical data -------------------------------------------------
dftot <- df_sums |>
  left_join(clinical, by = "sampleID") |>
  filter(EthnicityTot %in% c(ETH_DUTCH, ETH_SAS)) |>
  mutate(
    EthnicityTot = factor(EthnicityTot, levels = c(ETH_DUTCH, ETH_SAS)),
    timepoint    = factor(timepoint, levels = c("baseline", "follow-up"))
  ) |>
  droplevels()

# --- Wilcoxon tests per timepoint ---------------------------------------------
wilcox_res <- dftot |>
  group_by(timepoint) |>
  summarise(
    p_Mucin_DF = wilcox.test(log10_Mucin_DF ~ EthnicityTot)$p.value,
    p_GAG_DF   = wilcox.test(log10_GAG_DF   ~ EthnicityTot)$p.value,
    .groups = "drop"
  ) |>
  mutate(
    padj_Mucin_DF = p.adjust(p_Mucin_DF, method = "fdr"),
    padj_GAG_DF   = p.adjust(p_GAG_DF,   method = "fdr")
  )
print(wilcox_res)
write.csv2(wilcox_res,
           "results/4_functional_change/cayman/ratios/wilcoxon_ratios_by_ethnicity.csv",
           row.names = FALSE)

# --- LMMs: main effect of ethnicity + interaction with timepoint --------------
lmm_results <- list()
for (ratio_var in c("log10_Mucin_DF", "log10_GAG_DF")) {
  dftot$y <- dftot[[ratio_var]]
  tryCatch({
    mod <- lmer(y ~ EthnicityTot * timepoint + (1 | ID), data = dftot)
    res <- summary(mod)$coefficients
    main_row <- grep(paste0("^EthnicityTot", ETH_SAS, "$"), rownames(res))
    int_row  <- grep(paste0(ETH_SAS, ":timepointfollow-up"), rownames(res))
    lmm_results[[ratio_var]] <- tibble(
      ratio             = ratio_var,
      estimate_main     = res[main_row, 1],
      se_main           = res[main_row, 2],
      pval_main         = res[main_row, 5],
      estimate_interact = res[int_row,  1],
      se_interact       = res[int_row,  2],
      pval_interact     = res[int_row,  5]
    )
  }, error = function(e) message("LMM failed for ", ratio_var, ": ", e$message))
}
lmm_df <- bind_rows(lmm_results) |>
  mutate(
    padj_main     = p.adjust(pval_main,     method = "fdr"),
    padj_interact = p.adjust(pval_interact, method = "fdr")
  )
print(lmm_df)
write.csv2(lmm_df,
           "results/4_functional_change/cayman/ratios/lmm_ratio_ethnicity_timepoint.csv",
           row.names = FALSE)

# Unadjusted LMM in dietary subset (n=232): GAG/DF p=0.276, Mucin/DF p=0.579 — already non-significant without any dietary adjustment.

# --- Spearman correlations: baseline ratios vs baseline and follow-up outcomes -
cont_outcomes <- c("BMI", "WHR", "SBP", "DBP", "HbA1c", "Trig", "TC", "HDL", "LDL", "Fatperc")

df_baseline  <- dftot |> filter(timepoint == "baseline")
df_followup  <- dftot |> filter(timepoint == "follow-up")

# Join follow-up outcomes to baseline ratio data by ID
df_bl_fu <- df_baseline |>
  dplyr::select(ID, log10_Mucin_DF, log10_GAG_DF) |>
  left_join(
    df_followup |> dplyr::select(ID, all_of(cont_outcomes)),
    by = "ID"
  )

cor_pval <- function(x, y) {
  idx <- complete.cases(x, y)
  if (sum(idx) < 3) return(c(rho = NA_real_, pval = NA_real_))
  res <- suppressWarnings(cor.test(x[idx], y[idx], method = "spearman"))
  c(rho = unname(res$estimate), pval = res$p.value)
}

run_spearman <- function(data, tp_label) {
  expand.grid(
    ratio   = c("log10_Mucin_DF", "log10_GAG_DF"),
    outcome = cont_outcomes,
    stringsAsFactors = FALSE
  ) |>
    as_tibble() |>
    rowwise() |>
    mutate(
      timepoint = tp_label,
      n    = sum(complete.cases(data[[ratio]], data[[outcome]])),
      rho  = cor_pval(data[[ratio]], data[[outcome]])[["rho"]],
      pval = cor_pval(data[[ratio]], data[[outcome]])[["pval"]]
    ) |>
    ungroup()
}

spearman_bl <- run_spearman(df_baseline, "baseline")
spearman_fu <- run_spearman(df_bl_fu,   "follow-up")

spearman_res <- bind_rows(spearman_bl, spearman_fu) |>
  mutate(padj = p.adjust(pval, method = "fdr"))

print(spearman_res)
write.csv2(spearman_res,
           "results/4_functional_change/cayman/ratios/spearman_baseline_outcomes.csv",
           row.names = FALSE)

# --- ComplexHeatmap of Spearman correlations ----------------------------------
outcome_labels <- c(
  BMI = "BMI", WHR = "WHR", SBP = "SBP", DBP = "DBP",
  HbA1c = "HbA1c", Trig = "Triglycerides", TC = "Total Cholesterol",
  HDL = "HDL", LDL = "LDL", Fatperc = "Body fat %"
)
ratio_labels <- c(log10_Mucin_DF = "Mucin/DF", log10_GAG_DF = "GAG/DF")

make_matrices <- function(res, tp) {
  r <- res |> filter(timepoint == tp)
  cor_m <- r |>
    dplyr::select(ratio, outcome, rho) |>
    pivot_wider(names_from = ratio, values_from = rho) |>
    column_to_rownames("outcome") |>
    as.matrix()
  cor_m <- cor_m[cont_outcomes, ]
  rownames(cor_m) <- outcome_labels[rownames(cor_m)]
  colnames(cor_m) <- ratio_labels[colnames(cor_m)]

  # Heatmap stars show nominal (unadjusted) significance; nothing survives FDR here.
  pval_m <- r |>
    dplyr::select(ratio, outcome, pval) |>
    pivot_wider(names_from = ratio, values_from = pval) |>
    column_to_rownames("outcome") |>
    as.matrix()
  pval_m <- pval_m[cont_outcomes, ]
  rownames(pval_m) <- outcome_labels[rownames(pval_m)]
  colnames(pval_m) <- ratio_labels[colnames(pval_m)]

  list(cor = cor_m, pval = pval_m)
}

mat_bl <- make_matrices(spearman_res, "baseline")
mat_fu <- make_matrices(spearman_res, "follow-up")

col_fun <- colorRamp2(
  c(-0.3, 0, 0.3),
  c(pal_nejm()(6)[6], "white", pal_nejm()(3)[3])
)

make_heatmap <- function(cor_m, pval_m, title, show_row_names = TRUE) {
  Heatmap(
    cor_m,
    name              = "Spearman\nCorrelation",
    col               = col_fun,
    rect_gp           = gpar(col = "white", lwd = 2),
    na_col            = "grey95",
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_row_names    = show_row_names,
    show_column_names = TRUE,
    row_names_side    = "left",
    row_names_gp      = gpar(fontsize = 10),
    column_names_gp   = gpar(fontsize = 12, fontface = "bold"),
    column_names_rot  = 45,
    column_title      = title,
    column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
    show_heatmap_legend = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      pval <- pval_m[i, j]
      if (!is.na(pval)) {
        sig <- if (pval < 0.001) "***" else if (pval < 0.01) "**" else if (pval < 0.05) "*" else ""
        if (sig != "") grid.text(sig, x, y, gp = gpar(fontsize = 14), vjust = 0.75)
      }
    }
  )
}

ht_bl <- make_heatmap(mat_bl$cor, mat_bl$pval, "Baseline outcomes",  show_row_names = TRUE)
ht_fu <- make_heatmap(mat_fu$cor, mat_fu$pval, "Follow-up outcomes", show_row_names = FALSE)

# Give unique internal names to suppress ComplexHeatmap duplicate-name warning
ht_bl@name <- "bl"
ht_fu@name <- "fu"

lgd_cor <- Legend(
  col_fun = col_fun,
  title   = "Spearman\nCorrelation",
  at      = c(-0.3, -0.15, 0, 0.15, 0.3),
  labels  = c("-0.30", "-0.15", "0", "0.15", "0.30")
)
lgd_sig <- Legend(
  pch = c("*", "**", "***"), type = "points",
  labels = c("p < 0.05", "p < 0.01", "p < 0.001"),
  legend_gp = gpar(fontsize = 10)
)
lgd_packed <- packLegend(lgd_cor, lgd_sig, direction = "vertical", gap = unit(4, "mm"))

CairoPDF("results/4_functional_change/cayman/ratios/heatmap_spearman_outcomes.pdf",
         width = 7, height = 6)
draw(ht_bl + ht_fu, annotation_legend_list = list(lgd_packed),
     padding = unit(c(5, 30, 5, 5), "mm"))
dev.off()

# --- Spearman correlations: baseline ratios vs dietary intake -----------------
# Willett residual-adjusted macronutrients are pre-computed in dietarydata.R
# and stored in clinicaldata_wide.RDS as *_baseline_adj (scaled residuals from
# regressing each macronutrient on total energy intake across the full cohort).
helius_wide_diet <- readRDS("data/clinicaldata/clinicaldata_wide.RDS") |>
  dplyr::select(ID, TotalCalories_baseline,
                ends_with("_baseline_adj"))

df_baseline_diet <- df_baseline |>
  left_join(helius_wide_diet, by = "ID")

diet_vars <- c("TotalCalories_baseline",
               "Protein_baseline_adj", "Protein_animal_baseline_adj",
               "FattyAcids_baseline_adj", "SatFat_baseline_adj",
               "Carbohydrates_baseline_adj", "Fiber_baseline_adj",
               "Sodium_g_baseline_adj")
diet_vars <- diet_vars[diet_vars %in% names(df_baseline_diet)]

diet_labels <- c(
  TotalCalories_baseline        = "Total calories",
  Protein_baseline_adj          = "Protein",
  Protein_animal_baseline_adj   = "Animal protein",
  FattyAcids_baseline_adj       = "Fatty acids",
  SatFat_baseline_adj           = "Saturated fat",
  Carbohydrates_baseline_adj    = "Carbohydrates",
  Fiber_baseline_adj            = "Fiber",
  Sodium_g_baseline_adj         = "Sodium"
)

spearman_diet_all <- expand.grid(
  ratio   = c("log10_Mucin_DF", "log10_GAG_DF"),
  outcome = diet_vars,
  stringsAsFactors = FALSE
) |>
  as_tibble() |>
  rowwise() |>
  mutate(
    EthnicityTot = "All",
    n    = sum(complete.cases(df_baseline_diet[[ratio]], df_baseline_diet[[outcome]])),
    rho  = cor_pval(df_baseline_diet[[ratio]], df_baseline_diet[[outcome]])[["rho"]],
    pval = cor_pval(df_baseline_diet[[ratio]], df_baseline_diet[[outcome]])[["pval"]]
  ) |>
  ungroup()

spearman_diet_eth <- bind_rows(lapply(c(ETH_DUTCH, ETH_SAS), function(eth) {
  df_eth <- df_baseline_diet |> filter(EthnicityTot == eth)
  expand.grid(
    ratio   = c("log10_Mucin_DF", "log10_GAG_DF"),
    outcome = diet_vars,
    stringsAsFactors = FALSE
  ) |>
    as_tibble() |>
    rowwise() |>
    mutate(
      EthnicityTot = eth,
      n    = sum(complete.cases(df_eth[[ratio]], df_eth[[outcome]])),
      rho  = cor_pval(df_eth[[ratio]], df_eth[[outcome]])[["rho"]],
      pval = cor_pval(df_eth[[ratio]], df_eth[[outcome]])[["pval"]]
    ) |>
    ungroup()
}))

spearman_diet <- bind_rows(spearman_diet_all, spearman_diet_eth) |>
  group_by(EthnicityTot) |>
  mutate(padj = p.adjust(pval, method = "fdr")) |>
  ungroup()

print(spearman_diet)
write.csv2(spearman_diet,
           "results/4_functional_change/cayman/ratios/spearman_dietary_correlations.csv",
           row.names = FALSE)

make_diet_matrix <- function(res, eth_label) {
  r     <- res |> filter(EthnicityTot == eth_label)
  avail <- diet_vars[diet_vars %in% r$outcome]
  cor_m <- r |>
    dplyr::select(ratio, outcome, rho) |>
    pivot_wider(names_from = ratio, values_from = rho) |>
    column_to_rownames("outcome") |>
    as.matrix()
  cor_m <- cor_m[avail, , drop = FALSE]
  rownames(cor_m) <- diet_labels[rownames(cor_m)]
  colnames(cor_m) <- ratio_labels[colnames(cor_m)]

  # Heatmap stars show nominal (unadjusted) significance; nothing survives FDR here.
  pval_m <- r |>
    dplyr::select(ratio, outcome, pval) |>
    pivot_wider(names_from = ratio, values_from = pval) |>
    column_to_rownames("outcome") |>
    as.matrix()
  pval_m <- pval_m[avail, , drop = FALSE]
  rownames(pval_m) <- diet_labels[rownames(pval_m)]
  colnames(pval_m) <- ratio_labels[colnames(pval_m)]

  list(cor = cor_m, pval = pval_m)
}

mat_diet_all   <- make_diet_matrix(spearman_diet, "All")
mat_diet_dutch <- make_diet_matrix(spearman_diet, ETH_DUTCH)
mat_diet_sas   <- make_diet_matrix(spearman_diet, ETH_SAS)

ht_diet_all   <- make_heatmap(mat_diet_all$cor,   mat_diet_all$pval,   "All",     show_row_names = TRUE)
ht_diet_dutch <- make_heatmap(mat_diet_dutch$cor, mat_diet_dutch$pval, ETH_DUTCH, show_row_names = FALSE)
ht_diet_sas   <- make_heatmap(mat_diet_sas$cor,   mat_diet_sas$pval,   ETH_SAS,   show_row_names = FALSE)
ht_diet_all@name   <- "diet_all"
ht_diet_dutch@name <- "diet_dutch"
ht_diet_sas@name   <- "diet_sas"

CairoPDF("results/4_functional_change/cayman/ratios/heatmap_spearman_diet.pdf",
         width = 9, height = 7)
draw(ht_diet_all + ht_diet_dutch + ht_diet_sas,
     annotation_legend_list = list(lgd_packed),
     padding = unit(c(5, 30, 5, 5), "mm"))
dev.off()

# --- Stability: correlation between baseline and follow-up ratios -------------
df_wide_ratios <- dftot |>
  dplyr::select(ID, timepoint, log10_Mucin_DF, log10_GAG_DF) |>
  pivot_wider(names_from = timepoint, values_from = c(log10_Mucin_DF, log10_GAG_DF))

cat("\nBaseline vs follow-up Spearman correlations:\n")
for (ratio_var in c("log10_Mucin_DF", "log10_GAG_DF")) {
  bl_col <- paste0(ratio_var, "_baseline")
  fu_col <- paste0(ratio_var, "_follow-up")
  r <- cor(df_wide_ratios[[bl_col]], df_wide_ratios[[fu_col]],
           use = "pairwise.complete.obs", method = "spearman")
  p <- cor.test(df_wide_ratios[[bl_col]], df_wide_ratios[[fu_col]],
                method = "spearman")$p.value
  cat(ratio_var, ": rho =", round(r, 3), ", p =", formatC(p, format = "e", digits = 2), "\n")
}

# --- Prospective LMs: follow-up outcome ~ baseline ratio + baseline outcome ---
# Only for ratio-outcome pairs with q<0.05 in either heatmap panel.
# Tests whether baseline microbiome predicts future metabolic state beyond
# where the outcome already was at baseline.
sig_pairs <- spearman_res |>
  filter(padj < 0.05) |>
  distinct(ratio, outcome)

# Build a wide dataset: baseline ratios + baseline outcomes + follow-up outcomes
df_prosp <- df_baseline |>
  dplyr::select(ID, log10_Mucin_DF, log10_GAG_DF, all_of(cont_outcomes)) |>
  left_join(
    df_followup |> dplyr::select(ID, all_of(cont_outcomes)),
    by = "ID", suffix = c("_bl", "_fu")
  )

prosp_results <- list()
for (i in seq_len(nrow(sig_pairs))) {
  ratio_var <- sig_pairs$ratio[i]
  outcome   <- sig_pairs$outcome[i]
  bl_col    <- paste0(outcome, "_bl")
  fu_col    <- paste0(outcome, "_fu")
  tryCatch({
    mod <- lm(reformulate(c(ratio_var, bl_col), fu_col), data = df_prosp)
    res <- summary(mod)$coefficients
    ci  <- confint(mod)
    prosp_results[[paste(ratio_var, outcome)]] <- tibble(
      ratio    = ratio_var,
      outcome  = outcome,
      n        = nobs(mod),
      estimate = res[ratio_var, 1],
      se       = res[ratio_var, 2],
      conflow  = ci[ratio_var, 1],
      confhigh = ci[ratio_var, 2],
      pval     = res[ratio_var, 4]
    )
  }, error = function(e) message("LM failed for ", ratio_var, " ", outcome, ": ", e$message))
}
prosp_df <- bind_rows(prosp_results) |>
  mutate(padj = p.adjust(pval, method = "fdr"))
print(prosp_df)
write.csv2(prosp_df,
           "results/4_functional_change/cayman/ratios/lm_prospective_outcomes.csv",
           row.names = FALSE)

# --- Helper: p-value label for subtitle ---------------------------------------
plab <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                formatC(p, format = "f", digits = 3)))
}

# --- Helper: ratio violin with per-timepoint ethnicity comparison -------------
# Facets are timepoints; within each, Dutch vs South-Asian Surinamese are
# compared with an unpaired Wilcoxon test (all samples at that timepoint).
make_ratio_vln <- function(df, ratio_var, y_label, title_label, interact_pval, show_pval = TRUE) {
  pval_annot <- df %>%
    group_by(timepoint) %>%
    group_modify(~ {
      p <- tryCatch(wilcox.test(.x[[ratio_var]] ~ .x$EthnicityTot)$p.value, error = function(e) NA_real_)
      data.frame(p = p, y_pos = max(.x[[ratio_var]], na.rm = TRUE) + diff(range(.x[[ratio_var]], na.rm = TRUE)) * 0.08)
    }) %>%
    ungroup() %>%
    mutate(label = paste0("p=", plab(p)))

  ggplot(df, aes(x = EthnicityTot, y = .data[[ratio_var]], fill = EthnicityTot)) +
    geom_violin(colour = NA, aes(alpha = timepoint)) +
    geom_boxplot(fill = "white", width = 0.15, outlier.shape = NA) +
    { if (show_pval) geom_text(data = pval_annot,
                               aes(x = 1.5, y = y_pos, label = label),
                               inherit.aes = FALSE, size = 3) } +
    facet_wrap(~timepoint) +
    scale_x_discrete(labels = \(x) str_replace(x, "South-Asian Surinamese", "South-Asian\nSurinamese")) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
    # Extra headroom above the data so the Dutch-vs-SAS p-value label isn't clipped
    scale_y_continuous(expand = expansion(mult = c(0.05, if (show_pval) 0.22 else 0.05))) +
    theme_Publication() +
    labs(x = "", y = y_label, title = title_label,
         subtitle = paste0("Ethnicity × Timepoint p=", plab(interact_pval)))
}

# --- Plot: Mucin/DF ratio by ethnicity ----------------------------------------
p_mucin <- make_ratio_vln(
  dftot, "log10_Mucin_DF", "Mucin / DF (log10 ratio)", "Mucin-to-DF ratio",
  lmm_df$pval_interact[lmm_df$ratio == "log10_Mucin_DF"]
)

# --- Plot: GAG/DF ratio by ethnicity ------------------------------------------
p_gag <- make_ratio_vln(
  dftot, "log10_GAG_DF", "GAG / DF (log10 ratio)", "GAG-to-DF ratio",
  lmm_df$pval_interact[lmm_df$ratio == "log10_GAG_DF"]
)

# --- Assemble and save --------------------------------------------------------
fig_ratios <- ggarrange(p_mucin, p_gag, nrow = 1, ncol = 2, labels = c("A", "B"))

ggsave(
  "results/4_functional_change/cayman/ratios/cayman_ratios_ethnicity.pdf",
  fig_ratios, width = 10, height = 5
)
ggsave(
  "results/4_functional_change/cayman/ratios/cayman_ratios_ethnicity.png",
  fig_ratios, width = 10, height = 5, dpi = 300
)

# --- Ratio components: which substrate group drives the ratio differences? ----
# Mucin/DF and GAG/DF can differ by ethnicity because the numerator (Mucin,
# GAG) shifts, the denominator (DF) shifts, or both. Test each group on its
# own with the same models as the ratios. The six substrate classes are
# compositional, so besides raw log10 RPKM (not compositionally valid; a uniform
# offset in total CAZyme abundance moves every class) each class is also
# expressed relative to total CAZyme RPKM and as a CLR across the six-class
# vector. The ratios themselves are scale-invariant, so they are unaffected.
substrate_classes <- c("DF", "GAG", "Mucin", "Glycogen", "PG", "Other")
plot_classes      <- c("DF", "GAG", "Mucin")

rpkm_mat <- as.matrix(dftot[paste0(substrate_classes, "_rpkm")])
stopifnot(all(rpkm_mat > 0))   # CLR needs no zeros; no pseudocount used
colnames(rpkm_mat) <- substrate_classes
log_mat <- log(rpkm_mat)

comp_scales <- list(
  log10rpkm = list(label = "log10 RPKM",
                   mat   = log10(rpkm_mat)),
  log10rel  = list(label = "log10 relative to total CAZyme RPKM",
                   mat   = log10(rpkm_mat / rowSums(rpkm_mat))),
  clr       = list(label = "CLR",
                   mat   = log_mat - rowMeans(log_mat))
)

dftot_comp <- dftot
for (sc in names(comp_scales)) {
  m <- comp_scales[[sc]]$mat
  colnames(m) <- paste0(sc, "_", substrate_classes)
  dftot_comp <- bind_cols(dftot_comp, as_tibble(m))
}

# Ethnicity x timepoint LMM for one variable (SAS vs Dutch; sensitivity
# covariates optional)
fit_ethnicity_lmm <- function(df, var, covariates = NULL) {
  f   <- reformulate(c("EthnicityTot * timepoint", covariates, "(1 | ID)"), response = var)
  mod <- lmer(f, data = df)
  res <- summary(mod)$coefficients
  main_row <- grep(paste0("^EthnicityTot", ETH_SAS, "$"), rownames(res))
  int_row  <- grep(paste0(ETH_SAS, ":timepointfollow-up"), rownames(res))
  tibble(
    variable          = var,
    estimate_main     = res[main_row, 1],
    se_main           = res[main_row, 2],
    pval_main         = res[main_row, 5],
    estimate_interact = res[int_row,  1],
    se_interact       = res[int_row,  2],
    pval_interact     = res[int_row,  5]
  )
}

# Wilcoxon (Dutch vs SAS) per timepoint, FDR within timepoint
wilcox_by_timepoint <- function(df, vars) {
  df |>
    pivot_longer(all_of(vars), names_to = "variable", values_to = "value") |>
    group_by(variable, timepoint) |>
    summarise(
      p            = wilcox.test(value ~ EthnicityTot)$p.value,
      median_Dutch = median(value[EthnicityTot == ETH_DUTCH], na.rm = TRUE),
      median_SAS   = median(value[EthnicityTot == ETH_SAS],   na.rm = TRUE),
      .groups = "drop"
    ) |>
    group_by(timepoint) |>
    mutate(padj = p.adjust(p, method = "fdr")) |>
    ungroup()
}

lmm_by_variable <- function(df, vars, covariates = NULL) {
  map(vars, \(v) tryCatch(
    fit_ethnicity_lmm(df, v, covariates),
    error = function(e) { message("LMM failed for ", v, ": ", e$message); NULL }
  )) |>
    bind_rows() |>
    mutate(
      padj_main     = p.adjust(pval_main,     method = "fdr"),
      padj_interact = p.adjust(pval_interact, method = "fdr")
    )
}

comp_results <- list()
for (sc in names(comp_scales)) {
  vars   <- paste0(sc, "_", substrate_classes)
  wil    <- wilcox_by_timepoint(dftot_comp, vars)
  lmm_sc <- lmm_by_variable(dftot_comp, vars)
  cat("\n==== Components,", comp_scales[[sc]]$label, "====\n")
  print(wil)
  print(lmm_sc)
  write.csv2(wil,
             paste0("results/4_functional_change/cayman/ratios/wilcoxon_components_", sc, "_by_ethnicity.csv"),
             row.names = FALSE)
  write.csv2(lmm_sc,
             paste0("results/4_functional_change/cayman/ratios/lmm_components_", sc, "_ethnicity_timepoint.csv"),
             row.names = FALSE)
  comp_results[[sc]] <- lmm_sc

  comp_plots <- lapply(plot_classes, function(cl) {
    v <- paste0(sc, "_", cl)
    make_ratio_vln(
      dftot_comp, v,
      paste0(cl, " (", comp_scales[[sc]]$label, ")"),
      paste0(cl, " abundance"),
      lmm_sc$pval_interact[lmm_sc$variable == v]
    )
  })
  fig_components <- ggarrange(plotlist = comp_plots, nrow = 1, ncol = 3,
                              labels = c("A", "B", "C"))
  ggsave(paste0("results/4_functional_change/cayman/ratios/cayman_components_", sc, "_ethnicity.pdf"),
         fig_components, width = 15, height = 5)
  ggsave(paste0("results/4_functional_change/cayman/ratios/cayman_components_", sc, "_ethnicity.png"),
         fig_components, width = 15, height = 5, dpi = 300)
}

# Side-by-side: SAS-vs-Dutch main effect per class on each scale
comp_summary <- bind_rows(comp_results, .id = "scale") |>
  mutate(class = str_remove(variable, "^[a-z0-9]+_")) |>
  dplyr::select(scale, class, estimate_main, se_main, pval_main, padj_main)
print(comp_summary)
write.csv2(comp_summary,
           "results/4_functional_change/cayman/ratios/components_main_effect_all_scales.csv",
           row.names = FALSE)

# --- Secondary checks: CAZyme family richness and technical covariates --------
# Richness is families detected (n_families, sample_statistics.tsv). Technical
# variables are checked for ethnic differences to rule out the uniform offset
# in DF/GAG/Mucin being a depth, host-read or alignment artefact.
stats <- rio::import("data/shotgun/cayman_results/oct2026_results/sample_statistics.tsv") |>
  as_tibble() |>
  dplyr::rename(sampleID = sample) |>
  dplyr::select(sampleID, aligned_reads, filtered_reads, cazy_reads, pct_cazy_reads,
                richness = n_families)

# Per-sample read pairs entering and leaving host (human) removal, from MultiQC
# Bowtie2 summaries (nf-core/mag host removal step, as in qc_summary.R)
host_reads <- map_dfr(1:3, function(b) {
  y <- yaml::read_yaml(sprintf("data/shotgun/multiqc_data_%d/multiqc_bowtie2_bowtie2-1.yaml", b))
  map_dfr(str_subset(names(y), "^HELI"), \(nm) tibble(
    sampleID        = str_remove(nm, "_run[0-9]+$"),
    trimmed_pairs   = y[[nm]][["paired_total"]],
    post_host_pairs = y[[nm]][["paired_aligned_none"]]
  ))
}) |>
  mutate(human_frac = 1 - post_host_pairs / trimmed_pairs)
stopifnot(!anyDuplicated(host_reads$sampleID))

dftot_qc <- dftot_comp |>
  left_join(stats,      by = "sampleID") |>
  left_join(host_reads, by = "sampleID") |>
  mutate(
    log10_trimmed_pairs   = log10(trimmed_pairs),
    log10_post_host_pairs = log10(post_host_pairs),
    log10_human_frac      = log10(human_frac),
    log10_aligned_reads   = log10(aligned_reads),
    log10_cazy_reads      = log10(cazy_reads),
    aligned_per_pair      = aligned_reads / post_host_pairs,
    filtered_per_aligned  = filtered_reads / aligned_reads
  )
stopifnot(!anyNA(dftot_qc$richness), !anyNA(dftot_qc$post_host_pairs))

qc_labels <- c(
  richness              = "CAZyme family richness",
  log10_trimmed_pairs   = "Trimmed read pairs (log10)",
  log10_post_host_pairs = "Read pairs after host removal (log10)",
  log10_human_frac      = "Human read fraction (log10)",
  log10_aligned_reads   = "Aligned reads (log10)",
  pct_cazy_reads        = "CAZyme reads (% of filtered reads)",
  aligned_per_pair      = "Aligned reads per post-host pair",
  filtered_per_aligned  = "Filtered / aligned reads"
)
qc_vars <- names(qc_labels)

wilcox_qc <- wilcox_by_timepoint(dftot_qc, qc_vars)
lmm_qc    <- lmm_by_variable(dftot_qc, qc_vars)
print(wilcox_qc)
print(lmm_qc)
write.csv2(wilcox_qc,
           "results/4_functional_change/cayman/ratios/wilcoxon_richness_qc_by_ethnicity.csv",
           row.names = FALSE)
write.csv2(lmm_qc,
           "results/4_functional_change/cayman/ratios/lmm_richness_qc_ethnicity_timepoint.csv",
           row.names = FALSE)

# Richness scales with sequencing effort, so repeat with CAZyme read depth
# as a covariate
lmm_richness_adj <- lmm_by_variable(dftot_qc, "richness", covariates = "log10_cazy_reads")
print(lmm_richness_adj)
write.csv2(lmm_richness_adj,
           "results/4_functional_change/cayman/ratios/lmm_richness_adj_cazy_depth.csv",
           row.names = FALSE)

# Do the DF/GAG/Mucin ethnicity effects survive adjusting for depth and host reads?
adj_covariates <- c("log10_post_host_pairs", "log10_human_frac", "pct_cazy_reads")
lmm_comp_adj <- lmm_by_variable(dftot_qc,
                                paste0("clr_", plot_classes),
                                covariates = adj_covariates)
print(lmm_comp_adj)
write.csv2(lmm_comp_adj,
           "results/4_functional_change/cayman/ratios/lmm_components_clr_adj_technical.csv",
           row.names = FALSE)

qc_plots <- lapply(qc_vars, function(v) {
  make_ratio_vln(dftot_qc, v, qc_labels[[v]], qc_labels[[v]],
                 lmm_qc$pval_interact[lmm_qc$variable == v])
})
fig_qc <- ggarrange(plotlist = qc_plots, nrow = 2, ncol = 4, labels = LETTERS[seq_along(qc_vars)])
ggsave("results/4_functional_change/cayman/ratios/cayman_richness_qc_ethnicity.pdf",
       fig_qc, width = 20, height = 10)
ggsave("results/4_functional_change/cayman/ratios/cayman_richness_qc_ethnicity.png",
       fig_qc, width = 20, height = 10, dpi = 300)

cat("Done. Figures saved to results/4_functional_change/cayman/ratios/\n")
