## Cayman Mucin/DF and GAG/DF ratios by ethnicity
library(tidyverse)
library(readxl)
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

# --- Load data ----------------------------------------------------------------
anno <- read_excel("data/shotgun/cayman_results/mucin_df_gag_table.xlsx") |>
  dplyr::select(Family, FUNCTION_AT_DESTINATION_1)

df_raw <- rio::import("data/shotgun/cayman_results/families_cpm_table.tsv") |>
  dplyr::select(-HELIBA_103370, -HELIFU_103370)

clinical <- readRDS("data/clinicaldata/clinicaldata_long.RDS")


# --- Assign families to functional groups -------------------------------------
# A family with comma-separated annotations (e.g. "DF,Mucin") belongs to all
# listed groups.
anno_clean <- anno |>
  filter(!is.na(Family), !is.na(FUNCTION_AT_DESTINATION_1)) |>
  mutate(
    is_DF    = str_detect(FUNCTION_AT_DESTINATION_1, "\\bDF\\b"),
    is_Mucin = str_detect(FUNCTION_AT_DESTINATION_1, "\\bMucin\\b"),
    is_GAG   = str_detect(FUNCTION_AT_DESTINATION_1, "\\bGAG\\b")
  )

families_DF    <- anno_clean |> filter(is_DF)    |> pull(Family) |> unique()
families_Mucin <- anno_clean |> filter(is_Mucin) |> pull(Family) |> unique()
families_GAG   <- anno_clean |> filter(is_GAG)   |> pull(Family) |> unique()

cat("DF families (n =", length(families_DF), "):", paste(families_DF, collapse = ", "), "\n")
cat("Mucin families (n =", length(families_Mucin), "):", paste(families_Mucin, collapse = ", "), "\n")
cat("GAG families (n =", length(families_GAG), "):", paste(families_GAG, collapse = ", "), "\n")

# --- Compute per-sample group sums (CPM) --------------------------------------
all_families <- df_raw$family
df_mat <- df_raw |>
  column_to_rownames("family") |>
  as.matrix()

sum_group <- function(mat, families) {
  keep <- intersect(families, rownames(mat))
  if (length(keep) == 0) return(rep(0, ncol(mat)))
  if (length(keep) == 1) return(mat[keep, ])
  colSums(mat[keep, ])
}

df_sums <- tibble(
  sampleID = colnames(df_mat),
  CPM_DF    = sum_group(df_mat, families_DF),
  CPM_Mucin = sum_group(df_mat, families_Mucin),
  CPM_GAG   = sum_group(df_mat, families_GAG)
)

df_sums <- df_sums |>
  mutate(
    ratio_Mucin_DF = CPM_Mucin / CPM_DF,
    ratio_GAG_DF   = CPM_GAG   / CPM_DF,
    log10_Mucin_DF = log10(ratio_Mucin_DF),
    log10_GAG_DF   = log10(ratio_GAG_DF)
  )

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

  pval_m <- r |>
    dplyr::select(ratio, outcome, padj) |>
    pivot_wider(names_from = ratio, values_from = padj) |>
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
  labels = c("q < 0.05", "q < 0.01", "q < 0.001"),
  legend_gp = gpar(fontsize = 10)
)
lgd_packed <- packLegend(lgd_cor, lgd_sig, direction = "vertical", gap = unit(4, "mm"))

CairoPDF("results/4_functional_change/cayman/ratios/heatmap_spearman_outcomes.pdf",
         width = 7, height = 6)
draw(ht_bl + ht_fu, annotation_legend_list = list(lgd_packed),
     padding = unit(c(5, 30, 5, 5), "mm"))
dev.off()

# --- Spearman correlations: baseline ratios vs dietary intake -----------------
# Macronutrients are energy-adjusted using the Willett residual method: each
# macronutrient (g/day) is regressed on total energy intake (kcal/day) within
# the shotgun baseline subset, and the residuals are z-scored. This removes the
# confounding effect of overall energy intake, so the adjusted variable reflects
# macronutrient composition independent of how much a participant eats in total.
# Total calories is kept on its original scale as a separate predictor.
macro_vars_raw <- c("Protein", "Protein_animal", "FattyAcids", "MonoUnsatFat",
                    "PolyUnsatFat", "SatFat", "Carbohydrates", "Fiber", "Sodium_g")

df_baseline_diet <- df_baseline
for (mac in macro_vars_raw) {
  col_adj <- paste0(mac, "_adj")
  idx <- !is.na(df_baseline_diet[[mac]]) & !is.na(df_baseline_diet$TotalCalories)
  resid_vec <- rep(NA_real_, nrow(df_baseline_diet))
  if (sum(idx) > 2) {
    fit <- lm(df_baseline_diet[[mac]][idx] ~ df_baseline_diet$TotalCalories[idx])
    resid_vec[idx] <- residuals(fit)
  }
  df_baseline_diet[[col_adj]] <- as.numeric(scale(resid_vec))
}

diet_vars <- c("TotalCalories", paste0(macro_vars_raw, "_adj"))
diet_vars <- diet_vars[diet_vars %in% names(df_baseline_diet)]

diet_labels <- c(
  TotalCalories          = "Total calories",
  Protein_adj            = "Protein",
  Protein_animal_adj     = "Animal protein",
  FattyAcids_adj         = "Fatty acids",
  MonoUnsatFat_adj       = "Mono-unsat. fat",
  PolyUnsatFat_adj       = "Poly-unsat. fat",
  SatFat_adj             = "Saturated fat",
  Carbohydrates_adj      = "Carbohydrates",
  Fiber_adj              = "Fiber",
  Sodium_g_adj           = "Sodium"
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

  pval_m <- r |>
    dplyr::select(ratio, outcome, padj) |>
    pivot_wider(names_from = ratio, values_from = padj) |>
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
  if (is.na(p)) return("NA")
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  formatC(p, format = "f", digits = 3)
}

# --- Helper: ratio violin with per-facet paired Wilcoxon annotation -----------
make_ratio_vln <- function(df, ratio_var, y_label, title_label, interact_pval, show_pval = TRUE) {
  df_plot <- df %>%
    group_by(EthnicityTot, ID) %>%
    filter(n() == 2) %>%
    ungroup() %>%
    arrange(EthnicityTot, ID, timepoint)

  pval_annot <- df_plot %>%
    group_by(EthnicityTot) %>%
    group_modify(~ {
      bl <- .x[[ratio_var]][.x$timepoint == "baseline"]
      fu <- .x[[ratio_var]][.x$timepoint == "follow-up"]
      p  <- tryCatch(wilcox.test(bl, fu, paired = TRUE)$p.value, error = function(e) NA_real_)
      data.frame(p = p, y_pos = max(.x[[ratio_var]], na.rm = TRUE) + diff(range(.x[[ratio_var]], na.rm = TRUE)) * 0.08)
    }) %>%
    ungroup() %>%
    mutate(label = ifelse(p < 0.05,
                          paste0("p=", ifelse(p < 0.001,
                                              formatC(p, format = "e", digits = 2),
                                              formatC(p, format = "f", digits = 3))),
                          ""))

  ggplot(df_plot, aes(x = timepoint, y = .data[[ratio_var]], fill = EthnicityTot)) +
    geom_violin(colour = NA, aes(alpha = timepoint)) +
    geom_boxplot(fill = "white", width = 0.15, outlier.shape = NA) +
    { if (show_pval) geom_text(data = pval_annot,
                               aes(x = 1.5, y = y_pos, label = label),
                               inherit.aes = FALSE, size = 3) } +
    facet_wrap(~EthnicityTot) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
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

cat("Done. Figures saved to results/4_functional_change/cayman/ratios/\n")
