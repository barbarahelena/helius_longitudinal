## Alistipes putredinis microdiversity — multi-domain covariate screen
## Two analyses, mirroring the project's existing alpha-diversity /
## strain-sharing covariate forest plots:
##   A. Cross-sectional: baseline nucl_diversity ~ covariate (lm, like
##      6_strain_stability/strainsharing_covariates.R) on the full QC-passing
##      population (n ~117-120) — good power, but can't speak to change.
##   B. Longitudinal: nucl_diversity ~ covariate * timepoint + FUtime +
##      (1|subject_id) (lmer, like 2_cmb_microbiome/2_alphadiversity_cmb_16s.R),
##      extracting the covariate:timepoint interaction, restricted to RETAINED
##      strains only (same_strain == TRUE, n ~22-24) — this is the direct
##      "does this covariate predict the within-lineage diversity trajectory"
##      test, but is severely underpowered for binary covariates with few
##      exposed participants; treat as hypothesis-generating only.
## Reads results/3_species_change/5_instrain/instrain_microdiversity.csv,
## written by 3_instrain_microdiversity.R (run that first).
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

#### Libraries ####
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0),
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
}

#### Paths ####
microdiv_path <- "results/3_species_change/5_instrain/instrain_microdiversity.csv"
out_dir <- "results/3_species_change/5_instrain/covariates"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(microdiv_path))
    stop("Run 3_instrain_microdiversity.R first — ", microdiv_path, " not found.")

#### Load microdiversity data and attach clinical covariates ####
microdiv <- read.csv(microdiv_path, colClasses = c(subject_id = "character")) %>%
    mutate(ID = str_c("S", subject_id))

clin_long <- readRDS("data/clinicaldata/clinicaldata_long.RDS") %>%
    filter(timepoint == "baseline") %>%
    dplyr::select(ID, Sex, DM, HT_BPMed, Dyslipidemia, Statins, Metformin, AntiHT, PPI,
                  PsychoMed, Cortico, Smoking_current, Alcohol, ExerciseNorm)

clin_wide <- readRDS("data/clinicaldata/clinicaldata_wide.RDS")
diet_cols <- c("TotalCalories_baseline",
              "Protein_baseline_adj", "FattyAcids_baseline_adj",
              "Carbohydrates_baseline_adj", "Fiber_baseline_adj", "Sodium_g_baseline_adj")
diet <- clin_wide %>%
    dplyr::select(ID, all_of(diet_cols)) %>%
    rename_with(~ sub("_baseline_adj$", "_adj_bl", sub("_baseline$", "_bl", .)), -ID)

dat <- microdiv %>%
    left_join(clin_long, by = "ID") %>%
    left_join(diet, by = "ID")
cat("Microdiversity rows:", nrow(microdiv), " | after clinical join:", nrow(dat), "\n")

#### Covariate metadata: shared across both analyses ####
# Age/BMI/FUtime already on microdiv (from tip_meta_clades.RDS, see
# 3_instrain_microdiversity.R); EthnicityTot too, but that's covered
# separately (instrain_microdiversity_by_ethnicity.csv / the LMM in that
# script), so it's left out of this generic covariate screen.
var_meta <- data.frame(
    var   = c("Age", "Sex", "BMI", "Alcohol", "Smoking_current", "ExerciseNorm",
              "DM", "HT_BPMed", "Dyslipidemia",
              "Statins", "Metformin", "AntiHT", "PPI", "PsychoMed", "Cortico",
              "TotalCalories_bl", "Protein_adj_bl", "FattyAcids_adj_bl",
              "Carbohydrates_adj_bl", "Fiber_adj_bl", "Sodium_g_adj_bl"),
    label = c("Age (per SD)", "Sex (female)", "BMI (per SD)", "Alcohol use",
              "Current smoking", "Sufficient exercise",
              "Diabetes", "Hypertension", "Dyslipidemia",
              "Statins", "Metformin", "Antihypertensives", "PPI", "Psychotropics", "Corticosteroids",
              "Total calories (per SD)", "Protein (per SD)", "Fatty acids (per SD)",
              "Carbohydrates (per SD)", "Fiber (per SD)", "Sodium (per SD)"),
    group = c(rep("Risk factors", 6), rep("Disease", 3), rep("Medication", 6), rep("Diet", 6)),
    type  = c("continuous", "binary", "continuous", "binary", "binary", "binary",
              rep("binary", 9), rep("continuous", 6)),
    stringsAsFactors = FALSE
)

binary_vars    <- var_meta$var[var_meta$type == "binary"]
binary_vars_no <- binary_vars[binary_vars != "Sex"]

dat_scaled <- dat %>%
    mutate(
        across(all_of(var_meta$var[var_meta$type == "continuous"]), ~ as.numeric(scale(.))),
        across(all_of(intersect(binary_vars_no, names(dat))), ~ relevel(factor(.), ref = "No")),
        Sex = relevel(factor(Sex), ref = "Male")
    )

############################################################
#### A. Cross-sectional: baseline nucl_diversity ~ covariate (lm) ####
############################################################
# log(nucl_diversity_baseline), matching the log-outcome convention already
# used for the coverage-adjustment check in 3_instrain_microdiversity.R.
# Adjusted for FUtime (dropped when FUtime itself is the predictor) — there
# is no "same strain" restriction here: baseline diversity is a property of
# whichever population was present at baseline, regardless of what happened
# at follow-up, so the full QC-passing set is used for power.

extract_lm_cross <- function(data, predictor_var, label, group) {
    adj <- c("FUtime")
    adj <- adj[adj != predictor_var]
    rhs <- paste(c("predictor", adj), collapse = " + ")
    fml <- as.formula(paste("log(nucl_diversity_baseline) ~", rhs))

    df_use <- data %>%
        filter(!is.na(.data[[predictor_var]])) %>%
        mutate(predictor = .data[[predictor_var]])
    if (length(adj) > 0)
        df_use <- df_use %>% filter(complete.cases(dplyr::select(., all_of(adj))))
    if (is.factor(df_use$predictor) && min(table(df_use$predictor)) < 5) return(NULL)
    if (nrow(df_use) < 20) return(NULL)

    tryCatch({
        m  <- lm(fml, data = df_use)
        cf <- coef(summary(m))
        ci <- confint(m)
        data.frame(label = label, group = group,
                   estimate = cf[2, "Estimate"], conf.low = ci[2, 1], conf.high = ci[2, 2],
                   p.value = cf[2, "Pr(>|t|)"], n = nrow(df_use), stringsAsFactors = FALSE)
    }, error = function(e) NULL)
}

cross_effects <- purrr::map_dfr(seq_len(nrow(var_meta)), function(i)
    extract_lm_cross(dat_scaled, var_meta$var[i], var_meta$label[i], var_meta$group[i])
) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"),
           sig   = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR >= 0.05"),
           label = factor(label, levels = rev(var_meta$label)),
           group = factor(group, levels = c("Risk factors", "Disease", "Medication", "Diet")))

cat("\n=== A. Cross-sectional: covariates of BASELINE nucl_diversity ===\n")
print(cross_effects %>% arrange(p.value) %>% dplyr::select(label, group, estimate, p.value, p.adj, n))
write.csv(cross_effects, file.path(out_dir, "microdiversity_crosssectional_covariates.csv"), row.names = FALSE)

pl_cross <- ggplot(cross_effects, aes(x = estimate, y = label, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("FDR < 0.05" = "#F05C3BFF", "FDR >= 0.05" = "grey55"), name = NULL) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "log(nucl_diversity) at baseline: coefficient (95% CI)", y = NULL,
         title = "A. Baseline microdiversity") +
    theme_Publication() +
    theme(legend.position = "bottom", strip.text.y = element_blank(), strip.background.y = element_blank())

############################################################
#### B. Longitudinal: does the covariate predict the CHANGE over time? ####
############################################################
# Restricted to retained strains (same_strain == TRUE): this is the
# population where "change over time" is a within-lineage trajectory rather
# than a mix of persistence and replacement (see 3_instrain_microdiversity.R).
# n ~22-24 — most binary covariates will have too few exposed participants
# to fit; those are dropped by the guards below rather than reported with an
# unreliable estimate.

long_data <- dat_scaled %>%
    filter(same_strain) %>%
    dplyr::select(subject_id, nucl_diversity_baseline, nucl_diversity_followup, FUtime,
                  all_of(var_meta$var)) %>%
    pivot_longer(c(nucl_diversity_baseline, nucl_diversity_followup),
                 names_to = "timepoint", values_to = "nucl_diversity") %>%
    mutate(timepoint = factor(if_else(timepoint == "nucl_diversity_baseline", "baseline", "follow-up"),
                              levels = c("baseline", "follow-up")))
cat("\nRetained-strain participants available for the longitudinal screen:",
    n_distinct(long_data$subject_id), "\n")

extract_lmm_int <- function(data, predictor_var, label, group) {
    df_use <- data %>%
        filter(!is.na(.data[[predictor_var]])) %>%
        mutate(predictor = .data[[predictor_var]])
    if (is.factor(df_use$predictor) && min(table(df_use$predictor[df_use$timepoint == "baseline"])) < 5) return(NULL)
    if (n_distinct(df_use$subject_id) < 10) return(NULL)

    model <- tryCatch(
        lmer(log(nucl_diversity) ~ predictor * timepoint + FUtime + (1|subject_id), data = df_use),
        error = function(e) NULL)
    if (is.null(model)) return(NULL)
    cf <- coef(summary(model))
    int_row <- grep(":timepoint", rownames(cf), value = TRUE)[1]
    if (is.na(int_row)) return(NULL)
    ci <- tryCatch(confint(model, method = "Wald", parm = int_row), error = function(e) NULL)
    if (is.null(ci)) return(NULL)
    data.frame(label = label, group = group,
               estimate = cf[int_row, "Estimate"], conf.low = ci[1, 1], conf.high = ci[1, 2],
               p.value = cf[int_row, "Pr(>|t|)"], n = n_distinct(df_use$subject_id),
               stringsAsFactors = FALSE)
}

long_effects <- purrr::map_dfr(seq_len(nrow(var_meta)), function(i)
    extract_lmm_int(long_data, var_meta$var[i], var_meta$label[i], var_meta$group[i])
)

if (nrow(long_effects) > 0) {
    long_effects <- long_effects %>%
        mutate(p.adj = p.adjust(p.value, method = "BH"),
               sig   = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR >= 0.05"),
               label = factor(label, levels = rev(var_meta$label)),
               group = factor(group, levels = c("Risk factors", "Disease", "Medication", "Diet")))
} else {
    cat("\nNo covariate had enough retained-strain participants to fit the interaction model.\n")
}

cat("\n=== B. Longitudinal: covariates of the retained-strain nucl_diversity TRAJECTORY ===\n")
cat("(", nrow(var_meta) - nrow(long_effects), "of", nrow(var_meta),
    "covariates dropped for too few exposed participants in this n=",
    n_distinct(long_data$subject_id), "subset)\n")
if (nrow(long_effects) > 0)
    print(long_effects %>% arrange(p.value) %>% dplyr::select(label, group, estimate, p.value, p.adj, n))
write.csv(long_effects, file.path(out_dir, "microdiversity_longitudinal_covariates.csv"), row.names = FALSE)

if (nrow(long_effects) > 0) {
    pl_long <- ggplot(long_effects, aes(x = estimate, y = label, color = sig)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
        geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.2) +
        geom_point(size = 3) +
        scale_color_manual(values = c("FDR < 0.05" = "#F05C3BFF", "FDR >= 0.05" = "grey55"), name = NULL) +
        facet_grid(group ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
        labs(x = "nucl_diversity change: predictor x timepoint coefficient (95% CI)", y = NULL,
             title = sprintf("B. Change over time (retained strains, n=%d)", n_distinct(long_data$subject_id))) +
        theme_Publication() +
        theme(legend.position = "bottom", strip.text.y = element_blank(), strip.background.y = element_blank())

    (pl_combined <- ggarrange(pl_cross, pl_long, ncol = 2, labels = c("A", "B")))
    ggsave(file.path(out_dir, "microdiversity_covariates.pdf"), pl_combined, width = 13, height = 10)
} else {
    ggsave(file.path(out_dir, "microdiversity_covariates.pdf"), pl_cross, width = 7, height = 10)
}
cat("\nPlot saved to:", file.path(out_dir, "microdiversity_covariates.pdf"), "\n")
