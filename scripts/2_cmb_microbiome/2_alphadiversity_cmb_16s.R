## Figure 2 — Cardiometabolic disease and microbiome instability (16S)
## Alpha-diversity analyses
##
## Output (supplementary figure):
##   pl_shan_combined → supplementary_alphadiversity.pdf
##   (extended LMM forest plot + prevalence bar companion)

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(lme4)
library(afex)

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

#### Load data ####
df <- readRDS("data/16s/clin_alphadiversity.RDS") %>%
    dplyr::select(1, sampleID = ID, 4:5) %>%
    mutate(timepoint = case_when(
                        str_detect(sampleID, "HELIBA") ~ "baseline",
                        str_detect(sampleID, "HELIFU") ~ "follow-up"
                    ),
        ID = str_remove(str_remove(sampleID, "HELIFU_"), "HELIBA_"),
        ID = str_c("S", ID))

dfwide <- df %>% pivot_wider(., id_cols = "ID", names_from = "timepoint",
                             values_from = c(1,3,4)) %>%
    mutate(shannon_delta = `shannon_follow-up` - shannon_baseline,
           PD_delta = `PD_follow-up` - PD_baseline,
           richness_delta = `richness_follow-up` - richness_baseline) %>%
    dplyr::select(1, 8:10)
df2 <- left_join(df, dfwide, by = "ID")

helius <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
dftot <- left_join(df2 %>% filter(timepoint == "baseline"), helius %>% filter(timepoint == "baseline"),
                   by = c("ID", "timepoint", "sampleID")) %>%
    filter(EthnicityTot != "Other") %>% droplevels(.)

pairedids <- df2$sampleID[which(!is.na(df2$shannon_delta))]
dftot2 <- left_join(df, helius, by = c("ID", "timepoint", "sampleID")) %>%
            filter(sampleID %in% pairedids) %>%
            filter(! EthnicityTot %in% "Other" ) %>% droplevels(.)

#### Output folder ####
resultsfolder <- "results/2_cmb_microbiome"
dir.create(resultsfolder, showWarnings = FALSE)

#### Extended forest plot: multi-domain predictors of Shannon change ####

# Helper: extract disease x timepoint interaction from LMM
extract_lmm_int <- function(data, outcome_var, predictor_var, label, group = NULL) {
    df <- data %>%
        filter(!is.na(.data[[predictor_var]])) %>%
        mutate(
            outcome   = .data[[outcome_var]],
            predictor = .data[[predictor_var]],
            timepoint = factor(timepoint, levels = c("baseline", "follow-up"))
        )
    model   <- lmer(outcome ~ predictor * timepoint + FUtime + (1|ID), data = df)
    cf      <- coef(summary(model))
    ci      <- confint(model, method = "Wald")
    int_row_cf <- grep(":timepoint", rownames(cf), value = TRUE)[1]
    int_row_ci <- grep(":timepoint", rownames(ci), value = TRUE)[1]
    data.frame(
        label     = label,
        group     = group,
        estimate  = cf[int_row_cf, "Estimate"],
        conf.low  = ci[int_row_ci, 1],
        conf.high = ci[int_row_ci, 2],
        p.value   = cf[int_row_cf, "Pr(>|t|)"],
        stringsAsFactors = FALSE
    )
}

# Load wide clinical data and select baseline predictors
helius_wide <- readRDS("data/clinicaldata/clinicaldata_wide.RDS")

bl_vars <- c("DM", "HT_BPMed", "Dyslipidemia",
             "Statins", "Metformin", "AntiHT", "PPI", "GlucLowDrugs", "PsychoMed", "Cortico",
             "Smoking_current", "Alcohol", "ExerciseNorm",
             "Age", "BMI",
             "Protein", "FattyAcids", "Carbohydrates", "Fiber", "Sodium_g",
             "DiscrMean")

bl_cols <- paste0(bl_vars, "_baseline")
bl_cols <- bl_cols[bl_cols %in% names(helius_wide)]

helius_bl <- helius_wide %>%
    dplyr::select(ID, all_of(bl_cols)) %>%
    rename_with(~ sub("_baseline$", "_bl", .), all_of(bl_cols))

dftot2_ext <- dftot2 %>%
    left_join(helius_bl, by = "ID") %>%
    mutate(across(
        c(Age_bl, BMI_bl, Protein_bl, FattyAcids_bl, Carbohydrates_bl, Fiber_bl, Sodium_g_bl,
          DiscrMean_bl),
        ~ as.numeric(scale(.))
    )) %>%
    mutate(across(
        c(DM_bl, HT_BPMed_bl, Dyslipidemia_bl, Statins_bl, Metformin_bl, AntiHT_bl,
          PPI_bl, GlucLowDrugs_bl, PsychoMed_bl, Cortico_bl,
          Smoking_current_bl, Alcohol_bl, ExerciseNorm_bl),
        ~ relevel(factor(.), ref = "No"),
        .names = "{.col}"
    )) %>%
    mutate(Sex = relevel(factor(Sex), ref = "Male"))

shan_ext_effects <- bind_rows(
    # Risk factors
    extract_lmm_int(dftot2_ext, "shannon", "Age_bl",             "Age (per SD)",                      group = "Risk factors"),
    extract_lmm_int(dftot2_ext, "shannon", "Sex",                "Sex (female)",                      group = "Risk factors"),
    extract_lmm_int(dftot2_ext, "shannon", "BMI_bl",             "BMI (per SD)",                      group = "Risk factors"),
    extract_lmm_int(dftot2_ext, "shannon", "Smoking_current_bl", "Current smoking",                   group = "Risk factors"),
    extract_lmm_int(dftot2_ext, "shannon", "Alcohol_bl",         "Alcohol use",                       group = "Risk factors"),
    extract_lmm_int(dftot2_ext, "shannon", "ExerciseNorm_bl",    "Sufficient exercise",               group = "Risk factors"),
    extract_lmm_int(dftot2_ext, "shannon", "DiscrMean_bl",       "Perceived discrimination (per SD)", group = "Risk factors"),
    # Cardiometabolic disease
    extract_lmm_int(dftot2_ext, "shannon", "DM_bl",           "Diabetes",     group = "Disease"),
    extract_lmm_int(dftot2_ext, "shannon", "HT_BPMed_bl",     "Hypertension", group = "Disease"),
    extract_lmm_int(dftot2_ext, "shannon", "Dyslipidemia_bl", "Dyslipidemia", group = "Disease"),
    # Medications
    extract_lmm_int(dftot2_ext, "shannon", "Statins_bl",      "Statins",                group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "Metformin_bl",    "Metformin",              group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "AntiHT_bl",       "Antihypertensives",      group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "PPI_bl",          "PPI",                    group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "GlucLowDrugs_bl", "Glucose-lowering drugs", group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "PsychoMed_bl",    "Psychotropics",          group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "Cortico_bl",      "Corticosteroids",        group = "Medication"),
    # Diet
    extract_lmm_int(dftot2_ext, "shannon", "Protein_bl",       "Protein (per SD)",       group = "Diet"),
    extract_lmm_int(dftot2_ext, "shannon", "FattyAcids_bl",    "Fatty acids (per SD)",   group = "Diet"),
    extract_lmm_int(dftot2_ext, "shannon", "Carbohydrates_bl", "Carbohydrates (per SD)", group = "Diet"),
    extract_lmm_int(dftot2_ext, "shannon", "Fiber_bl",         "Fiber (per SD)",         group = "Diet"),
    extract_lmm_int(dftot2_ext, "shannon", "Sodium_g_bl",      "Sodium (per SD)",        group = "Diet")
) %>%
    mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        sig   = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR \u2265 0.05"),
        label = factor(label, levels = rev(unique(label))),
        group = factor(group, levels = c("Risk factors", "Disease", "Medication", "Diet"))
    )

(pl_shan_extended <- ggplot(shan_ext_effects, aes(x = estimate, y = label, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.2) +
    geom_point(size = 3) +
    scale_color_manual(
        values = c("FDR < 0.05" = "#e63946", "FDR \u2265 0.05" = "grey55"),
        name = NULL
    ) +
    facet_wrap(~ group, scales = "free_y", ncol = 1) +
    labs(x = "Shannon change: interaction estimate (\u00b1 95% CI)",
         y = NULL,
         title = "Predictors of Shannon diversity change") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_shannon_change_extended.pdf"), width = 6, height = 10)

# Companion bar panel
var_meta_shan <- data.frame(
    predictor = c(
        "Age_bl", "Sex", "BMI_bl",
        "Smoking_current_bl", "Alcohol_bl", "ExerciseNorm_bl", "DiscrMean_bl",
        "DM_bl", "HT_BPMed_bl", "Dyslipidemia_bl",
        "Statins_bl", "Metformin_bl", "AntiHT_bl", "PPI_bl",
        "GlucLowDrugs_bl", "PsychoMed_bl", "Cortico_bl",
        "Protein_bl", "FattyAcids_bl", "Carbohydrates_bl",
        "Fiber_bl", "Sodium_g_bl"
    ),
    label = c(
        "Age (per SD)", "Sex (female)", "BMI (per SD)",
        "Current smoking", "Alcohol use", "Sufficient exercise", "Perceived discrimination (per SD)",
        "Diabetes", "Hypertension", "Dyslipidemia",
        "Statins", "Metformin", "Antihypertensives", "PPI",
        "Glucose-lowering drugs", "Psychotropics", "Corticosteroids",
        "Protein (per SD)", "Fatty acids (per SD)", "Carbohydrates (per SD)",
        "Fiber (per SD)", "Sodium (per SD)"
    ),
    group = c(
        rep("Risk factors", 7), rep("Disease", 3), rep("Medication", 7), rep("Diet", 5)
    ),
    type = c(
        "continuous", "continuous", "continuous",
        rep("binary", 3), "continuous",
        rep("binary", 10),
        rep("continuous", 5)
    ),
    stringsAsFactors = FALSE
)

n_total_shan <- sum(dftot2_ext$timepoint == "baseline")

bar_data_shan <- purrr::map_dfr(seq_len(nrow(var_meta_shan)), function(i) {
    col <- var_meta_shan$predictor[i]
    if (!col %in% names(dftot2_ext)) return(NULL)
    x <- dftot2_ext[[col]]
    x <- x[dftot2_ext$timepoint == "baseline"]
    if (var_meta_shan$type[i] == "binary") {
        n_yes <- sum(x == "Yes", na.rm = TRUE)
        n_no  <- sum(x == "No",  na.rm = TRUE)
        n_tot <- n_yes + n_no
        data.frame(
            label    = var_meta_shan$label[i],
            group    = var_meta_shan$group[i],
            category = c("Yes", "No"),
            pct      = c(n_yes / n_tot * 100, n_no / n_tot * 100),
            n_label  = c(n_yes, NA_integer_)
        )
    } else {
        n_valid   <- sum(!is.na(x))
        pct_valid <- n_valid / n_total_shan * 100
        data.frame(
            label    = var_meta_shan$label[i],
            group    = var_meta_shan$group[i],
            category = "Non-missing",
            pct      = pct_valid,
            n_label  = n_valid
        )
    }
}) %>%
    mutate(
        label    = factor(label, levels = levels(shan_ext_effects$label)),
        group    = factor(group, levels = c("Risk factors", "Disease", "Medication", "Diet")),
        category = factor(category, levels = c("No", "Yes", "Non-missing"))
    )

pl_shan_bar <- ggplot(bar_data_shan, aes(x = pct, y = label, fill = category)) +
    geom_col(width = 0.6) +
    geom_text(
        data = bar_data_shan %>% filter(!is.na(n_label)),
        aes(label = paste0("n=", n_label), x = 101),
        hjust = 0, size = 2.8, color = "grey30"
    ) +
    scale_fill_manual(
        values = c("Yes" = "#4a90d9", "No" = "#d0d8e4", "Non-missing" = "#4a90d9"),
        guide = "none"
    ) +
    scale_x_continuous(limits = c(0, 140), breaks = c(0, 50, 100),
                       expand = expansion(mult = c(0, 0))) +
    facet_wrap(~ group, scales = "free_y", ncol = 1) +
    labs(x = "% Yes / % non-missing", y = NULL, title = " ") +
    theme_Publication() +
    theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position = "bottom"
    )

pl_shan_combined <- pl_shan_extended %>% aplot::insert_right(pl_shan_bar, width = 0.4)
ggsave(pl_shan_combined, filename = file.path(resultsfolder, "effectsize_shannon_change_extended_combined.pdf"),
       width = 9, height = 10)

#### Supplementary figure — Alpha-diversity × cardiometabolic disease ####

ggsave(pl_shan_combined, filename = file.path(resultsfolder, "supplementary_alphadiversity.pdf"),
       width = 9, height = 10, device = "pdf")
