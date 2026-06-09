## Figure 2 — Cardiometabolic disease and microbiome change (16S)
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

#### Ethnicity colour palette ####
eth_colors <- c(
    "Dutch"                  = "#709AE1FF",
    "South-Asian Surinamese" = "#FED439FF",
    "African Surinamese"     = "#8A9197FF",
    "Ghanaian"               = "#D2AF81FF",
    "Turkish"                = "#FD7446FF",
    "Moroccan"               = "#D5E4A2FF"
)

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
             "Statins", "Metformin", "AntiHT", "PPI", "PsychoMed", "Cortico",
             "Smoking_current", "Alcohol", "ExerciseNorm",
             "Age", "BMI", "TotalCalories",
             "Protein", "FattyAcids",
             "Carbohydrates", "Fiber", "Sodium_g")

bl_cols <- paste0(bl_vars, "_baseline")
bl_cols <- bl_cols[bl_cols %in% names(helius_wide)]

helius_bl <- helius_wide %>%
    dplyr::select(ID, all_of(bl_cols)) %>%
    rename_with(~ sub("_baseline$", "_bl", .), all_of(bl_cols))

nutrient_adj_vars <- c("Protein", "FattyAcids",
                       "Carbohydrates", "Fiber", "Sodium_g")
nutrient_adj_cols <- paste0(nutrient_adj_vars, "_baseline_adj")
nutrient_adj_cols <- nutrient_adj_cols[nutrient_adj_cols %in% names(helius_wide)]
helius_bl_adj <- helius_wide %>%
    dplyr::select(ID, all_of(nutrient_adj_cols)) %>%
    rename_with(~ sub("_baseline_adj$", "_adj_bl", .), all_of(nutrient_adj_cols))

dftot2_ext <- dftot2 %>%
    left_join(helius_bl, by = "ID") %>%
    left_join(helius_bl_adj, by = "ID") %>%
    mutate(across(
        c(Age_bl, BMI_bl, TotalCalories_bl),
        ~ as.numeric(scale(.))
    )) %>%
    mutate(across(
        c(DM_bl, HT_BPMed_bl, Dyslipidemia_bl, Statins_bl, Metformin_bl, AntiHT_bl,
          PPI_bl, PsychoMed_bl, Cortico_bl,
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
    # Cardiometabolic disease
    extract_lmm_int(dftot2_ext, "shannon", "DM_bl",           "Diabetes",     group = "Disease"),
    extract_lmm_int(dftot2_ext, "shannon", "HT_BPMed_bl",     "Hypertension", group = "Disease"),
    extract_lmm_int(dftot2_ext, "shannon", "Dyslipidemia_bl", "Dyslipidemia", group = "Disease"),
    # Medications
    extract_lmm_int(dftot2_ext, "shannon", "Statins_bl",      "Statins",                group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "Metformin_bl",    "Metformin",              group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "AntiHT_bl",       "Antihypertensives",      group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "PPI_bl",          "PPI",                    group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "PsychoMed_bl",    "Psychotropics",          group = "Medication"),
    extract_lmm_int(dftot2_ext, "shannon", "Cortico_bl",      "Corticosteroids",        group = "Medication"),
    # Diet (energy-adjusted via Willett residual method, pre-computed in dietarydata.R)
    extract_lmm_int(dftot2_ext, "shannon", "TotalCalories_bl",      "Total calories (per SD)",        group = "Diet"),
    extract_lmm_int(dftot2_ext, "shannon", "Protein_adj_bl",       "Protein (per SD)",               group = "Diet"),
    extract_lmm_int(dftot2_ext, "shannon", "FattyAcids_adj_bl",    "Fatty acids (per SD)",           group = "Diet"),
    extract_lmm_int(dftot2_ext, "shannon", "Carbohydrates_adj_bl", "Carbohydrates (per SD)",         group = "Diet"),
    extract_lmm_int(dftot2_ext, "shannon", "Fiber_adj_bl",         "Fiber (per SD)",                 group = "Diet"),
    extract_lmm_int(dftot2_ext, "shannon", "Sodium_g_adj_bl",      "Sodium (per SD)",                group = "Diet")
) %>%
    mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        sig   = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR \u2265 0.05"),
        label = factor(label, levels = rev(unique(label))),
        group = factor(group, levels = c("Risk factors", "Disease", "Medication", "Diet"))
    )

write.csv(shan_ext_effects, file.path(resultsfolder, "alphadiversity_shannon_main_effects.csv"), row.names = FALSE)

(pl_shan_extended <- ggplot(shan_ext_effects, aes(x = estimate, y = label, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.2) +
    geom_point(size = 3) +
    scale_color_manual(
        values = c("FDR < 0.05" = "#F05C3BFF", "FDR \u2265 0.05" = "grey55"),
        name = NULL
    ) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Shannon change: interaction estimate (\u00b1 95% CI)",
         y = NULL) +
    theme_Publication() +
    theme(legend.position = "bottom",
          strip.text.y = element_blank(),
          strip.background.y = element_blank()))
ggsave(file.path(resultsfolder, "effectsize_shannon_change_extended.pdf"), width = 6, height = 10)

# Companion bar panel
var_meta_shan <- data.frame(
    predictor = c(
        "Age_bl", "Sex", "BMI_bl",
        "Smoking_current_bl", "Alcohol_bl", "ExerciseNorm_bl",
        "DM_bl", "HT_BPMed_bl", "Dyslipidemia_bl",
        "Statins_bl", "Metformin_bl", "AntiHT_bl", "PPI_bl",
        "PsychoMed_bl", "Cortico_bl",
        "TotalCalories_bl",
        "Protein_adj_bl", "FattyAcids_adj_bl", "Carbohydrates_adj_bl",
        "Fiber_adj_bl", "Sodium_g_adj_bl"
    ),
    label = c(
        "Age (per SD)", "Sex (female)", "BMI (per SD)",
        "Current smoking", "Alcohol use", "Sufficient exercise",
        "Diabetes", "Hypertension", "Dyslipidemia",
        "Statins", "Metformin", "Antihypertensives", "PPI",
        "Psychotropics", "Corticosteroids",
        "Total calories (per SD)",
        "Protein (per SD)", "Fatty acids (per SD)", "Carbohydrates (per SD)",
        "Fiber (per SD)", "Sodium (per SD)"
    ),
    group = c(
        rep("Risk factors", 6), rep("Disease", 3), rep("Medication", 6), rep("Diet", 6)
    ),
    type = c(
        "continuous", "continuous", "continuous",
        rep("binary", 3),
        rep("binary", 9),
        rep("continuous", 6)
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
        values = c("Yes" = "#197EC0FF", "No" = "#d0d8e4", "Non-missing" = "#197EC0FF"),
        guide = "none"
    ) +
    scale_x_continuous(limits = c(0, 140), breaks = c(0, 50, 100),
                       expand = expansion(mult = c(0, 0))) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "% Yes / % non-missing", y = NULL) +
    theme_Publication() +
    theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position = "bottom"
    )

#### Ethnicity companion panel ####

extract_lmm_eth <- function(data, outcome_var, predictor_var, label, group) {
    eth_list <- setdiff(unique(data$EthnicityTot[!is.na(data$EthnicityTot)]), "Other")
    lapply(eth_list, function(eth) {
        sub <- data %>%
            filter(EthnicityTot == eth, !is.na(.data[[predictor_var]])) %>%
            mutate(
                outcome   = .data[[outcome_var]],
                predictor = .data[[predictor_var]],
                timepoint = factor(timepoint, levels = c("baseline", "follow-up"))
            )
        if (nrow(sub) < 10) return(NULL)
        if (is.factor(sub$predictor) && sum(sub$predictor == "Yes") < 5) return(NULL)
        model <- tryCatch(
            lmer(outcome ~ predictor * timepoint + FUtime + (1|ID), data = sub),
            error = function(e) NULL
        )
        if (is.null(model)) return(NULL)
        cf <- coef(summary(model))
        ci <- confint(model, method = "Wald")
        int_row_cf <- grep(":timepoint", rownames(cf), value = TRUE)[1]
        int_row_ci <- grep(":timepoint", rownames(ci), value = TRUE)[1]
        if (is.na(int_row_cf) || is.na(int_row_ci)) return(NULL)
        data.frame(
            label     = label,
            group     = group,
            ethnicity = eth,
            estimate  = cf[int_row_cf, "Estimate"],
            conf.low  = ci[int_row_ci, 1],
            conf.high = ci[int_row_ci, 2],
            p.value   = cf[int_row_cf, "Pr(>|t|)"],
            stringsAsFactors = FALSE
        )
    }) %>% bind_rows()
}

shan_eth_effects <- purrr::map_dfr(seq_len(nrow(var_meta_shan)), function(i) {
    extract_lmm_eth(
        dftot2_ext, "shannon",
        var_meta_shan$predictor[i],
        var_meta_shan$label[i],
        var_meta_shan$group[i]
    )
}) %>%
    mutate(
        p.adj     = p.adjust(p.value, method = "BH"),
        sig       = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR \u2265 0.05"),
        label     = factor(label, levels = levels(shan_ext_effects$label)),
        group     = factor(group, levels = c("Risk factors", "Disease", "Medication", "Diet")),
        ethnicity = factor(ethnicity, levels = names(eth_colors))
    )

write.csv(shan_eth_effects, file.path(resultsfolder, "alphadiversity_shannon_ethnicity_effects.csv"), row.names = FALSE)

pl_shan_eth <- ggplot(shan_eth_effects,
                      aes(x = estimate, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(fill = ethnicity), shape = 21, color = "black", size = 2.5, stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = eth_colors, name = NULL) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Estimate", y = NULL, title = "Baseline predictors of Shannon diversity change") +
    theme_Publication() +
    theme(
        legend.position = "bottom",
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        strip.text.y = element_blank(),
        strip.background.y = element_blank()
    )

pl_shan_combined <- pl_shan_extended + pl_shan_eth + pl_shan_bar +
    plot_layout(widths = c(0.6, 0.6, 0.4))
ggsave(pl_shan_combined, filename = file.path(resultsfolder, "effectsize_shannon_change_extended_combined.pdf"),
       width = 13, height = 10)

#### Supplementary figure — Alpha-diversity × cardiometabolic disease ####

ggsave(pl_shan_combined, filename = file.path(resultsfolder, "supplementary_alphadiversity.pdf"),
       width = 13, height = 10, device = "pdf")
