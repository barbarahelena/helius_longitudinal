## Figure 2 — Cardiometabolic disease and microbiome instability (16S)
## Beta-diversity (Bray-Curtis distance) analyses
##
## Panels (used in 3_assemble_figure.R):
##   A — pl_fig2_prev: disease prevalence by ethnicity (faceted bar, baseline + follow-up)
##   B — pl_bc_extended: multi-domain predictors of Bray-Curtis instability (forest plot)
##

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(patchwork)

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
helius    <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
heliusdist <- readRDS("data/16s/braydistance_delta.RDS") %>%
    dplyr::select(-any_of(names(helius)), ID) %>%
    left_join(helius %>% filter(timepoint == "baseline"), by = "ID")

#### Output folder ####
resultsfolder <- "results/2_cmb_microbiome"
dir.create(resultsfolder, showWarnings = FALSE)

#### Ethnicity colour palette (matches Figure 1 Simpsons factor-level order) ####
eth_colors <- c(
    "Dutch"                  = "#709AE1FF",
    "South-Asian Surinamese" = "#FED439FF",
    "African Surinamese"     = "#8A9197FF",
    "Ghanaian"               = "#D2AF81FF",
    "Turkish"                = "#FD7446FF",
    "Moroccan"               = "#D5E4A2FF"
)

#### Per-ethnicity Bray-Curtis effects — data used by ethnicity companion panel ####

dat_int <- heliusdist %>%
    filter(!is.na(DM), EthnicityTot != "Other") %>%
    mutate(DM           = relevel(factor(DM), ref = "No"),
           EthnicityTot = factor(EthnicityTot))

eth_effects_dm <- lapply(levels(dat_int$EthnicityTot), function(eth) {
    sub <- dat_int %>% filter(EthnicityTot == eth)
    if (sum(sub$DM == "Yes") < 5) return(NULL)
    m  <- lm(distance ~ DM, data = sub)
    cf <- coef(summary(m))
    ci <- confint(m)
    data.frame(
        ethnicity = eth,
        estimate  = cf[2, "Estimate"],
        conf.low  = ci[2, 1],
        conf.high = ci[2, 2],
        p.value   = cf[2, "Pr(>|t|)"],
        stringsAsFactors = FALSE
    )
}) %>%
    bind_rows() %>%
    mutate(sig       = ifelse(p.value < 0.05, "p < 0.05", "p \u2265 0.05"),
           ethnicity = forcats::fct_reorder(ethnicity, estimate))

dat_int_ht <- heliusdist %>%
    filter(!is.na(HT_BPMed), EthnicityTot != "Other") %>%
    mutate(HT_BPMed     = relevel(factor(HT_BPMed), ref = "No"),
           EthnicityTot = factor(EthnicityTot))

eth_effects_ht <- lapply(levels(dat_int_ht$EthnicityTot), function(eth) {
    sub <- dat_int_ht %>% filter(EthnicityTot == eth)
    if (sum(sub$HT_BPMed == "Yes") < 5) return(NULL)
    m  <- lm(distance ~ HT_BPMed, data = sub)
    cf <- coef(summary(m))
    ci <- confint(m)
    data.frame(
        ethnicity = eth,
        estimate  = cf[2, "Estimate"],
        conf.low  = ci[2, 1],
        conf.high = ci[2, 2],
        p.value   = cf[2, "Pr(>|t|)"],
        stringsAsFactors = FALSE
    )
}) %>%
    bind_rows() %>%
    mutate(sig       = ifelse(p.value < 0.05, "p < 0.05", "p \u2265 0.05"),
           ethnicity = forcats::fct_reorder(ethnicity, estimate))

dat_int_lld <- heliusdist %>%
    filter(!is.na(Dyslipidemia), EthnicityTot != "Other") %>%
    mutate(Dyslipidemia = relevel(factor(Dyslipidemia), ref = "No"),
           EthnicityTot = factor(EthnicityTot))

eth_effects_lld <- lapply(levels(dat_int_lld$EthnicityTot), function(eth) {
    sub <- dat_int_lld %>% filter(EthnicityTot == eth)
    if (sum(sub$Dyslipidemia == "Yes") < 5) return(NULL)
    m  <- lm(distance ~ Dyslipidemia, data = sub)
    cf <- coef(summary(m))
    ci <- confint(m)
    data.frame(
        ethnicity = eth,
        estimate  = cf[2, "Estimate"],
        conf.low  = ci[2, 1],
        conf.high = ci[2, 2],
        p.value   = cf[2, "Pr(>|t|)"],
        stringsAsFactors = FALSE
    )
}) %>%
    bind_rows() %>%
    mutate(sig       = ifelse(p.value < 0.05, "p < 0.05", "p \u2265 0.05"),
           ethnicity = forcats::fct_reorder(ethnicity, estimate))

#### Panel A — Disease prevalence by ethnicity: DM, HT, Dyslipidemia, MetSyn ####

paired_ids <- heliusdist$ID

prev_long <- bind_rows(
    helius %>%
        filter(ID %in% paired_ids, !is.na(DM), EthnicityTot != "Other") %>%
        group_by(EthnicityTot, timepoint) %>%
        summarise(prev = mean(DM == "Yes") * 100, .groups = "drop") %>%
        mutate(condition = "Diabetes"),
    helius %>%
        filter(ID %in% paired_ids, !is.na(HT_BPMed), EthnicityTot != "Other") %>%
        group_by(EthnicityTot, timepoint) %>%
        summarise(prev = mean(HT_BPMed == "Yes") * 100, .groups = "drop") %>%
        mutate(condition = "Hypertension"),
    helius %>%
        filter(ID %in% paired_ids, !is.na(Dyslipidemia), EthnicityTot != "Other") %>%
        group_by(EthnicityTot, timepoint) %>%
        summarise(prev = mean(Dyslipidemia == "Yes") * 100, .groups = "drop") %>%
        mutate(condition = "Dyslipidemia"),
    helius %>%
        filter(ID %in% paired_ids, !is.na(MetSyn), EthnicityTot != "Other") %>%
        group_by(EthnicityTot, timepoint) %>%
        summarise(prev = mean(MetSyn == "Yes") * 100, .groups = "drop") %>%
        mutate(condition = "Metabolic syndrome")
) %>%
    mutate(
        condition = factor(condition, levels = c("Diabetes", "Hypertension", "Dyslipidemia", "Metabolic syndrome")),
        timepoint = factor(timepoint, levels = c("follow-up", "baseline"))
    )

# Prevalence table with relative increase
prev_table <- prev_long %>%
    pivot_wider(names_from = timepoint, values_from = prev) %>%
    mutate(relative_increase_pct = ((`follow-up` - baseline) / baseline) * 100) %>%
    dplyr::select(condition, EthnicityTot, baseline, `follow-up`, relative_increase_pct) %>%
    rename(
        Condition        = condition,
        Ethnicity        = EthnicityTot,
        Baseline_prev    = baseline,
        Followup_prev    = `follow-up`,
        Relative_increase_pct = relative_increase_pct
    ) %>%
    arrange(Condition, Ethnicity)

write.csv(prev_table, file.path(resultsfolder, "prevalence_table.csv"), row.names = FALSE)

# Order ethnicities by baseline DM prevalence
eth_order <- prev_long %>%
    filter(condition == "Diabetes", timepoint == "baseline") %>%
    arrange(prev) %>%
    pull(EthnicityTot) %>%
    as.character()
prev_long <- prev_long %>% mutate(EthnicityTot = factor(EthnicityTot, levels = eth_order))

(pl_fig2_prev <- ggplot(prev_long %>%
                            mutate(EthnicityTot = fct_recode(EthnicityTot,
                                "South-Asian\nSurinamese" = "South-Asian Surinamese")),
                        aes(x = prev, y = EthnicityTot, fill = EthnicityTot, alpha = timepoint)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65,
             color = "black", linewidth = 0.3) +
    geom_text(
        data = prev_long %>% filter(condition == "Diabetes",
                                    EthnicityTot == head(eth_order, 1)),
        aes(x = prev + 0.5, y = EthnicityTot,
            label = ifelse(timepoint == "baseline", "Baseline", "Follow-up"),
            group = timepoint),
        position = position_dodge(width = 0.75),
        hjust = 0, size = 3, color = "grey30", inherit.aes = FALSE
    ) +
    scale_alpha_manual(values = c("baseline" = 0.35, "follow-up" = 0.85), guide = "none") +
    scale_fill_manual(values = setNames(eth_colors,
                                        gsub("South-Asian Surinamese", "South-Asian\nSurinamese",
                                             names(eth_colors))),
                      guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    facet_wrap(~condition, scales = "free_x", nrow = 1) +
    labs(x = "Prevalence (%)", y = NULL,
         title = "Disease prevalence") +
    theme_Publication() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = rel(1.0))))
ggsave(file.path(resultsfolder, "cmb_prevalence_combined.pdf"), width = 15, height = 5)

#### Panel B — Extended forest plot: multi-domain predictors of Bray-Curtis instability ####

# Helper: extract lm() effect (predictor in row 2 of coef table)
extract_lm_effect_ext <- function(data, predictor_var, label, group) {
    df <- data %>%
        filter(!is.na(.data[[predictor_var]])) %>%
        mutate(predictor = .data[[predictor_var]])
    model <- lm(distance ~ predictor + FUtime, data = df)
    cf <- coef(summary(model))
    ci <- confint(model)
    data.frame(
        label     = label,
        group     = group,
        estimate  = cf[2, "Estimate"],
        conf.low  = ci[2, 1],
        conf.high = ci[2, 2],
        p.value   = cf[2, "Pr(>|t|)"],
        stringsAsFactors = FALSE
    )
}

# Build analysis dataset from wide clinical data + Bray-Curtis distance
helius_wide <- readRDS("data/clinicaldata/clinicaldata_wide.RDS")

pcdiet <- readRDS("data/clinicaldata_long_pcdiet.RDS") %>%
    filter(timepoint == "baseline") %>%
    dplyr::select(ID, DietPC1, DietPC2)

heliusdist_ext <- helius_wide %>%
    left_join(heliusdist %>% dplyr::select(ID, distance), by = "ID") %>%
    left_join(pcdiet, by = "ID") %>%
    filter(!is.na(distance)) %>%
    mutate(across(
        c(Age_baseline, BMI_baseline, DiscrMean_baseline,
          Protein_baseline, FattyAcids_baseline,
          Carbohydrates_baseline, Fiber_baseline, Sodium_g_baseline,
          DietPC1, DietPC2),
        ~ as.numeric(scale(.))
    )) %>%
    mutate(across(
        c(DM_baseline, HT_BPMed_baseline, Dyslipidemia_baseline, MetSyn_baseline,
          Statins_baseline, Metformin_baseline, AntiHT_baseline, PPI_baseline,
          PsychoMed_baseline, Cortico_baseline,
          Smoking_current_baseline, Alcohol_baseline, ExerciseNorm_baseline),
        ~ relevel(factor(.), ref = "No"),
        .names = "{.col}"
    )) %>%
    mutate(Sex = relevel(factor(Sex), ref = "Male"))

bc_ext_effects <- bind_rows(
    # Risk factors
    extract_lm_effect_ext(heliusdist_ext, "Age_baseline",             "Age (per SD)",                      "Risk factors"),
    extract_lm_effect_ext(heliusdist_ext, "Sex",                      "Sex (female)",                      "Risk factors"),
    extract_lm_effect_ext(heliusdist_ext, "BMI_baseline",             "BMI (per SD)",                      "Risk factors"),
    extract_lm_effect_ext(heliusdist_ext, "Smoking_current_baseline", "Current smoking",                   "Risk factors"),
    extract_lm_effect_ext(heliusdist_ext, "Alcohol_baseline",         "Alcohol use",                       "Risk factors"),
    extract_lm_effect_ext(heliusdist_ext, "ExerciseNorm_baseline",    "Sufficient exercise",               "Risk factors"),
    extract_lm_effect_ext(heliusdist_ext, "DiscrMean_baseline",       "Perceived discrimination (per SD)", "Risk factors"),
    # Cardiometabolic disease
    extract_lm_effect_ext(heliusdist_ext, "DM_baseline",           "Diabetes",           "Disease"),
    extract_lm_effect_ext(heliusdist_ext, "HT_BPMed_baseline",     "Hypertension",       "Disease"),
    extract_lm_effect_ext(heliusdist_ext, "Dyslipidemia_baseline", "Dyslipidemia",       "Disease"),
    extract_lm_effect_ext(heliusdist_ext, "MetSyn_baseline",       "Metabolic syndrome", "Disease"),
    # Medications
    extract_lm_effect_ext(heliusdist_ext, "Statins_baseline",      "Statins",                "Medication"),
    extract_lm_effect_ext(heliusdist_ext, "Metformin_baseline",    "Metformin",              "Medication"),
    extract_lm_effect_ext(heliusdist_ext, "AntiHT_baseline",       "Antihypertensives",      "Medication"),
    extract_lm_effect_ext(heliusdist_ext, "PPI_baseline",          "PPI",                    "Medication"),
    extract_lm_effect_ext(heliusdist_ext, "PsychoMed_baseline",    "Psychotropics",          "Medication"),
    extract_lm_effect_ext(heliusdist_ext, "Cortico_baseline",      "Corticosteroids",        "Medication"),
    # Diet (energy-adjusted via Willett residual method, pre-computed in dietarydata.R)
    extract_lm_effect_ext(heliusdist_ext, "Protein_baseline",       "Protein (per SD)",       "Diet"),
    extract_lm_effect_ext(heliusdist_ext, "FattyAcids_baseline",    "Fatty acids (per SD)",   "Diet"),
    extract_lm_effect_ext(heliusdist_ext, "Carbohydrates_baseline", "Carbohydrates (per SD)", "Diet"),
    extract_lm_effect_ext(heliusdist_ext, "Fiber_baseline",         "Fiber (per SD)",         "Diet"),
    extract_lm_effect_ext(heliusdist_ext, "Sodium_g_baseline",      "Sodium (per SD)",        "Diet"),
    extract_lm_effect_ext(heliusdist_ext, "DietPC1",                "Diet PC1 (per SD)",      "Diet"),
    extract_lm_effect_ext(heliusdist_ext, "DietPC2",                "Diet PC2 (per SD)",      "Diet")
) %>%
    mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        sig   = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR \u2265 0.05"),
        label = factor(label, levels = rev(unique(label))),
        group = factor(group, levels = c("Risk factors", "Disease", "Medication", "Diet"))
    )

(pl_bc_extended <- ggplot(bc_ext_effects, aes(x = estimate, y = label, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.2) +
    geom_point(size = 3) +
    scale_color_manual(
        values = c("FDR < 0.05" = "#F05C3BFF", "FDR \u2265 0.05" = "grey55"),
        name = NULL
    ) +
    facet_wrap(~ group, scales = "free_y", ncol = 1) +
    labs(x = "Bray-Curtis instability: estimate (\u00b1 95% CI)",
         y = NULL) +
    guides(color = guide_legend(ncol = 1)) +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_braycurtis_extended.pdf"), width = 6, height = 10)

# Companion bar panel: prevalence (binary) and n non-missing (continuous)
var_meta_bc <- data.frame(
    predictor = c(
        "Age_baseline", "Sex", "BMI_baseline",
        "Smoking_current_baseline", "Alcohol_baseline", "ExerciseNorm_baseline", "DiscrMean_baseline",
        "DM_baseline", "HT_BPMed_baseline", "Dyslipidemia_baseline", "MetSyn_baseline",
        "Statins_baseline", "Metformin_baseline", "AntiHT_baseline", "PPI_baseline",
        "PsychoMed_baseline", "Cortico_baseline",
        "Protein_baseline", "FattyAcids_baseline", "Carbohydrates_baseline",
        "Fiber_baseline", "Sodium_g_baseline", "DietPC1", "DietPC2"
    ),
    label = c(
        "Age (per SD)", "Sex (female)", "BMI (per SD)",
        "Current smoking", "Alcohol use", "Sufficient exercise", "Perceived discrimination (per SD)",
        "Diabetes", "Hypertension", "Dyslipidemia", "Metabolic syndrome",
        "Statins", "Metformin", "Antihypertensives", "PPI",
        "Psychotropics", "Corticosteroids",
        "Protein (per SD)", "Fatty acids (per SD)", "Carbohydrates (per SD)",
        "Fiber (per SD)", "Sodium (per SD)", "Diet PC1 (per SD)", "Diet PC2 (per SD)"
    ),
    group = c(
        rep("Risk factors", 7), rep("Disease", 4), rep("Medication", 6), rep("Diet", 7)
    ),
    type = c(
        "continuous", "binary", "continuous",
        rep("binary", 3), "continuous",
        rep("binary", 10),
        rep("continuous", 7)
    ),
    stringsAsFactors = FALSE
)

n_total_bc <- nrow(heliusdist_ext)

bar_data_bc <- purrr::map_dfr(seq_len(nrow(var_meta_bc)), function(i) {
    col <- var_meta_bc$predictor[i]
    x <- heliusdist_ext[[col]]
    if (var_meta_bc$type[i] == "binary") {
        lvls  <- if (is.factor(x)) levels(x) else c("No", "Yes")
        n_no  <- sum(x == lvls[1], na.rm = TRUE)
        n_yes <- sum(x == lvls[2], na.rm = TRUE)
        n_tot <- n_yes + n_no
        data.frame(
            label    = var_meta_bc$label[i],
            group    = var_meta_bc$group[i],
            category = c("Yes", "No"),
            pct      = c(n_yes / n_tot * 100, n_no / n_tot * 100),
            n_label  = c(n_yes, NA)
        )
    } else {
        n_valid   <- sum(!is.na(x))
        pct_valid <- n_valid / n_total_bc * 100
        data.frame(
            label    = var_meta_bc$label[i],
            group    = var_meta_bc$group[i],
            category = "Non-missing",
            pct      = pct_valid,
            n_label  = n_valid
        )
    }
}) %>%
    mutate(
        label    = factor(label, levels = levels(bc_ext_effects$label)),
        group    = factor(group, levels = c("Risk factors", "Disease", "Medication", "Diet")),
        category = factor(category, levels = c("No", "Yes", "Non-missing"))
    )

pl_bc_bar <- ggplot(bar_data_bc, aes(x = pct, y = label, fill = category)) +
    geom_col(width = 0.6) +
    geom_text(
        data = bar_data_bc %>% filter(!is.na(n_label)),
        aes(label = paste0("n=", n_label), x = 101),
        hjust = 0, size = 2.8, color = "grey30"
    ) +
    scale_fill_manual(
        values = c("Yes" = "#197EC0FF", "No" = "#d0d8e4", "Non-missing" = "#197EC0FF"),
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

#### Ethnicity companion panel ####

extract_lm_eth <- function(data, predictor_var, label, group) {
    eth_list <- setdiff(unique(data$EthnicityTot[!is.na(data$EthnicityTot)]), "Other")
    lapply(eth_list, function(eth) {
        sub <- data %>%
            filter(EthnicityTot == eth, !is.na(.data[[predictor_var]])) %>%
            mutate(predictor = .data[[predictor_var]])
        if (nrow(sub) < 10) return(NULL)
        if (is.factor(sub$predictor) && "Yes" %in% levels(sub$predictor) && sum(sub$predictor == "Yes") < 5) return(NULL)
        m  <- lm(distance ~ predictor, data = sub)
        cf <- coef(summary(m))
        ci <- confint(m)
        data.frame(
            label     = label,
            group     = group,
            ethnicity = eth,
            estimate  = cf[2, "Estimate"],
            conf.low  = ci[2, 1],
            conf.high = ci[2, 2],
            p.value   = cf[2, "Pr(>|t|)"],
            stringsAsFactors = FALSE
        )
    }) %>% bind_rows()
}

bc_eth_effects <- purrr::map_dfr(seq_len(nrow(var_meta_bc)), function(i) {
    extract_lm_eth(
        heliusdist_ext,
        var_meta_bc$predictor[i],
        var_meta_bc$label[i],
        var_meta_bc$group[i]
    )
}) %>%
    mutate(
        p.adj     = p.adjust(p.value, method = "BH"),
        sig       = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR \u2265 0.05"),
        label     = factor(label, levels = levels(bc_ext_effects$label)),
        group     = factor(group, levels = c("Risk factors", "Disease", "Medication", "Diet")),
        ethnicity = factor(ethnicity, levels = names(eth_colors))
    )

# Per-variable interaction p-values (predictor × EthnicityTot F-test)
bc_int_pvals <- purrr::map_dfr(seq_len(nrow(var_meta_bc)), function(i) {
    col <- var_meta_bc$predictor[i]
    sub <- heliusdist_ext %>%
        filter(!is.na(.data[[col]]), EthnicityTot != "Other", !is.na(EthnicityTot)) %>%
        mutate(predictor    = .data[[col]],
               EthnicityTot = factor(EthnicityTot))
    if (nrow(sub) < 20) return(data.frame(label = var_meta_bc$label[i], p_int = NA_real_))
    m    <- lm(distance ~ predictor * EthnicityTot, data = sub)
    pval <- tryCatch(
        drop1(m, scope = ~predictor:EthnicityTot, test = "F")["predictor:EthnicityTot", "Pr(>F)"],
        error = function(e) NA_real_
    )
    data.frame(label = var_meta_bc$label[i], p_int = pval, stringsAsFactors = FALSE)
})

bc_eth_heatmap <- bc_eth_effects %>%
    left_join(bc_int_pvals, by = "label") %>%
    mutate(
        sig_label = case_when(p.adj < 0.001 ~ "***",
                              p.adj < 0.01  ~ "**",
                              p.adj < 0.05  ~ "*",
                              TRUE          ~ ""),
        label     = factor(label, levels = levels(bc_ext_effects$label)),
        ethnicity = factor(ethnicity, levels = names(eth_colors))
    )

est_lim <- quantile(abs(bc_eth_heatmap$estimate), 0.95, na.rm = TRUE)

(pl_bc_eth <- ggplot(bc_eth_heatmap,
                     aes(x = estimate, y = label, color = ethnicity)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(fill = ethnicity), shape = 21, color = "black", size = 2.5, stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = eth_colors, name = NULL) +
    facet_wrap(~ group, scales = "free_y", ncol = 1) +
    labs(x = "Estimate per ethnic group", y = NULL,
         title = "Baseline predictors of microbiota instability") +
    theme_Publication() +
    theme(
        legend.position = "bottom",
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank()
    ))

pl_bc_combined <- pl_bc_extended + pl_bc_eth + pl_bc_bar +
    plot_layout(widths = c(0.6, 0.6, 0.4))
ggsave(pl_bc_combined, filename = file.path(resultsfolder, "effectsize_braycurtis_extended_combined.pdf"),
       width = 12, height = 10)
