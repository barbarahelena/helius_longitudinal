## Sensitivity analysis — exclude baseline -> follow-up "starters" per predictor
##
## For each baseline disease/medication predictor of Bray-Curtis change that has
## a matched follow-up status, this re-fits the main-analysis model
##   distance ~ predictor + Age_baseline + Sex + BMI_baseline +
##              Metformin_baseline + PPI_baseline + FUtime
## (dropping Metformin/PPI from the covariate set when either is itself the
## predictor) after excluding participants who were "No" at baseline and
## "Yes" at follow-up for that specific variable ("starters"). Exclusion is
## done separately per variable — a starter for one variable stays in the
## analysis for all other variables. Participants who discontinued
## (Yes at baseline -> No at follow-up) are kept in the baseline "Yes" group,
## unchanged from the main analysis.
##
## Both a full-cohort model and ethnicity-stratified models (same formula,
## fit within each ethnic group) are re-run.
##
## Note: follow-up statin use was not collected in HELIUS (H1_Statines has no
## H2 equivalent, see scripts/0_run_workflows/3_data_cleaning/datacleaning.R).
## LLD (lipid-lowering drugs, H2_Antilipaemica) follow-up status is used as a
## proxy for statin starters; this is a broader drug class than statins alone.

## Libraries
library(tidyverse)

#### Output folder ####
resultsfolder <- "results/2_cmb_microbiome"
dir.create(resultsfolder, showWarnings = FALSE)

#### Load data (same construction as 1_betadiversity_cmb_16s.R) ####
heliusdist  <- readRDS("data/16s/braydistance_delta.RDS")
helius_wide <- readRDS("data/clinicaldata/clinicaldata_wide.RDS")

heliusdist_ext <- helius_wide %>%
    left_join(heliusdist %>% dplyr::select(ID, distance), by = "ID") %>%
    filter(!is.na(distance)) %>%
    mutate(across(
        c(DM_baseline, `DM_follow-up`,
          HT_BPMed_baseline, `HT_BPMed_follow-up`,
          Dyslipidemia_baseline, `Dyslipidemia_follow-up`,
          MetSyn_baseline, `MetSyn_follow-up`,
          Statins_baseline, LLD_baseline, `LLD_follow-up`,
          Metformin_baseline, `Metformin_follow-up`,
          AntiHT_baseline, `AntiHT_follow-up`,
          PPI_baseline, `PPI_follow-up`,
          Cortico_baseline, `Cortico_follow-up`,
          Smoking_current_baseline, `Smoking_current_follow-up`,
          Alcohol_baseline, `Alcohol_follow-up`),
        ~ relevel(factor(.), ref = "No"),
        .names = "{.col}"
    )) %>%
    mutate(Sex = relevel(factor(Sex), ref = "Male"))

#### Variables with matched baseline/follow-up status ####
# followup_col holds the transition variable used to define "starters"
# (No at baseline, Yes at followup_col); for Statins this is the LLD proxy.
var_list <- tribble(
    ~variable,               ~followup_col,              ~label,                ~group,        ~proxy_note,
    "DM_baseline",            "DM_follow-up",             "Diabetes",            "Disease",     NA_character_,
    "HT_BPMed_baseline",      "HT_BPMed_follow-up",       "Hypertension",        "Disease",     NA_character_,
    "Dyslipidemia_baseline",  "Dyslipidemia_follow-up",   "Dyslipidaemia",       "Disease",     NA_character_,
    "MetSyn_baseline",        "MetSyn_follow-up",         "Metabolic syndrome",  "Disease",     NA_character_,
    "Metformin_baseline",     "Metformin_follow-up",      "Metformin",           "Medication",  NA_character_,
    "AntiHT_baseline",        "AntiHT_follow-up",         "Antihypertensives",   "Medication",  NA_character_,
    "PPI_baseline",           "PPI_follow-up",            "PPI",                 "Medication",  NA_character_,
    "Cortico_baseline",       "Cortico_follow-up",        "Corticosteroids",     "Medication",  NA_character_,
    "Statins_baseline",       "LLD_follow-up",            "Statins",             "Medication",  "Follow-up statin use not collected; LLD (lipid-lowering drugs) follow-up used as proxy for starter status",
    "Smoking_current_baseline", "Smoking_current_follow-up", "Current smoking",  "Risk factor", NA_character_,
    "Alcohol_baseline",       "Alcohol_follow-up",        "Alcohol use",         "Risk factor", NA_character_
)

base_covars <- c("Age_baseline", "Sex", "BMI_baseline", "Metformin_baseline", "PPI_baseline", "FUtime")

#### Helper: fit distance ~ predictor + covariates, extract the Yes-vs-No effect ####
fit_predictor_effect <- function(df, covars) {
    if (nlevels(droplevels(df$predictor)) < 2 || sum(df$predictor == "Yes") < 5) return(NULL)
    f <- reformulate(c("predictor", covars), response = "distance")
    m <- lm(f, data = droplevels(df))
    cf <- coef(summary(m))
    if (!"predictorYes" %in% rownames(cf)) return(NULL)
    ci <- confint(m)
    tibble(
        n         = stats::nobs(m),
        estimate  = cf["predictorYes", "Estimate"],
        conf.low  = ci["predictorYes", 1],
        conf.high = ci["predictorYes", 2],
        p.value   = cf["predictorYes", "Pr(>|t|)"]
    )
}

#### Build the complete-case analysis frame for one variable, and its starter exclusion ####
build_variable_frame <- function(data, variable, followup_col) {
    covars <- setdiff(base_covars, variable)
    data %>%
        dplyr::select(ID, EthnicityTot, distance,
                      predictor = all_of(variable),
                      followup  = all_of(followup_col),
                      all_of(covars)) %>%
        filter(!is.na(predictor), !is.na(distance), if_all(all_of(covars), ~ !is.na(.))) %>%
        mutate(is_starter = predictor == "No" & !is.na(followup) & followup == "Yes")
}

#### Full-cohort main + sensitivity models, per variable ####
fullcohort_results <- pmap_dfr(var_list, function(variable, followup_col, label, group, proxy_note) {
    covars  <- setdiff(base_covars, variable)
    df_full <- build_variable_frame(heliusdist_ext, variable, followup_col)
    df_sens <- df_full %>% filter(!is_starter)

    main_fit <- fit_predictor_effect(df_full, covars)
    sens_fit <- fit_predictor_effect(df_sens, covars)

    tibble(
        variable         = variable,
        label            = label,
        group            = group,
        proxy_note       = proxy_note,
        n_excluded       = sum(df_full$is_starter),
        n_no_at_both     = sum(df_sens$predictor == "No"),
        n_yes_baseline   = sum(df_sens$predictor == "Yes"),
        main_n           = main_fit$n[1],
        main_estimate    = main_fit$estimate[1],
        main_conf.low    = main_fit$conf.low[1],
        main_conf.high   = main_fit$conf.high[1],
        main_p.value     = main_fit$p.value[1],
        sens_n           = sens_fit$n[1],
        sens_estimate    = sens_fit$estimate[1],
        sens_conf.low    = sens_fit$conf.low[1],
        sens_conf.high   = sens_fit$conf.high[1],
        sens_p.value     = sens_fit$p.value[1]
    )
}) %>%
    mutate(
        main_q.value = p.adjust(main_p.value, method = "BH"),
        sens_q.value = p.adjust(sens_p.value, method = "BH")
    ) %>%
    dplyr::select(variable, label, group, proxy_note, n_excluded, n_no_at_both, n_yes_baseline,
                  main_n, main_estimate, main_conf.low, main_conf.high, main_p.value, main_q.value,
                  sens_n, sens_estimate, sens_conf.low, sens_conf.high, sens_p.value, sens_q.value)

print(fullcohort_results, width = Inf)
write.csv(fullcohort_results,
          file.path(resultsfolder, "sensitivity_starters_fullcohort.csv"),
          row.names = FALSE)

#### Ethnicity-stratified main + sensitivity models, per variable ####
eth_levels <- setdiff(unique(heliusdist_ext$EthnicityTot[!is.na(heliusdist_ext$EthnicityTot)]), "Other")

ethnicity_results <- pmap_dfr(var_list, function(variable, followup_col, label, group, proxy_note) {
    covars  <- setdiff(base_covars, variable)
    df_full <- build_variable_frame(heliusdist_ext, variable, followup_col)
    df_sens <- df_full %>% filter(!is_starter)

    map_dfr(eth_levels, function(eth) {
        sub_full <- df_full %>% filter(EthnicityTot == eth)
        sub_sens <- df_sens %>% filter(EthnicityTot == eth)

        main_fit <- fit_predictor_effect(sub_full, covars)
        sens_fit <- fit_predictor_effect(sub_sens, covars)
        if (is.null(main_fit) && is.null(sens_fit)) return(NULL)

        tibble(
            variable        = variable,
            label           = label,
            group           = group,
            proxy_note      = proxy_note,
            ethnicity       = eth,
            n_excluded      = sum(sub_full$is_starter),
            n_no_at_both    = sum(sub_sens$predictor == "No"),
            n_yes_baseline  = sum(sub_sens$predictor == "Yes"),
            main_n          = if (!is.null(main_fit)) main_fit$n[1] else NA_integer_,
            main_estimate   = if (!is.null(main_fit)) main_fit$estimate[1] else NA_real_,
            main_conf.low   = if (!is.null(main_fit)) main_fit$conf.low[1] else NA_real_,
            main_conf.high  = if (!is.null(main_fit)) main_fit$conf.high[1] else NA_real_,
            main_p.value    = if (!is.null(main_fit)) main_fit$p.value[1] else NA_real_,
            sens_n          = if (!is.null(sens_fit)) sens_fit$n[1] else NA_integer_,
            sens_estimate   = if (!is.null(sens_fit)) sens_fit$estimate[1] else NA_real_,
            sens_conf.low   = if (!is.null(sens_fit)) sens_fit$conf.low[1] else NA_real_,
            sens_conf.high  = if (!is.null(sens_fit)) sens_fit$conf.high[1] else NA_real_,
            sens_p.value    = if (!is.null(sens_fit)) sens_fit$p.value[1] else NA_real_
        )
    })
}) %>%
    mutate(
        main_q.value = p.adjust(main_p.value, method = "BH"),
        sens_q.value = p.adjust(sens_p.value, method = "BH")
    ) %>%
    dplyr::select(variable, label, group, proxy_note, ethnicity, n_excluded, n_no_at_both, n_yes_baseline,
                  main_n, main_estimate, main_conf.low, main_conf.high, main_p.value, main_q.value,
                  sens_n, sens_estimate, sens_conf.low, sens_conf.high, sens_p.value, sens_q.value)

print(ethnicity_results, width = Inf)
write.csv(ethnicity_results,
          file.path(resultsfolder, "sensitivity_starters_ethnicity.csv"),
          row.names = FALSE)
