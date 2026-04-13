## Archive — Alpha-diversity exploratory plots (no longer in main/supplement figures)
##
## Sources the main alpha-diversity script to reuse data objects.

source("scripts/2_cmb_microbiome/2_alphadiversity_cmb_16s.R")

#### linearmixed helper (used below) ####

linearmixed <- function(data, var, groupvar){
    data <- data %>% mutate(var = {{ var }},
                            groupvar = {{ groupvar }},
                            timepoint = as.factor(timepoint)) %>%
                        filter(!is.na(MetSyn_new))

    model1 <- lmer(var ~ groupvar*timepoint + (1|ID), data = data)
    res <- summary(model1)
    print(res)
    pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres <- cbind(group1 = "baseline", group2 = "follow-up", pval)
    statres <- tibble::as_tibble(statres)
    statres$p.signif <- case_when(
        statres$pval < 0.05 ~paste0("*"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.001 ~paste0("***"),
        statres$pval > 0.05 ~paste0("")
    )
    statres <- statres %>% filter(p.signif != "")
    return(statres)
}

#### Propagate new-diagnosis status to follow-up rows ####

dftot2$DM_new[which(dftot2$timepoint == "follow-up")]     <- dftot2$DM_new[which(dftot2$ID == dftot2$ID & dftot2$timepoint == "baseline")]
dftot2$HT_new[which(dftot2$timepoint == "follow-up")]     <- dftot2$HT_new[which(dftot2$ID == dftot2$ID & dftot2$timepoint == "baseline")]
dftot2$MetSyn_new[which(dftot2$timepoint == "follow-up")] <- dftot2$MetSyn_new[which(dftot2$ID == dftot2$ID & dftot2$timepoint == "baseline")]

#### Delta Shannon by diabetes (shannon_dm.pdf) ####

ggplot(data = dftot %>% filter(!is.na(DM)), aes(x = DM, y = shannon_delta)) +
    geom_violin(colour = NA, aes(fill = DM)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Delta Shannon index", x = "Diabetes", title = "Diabetes") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0,
                       hide.ns = FALSE, label = "p.format") +
    theme_Publication()
ggsave(file.path(resultsfolder, "shannon_dm.pdf"), width = 4, height = 5)

shannon_dm <- linearmixed(dftot2, shannon, DM_new)

(pl1 <- dftot2 %>% filter(!is.na(DM_new)) %>% ggplot() +
        geom_line(aes(x = timepoint, y = shannon, group = ID), alpha = 0.1, color = "grey40") +
        geom_jitter(aes(x = timepoint, y = shannon, group = ID, color = DM_new),
                    alpha = 0.5, width = 0) +
        gghalves::geom_half_violin(aes(x = timepoint, y = shannon, fill = DM_new),
                                   side = c("l", "l", "r", "r"), nudge = 0.05) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "baseline"),
                                    aes(x = timepoint, y = shannon),
                                    nudge = 0.05, side = "l", fill = "white", color = "black",
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "follow-up"),
                                    aes(x = timepoint, y = shannon),
                                    nudge = 0.05, side = "r", fill = "white", color = "black",
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        stat_pvalue_manual(shannon_dm, label = "{p.signif}", y.position = max(dftot2$shannon), remove.bracket = FALSE) +
        scale_color_simpsons(guide = "none") +
        scale_fill_simpsons(guide = "none") +
        facet_wrap(~DM_new) +
        theme_Publication() +
        labs(x = "Timepoint", y = "Shannon index", title = "New diabetes diagnoses", color = ""))
ggsave(file.path(resultsfolder, "diabetes_deltashannon.pdf"), width = 4, height = 4)

#### Delta Shannon and cardiometabolic outcomes ####

ggplot(data = dftot %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = shannon_delta)) +
    geom_violin(colour = NA, aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Delta Shannon index", x = "Hypertension", title = "Hypertension") +
    theme_Publication()
ggsave(file.path(resultsfolder, "shannon_hypertension.pdf"), width = 4, height = 5)

ggplot(data = dftot %>% filter(!is.na(MetSyn)), aes(x = MetSyn, y = shannon_delta)) +
    geom_violin(colour = NA, aes(fill = MetSyn)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Delta Shannon index", x = "Metabolic syndrome", title = "Metabolic syndrome") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication()
ggsave(file.path(resultsfolder, "shannon_metsyn.pdf"), width = 4, height = 5)

ggplot(data = dftot %>% filter(!is.na(HbA1c_delta)), aes(x = shannon_delta, y = HbA1c_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta HbA1c", x = "Delta Shannon index", title = "Shannon and HbA1c change") +
    stat_cor() +
    theme_Publication()
ggsave(file.path(resultsfolder, "shannon_deltahba1c.pdf"), width = 7, height = 7)

ggplot(data = dftot %>% filter(!is.na(LDL_delta)), aes(x = shannon_delta, y = LDL_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta LDL", x = "Delta Shannon index", title = "Shannon and LDL change") +
    stat_cor() +
    theme_Publication()
ggsave(file.path(resultsfolder, "shannon_deltaldl.pdf"), width = 7, height = 7)

ggplot(data = dftot %>% filter(!is.na(BMI_delta)), aes(x = shannon_delta, y = BMI_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta BMI", x = "Delta Shannon index", title = "Shannon and BMI change") +
    stat_cor() +
    theme_Publication()
ggsave(file.path(resultsfolder, "shannon_deltabmi.pdf"), width = 7, height = 7)

#### LMM: new diagnoses and alpha diversity trajectories ####

df_means_ht <- dftot2 %>% filter(!is.na(HT_new)) %>%
    dplyr::select(ID, HT_new, timepoint, shannon, richness, PD) %>%
    group_by(timepoint, HT_new) %>%
    summarise(across(c(shannon, richness, PD),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)),
                     .names = "{.col}_{.fn}"),
              .groups = "drop_last")

df_means_ms <- dftot2 %>% filter(!is.na(MetSyn_new)) %>%
    dplyr::select(ID, MetSyn_new, timepoint, shannon, richness, PD) %>%
    group_by(timepoint, MetSyn_new) %>%
    summarise(across(c(shannon, richness, PD),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)),
                     .names = "{.col}_{.fn}"),
              .groups = "drop_last")

richness_dm <- linearmixed(dftot2, richness, DM_new)
pd_dm       <- linearmixed(dftot2, PD, DM_new)
shannon_ht  <- linearmixed(dftot2, shannon, HT_new)
richness_ht <- linearmixed(dftot2, richness, HT_new)
pd_ht       <- linearmixed(dftot2, PD, HT_new)
shannon_ms  <- linearmixed(dftot2, shannon, MetSyn_new)
richness_ms <- linearmixed(dftot2, richness, MetSyn_new)
pd_ms       <- linearmixed(dftot2, PD, MetSyn_new)

(pl1 <- dftot2 %>% filter(!is.na(DM_new)) %>% ggplot() +
        geom_line(aes(x = timepoint, y = richness, group = ID), alpha = 0.1, color = "grey40") +
        geom_jitter(aes(x = timepoint, y = richness, group = ID, color = DM_new),
                    alpha = 0.5, width = 0) +
        gghalves::geom_half_violin(aes(x = timepoint, y = richness, fill = DM_new),
                                   side = c("l", "l", "r", "r"), nudge = 0.05) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "baseline"),
                                    aes(x = timepoint, y = richness),
                                    nudge = 0.05, side = "l", fill = "white", color = "black",
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "follow-up"),
                                    aes(x = timepoint, y = richness),
                                    nudge = 0.05, side = "r", fill = "white", color = "black",
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        stat_pvalue_manual(richness_dm, label = "{p.signif}", y.position = max(dftot2$richness), remove.bracket = FALSE) +
        scale_color_simpsons(guide = "none") +
        scale_fill_simpsons(guide = "none") +
        facet_wrap(~DM_new) +
        theme_Publication() +
        labs(x = "Timepoint", y = "Richness", title = "New diabetes diagnoses", color = ""))
ggsave(file.path(resultsfolder, "diabetes_deltarichness.pdf"), width = 4, height = 4)

(pl2 <- dftot2 %>% filter(!is.na(HT_new)) %>% ggplot() +
        geom_line(aes(x = timepoint, y = shannon, group = ID), alpha = 0.1, color = "grey40") +
        geom_jitter(aes(x = timepoint, y = shannon, group = ID, color = HT_new),
                    alpha = 0.5, width = 0) +
        gghalves::geom_half_violin(aes(x = timepoint, y = shannon, fill = HT_new),
                                   side = c("l", "l", "r", "r"), nudge = 0.05) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "baseline"),
                                    aes(x = timepoint, y = shannon),
                                    nudge = 0.05, side = "l", fill = "white", color = "black",
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "follow-up"),
                                    aes(x = timepoint, y = shannon),
                                    nudge = 0.05, side = "r", fill = "white", color = "black",
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        stat_pvalue_manual(shannon_ht, label = "{p.signif}", y.position = max(dftot2$shannon), remove.bracket = FALSE) +
        scale_color_simpsons(guide = "none") +
        scale_fill_simpsons(guide = "none") +
        facet_wrap(~HT_new) +
        theme_Publication() +
        labs(x = "Timepoint", y = "Shannon index", title = "New hypertension diagnoses", color = ""))
ggsave(file.path(resultsfolder, "hypertension_deltashannon.pdf"), width = 4, height = 4)

(pl2 <- dftot2 %>% filter(!is.na(MetSyn_new)) %>% ggplot() +
        geom_line(aes(x = timepoint, y = shannon, group = ID), alpha = 0.1, color = "grey40") +
        geom_jitter(aes(x = timepoint, y = shannon, group = ID, color = MetSyn_new),
                    alpha = 0.5, width = 0) +
        gghalves::geom_half_violin(aes(x = timepoint, y = shannon, fill = MetSyn_new),
                                   side = c("l", "l", "r", "r"), nudge = 0.05) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "baseline"),
                                    aes(x = timepoint, y = shannon),
                                    nudge = 0.05, side = "l", fill = "white", color = "black",
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "follow-up"),
                                    aes(x = timepoint, y = shannon),
                                    nudge = 0.05, side = "r", fill = "white", color = "black",
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        stat_pvalue_manual(shannon_ht, label = "{p.signif}", y.position = max(dftot2$shannon), remove.bracket = FALSE) +
        scale_color_simpsons(guide = "none") +
        scale_fill_simpsons(guide = "none") +
        facet_wrap(~MetSyn_new) +
        theme_Publication() +
        labs(x = "Timepoint", y = "Shannon index", title = "New MetSyn diagnoses", color = ""))
ggsave(file.path(resultsfolder, "metsyn_deltashannon.pdf"), width = 4, height = 4)

#### pl_shan_A — Simple LMM effect sizes: DM, HTN, Dyslipidemia ####

dftot2_lmm <- dftot2 %>%
    group_by(ID) %>%
    mutate(
        DM_bl           = DM[timepoint == "baseline"][1],
        HT_BPMed_bl     = HT_BPMed[timepoint == "baseline"][1],
        Dyslipidemia_bl = Dyslipidemia[timepoint == "baseline"][1]
    ) %>%
    ungroup() %>%
    mutate(timepoint = factor(timepoint, levels = c("baseline", "follow-up")))

extract_lmm_interaction <- function(data, outcome, disease_bl, label) {
    df <- data %>%
        filter(!is.na(.data[[disease_bl]])) %>%
        mutate(
            outcome = .data[[outcome]],
            disease = relevel(factor(.data[[disease_bl]]), ref = "No")
        )
    model   <- lmer(outcome ~ disease * timepoint + (1|ID), data = df)
    cf      <- coef(summary(model))
    ci      <- confint(model, method = "Wald")
    int_row <- grep(":timepoint", rownames(cf), value = TRUE)[1]
    int_ci  <- grep(":timepoint", rownames(ci), value = TRUE)[1]
    data.frame(
        label     = label,
        estimate  = cf[int_row, "Estimate"],
        conf.low  = ci[int_ci, 1],
        conf.high = ci[int_ci, 2],
        p.value   = cf[int_row, "Pr(>|t|)"],
        stringsAsFactors = FALSE
    )
}

shan_lmm_effects_df <- bind_rows(
    extract_lmm_interaction(dftot2_lmm, "shannon", "DM_bl",           "Diabetes"),
    extract_lmm_interaction(dftot2_lmm, "shannon", "HT_BPMed_bl",     "Hypertension"),
    extract_lmm_interaction(dftot2_lmm, "shannon", "Dyslipidemia_bl", "Dyslipidemia")
) %>% mutate(
    sig   = ifelse(p.value < 0.05, "p < 0.05", "p \u2265 0.05"),
    label = factor(label, levels = c("Dyslipidemia", "Hypertension", "Diabetes"))
)

(pl_shan_A <- ggplot(shan_lmm_effects_df, aes(x = estimate, y = label, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.15) +
    geom_point(size = 4) +
    scale_color_manual(
        values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"),
        name = NULL
    ) +
    labs(x = "Shannon change (\u00b1 95% CI)", y = NULL, title = "Shannon diversity change") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_cmb_shannon_change.pdf"), width = 5, height = 4)

#### Per-ethnicity effect of DM on Shannon ####

dat_shan_dm <- dftot %>%
    filter(!is.na(DM), EthnicityTot != "Other") %>%
    mutate(DM           = relevel(factor(DM), ref = "No"),
           EthnicityTot = factor(EthnicityTot))

mod_shan_dm_int    <- lm(shannon ~ DM * EthnicityTot, data = dat_shan_dm)
pint_shan_dm       <- drop1(mod_shan_dm_int, scope = ~DM:EthnicityTot, test = "F")["DM:EthnicityTot", "Pr(>F)"]
pint_shan_dm_label <- paste0("Interaction p = ", format(round(pint_shan_dm, 3), nsmall = 3))

shan_effects_dm <- lapply(levels(dat_shan_dm$EthnicityTot), function(eth) {
    sub <- dat_shan_dm %>% filter(EthnicityTot == eth)
    if (sum(sub$DM == "Yes") < 5) return(NULL)
    m  <- lm(shannon ~ DM, data = sub)
    cf <- coef(summary(m)); ci <- confint(m)
    data.frame(ethnicity = eth, estimate = cf[2,"Estimate"], conf.low = ci[2,1],
               conf.high = ci[2,2], p.value = cf[2,"Pr(>|t|)"], stringsAsFactors = FALSE)
}) %>% bind_rows() %>%
    mutate(sig = ifelse(p.value < 0.05, "p < 0.05", "p \u2265 0.05"),
           ethnicity = forcats::fct_reorder(ethnicity, estimate))

(pl_shan_Aint_dm <- ggplot(shan_effects_dm, aes(x = estimate, y = ethnicity, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.15) +
    geom_point(size = 4) +
    scale_color_manual(values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"), name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8,
             label = pint_shan_dm_label, size = 3.2, fontface = "italic", color = "grey30") +
    labs(x = "Shannon (\u00b1 95% CI)", y = NULL, title = "Diabetes \u00d7 ethnicity") +
    theme_Publication() + theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_dm_shannon_byethnicity.pdf"), width = 5, height = 5)

#### Per-ethnicity effect of HT on Shannon ####

dat_shan_ht <- dftot %>%
    filter(!is.na(HT_BPMed), EthnicityTot != "Other") %>%
    mutate(HT_BPMed = relevel(factor(HT_BPMed), ref = "No"), EthnicityTot = factor(EthnicityTot))

mod_shan_ht_int    <- lm(shannon ~ HT_BPMed * EthnicityTot, data = dat_shan_ht)
pint_shan_ht       <- drop1(mod_shan_ht_int, scope = ~HT_BPMed:EthnicityTot, test = "F")["HT_BPMed:EthnicityTot", "Pr(>F)"]
pint_shan_ht_label <- paste0("Interaction p = ", format(round(pint_shan_ht, 3), nsmall = 3))

shan_effects_ht <- lapply(levels(dat_shan_ht$EthnicityTot), function(eth) {
    sub <- dat_shan_ht %>% filter(EthnicityTot == eth)
    if (sum(sub$HT_BPMed == "Yes") < 5) return(NULL)
    m  <- lm(shannon ~ HT_BPMed, data = sub)
    cf <- coef(summary(m)); ci <- confint(m)
    data.frame(ethnicity = eth, estimate = cf[2,"Estimate"], conf.low = ci[2,1],
               conf.high = ci[2,2], p.value = cf[2,"Pr(>|t|)"], stringsAsFactors = FALSE)
}) %>% bind_rows() %>%
    mutate(sig = ifelse(p.value < 0.05, "p < 0.05", "p \u2265 0.05"),
           ethnicity = forcats::fct_reorder(ethnicity, estimate))

(pl_shan_Aint_ht <- ggplot(shan_effects_ht, aes(x = estimate, y = ethnicity, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.15) +
    geom_point(size = 4) +
    scale_color_manual(values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"), name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8,
             label = pint_shan_ht_label, size = 3.2, fontface = "italic", color = "grey30") +
    labs(x = "Shannon (\u00b1 95% CI)", y = NULL, title = "Hypertension \u00d7 ethnicity") +
    theme_Publication() + theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_ht_shannon_byethnicity.pdf"), width = 5, height = 5)

#### Per-ethnicity effect of LLD on Shannon ####

dat_shan_lld <- dftot %>%
    filter(!is.na(Dyslipidemia), EthnicityTot != "Other") %>%
    mutate(Dyslipidemia = relevel(factor(Dyslipidemia), ref = "No"), EthnicityTot = factor(EthnicityTot))

mod_shan_lld_int    <- lm(shannon ~ Dyslipidemia * EthnicityTot, data = dat_shan_lld)
pint_shan_lld       <- drop1(mod_shan_lld_int, scope = ~Dyslipidemia:EthnicityTot, test = "F")["Dyslipidemia:EthnicityTot", "Pr(>F)"]
pint_shan_lld_label <- paste0("Interaction p = ", format(round(pint_shan_lld, 3), nsmall = 3))

shan_effects_lld <- lapply(levels(dat_shan_lld$EthnicityTot), function(eth) {
    sub <- dat_shan_lld %>% filter(EthnicityTot == eth)
    if (sum(sub$Dyslipidemia == "Yes") < 5) return(NULL)
    m  <- lm(shannon ~ Dyslipidemia, data = sub)
    cf <- coef(summary(m)); ci <- confint(m)
    data.frame(ethnicity = eth, estimate = cf[2,"Estimate"], conf.low = ci[2,1],
               conf.high = ci[2,2], p.value = cf[2,"Pr(>|t|)"], stringsAsFactors = FALSE)
}) %>% bind_rows() %>%
    mutate(sig = ifelse(p.value < 0.05, "p < 0.05", "p \u2265 0.05"),
           ethnicity = forcats::fct_reorder(ethnicity, estimate))

(pl_shan_Aint_lld <- ggplot(shan_effects_lld, aes(x = estimate, y = ethnicity, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.15) +
    geom_point(size = 4) +
    scale_color_manual(values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"), name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8,
             label = pint_shan_lld_label, size = 3.2, fontface = "italic", color = "grey30") +
    labs(x = "Shannon (\u00b1 95% CI)", y = NULL, title = "Dyslipidemia \u00d7 ethnicity") +
    theme_Publication() + theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_lld_shannon_byethnicity.pdf"), width = 5, height = 5)
