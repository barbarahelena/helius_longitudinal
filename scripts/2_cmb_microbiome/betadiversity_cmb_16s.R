## Figure 2 — Cardiometabolic disease and microbiome instability (16S)
## Beta-diversity (Bray-Curtis distance) analyses
##
## Figure panels (named pl_fig2_*):
##   2A: effect sizes (DM, HTN, MetSyn) — dot-whisker
##   2B: diabetes × ethnicity (Dutch + SAS focus + prevalence bar)
##   2D: new-onset diabetes, Dutch + SAS only
##
## Supplementary PDFs also produced:
##   distance_dm.pdf, braycurtis_deltahba1c.pdf
##   distance_ethnictiy_diabetes.pdf (all ethnicities)
##   clinicaloutcomes_bray.pdf, distance_hypertension.pdf, distance_metsyn.pdf
##   newdiabetes.pdf (all ethnicities)

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)

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
df <- readRDS("data/16s/archive/clin_betadiversity.RDS") %>% dplyr::select(1:2, sampleID = ID, 4:5)
helius <- readRDS("data/clinicaldata_long.RDS")
df <- left_join(df, helius, by = c("sampleID"))
ev_bray <- read.csv("results/1_longitudinal_change/ordination/expl_var_bray.csv", header = FALSE)
heliusdist <- readRDS("data/16s/archive/braydistance_delta.RDS")

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

#### Figure 2A — Effect sizes: DM, HTN, MetSyn on Bray-Curtis ####

# Supplementary individual plots
ggplot(data = heliusdist %>% filter(!is.na(DM)), aes(x = DM, y = distance)) +
    geom_violin(colour = NA, aes(fill = DM)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes (baseline)", title = "Diabetes") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication()
ggsave(file.path(resultsfolder, "distance_dm.pdf"), width = 4, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(HbA1c_delta)), aes(x = distance, y = HbA1c_delta)) +
    geom_jitter(color = "royalblue", alpha = 0.3) +
    geom_smooth(color = "black", method = "lm") +
    labs(y = "Delta HbA1c", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and HbA1c change") +
    stat_cor() +
    theme_Publication()
ggsave(file.path(resultsfolder, "braycurtis_deltahba1c.pdf"), width = 4.5, height = 5)

# Effect size dot-whisker: DM, HTN, MetSyn — uses base-R lm(), no extra dependencies
extract_lm_effect <- function(model, label) {
    cf <- coef(summary(model))
    ci <- confint(model)
    data.frame(
        label     = label,
        estimate  = cf[2, "Estimate"],
        conf.low  = ci[2, 1],
        conf.high = ci[2, 2],
        p.value   = cf[2, "Pr(>|t|)"],
        stringsAsFactors = FALSE
    )
}

mod_dm  <- lm(distance ~ DM,       data = heliusdist %>% filter(!is.na(DM))       %>% mutate(DM       = relevel(factor(DM),       ref = "No")))
mod_ht  <- lm(distance ~ HT_BPMed, data = heliusdist %>% filter(!is.na(HT_BPMed)) %>% mutate(HT_BPMed = relevel(factor(HT_BPMed), ref = "No")))
mod_lld <- lm(distance ~ LLD,      data = heliusdist %>% filter(!is.na(LLD))      %>% mutate(LLD      = relevel(factor(LLD),      ref = "No")))

effects_df <- bind_rows(
    extract_lm_effect(mod_dm,  "Diabetes"),
    extract_lm_effect(mod_ht,  "Hypertension"),
    extract_lm_effect(mod_lld, "Dyslipidemia")
) %>% mutate(
    sig   = ifelse(p.value < 0.05, "p < 0.05", "p \u2265 0.05"),
    label = factor(label, levels = c("Dyslipidemia", "Hypertension", "Diabetes"))
)

(pl_fig2_A <- ggplot(effects_df, aes(x = estimate, y = label, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8) +
    geom_point(size = 4) +
    scale_color_manual(
        values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"),
        name = NULL
    ) +
    labs(x = "Bray-Curtis (\u00b1 95% CI)",
         y = NULL,
         title = "Effect on microbiome instability") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_cmb_bray.pdf"), width = 5, height = 4)

#### Figure 2A (interaction) — Per-ethnicity effect of DM on Bray-Curtis ####

# Interaction model: DM × EthnicityTot
dat_int <- heliusdist %>%
    filter(!is.na(DM), EthnicityTot != "Other") %>%
    mutate(DM           = relevel(factor(DM), ref = "No"),
           EthnicityTot = factor(EthnicityTot))

mod_dm_int <- lm(distance ~ DM * EthnicityTot, data = dat_int)
pint_dm    <- drop1(mod_dm_int, scope = ~DM:EthnicityTot, test = "F")["DM:EthnicityTot", "Pr(>F)"]
pint_label <- paste0("Interaction p = ", format(round(pint_dm, 3), nsmall = 3))

# Per-ethnicity stratified LMs for clean per-group CIs
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

(pl_fig2_Aint <- ggplot(eth_effects_dm, aes(x = estimate, y = ethnicity, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8) +
    geom_point(size = 4) +
    scale_color_manual(
        values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"),
        name = NULL
    ) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8,
             label = pint_label, size = 3.2, fontface = "italic", color = "grey30") +
    labs(x = "Bray-Curtis (\u00b1 95% CI)",
         y = NULL,
         title = "Diabetes \u00d7 ethnicity") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_dm_byethnicity.pdf"), width = 5, height = 5)

#### Per-ethnicity effect of HT on Bray-Curtis ####

dat_int_ht <- heliusdist %>%
    filter(!is.na(HT_BPMed), EthnicityTot != "Other") %>%
    mutate(HT_BPMed     = relevel(factor(HT_BPMed), ref = "No"),
           EthnicityTot = factor(EthnicityTot))

mod_ht_int    <- lm(distance ~ HT_BPMed * EthnicityTot, data = dat_int_ht)
pint_ht       <- drop1(mod_ht_int, scope = ~HT_BPMed:EthnicityTot, test = "F")["HT_BPMed:EthnicityTot", "Pr(>F)"]
pint_ht_label <- paste0("Interaction p = ", format(round(pint_ht, 3), nsmall = 3))

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

(pl_fig2_Aint_ht <- ggplot(eth_effects_ht, aes(x = estimate, y = ethnicity, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8) +
    geom_point(size = 4) +
    scale_color_manual(
        values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"),
        name = NULL
    ) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8,
             label = pint_ht_label, size = 3.2, fontface = "italic", color = "grey30") +
    labs(x = "Bray-Curtis effect of hypertension (\u00b1 95% CI)",
         y = NULL,
         title = "Hypertension effect on\nmicrobiome instability by ethnicity") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_ht_byethnicity.pdf"), width = 5, height = 5)

#### Per-ethnicity effect of LLD on Bray-Curtis ####

dat_int_lld <- heliusdist %>%
    filter(!is.na(LLD), EthnicityTot != "Other") %>%
    mutate(LLD          = relevel(factor(LLD), ref = "No"),
           EthnicityTot = factor(EthnicityTot))

mod_lld_int    <- lm(distance ~ LLD * EthnicityTot, data = dat_int_lld)
pint_lld       <- drop1(mod_lld_int, scope = ~LLD:EthnicityTot, test = "F")["LLD:EthnicityTot", "Pr(>F)"]
pint_lld_label <- paste0("Interaction p = ", format(round(pint_lld, 3), nsmall = 3))

eth_effects_lld <- lapply(levels(dat_int_lld$EthnicityTot), function(eth) {
    sub <- dat_int_lld %>% filter(EthnicityTot == eth)
    if (sum(sub$LLD == "Yes") < 5) return(NULL)
    m  <- lm(distance ~ LLD, data = sub)
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

(pl_fig2_Aint_lld <- ggplot(eth_effects_lld, aes(x = estimate, y = ethnicity, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8) +
    geom_point(size = 4) +
    scale_color_manual(
        values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"),
        name = NULL
    ) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8,
             label = pint_lld_label, size = 3.2, fontface = "italic", color = "grey30") +
    labs(x = "Bray-Curtis (\u00b1 95% CI)",
         y = NULL,
         title = "Dyslipidemia \u00d7 ethnicity") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_lld_byethnicity.pdf"), width = 5, height = 5)

#### Figure 2B — Ethnicity × diabetes interaction ####

# Supplementary: all ethnicities
heliusdist %>% filter(!is.na(DM)) %>% group_by(EthnicityTot, DM) %>% summarise(count = length(DM), .groups = "drop_last")
ggplot(data = heliusdist %>% filter(!is.na(DM)), aes(x = DM, y = distance)) +
    geom_violin(colour = NA, aes(fill = DM)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes", title = "Diabetes") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.format") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave(file.path(resultsfolder, "distance_ethnictiy_diabetes.pdf"), width = 6, height = 11)

# Diabetes prevalence by ethnicity — left sub-panel of 2B
dm_prev <- heliusdist %>%
    filter(!is.na(DM), EthnicityTot != "Other") %>%
    group_by(EthnicityTot) %>%
    summarise(prev = mean(DM == "Yes") * 100, n = n(), .groups = "drop") %>%
    mutate(EthnicityTot = forcats::fct_reorder(EthnicityTot, prev))

pl_prev <- ggplot(dm_prev, aes(x = prev, y = EthnicityTot, fill = EthnicityTot)) +
    geom_col(alpha = 0.9, width = 0.6) +
    geom_text(aes(label = sprintf("n=%d", n)), hjust = -0.1, size = 2.8) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(x = "Diabetes prevalence (%)", y = NULL, title = "Diabetes prevalence\nby ethnicity") +
    theme_Publication() +
    theme(axis.text.y = element_text(size = rel(0.75)))

# Dutch + South-Asian Surinamese focus — right sub-panel of 2B
# (the two groups with highest diabetes burden and largest n)
pl_violin_sas <- ggplot(
    data = heliusdist %>% filter(!is.na(DM),
                                 EthnicityTot %in% c("Dutch", "South-Asian Surinamese")),
    aes(x = DM, y = distance)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot), alpha = 0.75) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Diabetes",
         title = "Diabetes \u00d7 ethnicity\n(Dutch vs. South-Asian Surinamese)") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0,
                       hide.ns = FALSE, label = "p.format") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()

(pl_fig2_B <- ggarrange(pl_prev, pl_violin_sas, ncol = 2, widths = c(1, 1.4)))
ggsave(file.path(resultsfolder, "distance_diabetes_ethnicity_focus.pdf"), width = 9, height = 5)

#### Figure 2: Hypertension compound panel ####

ht_prev <- heliusdist %>%
    filter(!is.na(HT_BPMed), EthnicityTot != "Other") %>%
    group_by(EthnicityTot) %>%
    summarise(prev = mean(HT_BPMed == "Yes") * 100, n = n(), .groups = "drop") %>%
    mutate(EthnicityTot = forcats::fct_reorder(EthnicityTot, prev))

pl_ht_prev <- ggplot(ht_prev, aes(x = prev, y = EthnicityTot, fill = EthnicityTot)) +
    geom_col(alpha = 0.9, width = 0.6) +
    geom_text(aes(label = sprintf("n=%d", n)), hjust = -0.1, size = 2.8) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(x = "Hypertension prevalence (%)", y = NULL, title = "Hypertension prevalence\nby ethnicity") +
    theme_Publication() +
    theme(axis.text.y = element_text(size = rel(0.75)))

pl_violin_sas_ht <- ggplot(
    data = heliusdist %>% filter(!is.na(HT_BPMed),
                                 EthnicityTot %in% c("Dutch", "South-Asian Surinamese")),
    aes(x = HT_BPMed, y = distance)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot), alpha = 0.75) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Hypertension",
         title = "Hypertension \u00d7 ethnicity\n(Dutch vs. South-Asian Surinamese)") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0,
                       hide.ns = FALSE, label = "p.format") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()

(pl_fig2_C <- ggarrange(pl_ht_prev, pl_violin_sas_ht, ncol = 2, widths = c(1, 1.4)))
ggsave(file.path(resultsfolder, "distance_hypertension_ethnicity_focus.pdf"), width = 9, height = 5)

#### Figure 2: Dyslipidemia compound panel ####

lld_prev <- heliusdist %>%
    filter(!is.na(LLD), EthnicityTot != "Other") %>%
    group_by(EthnicityTot) %>%
    summarise(prev = mean(LLD == "Yes") * 100, n = n(), .groups = "drop") %>%
    mutate(EthnicityTot = forcats::fct_reorder(EthnicityTot, prev))

pl_lld_prev <- ggplot(lld_prev, aes(x = prev, y = EthnicityTot, fill = EthnicityTot)) +
    geom_col(alpha = 0.9, width = 0.6) +
    geom_text(aes(label = sprintf("n=%d", n)), hjust = -0.1, size = 2.8) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(x = "Dyslipidemia prevalence (%)", y = NULL, title = "Dyslipidemia prevalence\nby ethnicity") +
    theme_Publication() +
    theme(axis.text.y = element_text(size = rel(0.75)))

pl_violin_sas_lld <- ggplot(
    data = heliusdist %>% filter(!is.na(LLD),
                                 EthnicityTot %in% c("Dutch", "South-Asian Surinamese")),
    aes(x = LLD, y = distance)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot), alpha = 0.75) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Dyslipidemia / LLD use",
         title = "Dyslipidemia \u00d7 ethnicity\n(Dutch vs. South-Asian Surinamese)") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0,
                       hide.ns = FALSE, label = "p.format") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()

(pl_fig2_D <- ggarrange(pl_lld_prev, pl_violin_sas_lld, ncol = 2, widths = c(1, 1.4)))
ggsave(file.path(resultsfolder, "distance_dyslipidemia_ethnicity_focus.pdf"), width = 9, height = 5)

#### Combined prevalence panel: DM, HT, LLD — baseline + follow-up ####

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
        filter(ID %in% paired_ids, !is.na(LLD), EthnicityTot != "Other") %>%
        group_by(EthnicityTot, timepoint) %>%
        summarise(prev = mean(LLD == "Yes") * 100, .groups = "drop") %>%
        mutate(condition = "Dyslipidemia")
) %>%
    mutate(
        condition = factor(condition, levels = c("Diabetes", "Hypertension", "Dyslipidemia")),
        timepoint = factor(timepoint, levels = c("follow-up", "baseline"))
    )

# Order ethnicities by baseline DM prevalence
eth_order <- prev_long %>%
    filter(condition == "Diabetes", timepoint == "baseline") %>%
    arrange(prev) %>%
    pull(EthnicityTot) %>%
    as.character()
prev_long <- prev_long %>% mutate(EthnicityTot = factor(EthnicityTot, levels = eth_order))

(pl_fig2_prev <- ggplot(prev_long,
                        aes(x = prev, y = EthnicityTot, fill = EthnicityTot, alpha = timepoint)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65,
             color = "black", linewidth = 0.3) +
    geom_text(
        data = prev_long %>% filter(condition == "Diabetes",
                                    EthnicityTot == tail(eth_order, 1)),
        aes(x = prev + 0.5, y = EthnicityTot,
            label = ifelse(timepoint == "baseline", "Baseline", "Follow-up"),
            group = timepoint),
        position = position_dodge(width = 0.75),
        hjust = 0, size = 3, color = "grey30", inherit.aes = FALSE
    ) +
    scale_alpha_manual(values = c("baseline" = 0.85, "follow-up" = 0.35), guide = "none") +
    scale_fill_manual(values = eth_colors, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    facet_wrap(~condition, scales = "free_x", nrow = 1) +
    labs(x = "Prevalence (%)", y = NULL,
         title = "CMB disease prevalence") +
    theme_Publication() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = rel(0.75))))
ggsave(file.path(resultsfolder, "cmb_prevalence_combined.pdf"), width = 12, height = 5)

#### Figure 2C — Hypertension & metabolic syndrome ####

# clinicaloutcomes_bray.pdf (PCoA panels for new CMB diagnoses)
(dmnewbray <- df %>% filter(!is.na(DM_new)) %>% filter(timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = DM_new, fill = DM_new), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = DM_new), size = 1, alpha = 0.5) +
        ggtitle("New diabetes") +
        xlab(paste0('PCo1 (', round(ev_bray$V1[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(ev_bray$V1[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        scale_alpha_manual(guide = "none") +
        theme_Publication() +
        labs(alpha = "") )

(htnewbray <- df %>% filter(!is.na(HT_new)) %>% filter(timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = HT_new, fill = HT_new), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = HT_new), size = 1, alpha = 0.5) +
        ggtitle("New hypertension") +
        xlab(paste0('PCo1 (', round(ev_bray$V1[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(ev_bray$V1[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        scale_alpha_manual(guide = "none") +
        theme_Publication() +
        labs(alpha = "") )

(metsynnewbray <- df %>% filter(!is.na(MetSyn_new)) %>% filter(timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = MetSyn_new, fill = MetSyn_new), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = MetSyn_new), size = 1, alpha = 0.5) +
        ggtitle("New metabolic syndrome") +
        xlab(paste0('PCo1 (', round(ev_bray$V1[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(ev_bray$V1[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        scale_alpha_manual(guide = "none") +
        theme_Publication() +
        labs(alpha = "") )

(lldnewbray <- df %>% filter(!is.na(LLD_new)) %>% filter(timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = LLD_new, fill = LLD_new), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = LLD_new), size = 1, alpha = 0.5) +
        ggtitle("New lipid lowering drug use") +
        xlab(paste0('PCo1 (', round(ev_bray$V1[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(ev_bray$V1[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        scale_alpha_manual(guide = "none") +
        theme_Publication() +
        labs(alpha = "") )

ggarrange(dmnewbray, htnewbray, metsynnewbray, lldnewbray, nrow = 1,
          labels = LETTERS[1:4])
ggsave(file.path(resultsfolder, "clinicaloutcomes_bray.pdf"), width = 18, height = 5)

# distance_hypertension.pdf
(pl_fig2_C1 <- ggplot(data = heliusdist %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = distance)) +
    geom_violin(colour = NA, aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension (baseline)", title = "Hypertension") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication())
ggsave(file.path(resultsfolder, "distance_hypertension.pdf"), width = 4, height = 5)

# distance_dyslipidemia.pdf
(pl_fig2_C2 <- ggplot(data = heliusdist %>% filter(!is.na(LLD)), aes(x = LLD, y = distance)) +
    geom_violin(colour = NA, aes(fill = LLD)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Dyslipidemia / LLD use (baseline)", title = "Dyslipidemia") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication())
ggsave(file.path(resultsfolder, "distance_dyslipidemia.pdf"), width = 4, height = 5)

#### Figure 2D — New-onset diabetes, Dutch + South-Asian Surinamese ####

# Supplementary: all ethnicities combined
heliusdist %>% filter(!is.na(DM_new)) %>%
    ggplot(aes(x = DM_new, y = distance, fill = DM_new)) +
    geom_violin(colour = NA) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "New diabetes", title = "New diabetes diagnosis") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0,
                       hide.ns = FALSE, label = "p.format") +
    theme_Publication()
ggsave(file.path(resultsfolder, "newdiabetes.pdf"), width = 4.5, height = 5)

# Main figure panel: restrict to Dutch + SAS where directional signal is strongest
(pl_fig2_D <- heliusdist %>%
     filter(!is.na(DM_new),
            EthnicityTot %in% c("Dutch", "South-Asian Surinamese")) %>%
     ggplot(aes(x = DM_new, y = distance, fill = EthnicityTot)) +
     geom_violin(colour = NA, alpha = 0.75) +
     geom_boxplot(fill = "white", width = 0.2) +
     scale_fill_manual(values = eth_colors, guide = "none") +
     labs(y = "Bray-Curtis dissimilarity over FU time", x = "New diabetes diagnosis",
          title = "New-onset diabetes\n(Dutch & South-Asian Surinamese)") +
     stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0,
                        hide.ns = FALSE, label = "p.format") +
     facet_wrap(~EthnicityTot) +
     theme_Publication())
ggsave(file.path(resultsfolder, "newdiabetes_sas_dutch.pdf"), width = 5.5, height = 5)
