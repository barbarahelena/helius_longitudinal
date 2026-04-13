## Archive — Beta-diversity exploratory plots (no longer in assembled Figure 2)
##
## These plots were produced during analysis but are not used in the main figure.
## This script sources the main beta-diversity script to reuse data objects.

source("scripts/2_cmb_microbiome/1_betadiversity_cmb_16s.R")

## Additional data needed for PCoA and new-onset panels
df      <- readRDS("data/16s/clin_bray.RDS") %>% dplyr::select(1:2, sampleID = ID, 4:5)
df      <- left_join(df, helius, by = c("sampleID"))
ev_bray <- read.csv("results/1_longitudinal_change/ordination/expl_var_bray.csv", header = FALSE)

#### Simple DM / HbA1c violin + scatter ####

ggplot(data = heliusdist %>% filter(!is.na(DM)), aes(x = DM, y = distance)) +
    geom_violin(colour = NA, aes(fill = DM)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Diabetes (baseline)", title = "Diabetes") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication()
ggsave(file.path(resultsfolder, "distance_dm.pdf"), width = 4, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(HbA1c_delta)), aes(x = distance, y = HbA1c_delta)) +
    geom_jitter(color = "royalblue", alpha = 0.3) +
    geom_smooth(color = "black", method = "lm") +
    labs(y = "Delta HbA1c", x = "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and HbA1c change") +
    stat_cor() +
    theme_Publication()
ggsave(file.path(resultsfolder, "braycurtis_deltahba1c.pdf"), width = 4.5, height = 5)

#### pl_fig2_A — Simple effect sizes: DM, HTN, Dyslipidemia ####

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

mod_dm  <- lm(distance ~ DM,          data = heliusdist %>% filter(!is.na(DM))          %>% mutate(DM          = relevel(factor(DM),          ref = "No")))
mod_ht  <- lm(distance ~ HT_BPMed,    data = heliusdist %>% filter(!is.na(HT_BPMed))    %>% mutate(HT_BPMed    = relevel(factor(HT_BPMed),    ref = "No")))
mod_lld <- lm(distance ~ Dyslipidemia, data = heliusdist %>% filter(!is.na(Dyslipidemia)) %>% mutate(Dyslipidemia = relevel(factor(Dyslipidemia), ref = "No")))

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
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.15) +
    geom_point(size = 4) +
    scale_color_manual(
        values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"),
        name = NULL
    ) +
    labs(x = "Bray-Curtis (\u00b1 95% CI)", y = NULL, title = "Microbiome instability") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_cmb_bray.pdf"), width = 5, height = 4)

#### pl_fig2_Aint — Per-ethnicity effect of DM ####

mod_dm_int <- lm(distance ~ DM * EthnicityTot, data = dat_int)
pint_dm    <- drop1(mod_dm_int, scope = ~DM:EthnicityTot, test = "F")["DM:EthnicityTot", "Pr(>F)"]
pint_label <- paste0("Interaction p = ", format(round(pint_dm, 3), nsmall = 3))

(pl_fig2_Aint <- ggplot(eth_effects_dm, aes(x = estimate, y = ethnicity, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.15) +
    geom_point(size = 4) +
    scale_color_manual(values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"), name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8,
             label = pint_label, size = 3.2, fontface = "italic", color = "grey30") +
    labs(x = "Bray-Curtis (\u00b1 95% CI)", y = NULL, title = "Diabetes \u00d7 ethnicity") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_dm_byethnicity.pdf"), width = 5, height = 5)

#### pl_fig2_Aint_ht — Per-ethnicity effect of HT ####

mod_ht_int    <- lm(distance ~ HT_BPMed * EthnicityTot, data = dat_int_ht)
pint_ht       <- drop1(mod_ht_int, scope = ~HT_BPMed:EthnicityTot, test = "F")["HT_BPMed:EthnicityTot", "Pr(>F)"]
pint_ht_label <- paste0("Interaction p = ", format(round(pint_ht, 3), nsmall = 3))

(pl_fig2_Aint_ht <- ggplot(eth_effects_ht, aes(x = estimate, y = ethnicity, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.15) +
    geom_point(size = 4) +
    scale_color_manual(values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"), name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8,
             label = pint_ht_label, size = 3.2, fontface = "italic", color = "grey30") +
    labs(x = "Bray-Curtis effect of hypertension (\u00b1 95% CI)", y = NULL,
         title = "Hypertension effect on\nmicrobiome instability by ethnicity") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_ht_byethnicity.pdf"), width = 5, height = 5)

#### pl_fig2_Aint_lld — Per-ethnicity effect of Dyslipidemia ####

mod_lld_int    <- lm(distance ~ Dyslipidemia * EthnicityTot, data = dat_int_lld)
pint_lld       <- drop1(mod_lld_int, scope = ~Dyslipidemia:EthnicityTot, test = "F")["Dyslipidemia:EthnicityTot", "Pr(>F)"]
pint_lld_label <- paste0("Interaction p = ", format(round(pint_lld, 3), nsmall = 3))

(pl_fig2_Aint_lld <- ggplot(eth_effects_lld, aes(x = estimate, y = ethnicity, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", linewidth = 0.8, width = 0.15) +
    geom_point(size = 4) +
    scale_color_manual(values = c("p < 0.05" = "#e63946", "p \u2265 0.05" = "grey55"), name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8,
             label = pint_lld_label, size = 3.2, fontface = "italic", color = "grey30") +
    labs(x = "Bray-Curtis (\u00b1 95% CI)", y = NULL, title = "Dyslipidemia \u00d7 ethnicity") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(file.path(resultsfolder, "effectsize_lld_byethnicity.pdf"), width = 5, height = 5)

#### pl_fig2_B — Diabetes compound panel (prevalence bar + Dutch/SAS violin) ####

heliusdist %>% filter(!is.na(DM)) %>%
    group_by(EthnicityTot, DM) %>% summarise(count = length(DM), .groups = "drop_last")

ggplot(data = heliusdist %>% filter(!is.na(DM)), aes(x = DM, y = distance)) +
    geom_violin(colour = NA, aes(fill = DM)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Diabetes", title = "Diabetes") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.format") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave(file.path(resultsfolder, "distance_ethnictiy_diabetes.pdf"), width = 6, height = 11)

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

#### pl_fig2_C — Hypertension compound panel ####

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

#### pl_fig2_D — Dyslipidemia compound panel ####

lld_prev <- heliusdist %>%
    filter(!is.na(Dyslipidemia), EthnicityTot != "Other") %>%
    group_by(EthnicityTot) %>%
    summarise(prev = mean(Dyslipidemia == "Yes") * 100, n = n(), .groups = "drop") %>%
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
    data = heliusdist %>% filter(!is.na(Dyslipidemia),
                                 EthnicityTot %in% c("Dutch", "South-Asian Surinamese")),
    aes(x = Dyslipidemia, y = distance)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot), alpha = 0.75) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Dyslipidemia",
         title = "Dyslipidemia \u00d7 ethnicity\n(Dutch vs. South-Asian Surinamese)") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0,
                       hide.ns = FALSE, label = "p.format") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()

(pl_fig2_D_lld <- ggarrange(pl_lld_prev, pl_violin_sas_lld, ncol = 2, widths = c(1, 1.4)))
ggsave(file.path(resultsfolder, "distance_dyslipidemia_ethnicity_focus.pdf"), width = 9, height = 5)

#### PCoA panels — new-onset CMB diagnoses ####

(dmnewbray <- df %>% filter(!is.na(DM_new), timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = DM_new, fill = DM_new), type = "norm", alpha = 0.1) +
        geom_point(aes(color = DM_new), size = 1, alpha = 0.5) +
        ggtitle("New diabetes") +
        xlab(paste0('PCo1 (', round(ev_bray$V1[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(ev_bray$V1[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        theme_Publication())

(htnewbray <- df %>% filter(!is.na(HT_new), timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = HT_new, fill = HT_new), type = "norm", alpha = 0.1) +
        geom_point(aes(color = HT_new), size = 1, alpha = 0.5) +
        ggtitle("New hypertension") +
        xlab(paste0('PCo1 (', round(ev_bray$V1[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(ev_bray$V1[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        theme_Publication())

(metsynnewbray <- df %>% filter(!is.na(MetSyn_new), timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = MetSyn_new, fill = MetSyn_new), type = "norm", alpha = 0.1) +
        geom_point(aes(color = MetSyn_new), size = 1, alpha = 0.5) +
        ggtitle("New metabolic syndrome") +
        xlab(paste0('PCo1 (', round(ev_bray$V1[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(ev_bray$V1[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        theme_Publication())

(lldnewbray <- df %>% filter(!is.na(LLD_new), timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = LLD_new, fill = LLD_new), type = "norm", alpha = 0.1) +
        geom_point(aes(color = LLD_new), size = 1, alpha = 0.5) +
        ggtitle("New lipid lowering drug use") +
        xlab(paste0('PCo1 (', round(ev_bray$V1[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(ev_bray$V1[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        theme_Publication())

ggarrange(dmnewbray, htnewbray, metsynnewbray, lldnewbray, nrow = 1, labels = LETTERS[1:4])
ggsave(file.path(resultsfolder, "clinicaloutcomes_bray.pdf"), width = 18, height = 5)

#### Simple HT / Dyslipidemia violins ####

(pl_fig2_C1 <- ggplot(data = heliusdist %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = distance)) +
    geom_violin(colour = NA, aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Hypertension (baseline)", title = "Hypertension") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication())
ggsave(file.path(resultsfolder, "distance_hypertension.pdf"), width = 4, height = 5)

(pl_fig2_C2 <- ggplot(data = heliusdist %>% filter(!is.na(Dyslipidemia)), aes(x = Dyslipidemia, y = distance)) +
    geom_violin(colour = NA, aes(fill = Dyslipidemia)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Dyslipidemia (baseline)", title = "Dyslipidemia") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication())
ggsave(file.path(resultsfolder, "distance_dyslipidemia.pdf"), width = 4, height = 5)

#### New-onset diabetes panels ####

heliusdist %>% filter(!is.na(DM_new)) %>%
    ggplot(aes(x = DM_new, y = distance, fill = DM_new)) +
    geom_violin(colour = NA) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "New diabetes", title = "New diabetes diagnosis") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0,
                       hide.ns = FALSE, label = "p.format") +
    theme_Publication()
ggsave(file.path(resultsfolder, "newdiabetes.pdf"), width = 4.5, height = 5)

(pl_fig2_D <- heliusdist %>%
     filter(!is.na(DM_new), EthnicityTot %in% c("Dutch", "South-Asian Surinamese")) %>%
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

#### pl_fig2_D_violin — Bray-Curtis by DM status, all ethnicities ####

eth_order_dm <- heliusdist %>%
    filter(!is.na(DM), EthnicityTot != "Other") %>%
    group_by(EthnicityTot, DM) %>%
    summarise(m = mean(distance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = DM, values_from = m) %>%
    mutate(diff = Yes - No) %>%
    arrange(desc(diff)) %>%
    pull(EthnicityTot) %>%
    as.character()

dat_violin_dm <- heliusdist %>%
    filter(!is.na(DM), EthnicityTot != "Other") %>%
    mutate(DM           = relevel(factor(DM), ref = "No"),
           EthnicityTot = factor(EthnicityTot, levels = eth_order_dm))

pvals_dm <- dat_violin_dm %>%
    group_by(EthnicityTot) %>%
    summarise(p.value = wilcox.test(distance ~ DM)$p.value, .groups = "drop") %>%
    mutate(group1     = "No", group2 = "Yes",
           label      = ifelse(p.value < 0.001, "p < 0.001", paste0("p = ", format(round(p.value, 3), nsmall = 3))),
           y.position = max(dat_violin_dm$distance, na.rm = TRUE) * 1.05)

(pl_fig2_D_violin <- ggplot(dat_violin_dm, aes(x = DM, y = distance, fill = EthnicityTot)) +
    geom_violin(aes(alpha = DM), colour = NA) +
    geom_boxplot(fill = "white", colour = "grey30", width = 0.2, outlier.shape = NA) +
    stat_pvalue_manual(pvals_dm %>% filter(p.value < 0.05), label = "label",
                       tip.length = 0, bracket.size = 0, color = "black", size = 3) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    scale_alpha_manual(values = c("No" = 0.35, "Yes" = 0.85), guide = "none") +
    facet_wrap(~EthnicityTot, nrow = 1) +
    labs(x = "Diabetes at baseline", y = "Bray-Curtis dissimilarity",
         title = "Microbiome instability by diabetes status") +
    theme_Publication() +
    theme(strip.text = element_text(size = rel(0.65), face = "bold")))
ggsave(file.path(resultsfolder, "violin_dm_byethnicity.pdf"), width = 14, height = 5)

#### pl_fig2_E_violin — Bray-Curtis by Dyslipidemia status, all ethnicities ####

eth_order_lld <- heliusdist %>%
    filter(!is.na(Dyslipidemia), EthnicityTot != "Other") %>%
    group_by(EthnicityTot, Dyslipidemia) %>%
    summarise(m = mean(distance, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Dyslipidemia, values_from = m) %>%
    mutate(diff = Yes - No) %>%
    arrange(desc(diff)) %>%
    pull(EthnicityTot) %>%
    as.character()

dat_violin_lld <- heliusdist %>%
    filter(!is.na(Dyslipidemia), EthnicityTot != "Other") %>%
    mutate(Dyslipidemia = relevel(factor(Dyslipidemia), ref = "No"),
           EthnicityTot = factor(EthnicityTot, levels = eth_order_lld))

pvals_lld <- dat_violin_lld %>%
    group_by(EthnicityTot) %>%
    summarise(p.value = wilcox.test(distance ~ Dyslipidemia)$p.value, .groups = "drop") %>%
    mutate(group1     = "No", group2 = "Yes",
           label      = ifelse(p.value < 0.001, "p < 0.001", paste0("p = ", format(round(p.value, 3), nsmall = 3))),
           y.position = max(dat_violin_lld$distance, na.rm = TRUE) * 1.05)

(pl_fig2_E_violin <- ggplot(dat_violin_lld, aes(x = Dyslipidemia, y = distance, fill = EthnicityTot)) +
    geom_violin(aes(alpha = Dyslipidemia), colour = NA) +
    geom_boxplot(fill = "white", colour = "grey30", width = 0.2, outlier.shape = NA) +
    stat_pvalue_manual(pvals_lld %>% filter(p.value < 0.05), label = "label",
                       tip.length = 0, bracket.size = 0, color = "black", size = 3) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    scale_alpha_manual(values = c("No" = 0.35, "Yes" = 0.85), guide = "none") +
    facet_wrap(~EthnicityTot, nrow = 1) +
    labs(x = "Dyslipidemia at baseline", y = "Bray-Curtis dissimilarity",
         title = "Microbiome instability by dyslipidemia status") +
    theme_Publication() +
    theme(strip.text = element_text(size = rel(0.65), face = "bold")))
ggsave(file.path(resultsfolder, "violin_lld_byethnicity.pdf"), width = 14, height = 5)
