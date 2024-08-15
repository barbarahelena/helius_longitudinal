## Calculate alpha diversity differences

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
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
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
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

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

helius <- readRDS("data/clinicaldata_long.RDS")
dftot <- left_join(df2 %>% filter(timepoint == "baseline"), helius %>% filter(timepoint == "baseline"), 
                   by = c("ID", "timepoint", "sampleID")) %>% 
    filter(EthnicityTot != "Other") %>% droplevels(.)

pairedids <- df2$sampleID[which(!is.na(df2$shannon_delta))]
dftot2 <- left_join(df, helius, by = c("ID", "timepoint", "sampleID")) %>% 
            filter(sampleID %in% pairedids) %>% 
            filter(! EthnicityTot %in% "Other" ) %>% droplevels(.)

betadiv <- readRDS("data/16s/braydistance_delta.RDS") %>% dplyr::select(1:4)
dftot3a <- left_join(dftot2 %>% filter(timepoint == "baseline"), 
                     betadiv %>% dplyr::select(sampleID = sampleID1, everything(.)))
dftot3b <- left_join(dftot2 %>% filter(timepoint == "follow-up"), 
                     betadiv %>% dplyr::select(sampleID = sampleID2, everything(.)))
dftot3 <- full_join(dftot3a, dftot3b) %>% droplevels(.) %>% 
    full_join(., df2)

#### Output folder ####
resultsfolder <- "results/alphadiversity"
dir.create(resultsfolder, showWarnings = FALSE)

#### Baseline and FU alpha diversity ####
comp <- rev(list(c("Dutch", "South-Asian Surinamese"), c("Dutch", "Ghanaian"), 
                 c("Dutch", "African Surinamese"),c("Dutch", "Turkish"), c("Dutch", "Moroccan")))
ggplot(data = dftot3 %>% filter(!is.na(EthnicityTot)), 
       aes(x = fct_reorder(EthnicityTot, .x = shannon, .fun = median), y = shannon)) +
    geom_violin(aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Shannon index", title = "Shannon index", x = "") +
    facet_wrap(~timepoint) +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "wilcox.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/alphadiversity/shannon_time_ethnicities.pdf", width = 8, height = 7)

ggplot(data = dftot3 %>% filter(!is.na(EthnicityTot)), 
       aes(x = timepoint, y = shannon)) +
    geom_violin(aes(fill = EthnicityTot, alpha = timepoint)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
    labs(y = "Shannon index", title = "Shannon index", x = "") +
    facet_wrap(~EthnicityTot) +
    stat_compare_means(tip.length = 0, hide.ns = TRUE, label.x = 1.5,
                       label = "p.signif", method = "wilcox.test") +
    theme_Publication()
ggsave("results/alphadiversity/shannon_ethnicities_time.pdf", width = 7, height = 7)

#### Delta alpha diversity ####
comp <- rev(list(c("Dutch", "South-Asian Surinamese"), c("Dutch", "Moroccan"), c("Dutch", "Turkish")))
ggplot(data = dftot %>% filter(!is.na(shannon_delta)), # selecting baseline samples
       aes(x = fct_reorder(EthnicityTot, .x = shannon_delta, .fun = median), y = shannon_delta)) +
    geom_violin(aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Shannon index difference", title = "Delta Shannon", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/alphadiversity/shannon_ethnicities.pdf", width = 6, height = 5)

comp <- list(c("Dutch", "South-Asian Surinamese"), c("Dutch", "African Surinamese"), c("Dutch", "Turkish"))
ggplot(data = dftot %>% filter(!is.na(richness_delta)), # selecting baseline samples (those have deltas available)
       aes(x = fct_reorder(EthnicityTot, .x = richness_delta, .fun = median), y = richness_delta)) +
    geom_violin(aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Difference in ASVs", title = "Delta richness", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/alphadiversity/richness_ethnicities.pdf", width = 6, height = 5)

comp <- list(c("Moroccan", "South-Asian Surinamese"), c("Dutch", "South-Asian Surinamese"), 
             c("Dutch", "African Surinamese"), c("Dutch", "Turkish"))
ggplot(data = dftot %>% filter(!is.na(PD_delta)), 
       aes(x = fct_reorder(EthnicityTot, .x = PD_delta, .fun = median), y = PD_delta)) +
    geom_violin(aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Difference in Faith's PD", title = "Delta Faith's PD", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/alphadiversity/faith_ethnicities.pdf", width = 6, height = 5)

#### Delta alpha diversity and outcomes ####
ggplot(data = dftot %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = shannon_delta)) +
    geom_violin(aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Delta Shannon index", x= "Hypertension", title = "Hypertension") +
    # stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
    #                    label = "p.signif") +
    theme_Publication()
ggsave("results/alphadiversity/shannon_hypertension.pdf", width = 4, height = 5)

ggplot(data = dftot %>% filter(!is.na(DM)), aes(x = DM, y = shannon_delta)) +
    geom_violin(aes(fill = DM)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Delta Shannon index", x= "Diabetes", title = "Diabetes") +
    # stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
    #                    label = "p.signif") +
    theme_Publication()
ggsave("results/alphadiversity/shannon_dm.pdf", width = 4, height = 5)

ggplot(data = dftot %>% filter(!is.na(MetSyn)), aes(x = MetSyn, y = shannon_delta)) +
    geom_violin(aes(fill = MetSyn)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Delta Shannon index", x= "Metabolic syndrome", title = "Metabolic syndrome") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication()
ggsave("results/alphadiversity/shannon_metsyn.pdf", width = 4, height = 5)

ggplot(data = dftot %>% filter(!is.na(HbA1c_delta)), aes(x = shannon_delta, y = HbA1c_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta HbA1c", x= "Delta Shannon index", title = "Shannon and HbA1c change") +
    stat_cor() +
    theme_Publication()
ggsave("results/alphadiversity/shannon_deltahba1c.pdf", width = 7, height = 7)

ggplot(data = dftot %>% filter(!is.na(LDL_delta)), aes(x = shannon_delta, y = LDL_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta LDL", x= "Delta Shannon index", title = "Shannon and LDL change") +
    stat_cor() +
    theme_Publication()
ggsave("results/alphadiversity/shannon_deltaldl.pdf", width = 7, height = 7)

ggplot(data = dftot %>% filter(!is.na(Age_delta)), aes(x = shannon_delta, y = Age_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta Age", x= "Delta Shannon index", title = "Shannon and age change") +
    stat_cor() +
    theme_Publication()
ggsave("results/alphadiversity/shannon_deltaage.pdf", width = 7, height = 7)

ggplot(data = dftot %>% filter(!is.na(BMI_delta)), aes(x = shannon_delta, y = BMI_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta BMI", x= "Delta Shannon index", title = "Shannon and BMI change") +
    stat_cor() +
    theme_Publication()
ggsave("results/alphadiversity/shannon_deltabmi.pdf", width = 7, height = 7)

## Linear mixed ##
dftot2$DM_new[which(dftot2$timepoint == "follow-up")] <- dftot2$DM_new[which(dftot2$ID == dftot2$ID & dftot2$timepoint == "baseline")]
dftot2$HT_new[which(dftot2$timepoint == "follow-up")] <- dftot2$HT_new[which(dftot2$ID == dftot2$ID & 
                                                                                 dftot2$timepoint == "baseline")]
dftot2$MetSyn_new[which(dftot2$timepoint == "follow-up")] <- dftot2$MetSyn_new[which(dftot2$ID == dftot2$ID & 
                                                                                         dftot2$timepoint == "baseline")]
df_means_dm <- dftot2 %>% filter(!is.na(DM_new)) %>% 
    select(ID, DM_new, timepoint, shannon, richness, PD) %>% 
    group_by(timepoint, DM_new) %>% 
    summarise(across(c(shannon, richness, PD),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ),
                     .names = "{.col}_{.fn}"))
df_means_ht <- dftot2 %>% filter(!is.na(HT_new)) %>% 
    select(ID, HT_new, timepoint, shannon, richness, PD) %>% 
    group_by(timepoint, HT_new) %>% 
    summarise(across(c(shannon, richness, PD),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ),
                     .names = "{.col}_{.fn}"))
df_means_ms <- dftot2 %>% filter(!is.na(MetSyn_new)) %>% 
    select(ID, MetSyn_new, timepoint, shannon, richness, PD) %>% 
    group_by(timepoint, MetSyn_new) %>% 
    summarise(across(c(shannon, richness, PD),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ),
                     .names = "{.col}_{.fn}"))

shannon_dm <- linearmixed(dftot2, shannon, DM_new)
richness_dm <- linearmixed(dftot2, richness, DM_new)
pd_dm <- linearmixed(dftot2, PD, DM_new)

shannon_ht <- linearmixed(dftot2, shannon, HT_new)
richness_ht <- linearmixed(dftot2, richness, HT_new)
pd_ht <- linearmixed(dftot2, PD, HT_new)

shannon_ms <- linearmixed(dftot2, shannon, MetSyn_new)
richness_ms <- linearmixed(dftot2, richness, MetSyn_new)
pd_ms <- linearmixed(dftot2, PD, MetSyn_new)

(pl1 <- dftot2 %>% filter(!is.na(DM_new)) %>%  ggplot() +
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
        labs(x = "Timepoint", y = "Shannon index", title = "New diabetes diagnoses",
             color = ""))
ggsave("results/alphadiversity/diabetes_deltashannon.pdf", width = 4, height = 4)

(pl1 <- dftot2 %>% filter(!is.na(DM_new)) %>%  ggplot() +
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
        labs(x = "Timepoint", y = "Richness", title = "New diabetes diagnoses",
             color = ""))
ggsave("results/alphadiversity/diabetes_deltarichness.pdf", width = 4, height = 4)

(pl2 <- dftot2 %>% filter(!is.na(HT_new)) %>%  ggplot() +
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
        labs(x = "Timepoint", y = "Shannon index", title = "New hypertension diagnoses",
             color = ""))
ggsave("results/alphadiversity/hypertension_deltashannon.pdf", width = 4, height = 4)

(pl2 <- dftot2 %>% filter(!is.na(MetSyn_new)) %>%  ggplot() +
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
        labs(x = "Timepoint", y = "Shannon index", title = "New MetSyn diagnoses",
             color = ""))
ggsave("results/alphadiversity/metsyn_deltashannon.pdf", width = 4, height = 4)

(pl2 <- dftot2 %>%  ggplot() +
        geom_line(aes(x = timepoint, y = shannon, group = ID), alpha = 0.1, color = "grey40") +
        geom_jitter(aes(x = timepoint, y = shannon, group = ID, color = EthnicityTot), 
                    alpha = 0.5, width = 0) +
        gghalves::geom_half_violin(aes(x = timepoint, y = shannon, fill = EthnicityTot),
                                   side = c(rep("l", nlevels(dftot2$EthnicityTot)),
                                                rep("r", nlevels(dftot2$EthnicityTot))), nudge = 0.05) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "baseline"), 
                                    aes(x = timepoint, y = shannon), 
                                    nudge = 0.05, side = "l", fill = "white", color = "black", 
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "follow-up"), 
                                    aes(x = timepoint, y = shannon), 
                                    nudge = 0.05, side = "r", fill = "white", color = "black", 
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        # stat_pvalue_manual(shannon_ht, label = "{p.signif}", y.position = max(dftot2$shannon), remove.bracket = FALSE) +
        # stat_compare_means(aes(x = timepoint, y = shannon)) +
        scale_color_simpsons(guide = "none") +
        scale_fill_simpsons(guide = "none") +
        facet_wrap(~EthnicityTot) +
        theme_Publication() +
        theme(strip.text = element_text(size = rel(0.6))) +
        labs(x = "Timepoint", y = "Shannon index", title = "Ethnicity and delta Shannon",
             color = ""))
ggsave("results/alphadiversity/ethnicity_deltashannon.pdf", width = 6, height = 7)

(pl2 <- dftot2 %>%  ggplot() +
        geom_line(aes(x = timepoint, y = richness, group = ID), alpha = 0.1, color = "grey40") +
        geom_jitter(aes(x = timepoint, y = richness, group = ID, color = EthnicityTot), 
                    alpha = 0.5, width = 0) +
        gghalves::geom_half_violin(aes(x = timepoint, y = richness, fill = EthnicityTot),
                                   side = c(rep("l", nlevels(dftot2$EthnicityTot)),
                                            rep("r", nlevels(dftot2$EthnicityTot))), nudge = 0.05) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "baseline"), 
                                    aes(x = timepoint, y = richness), 
                                    nudge = 0.05, side = "l", fill = "white", color = "black", 
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        gghalves::geom_half_boxplot(data = . %>% filter(timepoint == "follow-up"), 
                                    aes(x = timepoint, y = richness), 
                                    nudge = 0.05, side = "r", fill = "white", color = "black", 
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        # stat_pvalue_manual(shannon_ht, label = "{p.signif}", y.position = max(dftot2$shannon), remove.bracket = FALSE) +
        # stat_compare_means(aes(x = timepoint, y = shannon)) +
        scale_color_simpsons(guide = "none") +
        scale_fill_simpsons(guide = "none") +
        facet_wrap(~EthnicityTot) +
        theme_Publication() +
        theme(strip.text = element_text(size = rel(0.6))) +
        labs(x = "Timepoint", y = "Richness", title = "Ethnicity and delta richness",
             color = ""))
ggsave("results/alphadiversity/ethnicity_deltarichness.pdf", width = 6, height = 7)

#### Correlations alpha diversity and sample dissimilarity over FU ####
ggplot(data = dftot3 %>% filter(timepoint == "baseline"), aes(x = shannon, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Shannon index (baseline)", 
         title = "Alpha diversity and sample dissimilarity") +
    stat_cor() +
    theme_Publication()
ggsave("results/alphadiversity/braycurtis_shannonbaseline.pdf", width = 5, height = 5)

ggplot(data = dftot3 %>% filter(timepoint == "baseline"), aes(x = shannon_delta, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "loess") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Shannon index change", 
         title = "Shannon change and sample dissimilarity") +
    stat_cor() +
    theme_Publication()
ggsave("results/alphadiversity/braycurtis_shannonchange.pdf", width = 5, height = 5)

