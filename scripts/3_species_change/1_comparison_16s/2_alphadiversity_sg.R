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
library(rstatix)

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

#### Shotgun data ####
## Load data
df_new <- readRDS("data/clinicaldata/clinicaldata_long.RDS")

## Diversity metrics
# Shannon plots
sg <- readRDS("data/shotgun/shotgun_abundance.RDS")
sgtab <- as.data.frame(sg)
shannonsg <- vegan::diversity(sgtab, index = "shannon")
df_shansg <- data.frame(sampleID = names(shannonsg), shannon = shannonsg)
df_shansg$ID <- str_remove(str_remove(df_shansg$sampleID, "HELIBA_"), "HELIFU_")
df_shansg <- left_join(df_shansg, df_new)

## Species richness
specrichsg <- specnumber(sgtab)
dfspecsg <- data.frame(sampleID = names(specrichsg), richness = specrichsg)
dfspecsg$ID <- str_remove(str_remove(dfspecsg$sampleID, "HELIBA_"), "HELIFU_")
dfspecsg <- left_join(dfspecsg, df_shansg)

saveRDS(dfspecsg, "data/shotgun/clin_alphadiv_sg.RDS")

#### Load data ####
df <- readRDS("data/shotgun/clin_alphadiv_sg.RDS") %>% 
    dplyr::select(sampleID, ID, richness, shannon) %>% 
    mutate(timepoint = case_when(
                        str_detect(sampleID, "HELIBA") ~ "baseline",
                        str_detect(sampleID, "HELIFU") ~ "follow-up"
                    ),
           ID = str_c("S", ID))

dfwide <- df %>% pivot_wider(., id_cols = "ID", names_from = "timepoint",
                             values_from = c(1,3:4)) %>% 
    mutate(shannon_delta = `shannon_follow-up` - shannon_baseline,
           richness_delta = `richness_follow-up` - richness_baseline) %>% 
    dplyr::select(1, 8:9)
df2 <- left_join(df, dfwide, by = "ID") 

helius <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
dftot <- left_join(df2 %>% filter(timepoint == "baseline"), helius %>% filter(timepoint == "baseline"), 
                   by = c("ID", "timepoint", "sampleID"))

pairedids <- df2$sampleID[which(!is.na(df2$shannon_delta))]
dftot2 <- left_join(df, helius, by = c("ID", "timepoint", "sampleID")) %>% 
            filter(sampleID %in% pairedids) %>% 
            filter(! EthnicityTot %in% "Other" ) %>% droplevels(.)

betadiv <- readRDS("data/shotgun/braydistance_delta.RDS") %>% dplyr::select(1:4)
dftot3a <- left_join(dftot2 %>% filter(timepoint == "baseline"), 
                     betadiv %>% dplyr::select(sampleID = sampleID1, everything(.)))
dftot3b <- left_join(dftot2 %>% filter(timepoint == "follow-up"), 
                     betadiv %>% dplyr::select(sampleID = sampleID2, everything(.)))
dftot3 <- full_join(dftot3a, dftot3b) %>% droplevels(.) %>% 
    full_join(., df2)

saveRDS(dftot3, "data/shotgun/alphabetadiversity_shotgun.RDS")
df <- readRDS("data/shotgun/alphabetadiversity_shotgun.RDS")

#### Output folder ####
resultsfolder <- "results/3_species_change/1_comparison_16s/alphadiversity"
dir.create(resultsfolder, showWarnings = FALSE, recursive = TRUE)

#### Baseline and FU alpha diversity ####
comp <- list(c("Dutch", "South-Asian Surinamese"))
ggplot(data = dftot3 %>% filter(!is.na(EthnicityTot)), 
       aes(x = fct_reorder(EthnicityTot, .x = shannon, .fun = median), y = shannon)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Shannon index", title = "Shannon index", x = "") +
    facet_wrap(~timepoint) +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "wilcox.test") +
    theme_Publication()
ggsave("results/3_species_change/1_comparison_16s/alphadiversity/sg_shannon_time_ethnicities.pdf", width = 6, height = 5)

comp <- list(c("Dutch", "South-Asian Surinamese"))
ggplot(data = dftot3 %>% filter(!is.na(EthnicityTot)),
       aes(x = fct_reorder(EthnicityTot, .x = shannon, .fun = median), y = shannon)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Shannon index", title = "Shannon index", x = "") +
    facet_wrap(~timepoint) +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "wilcox.test") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/3_species_change/1_comparison_16s/alphadiversity/sg_shannon_ethnicities_time.pdf", width = 8, height = 5)

#### Delta alpha diversity ####
comp <- list(c("Dutch", "South-Asian Surinamese"))
ggplot(data = dftot %>% filter(!is.na(shannon_delta)), # selecting baseline samples
       aes(x = fct_reorder(EthnicityTot, .x = shannon_delta, .fun = median), y = shannon_delta)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Shannon index difference", title = "Delta Shannon", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication()
ggsave("results/3_species_change/1_comparison_16s/alphadiversity/sg_shannon_ethnicities.pdf", width = 4.5, height = 5)

comp <- list(c("Dutch", "South-Asian Surinamese"))
ggplot(data = dftot %>% filter(!is.na(richness_delta)), # selecting baseline samples (those have deltas available)
       aes(x = fct_reorder(EthnicityTot, .x = richness_delta, .fun = median), y = richness_delta)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Difference in ASVs", title = "Delta richness", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication()
ggsave("results/3_species_change/1_comparison_16s/alphadiversity/sg_richness_ethnicities.pdf", width = 4.5, height = 5)

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
        stat_compare_means(aes(x = timepoint, y = shannon), hide.ns = TRUE, label.x = 1.5, size = 4,
                           label = "p.signif", tip.length = 0) +
        scale_color_simpsons(guide = "none") +
        scale_fill_simpsons(guide = "none") +
        facet_wrap(~EthnicityTot) +
        theme_Publication() +
        theme(strip.text = element_text(size = rel(0.6))) +
        labs(x = "Timepoint", y = "Shannon index", title = "Ethnicity and Shannon",
             color = ""))
ggsave("results/3_species_change/1_comparison_16s/alphadiversity/sg_ethnicity_deltashannon.pdf", width = 6, height = 5)

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
        stat_compare_means(aes(x = timepoint, y = shannon), label = "p.signif", hide.ns = TRUE, label.x = 1.5,
                           size = 4) +
        scale_color_simpsons(guide = "none") +
        scale_fill_simpsons(guide = "none") +
        facet_wrap(~EthnicityTot) +
        theme_Publication() +
        theme(strip.text = element_text(size = rel(0.6))) +
        labs(x = "Timepoint", y = "Richness", title = "Ethnicity and richness",
             color = ""))
ggsave("results/3_species_change/1_comparison_16s/alphadiversity/sg_ethnicity_deltarichness.pdf", width = 6, height = 5)

#### Correlations alpha diversity and sample dissimilarity over FU ####
ggplot(data = dftot3 %>% filter(timepoint == "baseline"), aes(x = shannon, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Shannon index (baseline)", 
         title = "Alpha diversity and sample dissimilarity") +
    stat_cor() +
    theme_Publication()
ggsave("results/3_species_change/1_comparison_16s/alphadiversity/sg_braycurtis_shannonbaseline.pdf", width = 5, height = 5)

ggplot(data = dftot3 %>% filter(timepoint == "baseline"), aes(x = shannon_delta, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "loess", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Shannon index change", 
         title = "Shannon change and sample dissimilarity") +
    # stat_cor() +
    theme_Publication()
ggsave("results/3_species_change/1_comparison_16s/alphadiversity/sg_braycurtis_shannonchange.pdf", width = 5, height = 5)

#### Figure panels — Supplementary Figure 6 ####
eth_colors <- c(
    "Dutch"                  = "#709AE1FF",
    "South-Asian Surinamese" = "#FED439FF"
)

fmt_pval <- function(p) {
    ifelse(p < 0.0001, "p < 0.0001", paste0("p = ", formatC(p, format = "f", digits = 3)))
}

## Panel E — Shannon by ethnicity × timepoint
dftot3_fig <- dftot3 %>%
    filter(!is.na(EthnicityTot), EthnicityTot != "Other") %>%
    droplevels()

shannon_max_sg <- dftot3_fig %>%
    group_by(EthnicityTot) %>%
    summarize(y_max = max(shannon, na.rm = TRUE), .groups = "drop")

shannon_pvals_sg <- dftot3_fig %>%
    group_by(EthnicityTot) %>%
    rstatix::wilcox_test(shannon ~ timepoint) %>%
    ungroup() %>%
    filter(p < 0.05) %>%
    mutate(label = fmt_pval(p)) %>%
    left_join(shannon_max_sg, by = "EthnicityTot") %>%
    mutate(y.position = y_max * 1.08, xmin = 1, xmax = 2)

pl_sfig6_E <- ggplot(dftot3_fig, aes(x = timepoint, y = shannon)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot, alpha = timepoint)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
    labs(y = "Shannon index", title = "Shannon index", x = "") +
    facet_wrap(~EthnicityTot) +
    stat_pvalue_manual(shannon_pvals_sg, label = "label", tip.length = 0, label.size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_Publication() +
    theme(strip.text = element_text(size = 8))

## Panel F — Baseline Shannon vs Bray-Curtis dissimilarity
pl_sfig6_F <- dftot3_fig %>%
    filter(timepoint == "baseline", !is.na(distance)) %>%
    ggplot(aes(x = shannon, y = distance)) +
    geom_point(alpha = 0.35, size = 1.5, color = "#197EC0FF") +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.9,
                color = "black") +
    stat_cor(size = 5) +
    labs(x = "Baseline Shannon diversity",
         y = "Bray-Curtis dissimilarity over FU time",
         title = "Diversity-stability relationship") +
    theme_Publication()

