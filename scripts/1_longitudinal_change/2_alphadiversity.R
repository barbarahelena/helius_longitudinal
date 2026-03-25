## Calculate alpha diversity differences

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(gghalves)

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

#### Output folder ####
resultsfolder <- "results/1_longitudinal_change/alphadiversity"
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

# #### 16S data ####
# phydata <- readRDS("data/16s/phyloseq_withclinpaired.RDS")
# df_new <- rio:: import("data/clinicaldata_long.RDS") %>% mutate(ID = sampleID)
# tab <- as.data.frame(t(as(phydata@otu_table, 'matrix')))
# tab_matrix <- t(as(phydata@otu_table, 'matrix'))

# ## Diversity metrics
# # Shannon plots
# shannon <- vegan::diversity(tab, index = 'shannon')
# df_shan <- data.frame(ID = names(shannon), shannon = shannon)
# df_shan <- left_join(df_shan, df_new, by = "ID")

# ## Species richness
# specrich <- specnumber(tab)
# dfspec <- data.frame(ID = names(specrich), richness = specrich)
# dfspec <- left_join(dfspec, df_shan, by = "ID")

# ## Faith's PD
# faith <- picante::pd(samp = tab_matrix, tree = phydata@phy_tree)
# dffai <- as.data.frame(faith)
# dffai$ID <- rownames(faith)
# dffai <- left_join(dffai, dfspec, by = "ID")
# saveRDS(dffai, "data/16s/clin_alphadiversity.RDS")

#### Load data ####
df_raw <- readRDS("data/16s/clin_alphadiversity.RDS")

df <- df_raw %>%
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

betadiv <- readRDS("data/16s/archive/braydistance_delta.RDS") %>% dplyr::select(1:4)
dftot3a <- left_join(dftot2 %>% filter(timepoint == "baseline"), 
                     betadiv %>% dplyr::select(sampleID = sampleID1, everything(.)))
dftot3b <- left_join(dftot2 %>% filter(timepoint == "follow-up"), 
                     betadiv %>% dplyr::select(sampleID = sampleID2, everything(.)))
dftot3 <- full_join(dftot3a, dftot3b) %>% droplevels(.) %>% 
    full_join(., df2)

#### Overview: baseline vs follow-up ####
(plshan <- ggplot(data = df_raw, aes(x = timepoint, y = shannon, fill = timepoint)) +
    geom_violin(colour = NA) +
    scale_fill_manual(values = rev(pal_simpsons()(7)[c(1,7)]), guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    stat_compare_means(label.y = 5.5) +
    labs(title = "Shannon index", y = "Shannon index", x="") +
    theme_Publication())
# ggsave(plshan, filename = file.path(resultsfolder, "shannon.svg"), width = 4, height = 5)
ggsave(plshan, filename = file.path(resultsfolder, "shannon.pdf"), width = 4, height = 5)

(plrich <- ggplot(data = df_raw, aes(x = timepoint, y = SR, fill = timepoint)) +
    geom_violin(colour = NA) +
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() +
    scale_fill_manual(values = rev(pal_simpsons()(7)[c(1,7)]), guide = "none") +
    labs(title = "Species richness", y = "Number of species", x = "") +
    stat_compare_means(method = "wilcox.test"))
ggsave(plrich, filename = file.path(resultsfolder, "richness.pdf"), width = 4, height = 5)
# ggsave(plrich, filename = file.path(resultsfolder, "richness.svg"), width = 4, height = 5)

(plfaith <- ggplot(data = df_raw, aes(x = timepoint, y = PD, fill = timepoint)) +
    geom_violin(colour = NA) +
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() +
    scale_fill_manual(values = rev(pal_simpsons()(7)[c(1,7)]), guide = "none") +
    labs(title = "Faith's PD", y = "Faith's phylogenetic diversity", x = "") +
    stat_compare_means(method = "wilcox.test"))
ggsave(plfaith, filename = file.path(resultsfolder, "faiths.pdf"), width = 4, height = 5)
# ggsave(plfaith, filename = file.path(resultsfolder, "faiths.svg"), width = 4, height = 5)

pl_total <- ggarrange(plshan, plrich, plfaith, labels = c("A", "B", "C"), nrow = 1)
ggsave(pl_total, filename = file.path(resultsfolder, "alphadivplots.pdf"), width = 11, height = 5.5)
# ggsave(pl_total, filename = file.path(resultsfolder, "alphadivplots.svg"), width = 11, height = 5.5)

#### Baseline and FU alpha diversity ####
comp <- rev(list(c("Dutch", "South-Asian Surinamese"), c("Dutch", "Ghanaian"), 
                 c("Dutch", "African Surinamese"),c("Dutch", "Turkish"), c("Dutch", "Moroccan")))
ggplot(data = dftot3 %>% filter(!is.na(EthnicityTot)), 
       aes(x = fct_reorder(EthnicityTot, .x = shannon, .fun = median), y = shannon)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Shannon index", title = "Shannon index", x = "") +
    facet_wrap(~timepoint) +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "wilcox.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/1_longitudinal_change/alphadiversity/shannon_time_ethnicities.pdf", width = 8, height = 7)

ggplot(data = dftot3 %>% filter(!is.na(EthnicityTot)), 
       aes(x = timepoint, y = shannon)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot, alpha = timepoint)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
    labs(y = "Shannon index", title = "Shannon index", x = "") +
    facet_wrap(~EthnicityTot) +
    stat_compare_means(tip.length = 0, hide.ns = TRUE, label.x = 1.5,
                       label = "p.signif", method = "wilcox.test") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_Publication()
ggsave("results/1_longitudinal_change/alphadiversity/shannon_ethnicities_time.pdf", width = 7, height = 10)
pl_fig1_D <- last_plot()

#### Delta alpha diversity ####
comp <- rev(list(c("Dutch", "South-Asian Surinamese"), c("Dutch", "Moroccan"), c("Dutch", "Turkish")))
ggplot(data = dftot %>% filter(!is.na(shannon_delta)), # selecting baseline samples
       aes(x = fct_reorder(EthnicityTot, .x = shannon_delta, .fun = median), y = shannon_delta)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Shannon index difference", title = "Delta Shannon", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/1_longitudinal_change/alphadiversity/shannon_ethnicities.pdf", width = 6, height = 5)

comp <- list(c("Dutch", "South-Asian Surinamese"), c("Dutch", "African Surinamese"), c("Dutch", "Turkish"))
ggplot(data = dftot %>% filter(!is.na(richness_delta)), # selecting baseline samples (those have deltas available)
       aes(x = fct_reorder(EthnicityTot, .x = richness_delta, .fun = median), y = richness_delta)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Difference in ASVs", title = "Delta richness", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/1_longitudinal_change/alphadiversity/richness_ethnicities.pdf", width = 6, height = 5)

comp <- list(c("Moroccan", "South-Asian Surinamese"), c("Dutch", "South-Asian Surinamese"), 
             c("Dutch", "African Surinamese"), c("Dutch", "Turkish"))
ggplot(data = dftot %>% filter(!is.na(PD_delta)), 
       aes(x = fct_reorder(EthnicityTot, .x = PD_delta, .fun = median), y = PD_delta)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Difference in Faith's PD", title = "Delta Faith's PD", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/1_longitudinal_change/alphadiversity/faith_ethnicities.pdf", width = 6, height = 5)

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
        scale_color_manual(values = eth_colors, guide = "none") +
        scale_fill_manual(values = eth_colors, guide = "none") +
        facet_wrap(~EthnicityTot) +
        theme_Publication() +
        theme(strip.text = element_text(size = rel(0.6))) +
        labs(x = "Timepoint", y = "Shannon index", title = "Ethnicity and delta Shannon",
             color = ""))
ggsave("results/1_longitudinal_change/alphadiversity/ethnicity_deltashannon.pdf", width = 6, height = 7)

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
        scale_color_manual(values = eth_colors, guide = "none") +
        scale_fill_manual(values = eth_colors, guide = "none") +
        facet_wrap(~EthnicityTot) +
        theme_Publication() +
        theme(strip.text = element_text(size = rel(0.6))) +
        labs(x = "Timepoint", y = "Richness", title = "Ethnicity and delta richness",
             color = ""))
ggsave("results/1_longitudinal_change/alphadiversity/ethnicity_deltarichness.pdf", width = 6, height = 7)

#### Correlations alpha diversity and sample dissimilarity over FU ####

# Supplementary: all ethnicities
ggplot(data = dftot3 %>% filter(timepoint == "baseline"), aes(x = shannon, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Shannon index (baseline)",
         title = "Alpha diversity and sample dissimilarity") +
    stat_cor() +
    theme_Publication()
ggsave("results/1_longitudinal_change/alphadiversity/braycurtis_shannonbaseline.pdf", width = 5, height = 5)

# Figure 1 panel: all ethnicities, single regression
(pl_fig1_scatter <- dftot3 %>%
    filter(timepoint == "baseline",
           !is.na(distance),
           EthnicityTot != "Other") %>%
    ggplot(aes(x = shannon, y = distance)) +
    geom_point(alpha = 0.35, size = 1.5, color = "royalblue") +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.9, color = "black") +
    stat_cor(size = 3) +
    labs(x = "Baseline Shannon diversity",
         y = "Bray-Curtis dissimilarity over FU time",
         title = "Diversity-stability relationship") +
    theme_Publication())
ggsave("results/1_longitudinal_change/alphadiversity/braycurtis_shannonbaseline.pdf", width = 5, height = 5)

ggplot(data = dftot3 %>% filter(timepoint == "baseline"), aes(x = shannon_delta, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "loess") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "Shannon index change", 
         title = "Shannon change and sample dissimilarity") +
    stat_cor() +
    theme_Publication()
ggsave("results/1_longitudinal_change/alphadiversity/braycurtis_shannonchange.pdf", width = 5, height = 5)

