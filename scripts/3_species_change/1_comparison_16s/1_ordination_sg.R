## Calculate distances, plot PCoA and PCA

## Libraries
library(tidyverse)
library(vegan)
library(permute)
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

#### Output folder ####
resultsfolder <- "results/3_species_change/1_comparison_16s/ordination"
dir.create(resultsfolder, showWarnings = FALSE, recursive = TRUE)

#### Calculate distances ####
df_new <- readRDS("data/clinicaldata/clinicaldata_long.RDS")

expl_var_csv <- file.path(resultsfolder, "expl_var_bray_shotgun.csv")
if (!file.exists("data/shotgun/bray_shotgun.RDS") ||
    !file.exists("data/shotgun/clin_betadiversity_shotgun.RDS") ||
    !file.exists(expl_var_csv)) {
    shotdata <- readRDS("data/shotgun/shotgun_abundance.RDS")
    shotmat <- as.matrix(shotdata)

    print('Bray-Curtis distance total dataset')
    bray <- vegan::vegdist(shotmat, method = 'bray')
    saveRDS(bray, "data/shotgun/bray_shotgun.RDS")

    pcoord <- ape::pcoa(bray, correction = "cailliez")
    expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
    write_lines(as.character(expl_variance_bray), expl_var_csv)

    dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
    dbray <- as.data.frame(dbray)
    dbray$ID <- rownames(dbray)
    dbray <- left_join(dbray, df_new, by = 'ID') %>%
        dplyr::select(BrayPCo1 = `Axis.1`, BrayPCo2 = `Axis.2`, everything(.))
    saveRDS(dbray, "data/shotgun/clin_betadiversity_shotgun.RDS")
} else {
    print('Loading precomputed Bray-Curtis distance')
    bray <- readRDS("data/shotgun/bray_shotgun.RDS")
    expl_variance_bray <- as.numeric(readLines(expl_var_csv))
    dbray <- readRDS("data/shotgun/clin_betadiversity_shotgun.RDS")
}

#### Load plotting data ####
helius <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
df <- dbray %>% dplyr::select(1:2, sampleID = ID)
df <- left_join(df, helius, by = c("sampleID"))
dim(df)

#### Colour palette — timepoints ####
tp_colors <- setNames(pal_simpsons()(2), c("baseline", "follow-up"))

#### Colour palette — ethnicities ####
eth_colors <- c(
    "Dutch"                  = "#709AE1FF",
    "South-Asian Surinamese" = "#FED439FF"
)

#### PCoA per timepoint — ethnicity comparison ####
set.seed(1234)
permanova_tp <- setNames(lapply(c("baseline", "follow-up"), function(tp) {
    sub <- df %>% filter(timepoint == tp, !is.na(EthnicityTot), EthnicityTot != "Other")
    ids <- sub$sampleID
    braymat_tp <- as.matrix(bray)[ids, ids]
    if (nrow(sub) < 10 || length(unique(sub$EthnicityTot)) < 2) return(NULL)
    adonis2(as.dist(braymat_tp) ~ EthnicityTot, data = sub, by = "terms")
}), c("baseline", "follow-up"))
print(permanova_tp)

label_df_tp <- data.frame(
    timepoint = factor(c("baseline", "follow-up"), levels = c("baseline", "follow-up")),
    label = sapply(c("baseline", "follow-up"), function(tp) {
        res <- permanova_tp[[tp]]
        if (is.null(res)) return("")
        str_c("p = ", res$`Pr(>F)`[1], ", r² = ", format(round(res$R2[1], 3), nsmall = 3))
    })
)

df_tp_pcoa <- df %>%
    filter(!is.na(EthnicityTot), EthnicityTot != "Other") %>%
    mutate(timepoint = factor(timepoint, levels = c("baseline", "follow-up")))

(pl_pcoa_tp <- ggplot(df_tp_pcoa, aes(BrayPCo1, BrayPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = EthnicityTot, fill = EthnicityTot),
                 type = "norm", alpha = 0.1) +
    geom_point(aes(color = EthnicityTot), size = 0.8, alpha = 0.5) +
    geom_text(
        data = label_df_tp,
        aes(label = label), x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
        size = 2.5, color = "grey30", inherit.aes = FALSE
    ) +
    scale_color_manual(values = eth_colors, name = NULL) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    facet_wrap(~timepoint, scales = "free") +
    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
    labs(color = "", title = "PCoA Bray-Curtis — by timepoint") +
    theme_Publication() +
    theme(legend.position = "bottom"))
ggsave(pl_pcoa_tp, filename = file.path(resultsfolder, "PCoA_BrayCurtis_sg_timepoint.pdf"),
       device = "pdf", width = 10, height = 6)

## Distance between datapoints

braymat <- as.matrix(bray)
all_combinations <- t(combn(unique(rownames(braymat)), 2, simplify = TRUE))
data_long <- data.frame(
    sampleID1 = all_combinations[, 1],
    sampleID2 = all_combinations[, 2]
)
data_long <- data_long %>%
    filter(str_remove(sampleID1, "HELIBA_") == str_remove(sampleID2, "HELIFU_")) %>%
    mutate(ID = str_c("S", str_remove(sampleID1, "HELIBA_")))
for(a in 1:nrow(data_long)){
    distbray = braymat[paste0(data_long$sampleID1[a]), paste0(data_long$sampleID2[a])]
    data_long$distance[a] <- distbray
}
heliusdist <- inner_join(data_long, helius, by = "ID") %>% filter(timepoint == "baseline")
saveRDS(heliusdist, "data/shotgun/braydistance_delta.RDS")

#### Figure panels — Supplementary Figure 6 ####

## Panel A — Follow-up time by ethnicity
fu_data <- heliusdist %>%
    filter(!is.na(FUtime), !is.na(EthnicityTot), EthnicityTot != "Other") %>%
    mutate(EthnicityTot = fct_reorder(EthnicityTot, FUtime, median)) %>%
    droplevels()

pl_sfig6_A <- ggplot(fu_data, aes(x = EthnicityTot, y = FUtime)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Follow-up time (years)", title = "Follow-up time", x = "") +
    stat_compare_means(label = "p.format", tip.length = 0) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = 4:12) +
    theme_Publication() +
    coord_flip()
## Panel B — per-timepoint PCoA (computed above as pl_pcoa_tp)
pl_sfig6_B <- pl_pcoa_tp

## Panel C — Bray-Curtis dissimilarity vs follow-up time
pl_sfig6_C <- ggplot(
        heliusdist %>% filter(!is.na(FUtime), EthnicityTot != "Other"),
        aes(x = FUtime, y = distance)) +
    geom_jitter(color = "#197EC0FF", alpha = 0.3, position = position_jitter(seed = 1234, width = 0)) +
    geom_smooth(color = "black", method = "lm") +
    labs(y = "Bray-Curtis dissimilarity over FU time",
         x = "FU time (years)", title = "FU time and sample distance") +
    stat_cor() +
    theme_Publication()

## Panel D — Confounder-adjusted Bray-Curtis by ethnicity
heliusdist_adj <- heliusdist %>%
    filter(!is.na(EthnicityTot), EthnicityTot != "Other") %>%
    filter(!is.na(Age), !is.na(Sex), !is.na(BMI),
           !is.na(Metformin), !is.na(PPI), !is.na(AntiHT), !is.na(Statins), !is.na(AlcCons)) %>%
    droplevels()
lm_confounders <- lm(distance ~ Age + Sex + BMI + Metformin + PPI + FUtime,
                     data = heliusdist_adj)
heliusdist_adj <- heliusdist_adj %>%
    mutate(
        dist_adjusted = residuals(lm_confounders) + mean(distance, na.rm = TRUE),
        EthnicityTot  = fct_reorder(EthnicityTot, dist_adjusted, median)
    )
pl_sfig6_D <- ggplot(heliusdist_adj, aes(x = EthnicityTot, y = dist_adjusted)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Adjusted Bray-Curtis dissimilarity",
         title = "Distance baseline to follow-up", x = "") +
    stat_compare_means(aes(label = sprintf("p = %s", ..p.format..)),
                       tip.length = 0) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_Publication() +
    coord_flip()