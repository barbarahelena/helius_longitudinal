## Calculate distances, plot PCoA and PCA

## Libraries
library(tidyverse)
library(doParallel)
library(phyloseq)
library(mixOmics)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggsci)
registerDoParallel(parallel::detectCores() - 1)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family)
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
helius <- readRDS("data/clinicaldata_long.RDS")
df <- dbray %>% dplyr::select(1:2, sampleID = ID)
df <- left_join(df, helius, by = c("sampleID"))

#### Bray-Curtis distance ####
print('PERMANOVA..')
set.seed(1234)
dfanova <- df[match(attributes(bray)[["Labels"]], df$sampleID),]
all(dfanova$sampleID == attributes(bray)[["Labels"]]) # TRUE
dim(df)
res1 <- adonis2(bray ~ timepoint, data = df)
print(res1)

(braycurt <- df %>%
    ggplot(aes(BrayPCo1, BrayPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = timepoint, fill = timepoint), type = "norm",
                 alpha = 0.1) +
    geom_point(aes(color = timepoint), size = 1, alpha = 0.5) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
    scale_color_manual(values = pal_simpsons()(2)) +
    scale_fill_manual(values = pal_simpsons()(2), guide = "none") +
    theme_Publication() +
    labs(color = "") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res1$`Pr(>F)`, ", r2 = ",
                           format(round(res1$R2[1],3), nsmall = 3))
             ))
ggsave(braycurt, filename = file.path(resultsfolder, "PCoA_BrayCurtis_sg.pdf"), device = "pdf", width = 8, height = 8)


(dmnewbray <- df %>% filter(!is.na(DM_new)) %>% filter(timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = DM_new, fill = DM_new), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = DM_new), size = 1, alpha = 0.5) +
        ggtitle("New diabetes") +
        xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        theme_Publication() )

(htnewbray <- df %>% filter(!is.na(HT_new)) %>% filter(timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = HT_new, fill = HT_new), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = HT_new), size = 1, alpha = 0.5) +
        ggtitle("New hypertension") +
        xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        theme_Publication() )

(metsynnewbray <- df %>% filter(!is.na(MetSyn_new)) %>% filter(timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = MetSyn_new, fill = MetSyn_new), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = MetSyn_new), size = 1, alpha = 0.5) +
        ggtitle("New metabolic syndrome") +
        xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        theme_Publication() )

(lldnewbray <- df %>% filter(!is.na(LLD_new)) %>% filter(timepoint == "baseline") %>%
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = LLD_new, fill = LLD_new), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = LLD_new), size = 1, alpha = 0.5) +
        ggtitle("New lipid lowering drug use") +
        xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        theme_Publication() )

ggarrange(dmnewbray, htnewbray, metsynnewbray, lldnewbray, nrow = 1,
          labels = LETTERS[1:4])
ggsave(file.path(resultsfolder, "clinicaloutcomes_bray_sg.pdf"), width = 18, height = 5)

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

## Plots
ggplot(data = heliusdist %>% filter(!is.na(EthnicityTot)),
       aes(x = fct_reorder(EthnicityTot, .x = distance, .fun = median), y = distance)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis distance", title = "Distance baseline to follow-up", x = "") +
    stat_compare_means(tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave(file.path(resultsfolder, "sg_distance_ethnicities.pdf"), width = 6, height = 4)

ggplot(data = heliusdist %>% filter(!is.na(FUtime)),
       aes(x = fct_reorder(EthnicityTot, .x = FUtime, .fun = median), y = FUtime)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Follow-up time (years)", title = "Follow-up time", x = "") +
    stat_compare_means(tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "wilcox.test") +
    theme_Publication() +
    coord_flip()
ggsave(file.path(resultsfolder, "sg_futime_ethnicities.pdf"), width = 6, height = 4)

#### Categorical outcomes unstratified and stratified ####
ggplot(data = heliusdist %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = distance)) +
    geom_violin(colour = NA, aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension (baseline)", title = "Hypertension") +
    stat_compare_means(tip.length = 0, hide.ns = TRUE, label = "p.signif") +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_distance_hypertension.pdf"), width = 4, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = distance)) +
    geom_violin(colour = NA, aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension (baseline)", title = "Hypertension") +
    stat_compare_means(tip.length = 0, hide.ns = TRUE, label = "p.signif", method = "wilcox.test") +
    facet_wrap(~Ethnicity) +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_distance_hypertension_ethnicity.pdf"), width = 6, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(DM)), aes(x = DM, y = distance)) +
    geom_violin(colour = NA, aes(fill = DM)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes (baseline)", title = "Diabetes") +
    stat_compare_means(tip.length = 0, hide.ns = TRUE, label = "p.signif") +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_distance_dm.pdf"), width = 4, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = distance)) +
    geom_violin(colour = NA, aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes (baseline)", title = "Diabetes") +
    stat_compare_means(tip.length = 0, hide.ns = TRUE, label = "p.signif", method = "wilcox.test") +
    facet_wrap(~Ethnicity) +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_distance_diabetes_ethnicity.pdf"), width = 6, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(MetSyn)), aes(x = MetSyn, y = distance)) +
    geom_violin(colour = NA, aes(fill = MetSyn)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Metabolic syndrome (baseline)", title = "Metabolic syndrome") +
    stat_compare_means(tip.length = 0, hide.ns = TRUE, label = "p.signif") +
    theme_Publication()
ggsave(file.path(resultsfolder, "distance_metsyn.pdf"), width = 4, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(MetSyn)), aes(x = MetSyn, y = distance)) +
    geom_violin(colour = NA, aes(fill = MetSyn)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Metabolic syndrome (baseline)", title = "Metabolic syndrome") +
    stat_compare_means(tip.length = 0, hide.ns = TRUE, label = "p.signif", method = "wilcox.test") +
    facet_wrap(~Ethnicity) +
    theme_Publication()
ggsave(file.path(resultsfolder, "distance_metsyn_ethnicity.pdf"), width = 6, height = 5)

#### Continuous outcomes unstratified and stratified ####
ggplot(data = heliusdist %>% filter(!is.na(HbA1c_delta)), aes(x = distance, y = HbA1c_delta)) +
    geom_jitter(color = "royalblue", alpha = 0.5) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    labs(y = "Delta HbA1c", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and HbA1c change") +
    stat_cor() +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_braycurtis_deltahba1c.pdf"), width = 4.5, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(HbA1c_delta)), aes(x = distance, y = HbA1c_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta HbA1c", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and HbA1c change") +
    stat_cor() +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_braycurtis_deltahba1c_ethnicity.pdf"), width = 6, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(LDL_delta)), aes(x = distance, y = LDL_delta)) +
    geom_jitter(color = "royalblue", alpha = 0.5) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    labs(y = "Delta LDL", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and LDL change") +
    stat_cor(method = "spearman") +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_braycurtis_deltaldl.pdf"), width = 4.5, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(LDL_delta)), aes(x = distance, y = LDL_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta LDL", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and LDL change") +
    stat_cor(method = "spearman") +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_braycurtis_deltaldl_ethnicity.pdf"), width = 6, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(Age)), aes(x = Age, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Baseline age", title = "Bray-Curtis and baseline age") +
    stat_cor(method = "spearman") +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_braycurtis_baselineage.pdf"), width = 4.5, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(Age)), aes(x = Age, y = distance)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5, height = 0) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Baseline age", title = "Bray-Curtis and baseline age") +
    stat_cor(method = "spearman") +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_braycurtis_baseage_ethnicity.pdf"), width = 6, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(BMI_delta)), aes(x = BMI_delta, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.5, height = 0) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    labs(y = "Delta BMI", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and BMI change") +
    stat_cor() +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_braycurtis_deltabmi.pdf"), width = 4.5, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(BMI_delta)), aes(x = distance, y = BMI_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta BMI", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and BMI change") +
    stat_cor() +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_braycurtis_deltabmi_ethnicity.pdf"), width = 6, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(FUtime)), aes(x = FUtime, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "FU time (years)", title = "FU time and sample distance") +
    stat_cor() +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_braycurtis_futime.pdf"), width = 4.5, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(FUtime)), aes(x = FUtime, y = distance)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm", formula = y ~ x) +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "FU time (years)", title = "FU time and sample distance") +
    stat_cor(method = "pearson") +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_braycurtis_futime_ethnicity.pdf"), width = 6, height = 5)

#### New diagnoses ####
heliusdist %>% filter(!is.na(HT_BPMed)) %>% group_by(EthnicityTot, HT_BPMed) %>% summarise(count = length(HT_BPMed), .groups = "drop_last")
ggplot(data = heliusdist %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = distance)) +
    geom_violin(colour = NA, aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension", title = "Hypertension") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.format") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_distance_ethnicity_hypertension.pdf"), width = 6, height = 5)

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
ggsave(file.path(resultsfolder, "sg_distance_ethnictiy_diabetes.pdf"), width = 6, height = 5)

heliusdist %>% filter(!is.na(MetSyn)) %>% group_by(EthnicityTot, MetSyn) %>% summarise(count = length(MetSyn), .groups = "drop_last")
ggplot(data = heliusdist %>% filter(!is.na(MetSyn)), aes(x = MetSyn, y = distance)) +
    geom_violin(colour = NA, aes(fill = MetSyn)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "MetSyn", title = "MetSyn") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.format") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave(file.path(resultsfolder, "sg_distance_ethnicity_metsyn.pdf"), width = 6, height = 5)

(pl1 <- heliusdist %>% filter(!is.na(DM_new)) %>%
     ggplot(aes(x = DM_new, y = distance, fill = DM_new)) +
     geom_violin(colour = NA) +
     geom_boxplot(fill = "white", width = 0.2) +
     scale_fill_simpsons(guide = "none") +
     labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes", title = "New diabetes diagnosis") +
     stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                        label = "p.format") +
     theme_Publication())
ggsave(file.path(resultsfolder, "sg_newdiabetes.pdf"), width = 4.5, height = 5)

(pl1 <- heliusdist %>% filter(!is.na(HT_new)) %>%
        ggplot(aes(x = HT_new, y = distance, fill = HT_new)) +
        geom_violin(colour = NA) +
        geom_boxplot(fill = "white", width = 0.2) +
        scale_fill_simpsons(guide = "none") +
        labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension", title = "New hypertension diagnosis") +
        stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                           label = "p.format") +
        theme_Publication())
ggsave(file.path(resultsfolder, "sg_newhypertension.pdf"), width = 4.5, height = 5)

(pl1 <- heliusdist %>% filter(!is.na(MetSyn_new)) %>%
        ggplot(aes(x = MetSyn_new, y = distance, fill = MetSyn_new)) +
        geom_violin(colour = NA) +
        geom_boxplot(fill = "white", width = 0.2) +
        scale_fill_simpsons(guide = "none") +
        labs(y = "Bray-Curtis dissimilarity over FU time", x= "Metabolic syndrome", title = "New MetSyn diagnosis") +
        stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                           label = "p.format") +
        theme_Publication())
ggsave(file.path(resultsfolder, "sg_newmetsyn.pdf"), width = 4.5, height = 5)

(pl1 <- heliusdist %>% filter(!is.na(DM_new)) %>%
     ggplot(aes(x = DM_new, y = distance, fill = DM_new)) +
     geom_violin(colour = NA) +
     geom_boxplot(fill = "white", width = 0.2) +
     scale_fill_simpsons(guide = "none") +
     labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes", title = "New diabetes diagnosis") +
     stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                        label = "p.format") +
     facet_wrap(~EthnicityTot) +
     theme_Publication())
ggsave(file.path(resultsfolder, "sg_newdiabetes_ethnicity.pdf"), width = 6, height = 5)

(pl1 <- heliusdist %>% filter(!is.na(HT_new)) %>%
        ggplot(aes(x = HT_new, y = distance, fill = HT_new)) +
        geom_violin(colour = NA) +
        geom_boxplot(fill = "white", width = 0.2) +
        scale_fill_simpsons(guide = "none") +
        labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension", title = "New hypertension diagnosis") +
        stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                           label = "p.format") +
        facet_wrap(~EthnicityTot) +
        theme_Publication())
ggsave(file.path(resultsfolder, "sg_newhypertension_ethnicity.pdf"), width = 6, height = 5)
