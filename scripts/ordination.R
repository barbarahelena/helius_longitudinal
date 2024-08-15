## Calculate distances, plot PCoA and PCA

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

#### Load data ####
df <- readRDS("data/16s/clin_betadiversity.RDS") %>% dplyr::select(1:2, sampleID = ID, 4:5)
helius <- readRDS("data/clinicaldata_long.RDS")
df <- left_join(df, helius, by = c("sampleID"))
ev_bray <- read.csv("results/ordination/expl_var_bray.csv", header = FALSE)
ev_unifrac <- read.csv("results/ordination/expl_var_unifrac.csv", header = FALSE)
wunifrac <- readRDS("results/ordination/wunifrac.RDS")
bray <- readRDS("results/ordination/bray.RDS")

#### Output folder ####
resultsfolder <- "results/ordination"
dir.create(resultsfolder, showWarnings = FALSE)

#### Bray-Curtis distance ####
print('PERMANOVA..')
set.seed(1234)
# distance matrix and metadata must have the same sample order
dfanova <- df[match(attributes(bray)[["Labels"]], df$sampleID),]
all(dfanova$sampleID == attributes(bray)[["Labels"]]) # TRUE
dim(df)
res1 <- adonis2(bray ~ timepoint, data = df) # PERMANOVA
print(res1)

(braycurt <- df %>%
    ggplot(aes(BrayPCo1, BrayPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = DM_new, fill = DM_new), type = "norm",
                 alpha = 0.1) +
    geom_point(aes(color = timepoint), size = 1, alpha = 0.5) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(ev_bray$V1[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(ev_bray$V1[2], digits = 1),'%)')) +
    scale_color_manual(values = pal_lancet()(2)) +
    scale_fill_manual(values = pal_lancet()(2), guide = "none") +
    scale_alpha_manual(guide = "none") +
    theme_Publication() +
    labs(color = "", alpha = "") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res1$`Pr(>F)`, ", r2 = ", 
                           format(round(res1$R2[1],3), nsmall = 3))
             ))
# ggsave(braycurt, filename = "results/ordination/PCoA_BrayCurtis.pdf", device = "pdf", width = 8, height = 8)


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

#### Weighted UniFrac ####
dfanova <- df[match(attributes(wunifrac)[["Labels"]], df$sampleID),]
all(dfanova$sampleID == attributes(wunifrac)[["Labels"]]) # TRUE
dim(dfanova)
res2 <- adonis2(wunifrac ~ timepoint, data = dfanova) # PERMANOVA
print(res2)

(wunifracpl <- df %>% 
        ggplot(aes(UniFrac1, UniFrac2)) +
        stat_ellipse(geom = "polygon", aes(color = timepoint, fill = timepoint), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = timepoint), size = 1, alpha = 0.5) +
        ggtitle("PCoA Weighted UniFrac") +
        xlab(paste0('PCo1 (', round(ev_unifrac$V1[1]*100, digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(ev_unifrac$V1[2]*100, digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        scale_alpha_manual(guide = "none") +
        theme_Publication() +
        labs(color = "", alpha = "") +
        annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
                 label = str_c("PERMANOVA: p = ", res2$`Pr(>F)`, ", r2 = ", 
                               format(round(res2$R2[1],3), nsmall = 3))
        ))
# ggsave(wunifracpl, filename = "results/ordination/PCoA_WeightedUnifrac.pdf", device = "pdf", width = 8, height = 8)


## Distance between datapoints
braymat <- as.matrix(bray)
# Generate all possible combinations of IDs
all_combinations <- t(combn(unique(rownames(braymat)), 2, simplify = TRUE))
# Create a data frame with combinations and distances
data_long <- data.frame(
    sampleID1 = all_combinations[, 1],
    sampleID2 = all_combinations[, 2]
)
data_long <- data_long %>% 
    filter(str_remove(sampleID1, "HELIBA_") == str_remove(sampleID2, "HELIFU_")) %>% 
    mutate(
        ID = str_c("S", str_remove(sampleID1, "HELIBA_"))
    )
for(a in 1:nrow(data_long)){
    distbray = braymat[paste0(data_long$sampleID1[a]), paste0(data_long$sampleID2[a])]
    data_long$distance[a] <- distbray
}
heliusdist <- inner_join(data_long, helius, by = "ID") %>% filter(timepoint == "baseline")
saveRDS(heliusdist, "data/16s/braydistance_delta.RDS")

## Plots 
comp <- list(c("Dutch", "Moroccan"), c("South-Asian Surinamese", "Moroccan"))
ggplot(data = heliusdist %>% filter(!is.na(EthnicityTot)), 
       aes(x = fct_reorder(EthnicityTot, .x = distance, .fun = median), y = distance)) +
    geom_violin(aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis distance", title = "Distance baseline to follow-up", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/ordination/distance_ethnicities.pdf", width = 6, height = 5)

comp <- list(c("Ghanaian", "Moroccan"), c("Ghanaian", "Turkish"))
ggplot(data = heliusdist %>% filter(!is.na(FUtime)), 
       aes(x = fct_reorder(EthnicityTot, .x = FUtime, .fun = median), y = FUtime)) +
    geom_violin(aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Follow-up time (years)", title = "Follow-up time", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "wilcox.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/ordination/futime_ethnicities.pdf", width = 6, height = 5)

#### Categorical outcomes unstratified and stratified ####
ggplot(data = heliusdist %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = distance)) +
    geom_violin(aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension (baseline)", title = "Hypertension") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication()
ggsave("results/ordination/distance_hypertension.pdf", width = 4, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = distance)) +
    geom_violin(aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension (baseline)", title = "Hypertension") +
    # stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
    #                    label = "p.signif", method = "wilcox.test") +
    facet_wrap(~Ethnicity) +
    theme_Publication()
ggsave("results/ordination/distance_hypertension_ethnicity.pdf", width = 7, height = 7)

ggplot(data = heliusdist %>% filter(!is.na(DM)), aes(x = DM, y = distance)) +
    geom_violin(aes(fill = DM)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes (baseline)", title = "Diabetes") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication()
ggsave("results/ordination/distance_dm.pdf", width = 4, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = distance)) +
    geom_violin(aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes (baseline)", title = "Diabetes") +
    # stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
    #                    label = "p.signif", method = "wilcox.test") +
    facet_wrap(~Ethnicity) +
    theme_Publication()
ggsave("results/ordination/distance_diabetes_ethnicity.pdf", width = 7, height = 7)

ggplot(data = heliusdist %>% filter(!is.na(MetSyn)), aes(x = MetSyn, y = distance)) +
    geom_violin(aes(fill = MetSyn)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Metabolic syndrome (baseline)", title = "Metabolic syndrome") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif") +
    theme_Publication()
ggsave("results/ordination/distance_metsyn.pdf", width = 4, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(MetSyn)), aes(x = MetSyn, y = distance)) +
    geom_violin(aes(fill = MetSyn)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Metabolic syndrome (baseline)", title = "Metabolic syndrome") +
    # stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
    #                    label = "p.signif", method = "wilcox.test") +
    facet_wrap(~Ethnicity) +
    theme_Publication()
ggsave("results/ordination/distance_metsyn_ethnicity.pdf", width = 7, height = 7)


#### Continuous outcomes unstratified and stratified ####
ggplot(data = heliusdist %>% filter(!is.na(HbA1c_delta)), aes(x = distance, y = HbA1c_delta)) +
    geom_jitter(color = "royalblue", alpha = 0.3) +
    geom_smooth(color = "black", method = "lm") +
    labs(y = "Delta HbA1c", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and HbA1c change") +
    stat_cor() +
    theme_Publication()
ggsave("results/ordination/braycurtis_deltahba1c.pdf", width = 4.5, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(HbA1c_delta)), aes(x = distance, y = HbA1c_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta HbA1c", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and HbA1c change") +
    stat_cor() +
    theme_Publication()
ggsave("results/ordination/braycurtis_deltahba1c_ethnicity.pdf", width = 7, height = 7)

ggplot(data = heliusdist %>% filter(!is.na(LDL_delta)), aes(x = distance, y = LDL_delta)) +
    geom_jitter(color = "royalblue", alpha = 0.3) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Delta LDL", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and LDL change") +
    stat_cor(method = "spearman") +
    theme_Publication()
ggsave("results/ordination/braycurtis_deltaldl.pdf", width = 4.5, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(LDL_delta)), aes(x = distance, y = LDL_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta LDL", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and LDL change") +
    stat_cor(method = "spearman") +
    theme_Publication()
ggsave("results/ordination/braycurtis_deltaldl_ethnicity.pdf", width = 7, height = 7)

ggplot(data = heliusdist %>% filter(!is.na(Age_delta)), aes(x = distance, y = Age_delta)) +
    geom_jitter(color = "royalblue", alpha = 0.3, height = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Delta age", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and age change") +
    stat_cor(method = "spearman") +
    theme_Publication()
ggsave("results/ordination/braycurtis_deltaage.pdf", width = 4.5, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(Age_delta)), aes(x = distance, y = Age_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5, height = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta Age", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and age change") +
    stat_cor(method = "spearman") +
    theme_Publication()
ggsave("results/ordination/braycurtis_deltaage_ethnicity.pdf", width = 7, height = 7)

ggplot(data = heliusdist %>% filter(!is.na(BMI_delta)), aes(x = BMI_delta, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.3, height = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Delta BMI", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and BMI change") +
    stat_cor() +
    theme_Publication()
ggsave("results/ordination/braycurtis_deltabmi.pdf", width = 4.5, height = 5)

ggplot(data = heliusdist %>% filter(!is.na(BMI_delta)), aes(x = distance, y = BMI_delta)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Delta BMI", x= "Bray-Curtis dissimilarity over FU time", title = "Bray-Curtis and BMI change") +
    stat_cor() +
    theme_Publication()
ggsave("results/ordination/braycurtis_deltabmi_ethnicity.pdf", width = 7, height = 7)

ggplot(data = heliusdist %>% filter(!is.na(FUtime)), aes(x = FUtime, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "FU time (years)", title = "FU time and sample distance") +
    stat_cor() +
    theme_Publication()
ggsave("results/ordination/braycurtis_futime.pdf", width = 4, height = 4)

ggplot(data = heliusdist %>% filter(!is.na(FUtime)), aes(x = FUtime, y = distance)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "FU time (years)", title = "FU time and sample distance") +
    stat_cor(method = "pearson") +
    theme_Publication()
ggsave("results/ordination/braycurtis_futime_ethnicity.pdf", width = 7, height = 7)

# Differences between new diagnoses and controls stratified for EthnicityTot
heliusdist %>% filter(!is.na(HT_BPMed)) %>% group_by(EthnicityTot, HT_BPMed) %>% summarise(count = length(HT_BPMed))
ggplot(data = heliusdist %>% filter(!is.na(HT_BPMed)), aes(x = HT_BPMed, y = distance)) +
    geom_violin(aes(fill = HT_BPMed)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension", title = "Hypertension") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.format") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave("results/distance_ethnictiy_hypertension.pdf", width = 6, height = 11)

heliusdist %>% filter(!is.na(DM)) %>% group_by(EthnicityTot, DM) %>% summarise(count = length(DM))
ggplot(data = heliusdist %>% filter(!is.na(DM)), aes(x = DM, y = distance)) +
    geom_violin(aes(fill = DM)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes", title = "Diabetes") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.format") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave("results/distance_ethnictiy_diabetes.pdf", width = 6, height = 11)


heliusdist %>% filter(!is.na(MetSyn)) %>% group_by(EthnicityTot, MetSyn) %>% summarise(count = length(MetSyn))
ggplot(data = heliusdist %>% filter(!is.na(MetSyn)), aes(x = MetSyn, y = distance)) +
    geom_violin(aes(fill = MetSyn)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x= "MetSyn", title = "MetSyn") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.format") +
    facet_wrap(~EthnicityTot) +
    # geom_text(data = cts, aes(label=count), position=position_dodge(width=1.0)) +
    theme_Publication()
ggsave("results/distance_ethnictiy_metsyn.pdf", width = 6, height = 11)

#### New diagnoses ####
(pl1 <- heliusdist %>% filter(!is.na(DM_new)) %>%  
     ggplot(aes(x = DM_new, y = distance, fill = DM_new)) +
     geom_violin() +
     geom_boxplot(fill = "white", width = 0.2) +
     scale_fill_simpsons(guide = "none") +
     labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes", title = "New diabetes diagnosis") +
     stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                        label = "p.format") +
     theme_Publication())
ggsave("results/ordination/newdiabetes.pdf", width = 4.5, height = 5)

(pl1 <- heliusdist %>% filter(!is.na(HT_new)) %>%  
        ggplot(aes(x = HT_new, y = distance, fill = HT_new)) +
        geom_violin() +
        geom_boxplot(fill = "white", width = 0.2) +
        scale_fill_simpsons(guide = "none") +
        labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension", title = "New hypertension diagnosis") +
        stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                           label = "p.format") +
        theme_Publication())
ggsave("results/ordination/newhypertension.pdf", width = 4.5, height = 5)

(pl1 <- heliusdist %>% filter(!is.na(MetSyn_new)) %>%  
        ggplot(aes(x = MetSyn_new, y = distance, fill = MetSyn_new)) +
        geom_violin() +
        geom_boxplot(fill = "white", width = 0.2) +
        scale_fill_simpsons(guide = "none") +
        labs(y = "Bray-Curtis dissimilarity over FU time", x= "Metabolic syndrome", title = "New MetSyn diagnosis") +
        stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                           label = "p.format") +
        theme_Publication())
ggsave("results/ordination/newmetsyn.pdf", width = 4.5, height = 5)

(pl1 <- heliusdist %>% filter(!is.na(DM_new)) %>%  
     ggplot(aes(x = DM_new, y = distance, fill = DM_new)) +
     geom_violin() +
     geom_boxplot(fill = "white", width = 0.2) +
     scale_fill_simpsons(guide = "none") +
     labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes", title = "New diabetes diagnosis") +
     stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                        label = "p.format") +
     facet_wrap(~EthnicityTot) +
     theme_Publication())
ggsave("results/ordination/newdiabetes_ethnicity.pdf", width = 7, height = 7)

(pl1 <- heliusdist %>% filter(!is.na(HT_new)) %>%  
        ggplot(aes(x = HT_new, y = distance, fill = HT_new)) +
        geom_violin() +
        geom_boxplot(fill = "white", width = 0.2) +
        scale_fill_simpsons(guide = "none") +
        labs(y = "Bray-Curtis dissimilarity over FU time", x= "Hypertension", title = "New hypertension diagnosis") +
        stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                           label = "p.format") +
        facet_wrap(~EthnicityTot) +
        theme_Publication())
ggsave("results/ordination/newhypertension_ethnicity.pdf", width = 7, height = 7)
