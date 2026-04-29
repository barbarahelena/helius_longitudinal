## Bray-Curtis PCoA for CAZy (Cayman) gene family composition

## Libraries
library(tidyverse)
library(vegan)
library(ggsci)
library(ggpubr)

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family) +
     theme(plot.title = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA, fill = NA),
           plot.background = element_rect(colour = NA, fill = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold", size = rel(0.8)),
           axis.title.y = element_text(angle=90, vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(size = rel(0.7)),
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing  = unit(0, "cm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")))
}

#### Load data ####
df_raw <- rio::import("data/shotgun/cayman_results/families_cpm_table.tsv") |> select(-HELIBA_103370, -HELIFU_103370)
rownames(df_raw) <- df_raw$family
df_raw$family <- NULL
caymat <- t(as.matrix(df_raw))  # samples in rows, families in columns
caymat <- log10(caymat + 0.1)
clinical <- readRDS("data/clinicaldata/clinicaldata_long.RDS")

#### Output folder ####
resultsfolder <- "results/4_functional_change/cayman"
dir.create(resultsfolder, showWarnings = FALSE, recursive = TRUE)

#### Bray-Curtis distance ####
print('Bray-Curtis distance CAZy composition')
bray <- vegan::vegdist(caymat, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$sampleID <- rownames(dbray)
df <- left_join(dbray, clinical, by = 'sampleID') %>%
        select(BrayPCo1 = `Axis.1`, BrayPCo2 = `Axis.2`, everything(.))

#### PERMANOVA ####
set.seed(14)
dfanova <- df[match(attributes(bray)[["Labels"]], df$sampleID),]
all(dfanova$sampleID == attributes(bray)[["Labels"]]) # TRUE
res1 <- adonis2(bray ~ timepoint * EthnicityTot, data = dfanova, by = "term")
print(res1)

#### PCoA plot - timepoint ####
(braycurt <- df %>%
    ggplot(aes(BrayPCo1, BrayPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = fct_rev(timepoint), fill = fct_rev(timepoint)), type = "norm",
                 alpha = 0.1) +
    geom_point(aes(color = fct_rev(timepoint)), size = 1, alpha = 0.5) +
    ggtitle("PCoA Bray-Curtis distance (CAZy families)") +
    xlab(paste0('PCo1 (', round(expl_variance_bray[[1]], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance_bray[[2]], digits = 1),'%)')) +
    scale_color_manual(values = pal_simpsons()(2)) +
    scale_fill_manual(values = pal_simpsons()(2), guide = "none") +
    scale_alpha_manual(guide = "none") +
    theme_Publication() +
    labs(color = "", alpha = "") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res1$`Pr(>F)`[1], ", R2 = ",
                           format(round(res1$R2[1],3), nsmall = 3))
             ))
ggsave(braycurt, filename = "results/4_functional_change/cayman/PCoA_BrayCurtis_cayman.pdf", device = "pdf", width = 8, height = 8)

#### PCoA plot - ethnicity (baseline) ####
# Subset distance matrix for baseline samples with ethnicity data
bray_mat <- as.matrix(bray)
idx_bl <- which(dfanova$timepoint == "baseline" & !is.na(dfanova$EthnicityTot))
bray_bl <- as.dist(bray_mat[idx_bl, idx_bl])
res_eth <- adonis2(bray_bl ~ EthnicityTot, data = dfanova[idx_bl, ], by = "term")
print(res_eth)

(ethbray <- df %>% filter(!is.na(EthnicityTot)) %>% filter(timepoint == "baseline") %>%
    ggplot(aes(BrayPCo1, BrayPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = fct_rev(EthnicityTot), fill = fct_rev(EthnicityTot)), type = "norm",
                 alpha = 0.1) +
    geom_point(aes(color = fct_rev(EthnicityTot)), size = 1, alpha = 0.5) +
    ggtitle("PCoA Bray-Curtis (CAZy) - Baseline by Ethnicity") +
    xlab(paste0('PCo1 (', round(expl_variance_bray[[1]], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance_bray[[2]], digits = 1),'%)')) +
    scale_color_manual(values = pal_simpsons()(2)) +
    scale_fill_manual(values = pal_simpsons()(2), guide = "none") +
    theme_Publication() +
    labs(color = "", alpha = "") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res_eth$`Pr(>F)`[1], ", r2 = ",
                           format(round(res_eth$R2[1],3), nsmall = 3))
             ))
ggsave(ethbray, filename = "results/4_functional_change/cayman/PCoA_BrayCurtis_cayman_ethnicity.pdf", device = "pdf", width = 8, height = 8)

#### PCoA plot - ethnicity (FU) ####
# Subset distance matrix for FU samples with ethnicity data
bray_mat <- as.matrix(bray)
idx_bl <- which(dfanova$timepoint == "follow-up" & !is.na(dfanova$EthnicityTot))
bray_bl <- as.dist(bray_mat[idx_bl, idx_bl])
res_eth <- adonis2(bray_bl ~ EthnicityTot, data = dfanova[idx_bl, ], by = "term")
print(res_eth)

(ethbray <- df %>% filter(!is.na(EthnicityTot)) %>% filter(timepoint == "follow-up") %>%
    ggplot(aes(BrayPCo1, BrayPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = fct_rev(EthnicityTot), fill = fct_rev(EthnicityTot)), type = "norm",
                 alpha = 0.1) +
    geom_point(aes(color = fct_rev(EthnicityTot)), size = 1, alpha = 0.5) +
    ggtitle("PCoA Bray-Curtis (CAZy) - Follow-up by Ethnicity") +
    xlab(paste0('PCo1 (', round(expl_variance_bray[[1]], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance_bray[[2]], digits = 1),'%)')) +
    scale_color_manual(values = pal_simpsons()(2)) +
    scale_fill_manual(values = pal_simpsons()(2), guide = "none") +
    theme_Publication() +
    labs(color = "", alpha = "") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res_eth$`Pr(>F)`[1], ", r2 = ",
                           format(round(res_eth$R2[1],3), nsmall = 3))
             ))
ggsave(ethbray, filename = "results/4_functional_change/cayman/PCoA_BrayCurtis_cayman_ethnicity_FU.pdf", device = "pdf", width = 8, height = 8)

#### Bray-Curtis distance between baseline and follow-up per individual ####
braymat <- as.matrix(bray)
all_combinations <- t(combn(unique(rownames(braymat)), 2, simplify = TRUE))
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
heliusdist <- inner_join(data_long, clinical, by = "ID") %>% filter(timepoint == "baseline")
saveRDS(heliusdist, "data/shotgun/cayman_braydistance_delta.RDS")

#### Distance by ethnicity ####
ggplot(data = heliusdist %>% filter(!is.na(EthnicityTot)),
       aes(x = fct_reorder(fct_rev(EthnicityTot), .x = distance, .fun = median), y = distance)) +
    geom_violin(colour = NA, aes(fill = fct_rev(EthnicityTot))) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Bray-Curtis distance", title = "CAZy distance baseline to follow-up", x = "") +
    stat_compare_means(tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/4_functional_change/cayman/cayman_distance_ethnicities.pdf", width = 6, height = 4)

#### Correlation CAZy vs MB (shotgun) Bray-Curtis dissimilarity ----
mb_dist <- readRDS("data/shotgun/braydistance_delta.RDS") %>%
  dplyr::select(ID, mb_distance = distance)

cazy_mb_dist <- heliusdist %>%
  dplyr::select(ID, cazy_distance = distance, EthnicityTot) %>%
  inner_join(mb_dist, by = "ID")

# Spearman correlation
cor_res <- cor.test(cazy_mb_dist$cazy_distance, cazy_mb_dist$mb_distance, method = "spearman")
subtitle_cor <- paste0("Spearman rho = ", round(cor_res$estimate, 3),
                       ", p = ", formatC(cor_res$p.value, format = "e", digits = 2))

# Unstratified
(pl_cor <- ggplot(cazy_mb_dist, aes(x = mb_distance, y = cazy_distance)) +
  geom_point(color = "royalblue", alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", color = "black") +
  theme_Publication() +
  labs(x = "Microbiome Bray-Curtis dissimilarity (BL to FU)",
       y = "CAZy Bray-Curtis dissimilarity (BL to FU)",
       title = "Microbiome vs CAZy compositional change",
       subtitle = subtitle_cor))
ggsave("results/4_functional_change/cayman/cayman_vs_mb_braycurtis.pdf", pl_cor, width = 5.5, height = 5)

# Stratified by ethnicity
(pl_cor_eth <- ggplot(cazy_mb_dist %>% filter(!is.na(EthnicityTot)),
                      aes(x = mb_distance, y = cazy_distance)) +
  geom_point(aes(color = fct_rev(EthnicityTot)), alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_simpsons(guide = "none") +
  facet_wrap(~fct_rev(EthnicityTot)) +
  stat_cor(method = "spearman", size = 3) +
  theme_Publication() +
  labs(x = "Microbiome Bray-Curtis dissimilarity (BL to FU)",
       y = "CAZy Bray-Curtis dissimilarity (BL to FU)",
       title = "Microbiome vs CAZy compositional change"))
ggsave("results/4_functional_change/cayman/cayman_vs_mb_braycurtis_ethnicity.pdf", pl_cor_eth, width = 8, height = 6)

# Linear model: CAZy_BC ~ Microbiome_BC * Ethnicity
model_bc <- lm(cazy_distance ~ mb_distance * EthnicityTot, data = cazy_mb_dist %>% filter(!is.na(EthnicityTot)))
res_bc <- summary(model_bc)
print(res_bc)
ci_bc <- confint(model_bc)

# Extract all coefficients
bc_results <- data.frame(
  term = rownames(res_bc$coefficients),
  estimate = res_bc$coefficients[, 1],
  se = res_bc$coefficients[, 2],
  conflow = ci_bc[, 1],
  confhigh = ci_bc[, 2],
  pval = res_bc$coefficients[, 4]
) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv2(bc_results, "results/4_functional_change/cayman/lm_cazyBC_mbBC_ethnicity.csv", row.names = FALSE)

#### Canberra distance ####
print('Canberra distance CAZy composition')
canberra <- vegan::vegdist(caymat, method = 'canberra')
pcoord_can <- ape::pcoa(canberra, correction = "cailliez")
expl_variance_can <- pcoord_can$values$Rel_corr_eig * 100
dcan <- pcoord_can$vectors[, c('Axis.1', 'Axis.2')]
dcan <- as.data.frame(dcan)
dcan$sampleID <- rownames(dcan)
df_can <- left_join(dcan, clinical, by = 'sampleID') %>%
        select(CanberraPCo1 = `Axis.1`, CanberraPCo2 = `Axis.2`, everything(.))

#### PERMANOVA (Canberra) ####
set.seed(14)
dfcanova <- df_can[match(attributes(canberra)[["Labels"]], df_can$sampleID),]
all(dfcanova$sampleID == attributes(canberra)[["Labels"]]) # TRUE
res1_can <- adonis2(canberra ~ timepoint * EthnicityTot, data = dfcanova, by = "term")
print(res1_can)

#### PCoA plot - timepoint (Canberra) ####
(canberra_tp <- df_can %>%
    ggplot(aes(CanberraPCo1, CanberraPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = fct_rev(timepoint), fill = fct_rev(timepoint)), type = "norm",
                 alpha = 0.1) +
    geom_point(aes(color = fct_rev(timepoint)), size = 1, alpha = 0.5) +
    ggtitle("PCoA Canberra distance (CAZy families)") +
    xlab(paste0('PCo1 (', round(expl_variance_can[[1]], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance_can[[2]], digits = 1),'%)')) +
    scale_color_manual(values = pal_simpsons()(2)) +
    scale_fill_manual(values = pal_simpsons()(2), guide = "none") +
    scale_alpha_manual(guide = "none") +
    theme_Publication() +
    labs(color = "", alpha = "") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res1_can$`Pr(>F)`[1], ", R2 = ",
                           format(round(res1_can$R2[1],3), nsmall = 3))
             ))
ggsave(canberra_tp, filename = "results/4_functional_change/cayman/PCoA_Canberra_cayman.pdf", device = "pdf", width = 8, height = 8)

#### PCoA plot - ethnicity baseline (Canberra) ####
can_mat <- as.matrix(canberra)
idx_bl_can <- which(dfcanova$timepoint == "baseline" & !is.na(dfcanova$EthnicityTot))
can_bl <- as.dist(can_mat[idx_bl_can, idx_bl_can])
res_eth_can_bl <- adonis2(can_bl ~ EthnicityTot, data = dfcanova[idx_bl_can, ], by = "term")
print(res_eth_can_bl)

(ethcan_bl <- df_can %>% filter(!is.na(EthnicityTot)) %>% filter(timepoint == "baseline") %>%
    ggplot(aes(CanberraPCo1, CanberraPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = fct_rev(EthnicityTot), fill = fct_rev(EthnicityTot)), type = "norm",
                 alpha = 0.1) +
    geom_point(aes(color = fct_rev(EthnicityTot)), size = 1, alpha = 0.5) +
    ggtitle("PCoA Canberra (CAZy) - Baseline by Ethnicity") +
    xlab(paste0('PCo1 (', round(expl_variance_can[[1]], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance_can[[2]], digits = 1),'%)')) +
    scale_color_manual(values = pal_simpsons()(2)) +
    scale_fill_manual(values = pal_simpsons()(2), guide = "none") +
    theme_Publication() +
    labs(color = "", alpha = "") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res_eth_can_bl$`Pr(>F)`[1], ", r2 = ",
                           format(round(res_eth_can_bl$R2[1],3), nsmall = 3))
             ))
ggsave(ethcan_bl, filename = "results/4_functional_change/cayman/PCoA_Canberra_cayman_ethnicity.pdf", device = "pdf", width = 8, height = 8)

#### PCoA plot - ethnicity FU (Canberra) ####
idx_fu_can <- which(dfcanova$timepoint == "follow-up" & !is.na(dfcanova$EthnicityTot))
can_fu <- as.dist(can_mat[idx_fu_can, idx_fu_can])
res_eth_can_fu <- adonis2(can_fu ~ EthnicityTot, data = dfcanova[idx_fu_can, ], by = "term")
print(res_eth_can_fu)

(ethcan_fu <- df_can %>% filter(!is.na(EthnicityTot)) %>% filter(timepoint == "follow-up") %>%
    ggplot(aes(CanberraPCo1, CanberraPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = fct_rev(EthnicityTot), fill = fct_rev(EthnicityTot)), type = "norm",
                 alpha = 0.1) +
    geom_point(aes(color = fct_rev(EthnicityTot)), size = 1, alpha = 0.5) +
    ggtitle("PCoA Canberra (CAZy) - Follow-up by Ethnicity") +
    xlab(paste0('PCo1 (', round(expl_variance_can[[1]], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance_can[[2]], digits = 1),'%)')) +
    scale_color_manual(values = pal_simpsons()(2)) +
    scale_fill_manual(values = pal_simpsons()(2), guide = "none") +
    theme_Publication() +
    labs(color = "", alpha = "") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res_eth_can_fu$`Pr(>F)`[1], ", r2 = ",
                           format(round(res_eth_can_fu$R2[1],3), nsmall = 3))
             ))
ggsave(ethcan_fu, filename = "results/4_functional_change/cayman/PCoA_Canberra_cayman_ethnicity_FU.pdf", device = "pdf", width = 8, height = 8)

#### Canberra distance between baseline and follow-up per individual ####
canmat <- as.matrix(canberra)
all_combinations_can <- t(combn(unique(rownames(canmat)), 2, simplify = TRUE))
data_long_can <- data.frame(
    sampleID1 = all_combinations_can[, 1],
    sampleID2 = all_combinations_can[, 2]
)
data_long_can <- data_long_can %>%
    filter(str_remove(sampleID1, "HELIBA_") == str_remove(sampleID2, "HELIFU_")) %>%
    mutate(
        ID = str_c("S", str_remove(sampleID1, "HELIBA_"))
    )
for(a in 1:nrow(data_long_can)){
    distcan = canmat[paste0(data_long_can$sampleID1[a]), paste0(data_long_can$sampleID2[a])]
    data_long_can$distance[a] <- distcan
}
heliusdist_can <- inner_join(data_long_can, clinical, by = "ID") %>% filter(timepoint == "baseline")
saveRDS(heliusdist_can, "data/shotgun/cayman_canberradistance_delta.RDS")

#### Distance by ethnicity (Canberra) ####
ggplot(data = heliusdist_can %>% filter(!is.na(EthnicityTot)),
       aes(x = fct_reorder(fct_rev(EthnicityTot), .x = distance, .fun = median), y = distance)) +
    geom_violin(colour = NA, aes(fill = fct_rev(EthnicityTot))) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Canberra distance", title = "CAZy Canberra distance baseline to follow-up", x = "") +
    stat_compare_means(tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/4_functional_change/cayman/cayman_canberradistance_ethnicities.pdf", width = 6, height = 4)
