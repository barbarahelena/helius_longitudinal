## Calculate distances

## Libraries
library(tidyverse)
library(doParallel)
library(phyloseq)
library(mixOmics)
library(vegan)
registerDoParallel(24)

#### Load data ####
path <- '/projects/0/prjs0784/helius_longitudinal/'
phydata <- readRDS(file.path(path, "data/phyloseq_sampledata.RDS"))
df_new <- readRDS(file.path(path, "data/clinicaldata_long.RDS"))
tab <- as.data.frame(t(as(phydata@otu_table, 'matrix')))
tab_matrix <- t(as(phydata@otu_table, 'matrix'))

#### Bray-Curtis distance ####
print('Bray-Curtis distance total dataset')
bray <- vegan::vegdist(tab, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
write_lines(expl_variance_bray, "results/ordination/expl_var_bray.csv")
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, df_new, by = 'ID') %>%  # add metadata / covariates
        select(BrayPCo1 = `Axis.1`, BrayPCo2 = `Axis.2`, everything(.))

#### Weighted UniFrac ####
print('Weighted UniFrac')
wunifrac <- UniFrac(phydata, normalized = T, weighted = T, parallel = T)
pcoord <- ape::pcoa(wunifrac, correction = "cailliez")
expl_variance_unifrac <- pcoord$values$Rel_corr_eig * 100
write_lines(expl_variance_unifrac, "results/ordination/expl_var_unifrac.csv")
dfpc <- pcoord$vectors[, c('Axis.1', 'Axis.2')] # get PCoA coordinates
dfpc <- as.data.frame(dfpc)
dfpc$ID <- rownames(dfpc)
dfpc <- left_join(dfpc, dbray, by = 'ID') %>%  # add metadata / covariates
    select(UniFrac1 = `Axis.1`, UniFrac2 = `Axis.2`, everything(.))

#### CLR-transformed PCA ####
print('CLR-transformation and PCA')
pseudocount <- 1
# scale should be F when CLR is used
pc <- mixOmics::pca(tab_matrix + pseudocount, center = T, scale = F, logratio = 'CLR') 
expl_variance_clr <- pc$prop_expl_var$X * 100
write_lines(expl_variance_clr, "results/ordination/expl_var_clr.csv")
dfclr <- pc$x[, c('PC1', 'PC2')]
dfclr <- as.data.frame(dfclr)
dfclr$ID <- rownames(dfclr)
dfclr$timepoint <- left_join(dfclr, dfpc, by = 'ID') %>%  # add metadata / covariates
    select(CLR_PC1 = PC1, CLR_PC2 = PC2, everything(.))

saveRDS(dfclr, "data/16s/clin_betadiversity.RDS")
