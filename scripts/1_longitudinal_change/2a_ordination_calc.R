## Calculate distances

## Libraries
library(tidyverse)
library(doParallel)
library(phyloseq)
library(mixOmics)
library(vegan)
registerDoParallel(8)

#### Load data ####
phydata <- readRDS("/projects/prjs1784/helius16s/results/phyloseq/rarefied/phyloseq_rarefied.RDS")
df_new <- readRDS("/projects/prjs1784/helius16s/data/clinicaldata_long.RDS")
tab <- as.data.frame(t(as(phydata@otu_table, 'matrix')))
tab_matrix <- t(as(phydata@otu_table, 'matrix'))

#### Bray-Curtis distance ####
print('Bray-Curtis distance total dataset')
bray <- vegan::vegdist(tab, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
write_lines(expl_variance_bray, "expl_var_bray.csv")
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
write_lines(expl_variance_unifrac, "expl_var_unifrac.csv")
dfpc <- pcoord$vectors[, c('Axis.1', 'Axis.2')] # get PCoA coordinates
dfpc <- as.data.frame(dfpc)
dfpc$ID <- rownames(dfpc)
dfpc <- left_join(dfpc, dbray, by = 'ID') %>%  # add metadata / covariates
    select(UniFrac1 = `Axis.1`, UniFrac2 = `Axis.2`, everything(.))

saveRDS(dfpc, "clin_betadiversity.RDS")
