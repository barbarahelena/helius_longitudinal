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
shotdata <- readRDS("data/shotgun/shotgun_abundance.RDS")
df_new <- readRDS("data/clinicaldata_long.RDS")
shotmat <- as.matrix(shotdata)

#### Bray-Curtis distance ####
print('Bray-Curtis distance total dataset')
bray <- vegan::vegdist(shotmat, method = 'bray')
saveRDS(bray, "data/shotgun/bray_shotgun.RDS")
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
write_lines(expl_variance_bray, "results/ordination/expl_var_bray_shotgun.csv")
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, df_new, by = 'ID') %>%  # add metadata / covariates
        select(BrayPCo1 = `Axis.1`, BrayPCo2 = `Axis.2`, everything(.))

# #### Weighted UniFrac ####
# print('Weighted UniFrac')
# wunifrac <- UniFrac(shotmat, normalized = T, weighted = T, parallel = T)
# pcoord <- ape::pcoa(wunifrac, correction = "cailliez")
# expl_variance_unifrac <- pcoord$values$Rel_corr_eig * 100
# write_lines(expl_variance_unifrac, "results/ordination/expl_var_unifrac_shotgun.csv")
# dfpc <- pcoord$vectors[, c('Axis.1', 'Axis.2')] # get PCoA coordinates
# dfpc <- as.data.frame(dfpc)
# dfpc$ID <- rownames(dfpc)
# dfpc <- left_join(dfpc, dbray, by = 'ID') %>%  # add metadata / covariates
#     select(UniFrac1 = `Axis.1`, UniFrac2 = `Axis.2`, everything(.))

saveRDS(dbray, "data/shotgun/clin_betadiversity_shotgun.RDS")
