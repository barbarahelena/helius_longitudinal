#### Alpha diversity indices

## Libraries
library(tidyverse)
library(vegan)
library(phyloseq)

#### 16S data ####
phydata <- readRDS("data/phyloseq_sampledata.RDS")
df_new <- rio:: import("data/clinicaldata_long.RDS") %>% mutate(ID = sampleID)
tab <- as.data.frame(t(as(phydata@otu_table, 'matrix')))
tab_matrix <- t(as(phydata@otu_table, 'matrix'))

## Output folder
resultsfolder <- "results/alphadiversity"
dir.create(resultsfolder, showWarnings = FALSE)

## Diversity metrics
# Shannon plots
shannon <- vegan::diversity(tab, index = 'shannon')
df_shan <- data.frame(ID = names(shannon), shannon = shannon)
df_shan <- left_join(df_shan, df_new, by = "ID")

## Species richness
specrich <- specnumber(tab)
dfspec <- data.frame(ID = names(specrich), richness = specrich)
dfspec <- left_join(dfspec, df_shan, by = "ID")

## Faith's PD
faith <- picante::pd(samp = tab_matrix, tree = phydata@phy_tree)
dffai <- as.data.frame(faith)
dffai$ID <- rownames(faith)
dffai <- left_join(dffai, dfspec, by = "ID")

saveRDS(dffai, "data/16s/clin_alphadiversity.RDS")

#### Shotgun data ####
## Load data
df_new <- readRDS("data/clinicaldata_long.RDS")

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


