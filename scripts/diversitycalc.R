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
# ## Load data
# df_new <- rio:: import("data/clinicaldata_shotgun.RDS")
# 
# ## Diversity metrics
# # Shannon plots
# shannon <- rio::import("data/metaphlan/diversity/combined_table_shannon.tsv") 
# shannon <- shannon %>% select(ID = V1, shannon = diversity_shannon) %>% 
#     mutate(ID = str_remove(ID, "_T1"))
# df_shan <- left_join(shannon, df_new, by = "ID")
# 
# 
# ## Species richness
# richness <- rio::import("data/metaphlan/diversity/combined_table_richness.tsv")
# richness <- richness %>% select(ID = V1, richness = observed) %>% 
#     mutate(ID = str_remove(ID, "_T1"))
# dfspec <- left_join(richness, df_new, by = "ID")
#
# 

