## TBD
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
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
                strip.text = element_text(face="bold")
        ))
    
} 

# Data
bac <- rio::import("data/shotgun/bac_instrain_gtdbtk.tsv") |> 
  mutate(
    timepoint = case_when(timepoint == "HELIBA" ~ "baseline",
                          timepoint == "HELIFU" ~ "follow-up"),
    sampleID = case_when(
      str_detect(sample, "_HELIBA") ~ str_c("HELIBA_", str_remove(sample, "_HELIBA")),
      str_detect(sample, "_HELIFU") ~ str_c("HELIFU_", str_remove(sample, "_HELIFU"))
    )
  )
head(bac)
names(bac)
bac$completeness
hist(bac$completeness)
bac$taxonomy

odori <- bac |> filter(str_detect(taxonomy, "Odoribacter splanchnicus")) # 20 strains
alist <- bac |> filter(str_detect(taxonomy, "Alistipes putredinis")) # 227 strains
parab <- bac |> filter(str_detect(taxonomy, "Parabacteroides distasonis")) # 149 strains

df <- readRDS("data/clinicaldata_long.RDS")
dfodori <- left_join(odori, df, by = c("sampleID", "timepoint")) |> droplevels()
dim(dfodori)
dfalist <- left_join(alist, df, by = c("sampleID", "timepoint")) |> droplevels()
dim(dfalist)
dfparab <- left_join(parab, df, by = c("sampleID", "timepoint")) |> droplevels()
dim(dfparab)

summary(dfodori$EthnicityTot)
summary(dfalist$EthnicityTot)
summary(dfparab$EthnicityTot)
names(odori)

## SNV_distance
ggplot(data = dfodori, aes(x = EthnicityTot, y = SNV_distance_mean,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(dfodori$EthnicityTot, dfodori$timepoint)

ggplot(data = dfalist, aes(x = EthnicityTot, y = SNV_distance_mean,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(dfalist$EthnicityTot, dfalist$timepoint)

ggplot(data = dfparab, aes(x = EthnicityTot, y = SNV_distance_mean,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(dfparab$EthnicityTot, dfparab$timepoint)

## nucl_diversity_rarefied
ggplot(data = dfodori, aes(x = EthnicityTot, y = nucl_diversity_rarefied,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(dfodori$EthnicityTot, dfodori$timepoint)

ggplot(data = dfalist, aes(x = EthnicityTot, y = nucl_diversity_rarefied,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(dfalist$EthnicityTot, dfalist$timepoint)

ggplot(data = dfparab, aes(x = EthnicityTot, y = nucl_diversity_rarefied,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication()
table(dfparab$EthnicityTot, dfparab$timepoint)

## population_divergent_sites
ggplot(data = dfodori, aes(x = EthnicityTot, y = popANI_reference,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = dfalist, aes(x = EthnicityTot, y = popANI_reference,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = dfparab, aes(x = EthnicityTot, y = popANI_reference,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## population_divergent_sites
ggplot(data = dfodori, aes(x = EthnicityTot, y = iRep,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = dfalist, aes(x = EthnicityTot, y = iRep,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = dfparab, aes(x = EthnicityTot, y = iRep,
      fill = EthnicityTot)) +
  geom_violin(colour = NA) +
  geom_boxplot(fill = "white", width = 0.25) +
  geom_point() +
  scale_fill_jco(guide = "none") +
  stat_compare_means() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

names(bac)
## Maybe this should be LMMs instead of wilcox BA and FU
## Variables of interest: popANI_reference, nucl_diversity_rarefied,
## SNV_distance_mean, conANI_reference, population_divergent_sites
