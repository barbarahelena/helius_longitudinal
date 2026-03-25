## TBD
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(afex)

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

linearmixed <- function(data, var){
    data1 <- data %>% mutate(var = {{ var }})
    model1_v4 <- lmer(var ~ MetSyn*timenum + (1|ID), 
                      data = data1)
    res_v4 <- summary(model1_v4)
    print(res_v4)
    pval <- format(round(res_v4$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres <- cbind(group1 = "baseline", group2 = "follow-up", pval)
    statres <- tibble::as_tibble(statres)
    statres$p.signif <- case_when(
        statres$pval < 0.001 ~paste0("***"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.05 ~paste0("*"),
        statres$pval >= 0.05 ~paste0("")
    )
    statres <- statres %>% filter(p.signif != "")
    return(statres)
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

df <- readRDS("data/clinicaldata_long.RDS") |> 
  mutate(timenum = case_when(
    timepoint == "baseline" ~ 0,
    timepoint == "follow-up" ~ 1
  ),
   timenum = as.numeric(timenum))

dfodori <- left_join(odori, df, by = c("sampleID", "timepoint")) |> droplevels()
dim(dfodori)
dfalist <- left_join(alist, df, by = c("sampleID", "timepoint")) |> droplevels()
dim(dfalist)
dfparab <- left_join(parab, df, by = c("sampleID", "timepoint")) |> droplevels()

dim(dfparab)

(snvdist_parab_lm <- dfparab %>% linearmixed(SNV_distance_mean))
(snvdist_dfodori_lm <- dfodori %>% linearmixed(SNV_distance_mean))
(snvdist_dfalist_lm <- dfalist %>% linearmixed(SNV_distance_mean))

(snvdist_parab_lm <- dfparab %>% linearmixed(nucl_diversity_rarefied))
(snvdist_dfodori_lm <- dfodori %>% linearmixed(nucl_diversity_rarefied))
(snvdist_dfalist_lm <- dfalist %>% linearmixed(nucl_diversity_rarefied))

(snvdist_parab_lm <- dfparab %>% linearmixed(population_divergent_sites))
(snvdist_dfodori_lm <- dfodori %>% linearmixed(population_divergent_sites))
(snvdist_dfalist_lm <- dfalist %>% linearmixed(population_divergent_sites))

names(dfparab)
