# Gutsmash Longitudinal Analysis with Linear Mixed Models
library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggrepel)

# Theme
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
           strip.text = element_text(face="bold"),
           plot.caption = element_text(size = rel(0.5), face = "italic"),
           plot.subtitle = element_text(size=8, hjust = 0.5, face = "italic")))
}

# Data import
df <- rio::import("data/shotgun/gutsmash_results/population_pathways.tsv")
dim(df)
head(df)[1:5,1:5]
df <- df |> filter(! sample %in% c("HELIBA_103370", "HELIFU_103370"))
head(df)
rownames(df) <- df$sample
df$sample <- NULL
df <- as.matrix(df)
df_rel <- df / rowSums(df)
df_rel <- as.data.frame(df_rel)
dim(df_rel)
pathway_cols <- colnames(df_rel) # for later use
df_rel$sampleID <- rownames(df_rel)

clinical <- readRDS("data/clinicaldata_long.RDS")
df_clin  <- df_rel |> left_join(clinical, by = "sampleID") |> droplevels()
df_clin <- df_clin |> filter(!is.na(EthnicityTot))
table(df_clin$EthnicityTot) # sample per ethnicity

