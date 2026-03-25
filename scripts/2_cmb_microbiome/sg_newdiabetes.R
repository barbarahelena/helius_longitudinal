## Figure 2 — Cardiometabolic disease and microbiome instability (16S)
## Shotgun metagenomics — new diabetes Bray-Curtis distance
##
## Produces:
##   2E: sg_newdiabetes.pdf

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
    suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family)
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
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
}

#### Load data ####
df <- readRDS("data/shotgun/clin_betadiversity_shotgun.RDS") %>%
    dplyr::select(1:2, sampleID = ID)
helius <- readRDS("data/clinicaldata_long.RDS")
df <- left_join(df, helius, by = c("sampleID"))
bray <- readRDS("data/shotgun/bray_shotgun.RDS")

#### Output folder ####
resultsfolder <- "results/2_cmb_microbiome"
dir.create(resultsfolder, showWarnings = FALSE)

#### Compute pairwise Bray-Curtis distance (baseline → follow-up) ####
braymat <- as.matrix(bray)
all_combinations <- t(combn(unique(rownames(braymat)), 2, simplify = TRUE))
data_long <- data.frame(
    sampleID1 = all_combinations[, 1],
    sampleID2 = all_combinations[, 2]
)
data_long <- data_long %>%
    filter(str_remove(sampleID1, "HELIBA_") == str_remove(sampleID2, "HELIFU_")) %>%
    mutate(ID = str_c("S", str_remove(sampleID1, "HELIBA_")))
for(a in 1:nrow(data_long)){
    distbray = braymat[paste0(data_long$sampleID1[a]), paste0(data_long$sampleID2[a])]
    data_long$distance[a] <- distbray
}
heliusdist <- inner_join(data_long, helius, by = "ID") %>% filter(timepoint == "baseline")

#### Figure 2E — New-onset diabetes (Bray-Curtis, shotgun) ####

# sg_newdiabetes.pdf
(pl1 <- heliusdist %>% filter(!is.na(DM_new)) %>%
     ggplot(aes(x = DM_new, y = distance, fill = DM_new)) +
     geom_violin(colour = NA) +
     geom_boxplot(fill = "white", width = 0.2) +
     scale_fill_simpsons(guide = "none") +
     labs(y = "Bray-Curtis dissimilarity over FU time", x= "Diabetes", title = "New diabetes diagnosis") +
     stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                        label = "p.format") +
     theme_Publication())
ggsave(file.path(resultsfolder, "sg_newdiabetes.pdf"), width = 4.5, height = 5)
