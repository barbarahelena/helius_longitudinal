## Figure 2 — Cardiometabolic disease and microbiome instability (16S)
## Strain sharing — new diabetes
##
## Produces:
##   2E: dmnew.pdf

## libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggridges)

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

transfnum <- function(var) {
    var1 <- as.numeric(gsub(",", ".", gsub("\\.", "", var)))
    return(var1)
}

#### Output folder ####
resultsfolder <- "results/2_cmb_microbiome"
dir.create(resultsfolder, showWarnings = FALSE)

#### Data ####
df <- rio::import("data/shotgun/strainsharing_merged.csv")
thres <- rio::import("data/shotgun/thresholds_merged.csv") %>%
    mutate(
        across(c("n_markers", "n_samples", "aln_length", "avg_gap_prop",
                 "threshold_value", "max_youden", "false_positive_rate", "false_negative_rate"),
               transfnum))
thres$n_markers <- NULL # bug in pipeline: n_samples = n_markers, n_samples not extracted from info..
thres$n_markers <- thres$n_samples
thres$n_samples <- NULL

#### Calculation strain sharing metric ####
colnames(df) <- str_remove(colnames(df), "sharing_")
sharing_sgb <- apply(df[,3:ncol(df)], 2, function(x) sum(x, na.rm = TRUE))
df$sharing_sum <- apply(df[,3:ncol(df)], 1, function(x) sum(x, na.rm = TRUE))
df$strain_total <- apply(df[,3:ncol(df)], 1, function(x) sum(!is.na(x), na.rm = TRUE))
df$sharing_perc <- (df$sharing_sum / df$strain_total) * 100

#### Merge with clinical data ####
clin <- readRDS("data/clinicaldata_long.RDS")
dfsh <- df %>%
    filter(str_detect(sampleid_1, "HELIBA_")) %>%
    dplyr::select(sampleID = sampleid_1, sharing_perc)
dftot <- inner_join(clin, dfsh, by = "sampleID")

#### Figure 2E — New-onset diabetes (strain sharing) ####

# dmnew.pdf
dftot %>% group_by(DM_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc), .groups = "drop_last")
ggplot(data = dftot %>% filter(!is.na(DM_new)), aes(x = DM_new, y = sharing_perc)) +
    geom_violin(colour = NA, aes(fill = DM_new)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - new DM", x = "") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication()
ggsave(file.path(resultsfolder, "dmnew.pdf"), width = 4.5, height = 5)
