## Calculate distances, plot PCoA and PCA

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(doParallel)
registerDoParallel(6)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
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
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

#### Load data ####
df <- readRDS("data/16s/clin_betadiversity.RDS")

#### Output folder ####
resultsfolder <- "results/ordination"
dir.create(resultsfolder, showWarnings = FALSE)

#### Bray-Curtis distance: male-female ####
braycurt <- dbray %>% 
    ggplot(aes(PC1_BC, PC2_BC)) +
    geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    scale_color_manual(values = rev(pal_nejm()(2))) +
    theme_Publication() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    labs(color = "") +
    stat_ellipse(aes(color = Sex), type = "norm") + 
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res1$`Pr(>F)`, ", r2 = ", format(round(res1$R2[1],3), nsmall = 3)))
ggsave(braycurt, filename = "results/ordination/PCoA_BrayCurtis.pdf", device = "pdf", width = 8, height = 8)


#### Weighted UniFrac ####
unifracpl <- df %>% 
    ggplot(aes(PC1_unifrac, PC2_unifrac)) +
    geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
    xlab(paste0('PCo1 (', round(expl_variance[3], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[4], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_manual(values = rev(pal_nejm()(2))) +
    ggtitle('PCoA Weighted UniFrac') +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    labs(color = "") +
    stat_ellipse(aes(color = Sex), type = "norm")
ggsave(unifracpl, filename = "results/ordination/PCoA_WeightedUnifrac.pdf", device = "pdf", width = 8, height = 8)

#### CLR-transformed PCA ####
plclr <- df %>% 
    ggplot(aes(PC1_clr, PC2_clr)) +
    geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
    xlab(paste0('PC1 (', round(expl_variance[5], digits = 1),'%)')) + # PComp, not PCoord
    ylab(paste0('PC2 (', round(expl_variance[6], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_manual(values = rev(pal_nejm()(2))) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    ggtitle('PCA CLR-transformed') + 
    labs(color = "") +
    stat_ellipse(aes(color = Sex), type = "norm")
ggsave(plclr, filename = "results/ordination/PCA_CLR.pdf", device = "pdf", width = 8, height = 8)
