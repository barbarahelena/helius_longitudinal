#### Alpha diversity indices

## Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggsci)
library(phyloseq)

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

#### alpha diversity data ####
df <- readRDS("data/16s/clin_alphadiversity.RDS")
# counts <- sample_sums(phydata@otu_table)
# counts # samples should all sum up to 14932

## Output folder
resultsfolder <- "results/alphadiversity"
dir.create(resultsfolder, showWarnings = FALSE)

## Diversity metrics
# Shannon plots
(plshan <- ggplot(data = df, aes(x = timepoint, y = shannon, fill = timepoint)) +
    geom_violin() +
    scale_fill_manual(values = rev(pal_simpsons()(7)[c(1,7)]), guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    stat_compare_means(label.y = 5.5) +
    labs(title = "Shannon index", y = "Shannon index", x="") + 
    theme_Publication())
# ggsave(plshan, filename = "results/alphadiversity/shannon.svg", width = 4, height = 5)
# ggsave(plshan, filename = "results/alphadiversity/shannon.pdf", width = 4, height = 5)

# Richness plots
(plrich <- ggplot(data = dfspec, aes(x = timepoint, y = richness, fill = timepoint)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_manual(values = rev(pal_simpsons()(7)[c(1,7)]), guide = "none") + 
    labs(title = "Species richness", y = "Number of species", x = "") +
    stat_compare_means(method = "wilcox.test"))
ggsave(plrich, filename = "results/alphadiversity/richness.pdf", width = 4, height = 5)
ggsave(plrich, filename = "results/alphadiversity/richness.svg", width = 4, height = 5)

## Faith's PD
(plfaith <- ggplot(data = dffai, aes(x = timepoint, y = PD, fill = timepoint)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_manual(values = rev(pal_simpsons()(7)[c(1,7)]), guide = "none") + 
    labs(title = "Faith's PD", y = "Faith's phylogenetic diversity", x = "") +
    stat_compare_means(method = "wilcox.test"))
# ggsave(plfaith, filename = "results/alphadiversity/faiths.pdf", device = "pdf", width = 4, height = 5)
# ggsave(plfaith, filename = "results/alphadiversity/faiths.svg", device = "svg", width = 4, height = 5)

## Ggarrange
pl_total <- ggarrange(plshan, plrich, plfaith, labels = c("A", "B", "C"), nrow =1)
ggsave(pl_total, filename = "results/alphadiversity/alphadivplots.pdf", width = 11, height = 5.5)
ggsave(pl_total, filename = "results/alphadiversity/alphadivplots.svg", width = 11, height = 5.5)


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
# plshan <- ggplot(data = df_shan, aes(x = Sex, y = shannon, fill = Sex)) +
#     geom_violin() +
#     scale_fill_manual(values = rev(pal_nejm()(2)), guide = "none") +
#     geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
#     stat_compare_means(label.y = 5.5) +
#     labs(title = "Shannon index", y = "Shannon index", x="") + 
#     theme_Publication()
# ggsave(plshan, filename = "results/alphadiversity/shannon_shotgun.svg", width = 4, height = 5)
# ggsave(plshan, filename = "results/alphadiversity/shannon_shotgun.pdf", width = 4, height = 5)
# 
# 
# ## Species richness
# richness <- rio::import("data/metaphlan/diversity/combined_table_richness.tsv")
# richness <- richness %>% select(ID = V1, richness = observed) %>% 
#     mutate(ID = str_remove(ID, "_T1"))
# dfspec <- left_join(richness, df_new, by = "ID")
# 
# # Male-female
# plrich <- ggplot(data = dfspec, aes(x = Sex, y = richness, fill = Sex)) +
#     geom_violin()+
#     geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
#     theme_Publication() + 
#     scale_fill_manual(values = rev(pal_nejm()(2)), guide = "none") + 
#     labs(title = "Species richness", y = "Number of species", x = "") +
#     stat_compare_means(method = "wilcox.test")
# ggsave(plrich, filename = "results/alphadiversity/richness_shotgun_sex.pdf", width = 4, height = 5)
# ggsave(plrich, filename = "results/alphadiversity/richness_shotgun_sex.svg", width = 4, height = 5)
#
# 
# ## Ggarrange male-female
# pl_total <- ggarrange(plshan, plrich, labels = c("A", "B"), nrow =1)
# ggsave(pl_total, filename = "results/alphadiversity/alphadivplots_shotgun.pdf", width = 8, height = 5.5)
# ggsave(pl_total, filename = "results/alphadiversity/alphadivplots_shotgun.svg", width = 8, height = 5.5)
# 

