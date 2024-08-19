## Plotting strain sharing trees
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(tidytree)
library(ggtree)

## Open data
plot_strain_dm <- function(SGB, tax){
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
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "right",
                    legend.key.size= unit(0.5, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    plot.caption = element_text(size = rel(0.8), face = "italic")
            ))
        
    } 
    dir.create(str_c("results/strainsharing/trees/", SGB), showWarnings = FALSE)
    tree1 <- ape::read.tree(str_c("data/shotgun/trees/RAxML_bestTree.t__", SGB, ".StrainPhlAn4.tre"))
    df <- readRDS("data/clinicaldata_long.RDS") %>% dplyr::filter(sampleID %in% tree1$tip.label) %>% 
        filter(!is.na(DM))
    
    df$timedm <- case_when(is.na(df$DM_new) & df$timepoint == "baseline" ~ "DM - baseline",
                           is.na(df$DM_new) & df$timepoint == "follow-up" ~ "DM - follow-up",
                           !is.na(df$DM_new) & df$timepoint == "baseline" ~ "No DM baseline",
                           !is.na(df$DM_new) & df$timepoint == "follow-up" ~ "No DM follow-up")
    df$timedmyes <- case_when(str_detect(df$timedm, "No") ~ "No DM",
                              !str_detect(df$timedm, "No") ~ df$timedm) %>% as.factor() %>% 
        fct_relevel(., "No DM", after = 0L)
    df$timedmno <- case_when(str_detect(df$timedm, "No") ~ df$timedm,
                             str_detect(df$timedm, "DM") ~ "DM") %>% as.factor() %>% 
        fct_relevel(., "DM", after = 0L)
    
    groupInfo1 <- split(tree1$tip.label, df$timedm[match(df$sampleID, tree1$tip.label)])
    groupInfo2 <- split(tree1$tip.label, df$timedmyes[match(df$sampleID, tree1$tip.label)])
    groupInfo3 <- split(tree1$tip.label, df$timedmno[match(df$sampleID, tree1$tip.label)])
    
    tree1a <- groupOTU(tree1, groupInfo1)
    tree1b <- groupOTU(tree1, groupInfo2)
    tree1c <- groupOTU(tree1, groupInfo3)
    
    taxo <- tax$Species[which(tax$SGB == SGB)]
    
    pl2 <- ggtree(tree1b, aes(color = group), layout = "circular") +
        scale_color_manual(values = c(pal_simpsons()(8)[7:8], "lightgrey")) +
        labs(title = taxo, color = "") +
        theme_Publication()
    ggsave(pl2, filename = str_c("results/strainsharing/trees/", SGB, "/tree_", SGB, "_dm.pdf"),
           width = 10, height = 10)
    
    pl3 <- ggtree(tree1c, aes(color = group), layout = "circular") +
        scale_color_manual(values = c("lightgrey", pal_simpsons()(8)[7:8])) +
        labs(title = taxo, color = "") +
        theme_Publication()
    ggsave(pl3, filename = str_c("results/strainsharing/trees/", SGB, "/tree_", SGB, "_nodm.pdf"),
           width = 10, height = 10)
    
}

for(s in sgbs){
    plot_strain_dm(s, taxtrue)
}
