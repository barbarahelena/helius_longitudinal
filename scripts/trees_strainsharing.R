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
plot_strain <- function(SGB, tax){
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
    df <- readRDS("data/clinicaldata_long.RDS") %>% dplyr::filter(sampleID %in% tree1$tip.label)
    
    df$timeeth_factor <- case_when(df$EthnicityTot == "Dutch" & df$timepoint == "baseline" ~ "Dutch baseline",
                              df$EthnicityTot == "Dutch" & df$timepoint == "follow-up" ~ "Dutch follow-up",
                              df$EthnicityTot == "South-Asian Surinamese" & df$timepoint == "baseline" ~ "SAS baseline",
                              df$EthnicityTot == "South-Asian Surinamese" & df$timepoint == "follow-up" ~ "SAS follow-up")
    df$timedutch_factor <- case_when(str_detect(df$timeeth_factor, "SAS") ~ "SAS",
                                     !str_detect(df$timeeth_factor, "SAS") ~ df$timeeth_factor) %>% as.factor() %>% 
                                fct_relevel(., "SAS", after = 0L)
    df$timesas_factor <- case_when(!str_detect(df$timeeth_factor, "Dutch") ~ df$timeeth_factor,
                                   str_detect(df$timeeth_factor, "Dutch") ~ "Dutch") %>% as.factor() %>% 
                            fct_relevel(., "Dutch", after = 0L)
    
    groupInfo1 <- split(tree1$tip.label, df$timeeth_factor[match(df$sampleID, tree1$tip.label)])
    groupInfo2 <- split(tree1$tip.label, df$timedutch_factor[match(df$sampleID, tree1$tip.label)])
    groupInfo3 <- split(tree1$tip.label, df$timesas_factor[match(df$sampleID, tree1$tip.label)])
    
    tree1a <- groupOTU(tree1, groupInfo1)
    tree1b <- groupOTU(tree1, groupInfo2)
    tree1c <- groupOTU(tree1, groupInfo3)
    
    taxo <- tax$Species[which(tax$SGB == SGB)]
    
    pl1 <- ggtree(tree1a, aes(color = group), layout = "circular") +
                scale_color_manual(values = pal_simpsons()(8)[4:8]) +
                labs(title = taxo, color = "") +
                theme_Publication()
    ggsave(pl1, filename = str_c("results/strainsharing/trees/", SGB, "/tree_", SGB, "_allgroups.pdf"),
           width = 10, height = 10)
    
    pl2 <- ggtree(tree1b, aes(color = group), layout = "circular") +
                scale_color_manual(values = c(pal_simpsons()(8)[7:8], "lightgrey")) +
                labs(title = taxo, color = "") +
                theme_Publication()
    ggsave(pl2, filename = str_c("results/strainsharing/trees/", SGB, "/tree_", SGB, "_dutch.pdf"),
           width = 10, height = 10)
    
    pl3 <- ggtree(tree1c, aes(color = group), layout = "circular") +
                scale_color_manual(values = c("lightgrey", pal_simpsons()(8)[7:8])) +
                labs(title = taxo, color = "") +
                theme_Publication()
    ggsave(pl3, filename = str_c("results/strainsharing/trees/", SGB, "/tree_", SGB, "_sas.pdf"),
           width = 10, height = 10)
    
}

taxtrue <- readRDS("data/shotgun/sharing_tax.RDS")
sgbs <- levels(taxtrue$SGB)
# writeLines(sgbs, "data/shotgun/sgbs.txt")

for(s in sgbs){
    plot_strain(s, taxtrue)
}
