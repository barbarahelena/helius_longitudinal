## Calculate distances, plot PCoA and PCA

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(doParallel)
#registerDoParallel(6)

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
df <- readRDS("data/16s/clin_betadiversity.RDS") %>% dplyr::select(1:5)
helius <- readRDS("data/clinicaldata_delta.RDS")
df <- left_join(df, helius, by = "ID")
ev_bray <- read.csv("results/ordination/expl_var_bray.csv", header = FALSE)
ev_unifrac <- read.csv("results/ordination/expl_var_unifrac.csv", header = FALSE)
wunifrac <- readRDS("results/ordination/wunifrac.RDS")
bray <- readRDS("results/ordination/bray.RDS")

#### Output folder ####
resultsfolder <- "results/ordination"
dir.create(resultsfolder, showWarnings = FALSE)

#### Bray-Curtis distance ####
print('PERMANOVA..')
set.seed(1234)
# distance matrix and metadata must have the same sample order
dfanova <- df %>% slice(match(attributes(bray)[["Labels"]], ID)) 
all(dfanova$ID == attributes(bray)[["Labels"]]) # TRUE
dim(dfanova)
res1 <- adonis2(bray ~ timepoint, data = dfanova) # PERMANOVA
print(res1)

(braycurt <- df %>% 
    ggplot(aes(BrayPCo1, BrayPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = timepoint, fill = timepoint), type = "norm",
                 alpha = 0.1) +
    geom_point(aes(color = timepoint), size = 1, alpha = 0.5) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(ev_bray$V1[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(ev_bray$V1[2], digits = 1),'%)')) +
    scale_color_manual(values = pal_lancet()(2)) +
    scale_fill_manual(values = pal_lancet()(2), guide = "none") +
    scale_alpha_manual(guide = "none") +
    theme_Publication() +
    labs(color = "", alpha = "") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res1$`Pr(>F)`, ", r2 = ", 
                           format(round(res1$R2[1],3), nsmall = 3))
             ))
ggsave(braycurt, filename = "results/ordination/PCoA_BrayCurtis.pdf", device = "pdf", width = 8, height = 8)


#### Weighted UniFrac ####
dfanova <- df %>% slice(match(attributes(wunifrac)[["Labels"]], ID)) 
all(dfanova$ID == attributes(wunifrac)[["Labels"]]) # TRUE
dim(dfanova)
res2 <- adonis2(wunifrac ~ timepoint, data = dfanova) # PERMANOVA
print(res2)

(wunifracpl <- df %>% 
        ggplot(aes(UniFrac1, UniFrac2)) +
        stat_ellipse(geom = "polygon", aes(color = timepoint, fill = timepoint), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = timepoint), size = 1, alpha = 0.5) +
        ggtitle("PCoA Weighted UniFrac") +
        xlab(paste0('PCo1 (', round(ev_unifrac$V1[1]*100, digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(ev_unifrac$V1[2]*100, digits = 1),'%)')) +
        scale_color_manual(values = pal_lancet()(2)) +
        scale_fill_manual(values = pal_lancet()(2), guide = "none") +
        scale_alpha_manual(guide = "none") +
        theme_Publication() +
        labs(color = "", alpha = "") +
        annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
                 label = str_c("PERMANOVA: p = ", res2$`Pr(>F)`, ", r2 = ", 
                               format(round(res2$R2[1],3), nsmall = 3))
        ))
ggsave(wunifracpl, filename = "results/ordination/PCoA_WeightedUnifrac.pdf", device = "pdf", width = 8, height = 8)


## Distance between datapoints
braymat <- as.matrix(bray)
# Generate all possible combinations of IDs
all_combinations <- t(combn(unique(rownames(braymat)), 2, simplify = TRUE))
# Create a data frame with combinations and distances
data_long <- data.frame(
    ID1 = all_combinations[, 1],
    ID2 = all_combinations[, 2],
    Distance = braymat[all_combinations]
)
data_filt <- data_long %>%
    filter(str_remove(ID1, "HELIBA_") == str_remove(ID2, "HELIFU_")) %>% 
    mutate(ID = str_c("S", str_remove(ID1, "HELIBA_"))) %>% 
    select(-ID1, -ID2)

heliusdist <- inner_join(data_filt, helius, by = "ID")

## Plots 


