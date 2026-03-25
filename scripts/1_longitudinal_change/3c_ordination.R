## Calculate distances, plot PCoA and PCA

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
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
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

#### Output folder ####
resultsfolder <- "results/1_longitudinal_change/ordination"
dir.create(resultsfolder, showWarnings = FALSE)

#### Ethnicity colour palette ####
eth_colors <- c(
    "Dutch"                  = "#709AE1FF",
    "South-Asian Surinamese" = "#FED439FF",
    "African Surinamese"     = "#8A9197FF",
    "Ghanaian"               = "#D2AF81FF",
    "Turkish"                = "#FD7446FF",
    "Moroccan"               = "#D5E4A2FF"
)

#### Load data ####
df <- readRDS("data/16s/archive/clin_betadiversity.RDS") %>% dplyr::select(1:2, sampleID = ID, 4:5)
helius <- readRDS("data/clinicaldata_long.RDS")
df <- left_join(df, helius, by = c("sampleID"))
ev_bray <- read.csv("results/1_longitudinal_change/ordination/expl_var_bray.csv", header = FALSE)
bray <- readRDS("results/1_longitudinal_change/ordination/bray.RDS")

#### Bray-Curtis distance ####
print('PERMANOVA..')
set.seed(1234)
# distance matrix and metadata must have the same sample order
dfanova <- df[match(attributes(bray)[["Labels"]], df$sampleID),]
all(dfanova$sampleID == attributes(bray)[["Labels"]]) # TRUE
dim(df)
res1 <- adonis2(bray ~ timepoint, data = df) # PERMANOVA
print(res1)

# Figure 1B: clean timepoint-coloured PCoA
pl_fig1_B <- df %>%
    ggplot(aes(BrayPCo1, BrayPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = timepoint, fill = timepoint), type = "norm", alpha = 0.1) +
    geom_point(aes(color = timepoint), size = 1, alpha = 0.5) +
    xlab(paste0("PCo1 (", round(ev_bray$V1[1], 1), "%)")) +
    ylab(paste0("PCo2 (", round(ev_bray$V1[2], 1), "%)")) +
    scale_color_manual(values = pal_lancet()(2)) +
    scale_fill_manual(values = pal_lancet()(2), guide = "none") +
    labs(color = "", title = "Community composition shift") +
    theme_Publication() +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
             label = paste0("paste('PERMANOVA: R'^2*' = ",
                            format(round(res1$R2[1], 3), nsmall = 3),
                            ", p = ", res1$`Pr(>F)`[1], "')"),
             parse = TRUE)

ggsave(pl_fig1_B, filename = "results/1_longitudinal_change/ordination/PCoA_BrayCurtis.pdf", width = 8, height = 8)

## Distance between datapoints
braymat <- as.matrix(bray)
# Generate all possible combinations of IDs
all_combinations <- t(combn(unique(rownames(braymat)), 2, simplify = TRUE))
# Create a data frame with combinations and distances
data_long <- data.frame(
    sampleID1 = all_combinations[, 1],
    sampleID2 = all_combinations[, 2]
)
data_long <- data_long %>% 
    filter(str_remove(sampleID1, "HELIBA_") == str_remove(sampleID2, "HELIFU_")) %>% 
    mutate(
        ID = str_c("S", str_remove(sampleID1, "HELIBA_"))
    )
for(a in 1:nrow(data_long)){
    distbray = braymat[paste0(data_long$sampleID1[a]), paste0(data_long$sampleID2[a])]
    data_long$distance[a] <- distbray
}
heliusdist <- inner_join(data_long, helius, by = "ID") %>% filter(timepoint == "baseline")
saveRDS(heliusdist, "data/16s/braydistance_delta.RDS")

## Plots 
comp <- list(c("Dutch", "Moroccan"), c("South-Asian Surinamese", "Moroccan"))
ggplot(data = heliusdist %>% filter(!is.na(EthnicityTot) & EthnicityTot != "Other"), 
       aes(x = fct_reorder(EthnicityTot, .x = distance, .fun = median), y = distance)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Bray-Curtis distance", title = "Distance baseline to follow-up", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/1_longitudinal_change/ordination/distance_ethnicities.pdf", width = 6, height = 5)
pl_fig1_C <- last_plot()

comp <- list(c("Ghanaian", "Moroccan"), c("Ghanaian", "Turkish"))
ggplot(data = heliusdist %>% filter(!is.na(FUtime) & EthnicityTot != "Other"), 
       aes(x = fct_reorder(EthnicityTot, .x = FUtime, .fun = median), y = FUtime)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Follow-up time (years)", title = "Follow-up time", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "wilcox.test") +
    theme_Publication() +
    coord_flip()
ggsave("results/1_longitudinal_change/ordination/futime_ethnicities.pdf", width = 6, height = 5)
pl_fig1_futime <- last_plot()

ggplot(data = heliusdist %>% filter(!is.na(FUtime) & EthnicityTot != "Other"), aes(x = FUtime, y = distance)) +
    geom_jitter(color = "royalblue", alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_manual(values = eth_colors, guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "FU time (years)", title = "FU time and sample distance") +
    stat_cor() +
    theme_Publication()
ggsave("results/1_longitudinal_change/ordination/braycurtis_futime.pdf", width = 4, height = 4)
pl_fig1_A <- last_plot()

ggplot(data = heliusdist %>% filter(!is.na(FUtime) & EthnicityTot != "Other"), aes(x = FUtime, y = distance)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_manual(values = eth_colors, guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "FU time (years)", title = "FU time and sample distance") +
    stat_cor(method = "pearson") +
    theme_Publication()
ggsave("results/1_longitudinal_change/ordination/braycurtis_futime_ethnicity.pdf", width = 7, height = 7)
