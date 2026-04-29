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
df <- readRDS("data/16s/clin_betadiversity.RDS") 
df$sampleID <- NULL
df <- df |> dplyr::select(sampleID = ID, everything()) |> 
    dplyr::select(1:3)
heliusdf <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
df <- inner_join(df, heliusdf, by = "sampleID") |> dplyr::select(sampleID, ID, FUtime, Sex, EthnicityTot,
 BrayPCo1, BrayPCo2, timepoint)
ev_bray <- read.csv("data/16s/expl_var_bray.csv", header = FALSE)
bray <- readRDS("data/16s/bray.RDS")

#### Bray-Curtis distance ####
print('PERMANOVA..')
set.seed(1234)
# distance matrix and metadata must have the same sample order
dfanova <- df[match(attributes(bray)[["Labels"]], df$sampleID),]
all(dfanova$sampleID == attributes(bray)[["Labels"]]) # TRUE
dim(df)
res1 <- adonis2(bray ~ timepoint, data = dfanova) # PERMANOVA
print(res1)

# Figure 1B: clean timepoint-coloured PCoA
(pl_fig1_B <- df %>%
    ggplot(aes(BrayPCo1, BrayPCo2)) +
    stat_ellipse(geom = "polygon", aes(color = timepoint, fill = timepoint), type = "norm", alpha = 0.1) +
    geom_point(aes(color = timepoint), size = 1, alpha = 0.5) +
    xlab(paste0("PCo1 (", round(ev_bray$V1[1], 1), "%)")) +
    ylab(paste0("PCo2 (", round(ev_bray$V1[2], 1), "%)")) +
    scale_color_manual(values = c("#197EC0FF", "#F05C3BFF"), guide = guide_legend(position = "right")) +
    scale_fill_manual(values = c("#197EC0FF", "#F05C3BFF"), guide = "none") +
    labs(color = "", title = "Microbiota composition change") +
    theme_Publication() +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
             label = paste0("paste('PERMANOVA: R'^2*' = ",
             format(round(res1$R2[1], 3), nsmall = 3),
             ", p = ", res1$`Pr(>F)`[1], "')"), parse = TRUE, size = 4))

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
heliusdist <- inner_join(data_long, heliusdf, by = "ID") %>% filter(timepoint == "baseline")
saveRDS(heliusdist, "data/16s/braydistance_delta.RDS")

## Plots
fmt_pval <- function(p) {
    ifelse(p < 0.0001, "p < 0.0001", paste0("p = ", formatC(p, format = "f", digits = 3)))
}

comp <- list( c("Ghanaian", "Moroccan"), c("Ghanaian", "Turkish"), c("Ghanaian", "Dutch"))
fu_data <- heliusdist %>%
    filter(!is.na(FUtime) & EthnicityTot != "Other") %>%
    mutate(EthnicityTot = fct_reorder(EthnicityTot, FUtime, median))
(pl_fig1_A <- ggplot(data = fu_data,
       aes(x = EthnicityTot, y = FUtime)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Follow-up time (years)", title = "Follow-up time", x = "") +
    stat_compare_means(label = "p.format", tip.length = 0, comparisons = comp) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = 4:12) +
    theme_Publication() +
    coord_flip())
ggsave("results/1_longitudinal_change/ordination/futime_ethnicities.pdf", width = 6, height = 5)

(pl_fig1_C <- ggplot(data = heliusdist %>% filter(!is.na(FUtime) & EthnicityTot != "Other"), aes(x = FUtime, y = distance)) +
    geom_jitter(color = "#197EC0FF", alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_manual(values = eth_colors, guide = "none") +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "FU time (years)", title = "FU time and sample distance") +
    stat_cor() +
    theme_Publication())
ggsave("results/1_longitudinal_change/ordination/braycurtis_futime.pdf", width = 4, height = 4)

(ggplot(data = heliusdist %>% filter(!is.na(FUtime) & EthnicityTot != "Other"), aes(x = FUtime, y = distance)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.3, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_manual(values = eth_colors, guide = "none") +
    facet_wrap(~EthnicityTot) +
    labs(y = "Bray-Curtis dissimilarity over FU time", x = "FU time (years)", title = "FU time and sample distance") +
    stat_cor(method = "pearson") +
    theme_Publication())
ggsave("results/1_longitudinal_change/ordination/braycurtis_futime_ethnicity.pdf", width = 7, height = 7)

comp <- list(c("Dutch", "Moroccan"), c("South-Asian Surinamese", "Moroccan"), c("Turkish", "Moroccan"),
                    c("African Surinamese", "Moroccan"))
bc_data <- heliusdist %>%
    filter(!is.na(EthnicityTot) & EthnicityTot != "Other") %>%
    mutate(EthnicityTot = fct_reorder(EthnicityTot, distance, median))
(pl <- ggplot(data = bc_data,
       aes(x = EthnicityTot, y = distance)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(y = "Bray-Curtis distance", title = "Distance baseline to follow-up", x = "") +
    stat_compare_means(aes(label = sprintf("p = %s", ..p.format..)), tip.length = 0, comparisons = comp) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_Publication() +
    coord_flip())
ggsave("results/1_longitudinal_change/ordination/distance_ethnicities.pdf", width = 6, height = 5)

## Confounder-adjusted Bray-Curtis dissimilarity per ethnicity ####
heliusdist_adj <- heliusdist %>%
    filter(!is.na(EthnicityTot) & EthnicityTot != "Other") %>%
    filter(
        !is.na(Age) & !is.na(Sex) & !is.na(BMI) &
        !is.na(Metformin) & !is.na(PPI) & !is.na(AntiHT) & !is.na(Statins) &
        !is.na(DiscrMean_baseline) & !is.na(AlcCons)
    )

# Linear regression: distance ~ all confounders + ethnicity
lm_full <- lm(distance ~ Age + Sex + BMI + Metformin + PPI + FUtime + EthnicityTot,
              data = heliusdist_adj)
print(summary(lm_full))

# Confounder-only model: residuals + grand mean = adjusted dissimilarity
lm_confounders <- lm(distance ~ Age + Sex + BMI + Metformin + PPI + FUtime,
                     data = heliusdist_adj)

heliusdist_adj <- heliusdist_adj %>%
    mutate(
        dist_adjusted = residuals(lm_confounders) + mean(distance, na.rm = TRUE),
        EthnicityTot = fct_reorder(EthnicityTot, dist_adjusted, median)
    )

comp_adj <- list(c("South-Asian Surinamese", "Moroccan"),
 c("Dutch", "Moroccan"), c("Turkish", "Moroccan"), c("African Surinamese", "Moroccan"))
(pl_fig1_D <- ggplot(heliusdist_adj, aes(x = EthnicityTot, y = dist_adjusted)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(
        y = "Adjusted Bray-Curtis dissimilarity",
        title = "Distance baseline to follow-up",
        x = ""
    ) +
    stat_compare_means(
        aes(label = sprintf("p = %s", ..p.format..)),
        tip.length = 0,
        comparisons = comp_adj
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_Publication() +
    coord_flip())
ggsave("results/1_longitudinal_change/ordination/distance_ethnicities_adjusted.pdf", width = 6, height = 5)

## Confounder-adjusted Bray-Curtis dissimilarity per ethnicity ####

# Join with diet PCs (only available in clinicaldata_long_pcdiet.RDS)
pcdiet <- readRDS("data/clinicaldata_long_pcdiet.RDS") %>%
    filter(timepoint == "baseline") %>%
    dplyr::select(ID, DietPC1, DietPC2)

heliusdist_adj2 <- heliusdist %>%
    left_join(pcdiet, by = "ID") %>%
    filter(!is.na(EthnicityTot) & EthnicityTot != "Other") %>%
    filter(
        !is.na(Age) & !is.na(Sex) & !is.na(BMI) &
        !is.na(Metformin) & !is.na(PPI) & !is.na(AntiHT) & !is.na(Statins) &
        !is.na(DiscrMean_baseline) & !is.na(AlcCons) & !is.na(DietPC1) & !is.na(DietPC2)
    )

# Linear regression: distance ~ all confounders + ethnicity
lm_full <- lm(distance ~ Age + Sex + BMI + DietPC1 + DietPC2 +
                  Metformin + PPI + FUtime + EthnicityTot,
              data = heliusdist_adj2)
print(summary(lm_full))

# Confounder-only model: residuals + grand mean = adjusted dissimilarity
lm_confounders <- lm(distance ~ Age + Sex + BMI + DietPC1 + DietPC2 +
                         Metformin + PPI + FUtime,
                     data = heliusdist_adj2)

heliusdist_adj2 <- heliusdist_adj2 %>%
    mutate(
        dist_adjusted = residuals(lm_confounders) + mean(distance, na.rm = TRUE),
        EthnicityTot = fct_reorder(EthnicityTot, dist_adjusted, median)
    )

comp_adj <- list(c("South-Asian Surinamese", "Moroccan"), c("African Surinamese", "Moroccan"),
 c("Dutch", "Moroccan"), c("Turkish", "Moroccan"))
(pl_suppl_diet <- ggplot(heliusdist_adj2, aes(x = EthnicityTot, y = dist_adjusted)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(
        y = "Adjusted Bray-Curtis dissimilarity",
        title = "Confounder-adjusted\nincluding diet",
        x = ""
    ) +
    # stat_compare_means(
    #     aes(label = sprintf("p = %s", ..p.format..)),
    #     tip.length = 0,
    #     comparisons = comp_adj
    # ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_Publication() +
    coord_flip())
ggsave("results/1_longitudinal_change/ordination/suppl_distance_ethnicities_adjusted.pdf", width = 6, height = 5)

## Within-ethnicity spread of PCo1 and PCo2 at baseline and follow-up ####
df_eth <- df %>%
    filter(!is.na(EthnicityTot) & EthnicityTot != "Other") %>%
    mutate(
        EthnicityTot = factor(EthnicityTot, levels = names(eth_colors)),
        timepoint = factor(timepoint, levels = c("baseline", "follow-up"), labels = c("Baseline", "Follow-up"))
    )

# PCo1 per ethnicity, split by timepoint
(pl_pco1_eth <- ggplot(df_eth, aes(x = EthnicityTot, y = BrayPCo1)) +
    geom_violin(aes(fill = EthnicityTot), colour = NA, alpha = 0.8) +
    geom_boxplot(fill = "white", width = 0.2, outlier.size = 0.5) +
    facet_wrap(~timepoint) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(
        x = "",
        y = paste0("PCo1 (", round(ev_bray$V1[1], 1), "%)"),
        title = "Within-ethnicity spread – PCo1"
    ) +
    theme_Publication() +
    coord_flip())
ggsave("results/1_longitudinal_change/ordination/pco1_spread_per_ethnicity.pdf", width = 8, height = 5)

# PCo2 per ethnicity, split by timepoint
(pl_pco2_eth <- ggplot(df_eth, aes(x = EthnicityTot, y = BrayPCo2)) +
    geom_violin(aes(fill = EthnicityTot), colour = NA, alpha = 0.8) +
    geom_boxplot(fill = "white", width = 0.2, outlier.size = 0.5) +
    facet_wrap(~timepoint) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    labs(
        x = "",
        y = paste0("PCo2 (", round(ev_bray$V1[2], 1), "%)"),
        title = "Within-ethnicity spread – PCo2"
    ) +
    theme_Publication() +
    coord_flip())
ggsave("results/1_longitudinal_change/ordination/pco2_spread_per_ethnicity.pdf", width = 8, height = 5)

# Combined panel
(pl_pco_combined <- ggarrange(pl_pco1_eth, pl_pco2_eth, nrow = 2, labels = c("A", "B")))
ggsave("results/1_longitudinal_change/ordination/pco1_pco2_spread_per_ethnicity.pdf",
       pl_pco_combined, width = 8, height = 9)
