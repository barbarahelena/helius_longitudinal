#### HELIUS diet descriptives

## Libraries
library(dplyr)
library(ggsci)
library(ggplot2)
library(forcats)
library(ggpubr)
library(mixOmics)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
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

## Ethnicity colour palette
eth_colors <- c(
    "Dutch"                  = "#709AE1FF",
    "South-Asian Surinamese" = "#FED439FF",
    "African Surinamese"     = "#8A9197FF",
    "Ghanaian"               = "#D2AF81FF",
    "Turkish"                = "#FD7446FF",
    "Moroccan"               = "#D5E4A2FF"
)


## Output folder
resultsfolder1 <- "results/0_data_cleaning"
dir.create(resultsfolder1, showWarnings = FALSE)
resultsfolder <- "results/0_data_cleaning/diet"
dir.create(resultsfolder, showWarnings = FALSE)

## Load dataset
df <- readRDS("data/clinicaldata/clinicaldata_long.RDS")

## PCA diet
df_diet <- df %>% dplyr::select(ID, EthnicityTot, TotalCalories, Fiber, Protein, Protein_animal, FattyAcids,
                             Carbohydrates, Sodium_g) %>% 
    filter(!is.na(TotalCalories))
df_diet2 <- df_diet %>% dplyr::select(-ID, -EthnicityTot, -TotalCalories) %>% 
    mutate(across(everything(.), scale))
matdiet <- as.matrix(df_diet2)
tunediet <- tune.pca(matdiet, ncomp = 5, scale = TRUE)
# plot(tunediet)
pc <- mixOmics::pca(matdiet, ncomp = 2)
pcs <- as.data.frame(pc$variates$X)
pcs <- pcs %>% mutate(ID = df_diet$ID, EthnicityTot = df_diet$EthnicityTot)
expvar_diet <- pc$explained_variance[1:2]
loadings <- as.data.frame(pc$loadings$X)
loadings$Variables <- rownames(loadings)

(pcadiet <- pcs %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = EthnicityTot), size = 1, alpha = 1.0) +
        xlab(paste0('PC1 (', round(expvar_diet[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_diet[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = EthnicityTot, fill = EthnicityTot), linewidth = 1.0,
                     alpha = 0.1, type = "norm")+
        scale_color_manual(values = eth_colors) +
        scale_fill_manual(values = eth_colors, guide = "none") +
        labs(color = "", title = "PCA diet")+
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), 
                     arrow = arrow(length = unit(1/2, "picas")),
                     color = "black", linewidth = 0.9) +
        annotate("text", x = (loadings$PC1*13), y = (loadings$PC2*10),
                 label = loadings$Variables)
)
ggsave(pcadiet, filename = "results/0_data_cleaning/diet/PCA_diet_loading.pdf", width = 7, height = 7)

df <- left_join(df, pcs, by = c("ID", "EthnicityTot")) %>% 
    dplyr::select(everything(.), DietPC1=PC1, DietPC2=PC2)
saveRDS(df, "data/clinicaldata_long_pcdiet.RDS")

# Macronutrient groups
(pl1 <- ggplot(df_diet, aes(x=EthnicityTot, y=TotalCalories))+
    geom_violin(aes(fill=EthnicityTot), color = NA)+
    scale_fill_manual(values = eth_colors, guide = "none")+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'Kcal', title = "Total calories")+
    ggpubr::stat_compare_means(hide.ns = TRUE))

(pl2 <- ggplot(df_diet, aes(x=EthnicityTot, y=Fiber))+
    geom_violin(aes(fill=EthnicityTot), color = NA)+
    scale_fill_manual(values = eth_colors, guide = "none")+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'gram', title = "Fibers")+
    ggpubr::stat_compare_means(hide.ns = TRUE))

(pl3 <- ggplot(df_diet, aes(x=EthnicityTot, y=Protein))+
    geom_violin(aes(fill=EthnicityTot), color = NA)+
    scale_fill_manual(values = eth_colors, guide = "none")+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'gram', title = "Protein")+
    ggpubr::stat_compare_means(hide.ns = TRUE))

(pl4 <- ggplot(df_diet, aes(x=EthnicityTot, y=Protein_animal))+
    geom_violin(aes(fill=EthnicityTot), color = NA)+
    scale_fill_manual(values = eth_colors, guide = "none")+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'gram', title = "Animal protein")+
    ggpubr::stat_compare_means(hide.ns = TRUE))

(pl5 <- ggplot(df_diet, aes(x=EthnicityTot, y=FattyAcids))+
    geom_violin(aes(fill=EthnicityTot), color = NA)+
    scale_fill_manual(values = eth_colors, guide = "none")+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'gram', title = "Fatty acids")+
    ggpubr::stat_compare_means(hide.ns = TRUE))

(pl6 <- ggplot(df_diet, aes(x=EthnicityTot, y=Carbohydrates))+
    geom_violin(aes(fill=EthnicityTot), color = NA)+
    scale_fill_manual(values = eth_colors, guide = "none")+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'gram', title = "Carbohydrates")+
    ggpubr::stat_compare_means(hide.ns = TRUE))

(pl7 <- ggplot(df_diet, aes(x=EthnicityTot, y=Sodium_g))+
    geom_violin(aes(fill=EthnicityTot), color = NA)+
    scale_fill_manual(values = eth_colors, guide = "none")+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'gram', title = "Sodium")+
    ggpubr::stat_compare_means(hide.ns = TRUE))

(fig_macronutrients <- ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, pl7,
                                  ncol = 3, nrow = 3))
ggsave(fig_macronutrients, filename = file.path(resultsfolder, "macronutrients_by_ethnicity.pdf"),
       device = "pdf", width = 15, height = 15)

#### Energy-adjusted macronutrients (Willett residual method) ####
## For each macronutrient, regress on TotalCalories_baseline and store
## the scaled residuals as <var>_baseline_adj in clinicaldata_wide.RDS.
## These pre-adjusted variables are used by downstream analysis scripts
## instead of adjusting inline (which would introduce collinearity).

helius_wide_diet <- readRDS("data/clinicaldata/clinicaldata_wide.RDS")

macro_vars <- c("Protein", "FattyAcids", "Carbohydrates", "Fiber", "Sodium_g")

for (mac in macro_vars) {
    col     <- paste0(mac, "_baseline")
    col_adj <- paste0(mac, "_baseline_adj")
    if (col %in% names(helius_wide_diet) && "TotalCalories_baseline" %in% names(helius_wide_diet)) {
        complete_idx <- !is.na(helius_wide_diet[[col]]) & !is.na(helius_wide_diet[["TotalCalories_baseline"]])
        resid_vec <- rep(NA_real_, nrow(helius_wide_diet))
        fit <- lm(helius_wide_diet[[col]][complete_idx] ~ helius_wide_diet[["TotalCalories_baseline"]][complete_idx])
        resid_vec[complete_idx] <- residuals(fit)
        helius_wide_diet[[col_adj]] <- as.numeric(scale(resid_vec))
    }
}

saveRDS(helius_wide_diet, "data/clinicaldata/clinicaldata_wide.RDS")
