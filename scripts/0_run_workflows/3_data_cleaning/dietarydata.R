#### HELIUS diet descriptives

## Libraries
library(dplyr)
library(ggsci)
library(ggplot2)
library(forcats)
library(ggpubr)
library(mixOmics)
library(aplot)

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
   # "Ghanaian"               = "#D2AF81FF",
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
    filter(!is.na(TotalCalories)) |> droplevels()
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
        labs(color = "", title = "PCA diet") +
        theme(legend.position = "top") +
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), 
                     arrow = arrow(length = unit(1/2, "picas")),
                     color = "black", linewidth = 0.9) +
        annotate("text", x = (loadings$PC1*13), y = (loadings$PC2*10),
                 label = loadings$Variables)
)
ggsave(pcadiet, filename = "results/0_data_cleaning/diet/PCA_diet_loading.pdf", width = 7, height = 7)

(plright_diet <- ggplot(pcs, aes(x = EthnicityTot, y = PC2, fill = EthnicityTot)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    scale_x_discrete(expand = expansion(add = 0.3)) +
    theme_transparent())

(plbottom_diet <- ggplot(pcs, aes(x = fct_rev(EthnicityTot), y = PC1, fill = EthnicityTot)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    scale_x_discrete(expand = expansion(add = 0.3)) +
    theme_transparent() +
    coord_flip())

options("aplot_guides" = "keep")
ap_diet <- pcadiet %>%
    insert_bottom(plbottom_diet, height = 0.25) %>%
    insert_right(plright_diet, width = 0.25)
ggsave(ap_diet, filename = "results/0_data_cleaning/diet/PCA_diet_loading_box.pdf", device = "pdf", width = 9, height = 9)

df <- left_join(df, pcs, by = c("ID", "EthnicityTot")) %>% 
    dplyr::select(everything(.), DietPC1=PC1, DietPC2=PC2)
saveRDS(df, "data/clinicaldata_long_pcdiet.RDS")

# All pairwise ethnicity combinations for groups present in this dataset
eth_present <- intersect(names(eth_colors), unique(as.character(df_diet$EthnicityTot)))
eth_pairs <- combn(eth_present, 2, simplify = FALSE)

# Returns a stat_compare_means layer for significant pairs only, or NULL if none
sig_comparisons <- function(df, y_var, pairs, alpha = 0.05) {
    sig <- Filter(function(p) {
        g1 <- df[[y_var]][as.character(df$EthnicityTot) == p[1]]
        g2 <- df[[y_var]][as.character(df$EthnicityTot) == p[2]]
        tryCatch(wilcox.test(g1, g2)$p.value < alpha, error = function(e) FALSE)
    }, pairs)
    if (length(sig) == 0) return(NULL)
    ggpubr::stat_compare_means(comparisons = sig, label = "p.signif", tip.length = 0)
}
table(df_diet$EthnicityTot)
# Macronutrient groups - ethnicities ordered by median per nutrient
(pl1 <- ggplot(df_diet, aes(x = fct_reorder(EthnicityTot, TotalCalories, median, na.rm = TRUE), y = TotalCalories)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'Kcal', title = "Total calories") +
    sig_comparisons(df_diet, "TotalCalories", eth_pairs))
ggsave(pl1, filename = file.path(resultsfolder, "violin_TotalCalories.pdf"), device = "pdf", width = 6, height = 6)

(pl2 <- ggplot(df_diet, aes(x = fct_reorder(EthnicityTot, Fiber, median, na.rm = TRUE), y = Fiber)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Fibers") +
    sig_comparisons(df_diet, "Fiber", eth_pairs))
ggsave(pl2, filename = file.path(resultsfolder, "violin_Fiber.pdf"), device = "pdf", width = 6, height = 6)

(pl3 <- ggplot(df_diet, aes(x = fct_reorder(EthnicityTot, Protein, median, na.rm = TRUE), y = Protein)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Protein") +
    sig_comparisons(df_diet, "Protein", eth_pairs))
ggsave(pl3, filename = file.path(resultsfolder, "violin_Protein.pdf"), device = "pdf", width = 6, height = 6)

(pl4 <- ggplot(df_diet, aes(x = fct_reorder(EthnicityTot, Protein_animal, median, na.rm = TRUE), y = Protein_animal)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Animal protein") +
    sig_comparisons(df_diet, "Protein_animal", eth_pairs))
ggsave(pl4, filename = file.path(resultsfolder, "violin_Protein_animal.pdf"), device = "pdf", width = 6, height = 6)

(pl5 <- ggplot(df_diet, aes(x = fct_reorder(EthnicityTot, FattyAcids, median, na.rm = TRUE), y = FattyAcids)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Fatty acids") +
    sig_comparisons(df_diet, "FattyAcids", eth_pairs))
ggsave(pl5, filename = file.path(resultsfolder, "violin_FattyAcids.pdf"), device = "pdf", width = 6, height = 6)

(pl6 <- ggplot(df_diet, aes(x = fct_reorder(EthnicityTot, Carbohydrates, median, na.rm = TRUE), y = Carbohydrates)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Carbohydrates") +
    sig_comparisons(df_diet, "Carbohydrates", eth_pairs))
ggsave(pl6, filename = file.path(resultsfolder, "violin_Carbohydrates.pdf"), device = "pdf", width = 6, height = 6)

(pl7 <- ggplot(df_diet, aes(x = fct_reorder(EthnicityTot, Sodium_g, median, na.rm = TRUE), y = Sodium_g)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Sodium") +
    sig_comparisons(df_diet, "Sodium_g", eth_pairs))
ggsave(pl7, filename = file.path(resultsfolder, "violin_Sodium.pdf"), device = "pdf", width = 6, height = 6)

(fig_macronutrients <- ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, pl7,
                                  ncol = 3, nrow = 3))
ggsave(fig_macronutrients, filename = file.path(resultsfolder, "macronutrients_by_ethnicity.pdf"),
       device = "pdf", width = 15, height = 18)

#### Calorie-normalized macronutrients ####
df_diet_norm <- df_diet %>%
    filter(!is.na(TotalCalories) & TotalCalories > 0) %>%
    mutate(
        Protein_per1000        = Protein / TotalCalories * 1000,
        Protein_animal_per1000 = Protein_animal / TotalCalories * 1000,
        FattyAcids_per1000     = FattyAcids / TotalCalories * 1000,
        Carbohydrates_per1000  = Carbohydrates / TotalCalories * 1000,
        Fiber_per1000          = Fiber / TotalCalories * 1000,
        Sodium_per1000         = Sodium_g / TotalCalories * 1000
    )

resultsfolder_norm <- file.path(resultsfolder, "calorie_normalized")
dir.create(resultsfolder_norm, showWarnings = FALSE)

(pln1 <- ggplot(df_diet_norm, aes(x = fct_reorder(EthnicityTot, Protein_per1000, median, na.rm = TRUE), y = Protein_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Protein (per 1000 kcal)") +
    sig_comparisons(df_diet_norm, "Protein_per1000", eth_pairs))
ggsave(pln1, filename = file.path(resultsfolder_norm, "violin_norm_Protein.pdf"), device = "pdf", width = 6, height = 6)

(pln2 <- ggplot(df_diet_norm, aes(x = fct_reorder(EthnicityTot, Protein_animal_per1000, median, na.rm = TRUE), y = Protein_animal_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Animal protein (per 1000 kcal)") +
    sig_comparisons(df_diet_norm, "Protein_animal_per1000", eth_pairs))
ggsave(pln2, filename = file.path(resultsfolder_norm, "violin_norm_Protein_animal.pdf"), device = "pdf", width = 6, height = 6)

(pln3 <- ggplot(df_diet_norm, aes(x = fct_reorder(EthnicityTot, FattyAcids_per1000, median, na.rm = TRUE), y = FattyAcids_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Fatty acids (per 1000 kcal)") +
    sig_comparisons(df_diet_norm, "FattyAcids_per1000", eth_pairs))
ggsave(pln3, filename = file.path(resultsfolder_norm, "violin_norm_FattyAcids.pdf"), device = "pdf", width = 6, height = 6)

(pln4 <- ggplot(df_diet_norm, aes(x = fct_reorder(EthnicityTot, Carbohydrates_per1000, median, na.rm = TRUE), y = Carbohydrates_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Carbohydrates (per 1000 kcal)") +
    sig_comparisons(df_diet_norm, "Carbohydrates_per1000", eth_pairs))
ggsave(pln4, filename = file.path(resultsfolder_norm, "violin_norm_Carbohydrates.pdf"), device = "pdf", width = 6, height = 6)

(pln5 <- ggplot(df_diet_norm, aes(x = fct_reorder(EthnicityTot, Fiber_per1000, median, na.rm = TRUE), y = Fiber_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Fiber (per 1000 kcal)") +
    sig_comparisons(df_diet_norm, "Fiber_per1000", eth_pairs))
ggsave(pln5, filename = file.path(resultsfolder_norm, "violin_norm_Fiber.pdf"), device = "pdf", width = 6, height = 6)

(pln6 <- ggplot(df_diet_norm, aes(x = fct_reorder(EthnicityTot, Sodium_per1000, median, na.rm = TRUE), y = Sodium_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Sodium (per 1000 kcal)") +
    sig_comparisons(df_diet_norm, "Sodium_per1000", eth_pairs))
ggsave(pln6, filename = file.path(resultsfolder_norm, "violin_norm_Sodium.pdf"), device = "pdf", width = 6, height = 6)

(fig_macronutrients_norm <- ggarrange(pln1, pln2, pln3, pln4, pln5, pln6,
                                       ncol = 3, nrow = 2))
ggsave(fig_macronutrients_norm, filename = file.path(resultsfolder_norm, "macronutrients_norm_by_ethnicity.pdf"),
       device = "pdf", width = 15, height = 12)

#### Shotgun subset - macronutrients ####
shotids <- read.csv('data/shotgun/shotgunseq_ids.csv') %>%
    dplyr::select(ID = x) %>%
    mutate(ID = str_c("S", ID))

df_diet_sg <- df_diet %>% filter(ID %in% shotids$ID) %>% droplevels()

eth_present_sg <- intersect(names(eth_colors), unique(as.character(df_diet_sg$EthnicityTot)))
eth_pairs_sg   <- combn(eth_present_sg, 2, simplify = FALSE)

resultsfolder_sg <- file.path(resultsfolder, "shotgun_subset")
dir.create(resultsfolder_sg, showWarnings = FALSE)

(sg_pl1 <- ggplot(df_diet_sg, aes(x = fct_reorder(EthnicityTot, TotalCalories, median, na.rm = TRUE), y = TotalCalories)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'Kcal', title = "Total calories") +
    sig_comparisons(df_diet_sg, "TotalCalories", eth_pairs_sg))
ggsave(sg_pl1, filename = file.path(resultsfolder_sg, "violin_TotalCalories.pdf"), device = "pdf", width = 6, height = 6)

(sg_pl2 <- ggplot(df_diet_sg, aes(x = fct_reorder(EthnicityTot, Fiber, median, na.rm = TRUE), y = Fiber)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Fibers") +
    sig_comparisons(df_diet_sg, "Fiber", eth_pairs_sg))
ggsave(sg_pl2, filename = file.path(resultsfolder_sg, "violin_Fiber.pdf"), device = "pdf", width = 6, height = 6)

(sg_pl3 <- ggplot(df_diet_sg, aes(x = fct_reorder(EthnicityTot, Protein, median, na.rm = TRUE), y = Protein)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Protein") +
    sig_comparisons(df_diet_sg, "Protein", eth_pairs_sg))
ggsave(sg_pl3, filename = file.path(resultsfolder_sg, "violin_Protein.pdf"), device = "pdf", width = 6, height = 6)

(sg_pl4 <- ggplot(df_diet_sg, aes(x = fct_reorder(EthnicityTot, Protein_animal, median, na.rm = TRUE), y = Protein_animal)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Animal protein") +
    sig_comparisons(df_diet_sg, "Protein_animal", eth_pairs_sg))
ggsave(sg_pl4, filename = file.path(resultsfolder_sg, "violin_Protein_animal.pdf"), device = "pdf", width = 6, height = 6)

(sg_pl5 <- ggplot(df_diet_sg, aes(x = fct_reorder(EthnicityTot, FattyAcids, median, na.rm = TRUE), y = FattyAcids)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Fatty acids") +
    sig_comparisons(df_diet_sg, "FattyAcids", eth_pairs_sg))
ggsave(sg_pl5, filename = file.path(resultsfolder_sg, "violin_FattyAcids.pdf"), device = "pdf", width = 6, height = 6)

(sg_pl6 <- ggplot(df_diet_sg, aes(x = fct_reorder(EthnicityTot, Carbohydrates, median, na.rm = TRUE), y = Carbohydrates)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Carbohydrates") +
    sig_comparisons(df_diet_sg, "Carbohydrates", eth_pairs_sg))
ggsave(sg_pl6, filename = file.path(resultsfolder_sg, "violin_Carbohydrates.pdf"), device = "pdf", width = 6, height = 6)

(sg_pl7 <- ggplot(df_diet_sg, aes(x = fct_reorder(EthnicityTot, Sodium_g, median, na.rm = TRUE), y = Sodium_g)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'gram', title = "Sodium") +
    sig_comparisons(df_diet_sg, "Sodium_g", eth_pairs_sg))
ggsave(sg_pl7, filename = file.path(resultsfolder_sg, "violin_Sodium.pdf"), device = "pdf", width = 6, height = 6)

(fig_macronutrients_sg <- ggarrange(sg_pl1, sg_pl2, sg_pl3, sg_pl4, sg_pl5, sg_pl6, sg_pl7,
                                     ncol = 3, nrow = 3))
ggsave(fig_macronutrients_sg, filename = file.path(resultsfolder_sg, "macronutrients_by_ethnicity.pdf"),
       device = "pdf", width = 15, height = 18)

#### Shotgun subset - calorie-normalized macronutrients ####
df_diet_norm_sg <- df_diet_norm %>% filter(ID %in% shotids$ID) %>% droplevels()

eth_present_norm_sg <- intersect(names(eth_colors), unique(as.character(df_diet_norm_sg$EthnicityTot)))
eth_pairs_norm_sg   <- combn(eth_present_norm_sg, 2, simplify = FALSE)

resultsfolder_norm_sg <- file.path(resultsfolder_sg, "calorie_normalized")
dir.create(resultsfolder_norm_sg, showWarnings = FALSE)

(sg_pln1 <- ggplot(df_diet_norm_sg, aes(x = fct_reorder(EthnicityTot, Protein_per1000, median, na.rm = TRUE), y = Protein_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Protein (per 1000 kcal)") +
    sig_comparisons(df_diet_norm_sg, "Protein_per1000", eth_pairs_norm_sg))
ggsave(sg_pln1, filename = file.path(resultsfolder_norm_sg, "violin_norm_Protein.pdf"), device = "pdf", width = 6, height = 6)

(sg_pln2 <- ggplot(df_diet_norm_sg, aes(x = fct_reorder(EthnicityTot, Protein_animal_per1000, median, na.rm = TRUE), y = Protein_animal_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Animal protein (per 1000 kcal)") +
    sig_comparisons(df_diet_norm_sg, "Protein_animal_per1000", eth_pairs_norm_sg))
ggsave(sg_pln2, filename = file.path(resultsfolder_norm_sg, "violin_norm_Protein_animal.pdf"), device = "pdf", width = 6, height = 6)

(sg_pln3 <- ggplot(df_diet_norm_sg, aes(x = fct_reorder(EthnicityTot, FattyAcids_per1000, median, na.rm = TRUE), y = FattyAcids_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Fatty acids (per 1000 kcal)") +
    sig_comparisons(df_diet_norm_sg, "FattyAcids_per1000", eth_pairs_norm_sg))
ggsave(sg_pln3, filename = file.path(resultsfolder_norm_sg, "violin_norm_FattyAcids.pdf"), device = "pdf", width = 6, height = 6)

(sg_pln4 <- ggplot(df_diet_norm_sg, aes(x = fct_reorder(EthnicityTot, Carbohydrates_per1000, median, na.rm = TRUE), y = Carbohydrates_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Carbohydrates (per 1000 kcal)") +
    sig_comparisons(df_diet_norm_sg, "Carbohydrates_per1000", eth_pairs_norm_sg))
ggsave(sg_pln4, filename = file.path(resultsfolder_norm_sg, "violin_norm_Carbohydrates.pdf"), device = "pdf", width = 6, height = 6)

(sg_pln5 <- ggplot(df_diet_norm_sg, aes(x = fct_reorder(EthnicityTot, Fiber_per1000, median, na.rm = TRUE), y = Fiber_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Fiber (per 1000 kcal)") +
    sig_comparisons(df_diet_norm_sg, "Fiber_per1000", eth_pairs_norm_sg))
ggsave(sg_pln5, filename = file.path(resultsfolder_norm_sg, "violin_norm_Fiber.pdf"), device = "pdf", width = 6, height = 6)

(sg_pln6 <- ggplot(df_diet_norm_sg, aes(x = fct_reorder(EthnicityTot, Sodium_per1000, median, na.rm = TRUE), y = Sodium_per1000)) +
    geom_violin(aes(fill = EthnicityTot), color = NA) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = 'g / 1000 kcal', title = "Sodium (per 1000 kcal)") +
    sig_comparisons(df_diet_norm_sg, "Sodium_per1000", eth_pairs_norm_sg))
ggsave(sg_pln6, filename = file.path(resultsfolder_norm_sg, "violin_norm_Sodium.pdf"), device = "pdf", width = 6, height = 6)

(fig_macronutrients_norm_sg <- ggarrange(sg_pln1, sg_pln2, sg_pln3, sg_pln4, sg_pln5, sg_pln6,
                                          ncol = 3, nrow = 2))
ggsave(fig_macronutrients_norm_sg, filename = file.path(resultsfolder_norm_sg, "macronutrients_norm_by_ethnicity.pdf"),
       device = "pdf", width = 15, height = 12)

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
