## Compositional plots and diversity metrics

## Libraries
library(tidyverse)
library(ggpubr)
library(ggsci)
library(mixOmics)
library(aplot)
library(ggplotify)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
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
                strip.text = element_text(face="bold")
        ))
}

## Data
ba <- rio::import("data/240605_BileAcids_HELIUS.xlsx")
ba$IDlong <- case_when(
    ba$timepoint == "helius covid" ~ str_c("HELICOV_", ba$ID),
    .default = str_c("HELIFU_", ba$ID)
)
ba$ID <- str_c("S", ba$ID)
ba$timepoint <- case_when(
    ba$timepoint == "helius covid" ~ "COV",
    .default = "H2"
)
ba <- ba %>% dplyr::select(!any_of(c("TLCA", "TUDCA")))
bileacids <- names(ba)[6:ncol(ba)-2]
primba <- c("CA", "DCA")
secba <- bileacids[which(!bileacids %in% primba)]
clin <- readRDS("data/clinicaldata_wide.RDS")
names(clin)
head(ba)
clinfull <- left_join(ba, clin, by = "ID") %>% filter(!is.na(timepoint))
clinfull$ID <- as.factor(clinfull$ID)

plist <- list()
for(a in bileacids){
    clinfull$var <- clinfull[[a]]
    clinfull$var <- case_when(is.na(clinfull$var) ~ 0, .default = clinfull$var)
    clinfull$var <- as.numeric(clinfull$var)
    ps <- min(clinfull$var[which(clinfull$var != 0)]) /2
    pl <- ggplot(data = clinfull, aes(x = timepoint, y = log10(var+ps))) +
        geom_line(aes(group = ID), alpha = 0.5, color = "grey40") +
        geom_jitter(aes(group = ID, color = Ethnicity), alpha = 0.5, width = 0) +
        gghalves::geom_half_violin(data = clinfull %>% filter(timepoint == "COV"), 
                                    aes(x = timepoint, y = log10(var+ps), fill = Ethnicity), 
                                    nudge = 0.05, side = "l", color = "black", 
                                    width = 0.25) +
        gghalves::geom_half_violin(data = clinfull %>% filter(timepoint == "H2"), 
                                    aes(x = timepoint, y = log10(var+ps), fill = Ethnicity), 
                                    nudge = 0.05, side = "r", color = "black", 
                                    width = 0.25) +
        gghalves::geom_half_boxplot(data = clinfull %>% filter(timepoint == "COV"), 
                                    aes(x = timepoint, y = log10(var+ps)), 
                                    nudge = 0.05, side = "l", fill = "white", color = "black", 
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        gghalves::geom_half_boxplot(data = clinfull %>% filter(timepoint == "H2"), 
                                    aes(x = timepoint, y = log10(var+ps)), 
                                    nudge = 0.05, side = "r", fill = "white", color = "black", 
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        # stat_pvalue_manual(res_lmm, y.position = metmax, label = "{sig}",
        #                    remove.bracket = FALSE) +
        stat_compare_means(aes(x = timepoint, y = log(var)), 
                           label = "p.signif", tip.length = 0, paired = TRUE,
                           comparisons = list(c("COV", "H2")))+
        scale_color_simpsons(guide = "none") +
        scale_fill_simpsons(guide = "none") +
        # coord_cartesian(ylim = c(0,metmax)) +
        theme_Publication() +
        labs(x = "Timepoint", y = paste0("log10(", a, ")"), title = paste0(a),
             color = "") +
        facet_wrap(~Ethnicity)
    # print(pl)
    plist[[a]] <- pl
}
ggarrange(plotlist = plist, labels = LETTERS[1:length(plist)],
          nrow = 4, ncol = 4)
ggsave("results/bileacidconc_time_ethn.pdf", width = 22, height = 18)

plist2 <- list()
for(a in bileacids){
    clinfull$var <- clinfull[[a]]
    clinfull$var <- case_when(is.na(clinfull$var) ~ 0, .default = clinfull$var)
    clinfull$var <- as.numeric(clinfull$var)
    ps <- min(clinfull$var[which(clinfull$var != 0)]) /2
    pl <- ggplot(data = clinfull, aes(x = Ethnicity, y = log10(var+ps))) +
        geom_violin(aes(fill = Ethnicity)) +
        geom_boxplot(width = 0.2, fill = "white") +
        stat_compare_means(aes(x = Ethnicity, y = log(var)), 
                           label = "p.signif", tip.length = 0, paired = FALSE,
                           comparisons = list(c("Dutch", "Surinamese")))+
        scale_color_simpsons(guide = "none") +
        scale_fill_simpsons(guide = "none") +
        # coord_cartesian(ylim = c(0,metmax)) +
        theme_Publication() +
        labs(x = "Ethnicity", y = paste0("log10(", a, ")"), title = paste0(a),
             color = "") +
        facet_wrap(~timepoint)
    # print(pl)
    plist2[[a]] <- pl
}
ggarrange(plotlist = plist2, labels = LETTERS[1:length(plist2)],
          nrow = 4, ncol = 4)
ggsave("results/bileacidconc_timepoints.pdf", width = 22, height = 22)

## PCA diet
df_ba <- clinfull %>% dplyr::select(ID, Ethnicity, timepoint, 4:19)
df_ba2 <- df_ba %>% dplyr::select(-ID, -Ethnicity, -timepoint) %>% 
    mutate(across(everything(.), ~log10(.x+0.001))) %>% 
    mutate(across(everything(.), scale))
matba <- as.matrix(df_ba2)
tunediet <- mixOmics::tune.pca(matba, ncomp = 5, scale = TRUE)
plot(tunediet)
pc <- mixOmics::pca(matba, ncomp = 2)
pcs <- as.data.frame(pc$variates$X)
pcs <- pcs %>% mutate(ID = df_ba$ID, timepoint = df_ba$timepoint, Ethnicity = df_ba$Ethnicity)
expvar_ba <- pc$prop_expl_var$X[1:2]
loadings <- as.data.frame(pc$loadings$X)
loadings$Variables <- rownames(loadings)
pc2highest <- loadings %>% arrange(-PC2) %>% slice(c(1:2, 16)) %>% dplyr::select(Variables)
pc1highest <- loadings %>% arrange(-PC1) %>% slice(1:3) %>% dplyr::select(Variables)
pchighest <- c(pc2highest$Variables, pc1highest$Variables)
loadings <- loadings %>% filter(Variables %in% pchighest)
xco <- loadings$PC1
yco <- loadings$PC2
lab <- loadings$Variables

(pcaba <- ggplot(data = pcs, aes(PC1, PC2)) +
            geom_point(aes(color = timepoint), size = 1, alpha = 1.0) +
            xlab(paste0('PC1 (', round(expvar_ba[1]*100, digits = 1),'%)')) +
            ylab(paste0('PC2 (', round(expvar_ba[2]*100, digits = 1),'%)')) +
            theme_Publication() +
            stat_ellipse(geom = "polygon", aes(color = timepoint, fill = timepoint), linewidth = 1.0,
                         alpha = 0.1, type = "norm") +
            scale_color_simpsons() +
            scale_fill_simpsons(guide = "none") +
            # labs(color = "", title = "PCA bile acids")+
            geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*10), yend = (PC2*10)),
                         arrow = arrow(length = unit(1/2, "picas")),
                         color = "black", linewidth = 0.9) +
            facet_wrap(~Ethnicity) + 
            geom_text(data = loadings, aes(x = PC1*12, y = PC2*12, label = Variables))
)

(pcaba <- ggplot(data = pcs, aes(PC1, PC2)) +
        geom_point(aes(color = timepoint), size = 1, alpha = 1.0) +
        xlab(paste0('PC1 (', round(expvar_ba[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_ba[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = timepoint, fill = timepoint), linewidth = 1.0, alpha = 0.1, type = "norm") +
        scale_color_simpsons() +
        scale_fill_simpsons(guide = "none") +
        # labs(color = "", title = "PCA bile acids")+
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*10), yend = (PC2*10)),
                     arrow = arrow(length = unit(1/2, "picas")),
                     color = "black", linewidth = 0.9) +
        geom_text(data = loadings, aes(x = PC1*12, y = PC2*12, label = Variables))
)

df_bafilt <- df_ba %>% filter(timepoint != "COV")
df_ba3 <- df_bafilt %>% 
    dplyr::select(-ID, -Ethnicity, -timepoint) %>% 
    mutate(across(everything(.), ~log10(.x+0.001))) %>% 
    mutate(across(everything(.), scale))
matba <- as.matrix(df_ba3)
tunediet <- mixOmics::tune.pca(matba, ncomp = 4, scale = TRUE)
plot(tunediet)
pc <- mixOmics::pca(matba, ncomp = 2)
pcs <- as.data.frame(pc$variates$X)
pcs <- pcs %>% mutate(ID = df_bafilt$ID, timepoint = df_bafilt$timepoint, 
                      Ethnicity = df_bafilt$Ethnicity, 
                      outlier = case_when(
                          PC2 > 3 | PC2 < -4 | PC1 > 5 ~ ID,
                          .default = NA
                      ))
expvar_ba <- pc$prop_expl_var$X[1:2]
loadings <- as.data.frame(pc$loadings$X)
loadings$Variables <- rownames(loadings)
pc2highest <- loadings %>% arrange(-PC2) %>% slice(c(1, 15:16)) %>% 
    dplyr::select(Variables)
pc1highest <- loadings %>% arrange(-PC1) %>% slice(1:3) %>% 
    dplyr::select(Variables)
pchighest <- c(pc2highest$Variables, pc1highest$Variables)
loadings <- loadings %>% filter(Variables %in% pchighest)

(pcaba <- ggplot(data = pcs, aes(PC1, PC2)) +
        geom_point(aes(color = Ethnicity), size = 1, alpha = 1.0) +
        xlab(paste0('PC1 (', round(expvar_ba[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_ba[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = Ethnicity, fill = Ethnicity), 
                     linewidth = 1.0, alpha = 0.1, type = "norm") +
        scale_color_simpsons() +
        scale_fill_simpsons(guide = "none") +
        # labs(color = "", title = "PCA bile acids")+
        geom_label(aes(label = outlier), nudge_x = 0, nudge_y = 0.8) +
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*6), yend = (PC2*6)),
                     arrow = arrow(length = unit(1/2, "picas")),
                     color = "black", linewidth = 0.9) +
        geom_text(data = loadings, aes(x = PC1*8, y = PC2*8, label = Variables))
)

set.seed(124)
pcs2 <- pcs %>% filter(!is.na(outlier) | ID %in% sample(ID, size = 7))
(pcaba2 <- ggplot(data = pcs2, aes(PC1, PC2)) +
        geom_point(aes(color = Ethnicity), size = 2, alpha = 1.0) +
        xlab(paste0('PC1 (', round(expvar_ba[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_ba[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        scale_color_manual(values = pal_simpsons()(5)[c(1:3,5)]) +
        scale_fill_simpsons(guide = "none") +
        # labs(color = "", title = "PCA bile acids")+
        geom_label(aes(label = ID, color = !is.na(outlier)), 
                   nudge_x = -0.1, nudge_y = 0.8)
)

df_ba2 <- clinfull %>% dplyr::select(ID, Ethnicity, timepoint, 4:19) %>% 
    filter(ID %in% pcs2$ID & timepoint == "H2") %>% 
    replace(is.na(.), 0) %>%
    mutate(outlier = pcs2$outlier[match(pcs2$ID, ID)],
            outlier = case_when(!is.na(outlier) ~ TRUE, .default = FALSE),
           # sumba = rowSums(across(where(is.numeric))),
           # across(c(4:19), ~.x / sumba)
           )
rowSums(df_ba2[,4:19])
df_balong <- df_ba2 %>% pivot_longer(., cols = 4:19, names_to = "bileacid", 
                                     values_to = "conc") %>% 
    arrange(outlier) %>% 
    mutate(ID = fct_inorder(ID),
           bileacid = fct_reorder(bileacid, conc))

bileplot <- ggplot(df_balong, aes(fill=bileacid, y=conc, x=ID)) + 
                geom_bar(position="fill", stat="identity") +
                scale_fill_simpsons() +
                theme_Publication() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      legend.position = "right")
plbottom <- ggplot(df_balong, aes(x = ID, y = 1, fill = outlier)) +
        geom_tile() +
        scale_fill_manual(values = c("grey70", pal_simpsons()(5)[5])) +
        theme_transparent()
approp <- bileplot %>% insert_top(plbottom, height=.1)    
approp

bileplot <- ggplot(df_balong, aes(fill=bileacid, y=conc, x=ID)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_simpsons() +
    labs(y = "concentration") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
plbottom <- ggplot(df_balong, aes(x = ID, y = 1, fill = outlier)) +
        geom_tile() +
    scale_fill_manual(values = c("grey70", pal_simpsons()(5)[5])) +
        theme_transparent()
apstack <- bileplot %>% insert_top(plbottom, height=.1)    
apstack
ggarrange(approp, apstack)
ggarrange(pcaba, pcaba2, as.ggplot(approp), as.ggplot(apstack),
          nrow = 2, ncol = 2, 
          labels = LETTERS[1:4])
ggsave("results/bileacids_H2_outliers.pdf", width = 14, height = 10)
