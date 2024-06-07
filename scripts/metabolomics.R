## Cleaning metabolomics data

library(rio)
library(haven)
library(tidyverse)
library(tableone)
library(ggsci)

theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                #family = 'Helvetica'
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.line.y = element_line(colour="black"),
                axis.ticks.x = element_line(),
                axis.ticks.y = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(5,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(face = "italic", size=rel(0.6))
        ))
} 

## Data
helius <- readRDS('data/clinicaldata_long.RDS')
infomet <- rio::import('data/metabolomics/Info_HELIUS_metabolomics.xlsx')
colnames(infomet)[1] <- "MetabolonID"
met <- rio::import('data/metabolomics/HELIUS_metabolomics_abundance.xlsx')
colnames(met)[1] <- "sampleID"
colnames(met)[2:ncol(met)] <- infomet$CHEMICAL_NAME[match(colnames(met)[2:ncol(met)],infomet$MetabolonID)]
meta <- rio::import("data/metabolomics/HELIUS_metabolomics_metadata.xlsx")
colnames(meta)[1] <- "sampleID"
head(meta$sampleID)

meta <- meta %>% mutate(
    subjectID = str_replace_all(CLIENT_IDENTIFIER, "Helius ", "HELIUS_"),
    subjectID = case_when(
        str_detect(CLIENT_IDENTIFIER, "Covid ") ~ str_replace(subjectID, "HELIUS", "HELICOV"),
        .default = str_replace(subjectID, "HELIUS", "HELIFU")
    ),
    subjectID = str_replace_all(subjectID, "Covid ", ""),
    timepoint = case_when(
        str_detect(subjectID, "HELICOV") ~ "baseline",
        str_detect(subjectID, "HELIFU") ~ "follow-up"
    ),
    timepoint = as.factor(timepoint),
    across(c("NEG", "POLAR", "POS.EARLY", "POS.LATE"), as.factor)
)

write.csv2(meta$subjectID, 'data/metabolomics/ids_metabolomics.csv')

summary(meta$NEG)
summary(meta$POLAR)
summary(meta$POS.EARLY)
summary(meta$POS.LATE)

dim(met)
met$subjectID <- meta$subjectID[match(met$sampleID, meta$sampleID)]
rownames(met) <- met$subjectID
met$subjectID <- NULL
met$sampleID <- NULL
metmatrix <- as.matrix(met)
constants <- apply(metmatrix, 2, var)
any(constants == 0)
print('In total there are this number of constants in the data:')
length(constants[which(constants == 0)])
nameconstants <- names(constants[which(constants == 0)])
write.csv(nameconstants, "results/constant_metabolites.csv")
const <- metmatrix[, constants == 0]
metmatrix <- metmatrix[,constants != 0]
dim(metmatrix) # so 75 of 1468 lost

infomet <- infomet %>% filter(!CHEMICAL_NAME %in% nameconstants)

namespolar <- infomet$CHEMICAL_NAME[which(infomet$PLATFORM == "Polar")]
namesposearly <- infomet$CHEMICAL_NAME[which(infomet$PLATFORM == "Pos Early")]
namesposlate <- infomet$CHEMICAL_NAME[which(infomet$PLATFORM == "Pos Late")]
namesneg <- infomet$CHEMICAL_NAME[which(infomet$PLATFORM == "Neg")]

length(namespolar)
length(namesposlate)
length(namesposearly)
length(namesneg)

metpolar <- metmatrix[,namespolar]
metposearly <- metmatrix[,namesposearly]
metposlate <- metmatrix[,namesposlate]
metneg <- metmatrix[,namesneg]

# PCA metabolites and timepoints
pcatimepoints <- function(matrix, helius, plotname){
    pca <- prcomp(matrix, scale. = T)
    df <- as.data.frame(pca$x[, 1:2])
    expvar1 <- format(round(summary(pca)['importance'][[1]][2,1]*100, digits = 1), nsmall = 1)
    expvar2 <- format(round(summary(pca)['importance'][[1]][2,2]*100, digits = 1), nsmall = 1)
    head(df)
    df$ID <- rownames(df)
    df$timepoint <- case_when(
        str_detect(df$ID, "HELICOV") ~ "COVID",
        str_detect(df$ID, "HELIFU") ~ "follow-up"
    )
    
    df$posearly <- meta$POS.EARLY[match(df$ID, meta$subjectID)]
    df$poslate <- meta$POS.LATE[match(df$ID, meta$subjectID)]
    df$neg <- meta$NEG[match(df$ID, meta$subjectID)]
    df$polar <- meta$POLAR[match(df$ID, meta$subjectID)]
    
    (pl <- ggplot(df, aes(x=PC1, y=PC2, color=timepoint)) +
            geom_point()+
            ggtitle(str_c('PCA metabolites - ', plotname)) +
            scale_color_lancet() +
            theme_minimal() +
            stat_ellipse() +
            labs(x=str_c('PC1 ', expvar1, '%'), y=str_c('PC2 ', expvar2, '%')) +
            theme_Publication())
    ggsave(pl, filename = str_c('results/pca_', str_replace(plotname, " ", ""), '.pdf'), width =5, height = 5)
    
    platforms <- infomet$PLATFORM[match(colnames(matrix), infomet$CHEMICAL_NAME)]
    platforms <- platforms[!is.na(platforms)]
    
    if(all(platforms == "Pos Early")){
        print("pos early")
        (pl2 <- ggplot(df, aes(x=PC1, y=PC2, color=posearly)) +
             geom_point()+
             ggtitle(str_c('PCA metabolites (selection) - positive early')) +
             scale_color_lancet() +
             theme_minimal() +
             stat_ellipse() +
             labs(x=str_c('PC1 ', expvar1, '%'), y=str_c('PC2 ', expvar2, '%')) +
             theme_Publication())
        ggsave(pl2, filename = str_c('results/pca_onlyposearly.pdf'), width = 7, height = 5)
    } else if(all(platforms == "Pos Late")){
        print("pos late")
        (pl2 <- ggplot(df, aes(x=PC1, y=PC2, color=poslate)) +
             geom_point()+
             ggtitle(str_c('PCA metabolites (selection) - positive late')) +
             scale_color_lancet() +
             theme_minimal() +
             stat_ellipse() +
             labs(x=str_c('PC1 ', expvar1, '%'), y=str_c('PC2 ', expvar2, '%')) +
             theme_Publication())
        ggsave(pl2, filename = str_c('results/pca_onlyposlate.pdf'), width = 7, height = 5)
    }  else if(all(platforms == "Neg")){
        print("neg")
        (pl2 <- ggplot(df, aes(x=PC1, y=PC2, color=neg)) +
             geom_point()+
             ggtitle(str_c('PCA metabolites (selection) - negative')) +
             scale_color_lancet() +
             theme_minimal() +
             stat_ellipse() +
             labs(x=str_c('PC1 ', expvar1, '%'), y=str_c('PC2 ', expvar2, '%')) +
             theme_Publication())
        ggsave(pl2, filename = str_c('results/pca_onlyneg.pdf'), width = 7, height = 5)
    } else if(all(platforms == "Polar")){
        print("polar")
        (pl2 <- ggplot(df, aes(x=PC1, y=PC2, color=polar)) +
             geom_point()+
             ggtitle(str_c('PCA metabolites (selection) - polar')) +
             scale_color_lancet() +
             theme_minimal() +
             stat_ellipse() +
             labs(x=str_c('PC1 ', expvar1, '%'), y=str_c('PC2 ', expvar2, '%')) +
             theme_Publication())
        ggsave(pl2, filename = str_c('results/pca_onlypolar.pdf'), width = 7, height = 5)
        return(pl2)
    }
}

pcatimepoints(metneg, helius, 'platform negative')
pcatimepoints(metposearly, helius, 'platform positive early')
pcatimepoints(metposlate, helius, 'platform positive late')
pcatimepoints(metpolar, helius, 'platform polar')
pcatimepoints(metmatrix, helius, 'all')

######################## Metabolite info ########################

summary(infomet$SUPER_PATHWAY)
infomet <- infomet %>% select(MetabolonID, sup = SUPER_PATHWAY, sub = SUB_PATHWAY,
                              metname = PLOT_NAME, platform = PLATFORM) %>% 
    mutate(sup = case_when(
        is.na(sup) ~ "Unknown pathway",
        .default = sup
            )
)
            

colpal <- c("Lipid"="#00468BFF", "Amino Acid"="#ED0000FF", "Nucleotide"="#42B540FF",
            "Cofactors and Vitamins"="#0099B4FF", "Carbohydrate"="#925E9FFF", "Peptide"="#FDAF91FF", "Partially Characterized Molecules"="#AD002AFF", "Energy"="gold", "Xenobiotics"="#1B1919FF",
            "Unknown pathway" = "darkgrey")

infomet %>% 
    group_by(sup) %>% 
    summarize(count=n()) %>% 
    arrange(desc(count)) %>% 
    ggplot(.) +
        geom_histogram(aes(x=reorder(sup, count), y=count, fill = sup), stat = "identity") +
        scale_fill_manual(values = colpal) +
        labs(title="Number of metabolites per group", x = element_blank(), y='Count') +
        coord_flip() +
        theme_Publication() +
        theme(legend.position = "none")
ggsave("results/info_metabolites.pdf", device = "pdf", width = 8, height = 6)
ggsave("results/info_metabolites.svg", device = "svg", width = 8, height = 6)

infomet %>% 
    group_by(sup, platform) %>% 
    summarize(count=n()) %>% 
    arrange(desc(count)) %>% 
    ggplot(.) +
    geom_histogram(aes(x=reorder(sup, count), y=count, fill = sup), stat = "identity") +
    scale_fill_manual(values = colpal) +
    labs(title="Number of metabolites per group", x = element_blank(), y='Count') +
    coord_flip() +
    facet_wrap(~platform) +
    theme_Publication() +
    theme(legend.position = "none")
ggsave("results/platform_info_metabolites.pdf", device = "pdf", width = 8, height = 6)
ggsave("results/platform_info_metabolites.svg", device = "svg", width = 8, height = 6)

(pl_noxeno <- infomet %>% 
        group_by(sup) %>% 
        summarize(count=n()) %>% 
        arrange(desc(count)) %>% 
        filter(sup!="Xenobiotics") %>% 
        ggplot(.) +
        geom_histogram(aes(x=reorder(sup, count), y=count, fill = sup), stat = "identity") +
        scale_fill_manual(values = colpal) +
        labs(title="Groups of metabolites", x = element_blank(), y='Count') +
        coord_flip() +
        theme_Publication() +
        theme(legend.position = "none", axis.line.y = element_blank(),
              axis.ticks.y = element_blank()))
ggsave("metabolites_noxeno.pdf", device = "pdf", width = 8, height = 6)
ggsave("metabolites_noxeno.svg", device = "svg", width = 8, height = 6)

