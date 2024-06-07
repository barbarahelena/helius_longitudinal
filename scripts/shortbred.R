## ShortBRED analyses

## Libraries
library(tidyverse)
library(ggsci)
library(mixOmics)

## Theme
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

## Data
df <- read.csv2("data/shotgun/merged_shortbred.csv", row.names = "X")
dim(df)
df <- df %>% mutate_all(as.numeric)

## Filter out too low prevalent genes
sums <- colSums(df[2:ncol(df)])
sumsdf <- as.data.frame(sums) %>% arrange(-sums)
ggplot(sumsdf, aes(x = 1:nrow(sumsdf), y = sums)) +
    geom_line(color = "firebrick") +
    theme_Publication() +
    scale_y_log10() +
    labs(x = "Antibiotic resistance genes", y = "total sum counts", title = "Prevalence of genes")
ggsave("results/ARG_prevalence.pdf", width = 5, height = 5)
prevalent <- names(sums[which(sums != 0)])
df <- df %>% dplyr::select(all_of(prevalent))
dim(df) # 350 genes left

## Clean the IDs
extract_info <- function(name) {
    # Extract the organism name
    gene_name <- sub("^gb\\.[^.]+\\.ARO_\\d+\\.", "", name)
    accession_id <- str_match(name, "(ARO_\\d+)")[,2]
    new_name <- paste0(gene_name, " (", accession_id, ")")
    return(new_name)
}
names <- extract_info(colnames(df))
any(duplicated(names))

for (i in which(duplicated(names))) {
    gene_name <- sub("^gb\\.[^.]+\\.ARO_\\d+\\.", "", colnames(df)[i])
    accession_id <- str_match(colnames(df)[i], "(ARO_\\d+)")[,2]
    db_id <- sub("^gb\\.([A-Za-z0-9_]+)\\..*", "\\1", colnames(df)[i])
    names[i] <- paste0(gene_name, " (", accession_id, ")", " (", db_id, ")")
}
any(duplicated(names))

names_with_two_ARO <- names[grep("ARO.*ARO", names)]
print(paste0(length(names_with_two_ARO), " names with more than one ARO number"))
combined_names <- c()
for (i in seq_along(names_with_two_ARO)) {
    first_part <- sub("_gb.*", "", names_with_two_ARO[i])
    print(first_part)
    # Extract ARO numbers
    aro_numbers <- regmatches(names_with_two_ARO[i], gregexpr("ARO_\\d+", names_with_two_ARO[i]))[[1]]
    print(aro_numbers)
    # Combine the first part with ARO numbers
    combined_names[i] <- paste0(first_part, " (", paste(aro_numbers, collapse = ", "), ")")
}
names[grep("ARO.*ARO", names)] <- combined_names
print(names[grep("ARO.*ARO", names)])
colnames(df) <- names

## Make baseline set
rownames(df)
which(duplicated(colnames(df)))
dfbaseline <- df %>% filter(str_detect(rownames(.), "HELIBA"))
dim(dfbaseline)
rownames(dfbaseline) <- str_c("S", str_remove(rownames(dfbaseline), "HELIBA_"))
sums <- colSums(dfbaseline[2:ncol(dfbaseline)])
prevalent <- names(sums[which(sums != 0)])
abundant <- names(sums[which(sums > 5)])
dfbaseline_prev <- dfbaseline %>% dplyr::select(all_of(prevalent))
dim(dfbaseline_prev) # 287 genes left
dfbaseline_ab <- dfbaseline %>% dplyr::select(all_of(abundant)) %>% 
        mutate_all(as.integer)
dim(dfbaseline_ab) # 188 genes
matbase <- t(as.matrix(dfbaseline_ab))
clindf <- readRDS("data/clinicaldata_wide.RDS")
clindf <- clindf %>% filter(ID %in% rownames(dfbaseline_ab)) %>% droplevels(.)
all(clindf$ID == colnames(matbase)) # TRUE

# PCA 
clindf <- readRDS("data/clinicaldata_wide.RDS")
clindf <- clindf %>% filter(ID %in% rownames(dfbaseline_prev)) %>% droplevels(.)
all(clindf$ID == rownames(matbase)) # TRUE
matbase <- log10(matbase+1)
tunebase <- tune.pca(matbase, ncomp = 5, scale = TRUE)
plot(tunebase)
pc <- mixOmics::pca(matbase, ncomp = 2, scale = TRUE)
pcs <- as.data.frame(pc$variates$X)
pcs <- pcs %>% mutate(ID = clindf$ID, Ethnicity = clindf$Ethnicity, Age = clindf$Age_baseline,
                      Diabetes = clindf$DM_baseline)
expvar_base <- pc$prop_expl_var$X[1:2]
loadings <- as.data.frame(pc$loadings$X)
loadings$Variables <- rownames(loadings)
(pcadiet <- pcs %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = Ethnicity), size = 1, alpha = 1.0) +
        xlab(paste0('PC1 (', round(expvar_delta[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_delta[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = Ethnicity, fill = Ethnicity), linewidth = 1.0,
                     alpha = 0.1, type = "norm")+
        scale_color_manual(values = pal_jco()(2)) +
        scale_fill_manual(values = pal_jco()(2), guide = "none") +
        labs(color = "", title = "PCA baseline ARG")
)

(pcadiet <- pcs %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = Diabetes), size = 1, alpha = 1.0) +
        xlab(paste0('PC1 (', round(expvar_delta[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_delta[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = Diabetes, fill = Diabetes), linewidth = 1.0,
                     alpha = 0.1, type = "norm")+
        scale_color_manual(values = pal_jco()(2)) +
        scale_fill_manual(values = pal_jco()(2), guide = "none") +
        labs(color = "", title = "PCA baseline ARG"))

# Maybe some of the participants don't have any ARGs
argsubj <- colSums(matbase)
summary(argsubj == 0) # only 1 participant has none

# Which genes are most prevalent
sumsdf <- as.data.frame(sums) %>% arrange(-sums)
head(sumsdf)
mostprev <- sumsdf %>% filter(sums > 1000)
baseprev <- dfbaseline %>% dplyr::select(all_of(rownames(mostprev))) %>% 
    mutate(ID = rownames(.)) %>% 
    right_join(., clindf, by = "ID")
baseprevlong <- pivot_longer(baseprev, cols = 1:46, names_to = "ARG", 
                             values_to = "counts")
ggplot(baseprevlong, aes(x = Ethnicity, y = log10(counts +1))) +
    geom_violin(aes(fill = Ethnicity)) +
    scale_fill_jco()+
    geom_boxplot(fill = "white", width = 0.25) +
    ggpubr::stat_compare_means() +
    facet_wrap(~ARG) +
    theme_Publication() +
    labs(x = "", y = "counts (log10 scale)")
ggsave("results/ARG_baseline.pdf", width = 20, height = 20)

## Make baseline + follow-up set for Ethnicity analyses
rownames(df)
which(duplicated(colnames(df)))
sums <- colSums(df[1:ncol(df)])
length(sums)
abundant <- names(sums[which(sums > 300)])
prevalent <- names(df[,colSums(df > 0) > 100])
abprev <- abundant[abundant %in% prevalent]
length(abprev)
dftotal <- df %>% dplyr::select(all_of(abprev)) %>% 
    mutate_all(as.integer)
dim(dftotal) #  genes left
mattot <- t(as.matrix(dftotal))
nocounts <- colnames(mattot)[colSums(mattot) == 0]
removeid <- c("HELIBA_103370", "HELIFU_103370")
mattot <- mattot[,!colnames(mattot) %in% removeid]
clindf <- readRDS("data/clinicaldata_long.RDS")
clindf <- clindf %>% filter(sampleID %in% colnames(mattot)) %>% droplevels(.) %>% 
    mutate(timepoint = fct_recode(timepoint, "followup" = "follow-up"))
clindf <- clindf[order(clindf$sampleID, colnames(mattot)),]
all(clindf$sampleID == colnames(mattot)) # TRUE

# DESeq2
library(DESeq2)
# first + 1 for log transformation
mattot <- mattot + 1
deseq2Data <- DESeqDataSetFromMatrix(countData=mattot, colData=clindf, 
                                     design= ~ timepoint + timepoint:Ethnicity)
deseq2Data <- DESeq(deseq2Data, fitType = "mean")
deseq2Results <- results(deseq2Data, 
                         name = "timepointfollowup.EthnicitySurinamese",
                         alpha = 0.05)
deseq2Results$ENSEMBL <- rownames(deseq2Results)
summary(deseq2Results)
deseq2Results[which(deseq2Results$padj < 0.05),]
res <- as.data.frame(deseq2Results)
write.csv2(res, "results/ARG_deseq_baseline_ethnicity.csv")
head(res)
res %>% filter(padj < 0.05) %>% nrow(.)


# counts <- as.data.frame(counts(estimateSizeFactors(deseq2Data), normalized=TRUE))
# counts$ENSEMBL <- rownames(counts)