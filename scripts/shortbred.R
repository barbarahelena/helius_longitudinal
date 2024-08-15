## ShortBRED analyses

## Libraries
library(tidyverse)
library(ggsci)
library(mixOmics)
library(DESeq2)
library(ggrepel)
library(lme4)
library(afex)
library(ggpubr)

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

boxplots <- function(data, var){
    data$dep <- data[[var]]
    pl <- ggplot(data, aes(x = timepoint, y = dep+1)) +
        scale_y_log10() +
        geom_violin(aes(fill = timepoint)) +
        geom_boxplot(fill = "white", width = 0.2) +
        ggpubr::stat_compare_means() +
        ggsci::scale_fill_simpsons(guide = "none") +
        labs(x = "", y = str_c("log10(counts + 1)"), title = var) +
        theme_Publication()
}

## Data
df <- read.csv2("data/shotgun/merged_shortbred.csv", row.names = "X")
dim(df)
df <- df %>% mutate_all(as.numeric)

## Filter out too low prevalent genes
sums <- colSums(df[2:ncol(df)])
sumsdf <- as.data.frame(sums) %>% arrange(-sums) %>% filter(sums != 0)
ggplot(sumsdf, aes(x = 1:nrow(sumsdf), y = sums)) +
    geom_line(color = "firebrick", size = 1.0) +
    theme_Publication() +
    scale_y_log10() +
    labs(x = "Antibiotic resistance genes", y = "total sum counts", title = "Prevalence of genes")
ggsave("results/shortbred/ARG_prevalence.pdf", width = 5, height = 5)
prevalent <- names(sums[which(sums != 0)])
df <- df %>% dplyr::select(all_of(prevalent))
dim(df) # 350 genes left

sums <- rowSums(df)
sumsdf <- as.data.frame(sums) %>% arrange(-sums) %>% filter(sums != 0)
ggplot(sumsdf, aes(x = 1:nrow(sumsdf), y = sums)) +
    geom_line(color = "royalblue", size = 1.0) +
    theme_Publication() +
    scale_y_log10() +
    labs(x = "Subjects", y = "total sum counts", title = "Subjects")
ggsave("results/shortbred/subjects_prevalence.pdf", width = 5, height = 5)
prevalent <- names(sums[which(sums != 0)])
rownames(df)[which(!rownames(df) %in% prevalent)]
df <- df %>% dplyr::filter(!rownames(.) %in% c("HELIBA_103370", "HELIFU_103370"))
dim(df) # 950 subjects left

presentabsent <- as.matrix(df)
presentabsent <- ifelse(presentabsent != 0, 1, 0)
sums <- rowSums(presentabsent)
median(sums)
mean(sums)
sumsdf <- as.data.frame(sums) 
ggplot(sumsdf, aes(x = sums)) + 
    geom_histogram(fill = "royalblue3", alpha = 0.7) +
    labs(x = "total number of ARG", y = "number of subjects", 
         title = "Number of genes per subject") +
    theme_Publication()
ggsave("results/shortbred/ARG_subjects.pdf", width = 5, height = 5)

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
dfmat <- as.matrix(dfbaseline)
th <- apply(dfmat, 2, function(x) sum(x > 1) > (0.05*length(x)))
matbase <- t(dfmat[,th]) # 107 genes left
mode(matbase) <- "integer"
clindf <- readRDS("data/clinicaldata_wide.RDS")
clindf <- clindf %>% filter(ID %in% rownames(dfmat)) %>% droplevels(.)
all(clindf$ID == colnames(matbase)) # TRUE

# Which genes are most prevalent
# sumsdf <- as.data.frame(sums) %>% arrange(-sums)
# head(sumsdf)
# mostprev <- sumsdf %>% filter(sums > 1000)
# baseprev <- dfbaseline %>% dplyr::select(all_of(rownames(mostprev))) %>% 
#     mutate(ID = rownames(.)) %>% 
#     right_join(., clindf, by = "ID")
# baseprevlong <- pivot_longer(baseprev, cols = 1:46, names_to = "ARG", 
#                              values_to = "counts")
# ggplot(baseprevlong, aes(x = Ethnicity, y = log10(counts +1))) +
#     geom_violin(aes(fill = Ethnicity)) +
#     scale_fill_jco()+
#     geom_boxplot(fill = "white", width = 0.25) +
#     ggpubr::stat_compare_means() +
#     facet_wrap(~ARG) +
#     theme_Publication() +
#     labs(x = "", y = "counts (log10 scale)")
# ggsave("results/highestARG_baseline.pdf", width = 20, height = 20)

#### DESeq2 baseline ####
# first + 1 for log transformation
matbase <- matbase + 1
levels(clindf$Ethnicity) # Dutch is ref (alphabetically)
deseq2Data <- DESeqDataSetFromMatrix(countData=matbase, colData=clindf, 
                                     design= ~ Ethnicity)
deseq2Data <- DESeq(deseq2Data, fitType = "mean")
deseq2Results <- results(deseq2Data, alpha = 0.05)
deseq2Results$Gene.name <- rownames(deseq2Results)
summary(deseq2Results)
deseq2Results[which(deseq2Results$padj < 0.05),]
res <- as.data.frame(deseq2Results)
res <- res %>% mutate(
    group = case_when(
        log2FoldChange < 0 ~ paste0("lower in Surinamese"),
        log2FoldChange >= 0 ~ paste0("higher in Surinamese")
    ),
    sig = case_when(padj < 0.05 ~ paste0("sig"), 
                    padj > 0.05 ~ paste0("insig")),
    sigdir = case_when(
        sig == "sig" & group == "higher in Surinamese" ~ paste0("up"),
        sig == "sig" & group == "lower in Surinamese" ~ paste0("down"),
        sig == "insig" ~ paste0("no")
    ),
    sigdir = as.factor(sigdir),
    log2eq = ifelse(log2FoldChange < 0, log2FoldChange * -1, log2FoldChange)
) %>% 
    # filter(padj < 0.05) %>% 
    arrange(log2FoldChange) %>% 
    mutate(Gene.name = fct_inorder(Gene.name)) %>% droplevels(.)
highestfold <- res %>% arrange(-log2FoldChange) %>% dplyr::slice(1:5)
lowestfold <- res %>% arrange(log2FoldChange) %>% dplyr::slice(1:5)
highestsig <- res %>% arrange(-padj) %>% dplyr::slice(1:10)
listint <- unique(c(highestfold$Gene.name, lowestfold$Gene.name, highestsig$Gene.name)) %>% 
    droplevels(.)
res <- res %>% mutate(labels = case_when(
    Gene.name %in% listint ~ Gene.name,
    .default = ""
))
write.csv2(res, "results/shortbred/ARG_deseq_baseline_ethnicity.csv")
head(res)
tail(res)
res %>% filter(padj < 0.05) %>% nrow(.)

ggplot(res %>% filter(padj < 0.05), aes(x = Gene.name, y = log2FoldChange, fill = group)) +
    theme_Publication() +
    geom_bar(stat = "identity", alpha = 0.8) +
    ggsci::scale_fill_simpsons() +
    coord_flip() +
    labs(y = 'Log 2 fold difference in Surinamese compared to Dutch', x='', 
         title = 'Differential gene expression baseline',
         fill = '') +
    theme(axis.text.x = element_text(size=9)) + 
    theme(axis.text.y = element_text(size=8)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')
ggsave("results/shortbred/DGE_baseline_ethnicity.pdf", height = 8, width = 10, device = "pdf")

res2 <- res %>% filter(padj < 0.05) %>% arrange(-padj)
res2 <- res2[1:20,]

ggplot(res2, aes(x = Gene.name, y = log2FoldChange, fill = group)) +
    theme_Publication() +
    geom_bar(stat = "identity", alpha = 0.8) +
    ggsci::scale_fill_simpsons() +
    coord_flip() +
    labs(y = 'Log 2 fold difference in Surinamese compared to Dutch', x='', 
         title = 'Differential gene expression baseline',
         fill = '') +
    theme(axis.text.x = element_text(size=9)) + 
    theme(axis.text.y = element_text(size=10)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')
ggsave("results/shortbred/DGE_baseline_ethnicity_top20sig.pdf", height = 8, width = 10, device = "pdf")

set.seed(1234)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = group, label = labels)) +
    theme_Publication() +
    geom_hline(aes(yintercept = -log10(0.1)), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 1.1, min.segment.length = 0,
                    point.padding = 0, color = "black",
                    fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey50", 
                    force = 1, max.overlaps = 10) +
    scale_color_manual(values = c(ggsci::pal_simpsons()(2)), guide = "none") +
    labs(x = "Log2 fold difference in Surinamese compared to Dutch",
         y = "-log10(p-value)") 
ggsave("results/shortbred/volcanoplot_baseline.pdf", width = 7, height = 7, device = "pdf")    

## Make baseline + follow-up set for Ethnicity analyses
rownames(df)
which(duplicated(colnames(df)))
dfmat <- as.matrix(df)
th <- apply(dfmat, 2, function(x) sum(x > 1) > (0.05*length(x)))
mattot <- t(dfmat[,th]) # 109 genes left
mode(mattot) <- "integer"
clindf <- readRDS("data/clinicaldata_long.RDS")
clindf <- clindf %>% filter(sampleID %in% colnames(mattot)) %>% droplevels(.) %>% 
    mutate(timepoint = fct_recode(timepoint, "followup" = "follow-up"))
clindf <- clindf[order(clindf$sampleID, colnames(mattot)),]
all(clindf$sampleID == colnames(mattot)) # TRUE

# DESeq2
# first + 1 for log transformation
mattot <- mattot + 1
deseq2Data <- DESeqDataSetFromMatrix(countData=mattot, colData=clindf, 
                                     design= ~ timepoint)
deseq2Data <- DESeq(deseq2Data, fitType = "mean")
deseq2Results <- results(deseq2Data, alpha = 0.05)
deseq2Results$Gene.name <- rownames(deseq2Results)
summary(deseq2Results)
deseq2Results[which(deseq2Results$padj < 0.05),]
res <- as.data.frame(deseq2Results)
res <- res %>% mutate(
    group = case_when(
        log2FoldChange < 0 ~ paste0("lower at follow-up"),
        log2FoldChange >= 0 ~ paste0("higher at follow-up")
    ),
    sig = case_when(padj < 0.05 ~ paste0("sig"), 
                    padj > 0.05 ~ paste0("insig")),
    sigdir = case_when(
        sig == "sig" & group == "higher at follow-up" ~ paste0("up"),
        sig == "sig" & group == "lower at follow-up" ~ paste0("down"),
        sig == "insig" ~ paste0("no")
    ),
    sigdir = as.factor(sigdir)
) %>% 
    # filter(padj < 0.05) %>% 
    arrange(log2FoldChange) %>% 
    mutate(Gene.name = fct_inorder(Gene.name)) %>% droplevels(.)
highestfold <- res %>% arrange(-log2FoldChange) %>% dplyr::slice(1:5)
lowestfold <- res %>% arrange(log2FoldChange) %>% dplyr::slice(1:5)
highestsig <- res %>% arrange(-padj) %>% dplyr::slice(1:10)
listint <- unique(c(highestfold$Gene.name, lowestfold$Gene.name, highestsig$Gene.name)) %>% 
    droplevels(.)
res <- res %>% mutate(labels = case_when(
    Gene.name %in% listint ~ Gene.name,
    .default = ""
))
gene_lmm <- res %>% filter(labels != "") %>% dplyr::select(Gene.name)
write.csv2(res, "results/shortbred/ARG_deseq_FU.csv")
head(res)
tail(res)
res %>% filter(padj < 0.05) %>% nrow(.)

ggplot(res %>% filter(padj < 0.05), aes(x = Gene.name, y = log2FoldChange, fill = group)) +
    theme_Publication() +
    geom_bar(stat = "identity", alpha = 0.8) +
    ggsci::scale_fill_simpsons() +
    coord_flip() +
    labs(y = 'Log 2 fold difference at follow-up', x='', 
         title = 'Differential gene expression',
         fill = '') +
    theme(axis.text.x = element_text(size=9)) + 
    theme(axis.text.y = element_text(size=8)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')
ggsave("results/shortbred/DGE_followup.pdf", height = 14, width = 10, device = "pdf")

res2 <- res %>% filter(padj < 0.05) %>% arrange(-padj)
res2 <- res2[1:20,]

ggplot(res2, aes(x = Gene.name, y = log2FoldChange, fill = group)) +
    theme_Publication() +
    geom_bar(stat = "identity", alpha = 0.8) +
    ggsci::scale_fill_simpsons() +
    coord_flip() +
    labs(y = 'Log 2 fold difference at follow-up', x='', 
         title = 'Differential gene expression',
         fill = '') +
    theme(axis.text.x = element_text(size=9)) + 
    theme(axis.text.y = element_text(size=10)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')
ggsave("results/shortbred/DGE_followup_top20sig.pdf", height = 8, width = 10, device = "pdf")

set.seed(1234)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = group, label = labels)) +
    theme_Publication() +
    geom_hline(aes(yintercept = -log10(0.1)), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 1.1, min.segment.length = 0,
                    point.padding = 0, color = "black",
                    fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey50", 
                    force = 1, max.overlaps = 50) +
    scale_color_manual(values = c(ggsci::pal_simpsons()(2)), guide = "none") +
    labs(x = "Log2 fold change",
         y = "-log10(p-value)") 
ggsave("results/shortbred/volcanoplot_followup.pdf", width = 7, height = 7, device = "pdf")


## Dutch subjects change over time
clindutch <- clindf %>% filter(Ethnicity == "Dutch")
matdutch <- mattot[,colnames(mattot) %in% clindutch$sampleID]
deseq2Data <- DESeqDataSetFromMatrix(countData=matdutch, colData=clindutch, 
                                     design= ~ timepoint)
deseq2Data <- DESeq(deseq2Data, fitType = "mean")
deseq2Results <- results(deseq2Data, alpha = 0.05)
deseq2Results$Gene.name <- rownames(deseq2Results)
summary(deseq2Results)
deseq2Results[which(deseq2Results$padj < 0.05),]
res <- as.data.frame(deseq2Results)
res <- res %>% mutate(
    group = case_when(
        log2FoldChange < 0 ~ paste0("lower at follow-up"),
        log2FoldChange >= 0 ~ paste0("higher at follow-up")
    ),
    sig = case_when(padj < 0.05 ~ paste0("sig"), 
                    padj > 0.05 ~ paste0("insig")),
    sigdir = case_when(
        sig == "sig" & group == "higher at follow-up" ~ paste0("up"),
        sig == "sig" & group == "lower at follow-up" ~ paste0("down"),
        sig == "insig" ~ paste0("no")
    ),
    sigdir = as.factor(sigdir),
    log2eq = ifelse(log2FoldChange < 0, log2FoldChange * -1, log2FoldChange)
) %>% 
    arrange(log2FoldChange) %>% 
    mutate(Gene.name = fct_inorder(Gene.name)) %>% droplevels(.)
highestfold <- res %>% arrange(-log2FoldChange) %>% dplyr::slice(1:5)
lowestfold <- res %>% arrange(log2FoldChange) %>% dplyr::slice(1:5)
highestsig <- res %>% arrange(-padj) %>% dplyr::slice(1:10)
listint <- unique(c(highestfold$Gene.name, lowestfold$Gene.name, highestsig$Gene.name)) %>% 
    droplevels(.)
res <- res %>% mutate(labels = case_when(
    Gene.name %in% listint ~ Gene.name,
    .default = ""
))
gene_lmm_dutch <- res %>% filter(labels != "") %>% dplyr::select(Gene.name)
write.csv2(res, "results/shortbred/ARG_deseq_FU_Dutch.csv")
head(res)
tail(res)
res %>% filter(padj < 0.05) %>% nrow(.)

ggplot(res %>% filter(padj < 0.05), aes(x = Gene.name, y = log2FoldChange, fill = group)) +
    theme_Publication() +
    geom_bar(stat = "identity", alpha = 0.8) +
    ggsci::scale_fill_simpsons() +
    coord_flip() +
    labs(y = 'Log 2 fold difference at follow-up', x='', 
         title = 'Differential gene expression - Dutch',
         fill = '') +
    theme(axis.text.x = element_text(size=9)) + 
    theme(axis.text.y = element_text(size=8)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')
ggsave("results/shortbred/DGE_followup_dutch.pdf", height = 8, width = 10, device = "pdf")

res2 <- res %>% filter(padj < 0.05) %>% arrange(-padj)
res2 <- res2[1:20,]

ggplot(res2, aes(x = Gene.name, y = log2FoldChange, fill = group)) +
    theme_Publication() +
    geom_bar(stat = "identity", alpha = 0.8) +
    ggsci::scale_fill_simpsons() +
    coord_flip() +
    labs(y = 'Log 2 fold difference at follow-up', x='', 
         title = 'Differential gene expression - Dutch',
         fill = '') +
    theme(axis.text.x = element_text(size=9)) + 
    theme(axis.text.y = element_text(size=10)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')
ggsave("results/shortbred/DGE_followupdutch_top20sig.pdf", height = 8, width = 10, device = "pdf")

set.seed(1234)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = group, label = labels)) +
    theme_Publication() +
    geom_hline(aes(yintercept = -log10(0.1)), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 1.1, min.segment.length = 0,
                    point.padding = 0, color = "black",
                    fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey50", 
                    force = 1, max.overlaps = 50) +
    scale_color_manual(values = c(ggsci::pal_simpsons()(2)), guide = "none") +
    labs(x = "Log2 fold change",
         y = "-log10(p-value)", title = "Change over time - Dutch") 
ggsave("results/shortbred/volcanoplot_followupdutch.pdf", width = 7, height = 7, device = "pdf")

dfdutch <- as.data.frame(t(matdutch))
dfdutch$sampleID <- rownames(dfdutch)
geneplot <- res %>% arrange(padj) %>% dplyr::slice(1:20) %>% dplyr::select(Gene.name)
geneplot <- geneplot$Gene.name
dfplots <- left_join(clindutch, dfdutch, by = "sampleID") %>% dplyr::select(all_of(geneplot), timepoint)

plist <- list()
for(var in geneplot){
    pl <- boxplots(dfplots, var)
    plist[[var]] <- pl
}

ggarrange(plotlist = plist, labels = LETTERS[1:length(geneplot)], common.legend = TRUE)
ggsave("results/shortbred/differences_dutch.pdf", height = 14, width = 18)


## SAS subjects change over time
clinsas <- clindf %>% filter(Ethnicity == "Surinamese")
matsas <- mattot[,colnames(mattot) %in% clinsas$sampleID]
deseq2Data <- DESeqDataSetFromMatrix(countData=matsas, colData=clinsas, 
                                     design= ~ timepoint)
deseq2Data <- DESeq(deseq2Data, fitType = "mean")
deseq2Results <- results(deseq2Data, alpha = 0.05)
deseq2Results$Gene.name <- rownames(deseq2Results)
summary(deseq2Results)
deseq2Results[which(deseq2Results$padj < 0.05),]
res <- as.data.frame(deseq2Results)
res <- res %>% mutate(
    group = case_when(
        log2FoldChange < 0 ~ paste0("lower at follow-up"),
        log2FoldChange >= 0 ~ paste0("higher at follow-up")
    ),
    sig = case_when(padj < 0.05 ~ paste0("sig"), 
                    padj > 0.05 ~ paste0("insig")),
    sigdir = case_when(
        sig == "sig" & group == "higher at follow-up" ~ paste0("up"),
        sig == "sig" & group == "lower at follow-up" ~ paste0("down"),
        sig == "insig" ~ paste0("no")
    ),
    sigdir = as.factor(sigdir),
    log2eq = ifelse(log2FoldChange < 0, log2FoldChange * -1, log2FoldChange)
) %>% 
    arrange(log2FoldChange) %>% 
    mutate(Gene.name = fct_inorder(Gene.name)) %>% droplevels(.)
highestfold <- res %>% arrange(-log2FoldChange) %>% dplyr::slice(1:5)
lowestfold <- res %>% arrange(log2FoldChange) %>% dplyr::slice(1:5)
highestsig <- res %>% arrange(-padj) %>% dplyr::slice(1:10)
listint <- unique(c(highestfold$Gene.name, lowestfold$Gene.name, highestsig$Gene.name)) %>% 
    droplevels(.)
res <- res %>% mutate(labels = case_when(
    Gene.name %in% listint ~ Gene.name,
    .default = ""
))
gene_lmm_sas <- res %>% filter(labels != "") %>% dplyr::select(Gene.name)
write.csv2(res, "results/shortbred/ARG_deseq_FU_SAS.csv")
head(res)
tail(res)
res %>% filter(padj < 0.05) %>% nrow(.)

ggplot(res %>% filter(padj < 0.05), aes(x = Gene.name, y = log2FoldChange, fill = group)) +
    theme_Publication() +
    geom_bar(stat = "identity", alpha = 0.8) +
    ggsci::scale_fill_simpsons() +
    coord_flip() +
    labs(y = 'Log 2 fold difference at follow-up', x='', 
         title = 'Differential gene expression - Surinamese',
         fill = '') +
    theme(axis.text.x = element_text(size=9)) + 
    theme(axis.text.y = element_text(size=8)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')
ggsave("results/shortbred/DGE_followup_sas.pdf", height = 8, width = 10, device = "pdf")

res2 <- res %>% filter(padj < 0.05) %>% arrange(-padj)
res2 <- res2[1:20,]

ggplot(res2, aes(x = Gene.name, y = log2FoldChange, fill = group)) +
    theme_Publication() +
    geom_bar(stat = "identity", alpha = 0.8) +
    ggsci::scale_fill_simpsons() +
    coord_flip() +
    labs(y = 'Log 2 fold difference at follow-up', x='', 
         title = 'Differential gene expression - Surinamese',
         fill = '') +
    theme(axis.text.x = element_text(size=9)) + 
    theme(axis.text.y = element_text(size=10)) +
    theme(legend.key.size= unit(0.5, "cm")) +
    theme(legend.position = 'none', legend.justification = 'center')
ggsave("results/shortbred/DGE_followupsas_top20sig.pdf", height = 8, width = 10, device = "pdf")

set.seed(1234)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = group, label = labels)) +
    theme_Publication() +
    geom_hline(aes(yintercept = -log10(0.1)), color = "darkgrey", linetype = "dashed") +
    geom_point(alpha = 0.6) +
    geom_text_repel(size = 3, seed = 24, box.padding = 1.1, min.segment.length = 0,
                    point.padding = 0, color = "black",
                    fontface = "bold", force_pull = 0,
                    nudge_x = 0.05, nudge_y = 0.5, segment.color = "grey50", 
                    force = 1, max.overlaps = 50) +
    scale_color_manual(values = c(ggsci::pal_simpsons()(2)), guide = "none") +
    labs(x = "Log2 fold change",
         y = "-log10(p-value)", title = "Change over time - Surinamese") 
ggsave("results/shortbred/volcanoplot_followupsas.pdf", width = 7, height = 7, device = "pdf")


dfsas <- as.data.frame(t(matsas))
dfsas$sampleID <- rownames(dfsas)
geneplot <- res %>% arrange(padj) %>% dplyr::slice(1:20) %>% dplyr::select(Gene.name)
geneplot <- geneplot$Gene.name
dfplots <- left_join(clinsas, dfsas, by = "sampleID") %>% dplyr::select(all_of(geneplot), timepoint)

plist <- list()
for(var in geneplot){
    pl <- boxplots(dfplots, var)
    plist[[var]] <- pl
}

ggarrange(plotlist = plist, labels = LETTERS[1:length(geneplot)], common.legend = TRUE)
ggsave("results/shortbred/differences_sas.pdf", height = 14, width = 18)


#### LMM ####
linearmixed_arg <- function(data, var){
    data$var <- log10(data[[var]]+1)
    model1_v4 <- lmer(var ~ Ethnicity*timepoint + (1|ID),
                      data = data)
    res_v4 <- summary(model1_v4)
    pval <- format(round(res_v4$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres <- cbind(gene = var, group1 = 0, group2 = 1, pval)
    statres <- tibble::as_tibble(statres)
    statres$p_signif <- case_when(
        statres$pval < 0.05 ~paste0("*"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.001 ~paste0("***"),
        .default = ""
    )
    return(statres)
}

plot_lmm <- function(data, var, lmm){
    data$var <- log10(data[[var]]+1)
    df_means <- data %>% 
        group_by(Ethnicity, timepoint) %>% 
        summarise(mean_var=mean(var, na.rm = TRUE), sd_var = sd(var, na.rm = TRUE),
                  n_var = length(var))
    pl <- ggplot() +
         geom_line(data = df_means, aes(x = timepoint, y = mean_var, 
                color = Ethnicity, group = Ethnicity), alpha = 1, linewidth = 0.75) +
         geom_line(data = data, aes(x = timepoint, y = var, color = Ethnicity, group = ID), alpha = 0.1) +
         geom_errorbar(data = df_means,
                       aes(ymin = mean_var - (sd_var/sqrt(n_var)),
                           ymax = mean_var + (sd_var/sqrt(n_var)),
                           x = timepoint,
                           color = Ethnicity), width=0.05, linewidth = 0.75) +
         ggpubr::stat_pvalue_manual(lmm, y.position = max(data$var), label = "p_signif", 
                                    tip.length = 0, bracket.shorten = 0.1, size = 5) +
         scale_color_manual(values = rev(pal_simpsons()(2))) + 
         # scale_y_continuous(limits = c(0,220), breaks = seq(from = 0, to = 220, by = 20)) +
         theme_Publication() +
         labs(x = "Timepoint", y = "log10(cpm+1)", title = var,
              color = "")
    print(pl)
}

names(df)
df2 <- df %>% dplyr::select(all_of(gene_lmm$Gene.name))
df2$sampleID <- rownames(df2)
dftot <- left_join(df2, clindf, by = "sampleID") %>% 
    mutate(timepoint = case_when(
        timepoint == "baseline" ~ 0,
        timepoint == "followup" ~ 1
    ),
    timepoint = as.numeric(timepoint),
    ID = as.factor(ID)
    )
names(dftot)

names(dftot)
statres <- c()
for(a in colnames(dftot)[1:ncol(df2)-1]){
    print(a)
    statresline <- linearmixed_arg(dftot, a)
    statres <- rbind(statres, statresline)
}

names(df)
df2 <- df %>% dplyr::select(all_of(gene_lmm_dutch$Gene.name))
df2$sampleID <- rownames(df2)
dftot <- left_join(df2, clindf, by = "sampleID") %>% 
    mutate(timepoint = case_when(
        timepoint == "baseline" ~ 0,
        timepoint == "followup" ~ 1
    ),
    timepoint = as.numeric(timepoint),
    ID = as.factor(ID)
    )
names(dftot)

names(dftot)
statresdutch <- c()
for(a in colnames(dftot)[1:ncol(df2)-1]){
    print(a)
    statresline <- linearmixed_arg(dftot, a)
    statresdutch <- rbind(statresdutch, statresline)
}

names(df)
df2 <- df %>% dplyr::select(all_of(gene_lmm_sas$Gene.name))
df2$sampleID <- rownames(df2)
dftot <- left_join(df2, clindf, by = "sampleID") %>% 
    mutate(timepoint = case_when(
        timepoint == "baseline" ~ 0,
        timepoint == "followup" ~ 1
    ),
    timepoint = as.numeric(timepoint),
    ID = as.factor(ID)
    )
names(dftot)

names(dftot)
statressas <- c()
for(a in colnames(dftot)[1:ncol(df2)-1]){
    print(a)
    statresline <- linearmixed_arg(dftot, a)
    statressas <- rbind(statressas, statresline)
}

statsig <- statressas %>% filter(pval < 0.05) %>% 
    mutate(group1 = "baseline",
           group2 = "followup")
statsig$gene

dfsig <- df %>% dplyr::select(all_of(statsig$gene))
dfsig$sampleID <- rownames(dfsig)
dfsig <- left_join(dfsig, clindf, by = "sampleID")
for(a in colnames(dfsig)[1:nrow(statsig)-1]){
    print(a)
    plot_lmm(data = dfsig, var = a, lmm = statsig)
    statressas <- rbind(statressas, statresline)
}

