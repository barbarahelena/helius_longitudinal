# Prediction models new diagnoses: XGBoost input

library(dplyr)
library(phyloseq)
library(ggpubr)
rm(list=ls())

# make data for machine learning XGB classification models

# writes input data files for XGB models as tab-delimited 
# subject ids and feature ids are written as separate tab-delimited files
# write X data / predictors
write_data <- function(x, data_path){
    x <- as.matrix(x)
    if(any(is.na(x))){
        cat('There are missing values in the input data!\n')
    }
    write.table(x, file.path(data_path, 'X_data.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(colnames(x), file.path(data_path,'feat_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(rownames(x), file.path(data_path,'subject_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
}

# write y / predicted outcome
write_y <- function(x, name_y, data_path){
    if(missing(name_y)){
        cat('\n\nYou need to provide a name for the y data file!\n')
    }
    if(!name_y %in% c('y_binary.txt', 'y_reg.txt')){
        cat('\nThe file name is not compatible with XGBeast!\n' )
    }
    if(any(is.na(x))){
        cat('\nThere are missing values in the outcome data!\n')
    }
    write.table(x, file = file.path(data_path, name_y), row.names = F, col.names = F, sep = '\t', quote = F)
}

## Screen with wilcoxon-tests
screen_wilcox <- function(clindf, var, mbdf){
    clindf$var <- as.factor(clindf[[var]])
    clindf$ID <- clindf$sampleID
    mbdf$ID <- rownames(mbdf) 
    tot <- left_join(clindf, mbdf, by = "ID")
    res <- c()
    for(a in colnames(mbdf)[!str_detect(colnames(mbdf), "ID")]){
        tot$dep <- log10(tot[[a]]+1)
        test <- wilcox.test(dep ~ var, data = tot)
        rowtest <- cbind(a, test$p.value)
        colnames(rowtest) <- c("ASV", "pvalue")
        res <- rbind(res, rowtest)
    }
    res <- as.data.frame(res)
    res <- res %>% mutate(group = var,
                          pvalue = as.numeric(pvalue),
                          padj = p.adjust(pvalue, method = "BH")) %>% 
        arrange(pvalue)
    return(res)
}
plot_sig <- function(tab, clindf, mbdf){
    tab <- tab %>% filter(padj < 0.05)
    if(nrow(tab) > 20){ tab <- tab %>% arrange(padj) %>% slice(1:20)}
    clindf$var <- as.factor(clindf[[tab$group[1]]])
    clindf$ID <- clindf$sampleID
    mbdf$ID <- rownames(mbdf) 
    tot <- left_join(clindf, mbdf, by = "ID") %>% dplyr::select(var, all_of(tab$ASV)) %>% 
        filter(!is.na(var))
    plist <- c()
    for(a in colnames(tot)[2:ncol(tot)]){
        tot$dep <- tot[[a]]
        pl <- ggplot(tot, aes(x = var, y = dep+1)) +
                    scale_y_log10() +
                    geom_violin(aes(fill = var)) +
                    geom_boxplot(fill = "white", width = 0.2) +
                    stat_compare_means() +
                    ggsci::scale_fill_simpsons(guide = "none") +
                    labs(x = tab$group[1], y = str_c(a, " + 1"), title = a)
                    theme_light()
        plist[[a]] <- pl
    }
    print(ggarrange(plotlist = plist, labels = LETTERS[1:nrow(tab)], common.legend = TRUE))
}


## Open dataframe
mbc <- readRDS("data/16s/phyloseq_paired16s.RDS")
mbc <- as.data.frame(t(as(mbc@otu_table, "matrix")))
rownames(mbc)
df <- readRDS('data/clinicaldata_long.RDS') %>% filter(sampleID %in% rownames(mbc))
df <- df %>% filter(EthnicityTot %in% c("Dutch", "South-Asian Surinamese")) %>% mutate(timepoint = case_when(
    timepoint == "baseline" ~ 0,
    timepoint == "follow-up" ~1),
    EthnicityTot = case_when(
        EthnicityTot == "Dutch" ~ 0,
        EthnicityTot == "South-Asian Surinamese" ~ 1
    )
)
mbc <- mbc[rownames(mbc) %in% df$sampleID,]
dim(mbc) # 2080 subjects

dutch <- df %>% filter(Ethnicity == "Dutch") # 1284 subjects
dutchids <- dutch$sampleID
sas <- df %>% filter(Ethnicity == "Surinamese") # 796 subjects
sasids <- sas$sampleID
ba <- df %>% filter(str_detect(sampleID, "HELIBA"))
fu <- df %>% filter(str_detect(sampleID, "HELIFU"))

## Baseline - follow-up difference
otu <- mbc[which(rownames(mbc) %in% df$sampleID),]
mb1 <- otu[str_detect(rownames(otu), "HELIBA"),]
tk1 <- apply(mb1, 2, function(x) sum(x > 5) > (0.30*length(x)))
mb2 <- otu[str_detect(rownames(otu), "HELIFU"),]
tk2 <- apply(mb2, 2, function(x) sum(x > 5) > (0.30*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
summary(tk)
mbdf <- otu %>% select(all_of(names(tk[which(tk == TRUE)])))
clindf <- df %>% filter(sampleID %in% rownames(mbdf))
clindf <- clindf[match(rownames(mbc), clindf$sampleID),]
all(clindf$sampleID == rownames(mbc)) # TRUE
clindf$sampleID; rownames(mbc)

tab <- screen_wilcox(clindf, "timepoint", mbdf)
head(tab)
tab %>% filter(padj < 0.05)
write.csv2(tab, file = "results/wilcoxon_timepoint_16s.csv")

path <- 'timepoint_16s'
dir.create(path, showWarnings = FALSE)
dir.create("timepoint_16s/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$timepoint)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Dutch
otu <- mbc[which(rownames(mbc) %in% dutch$sampleID),]
mb1 <- otu[str_detect(rownames(otu), "HELIBA"),]
tk1 <- apply(mb1, 2, function(x) sum(x > 5) > (0.30*length(x)))
mb2 <- otu[str_detect(rownames(otu), "HELIFU"),]
tk2 <- apply(mb2, 2, function(x) sum(x > 5) > (0.30*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% dplyr::select(all_of(names(tk[which(tk == TRUE)])))
clindf <- dutch %>% filter(sampleID %in% rownames(mbdf))
clindf <- clindf[match(rownames(mbdf), clindf$sampleID),]

all(clindf$sampleID == rownames(mbdf)) # TRUE
clindf$sampleID; rownames(mbdf)

tabdutch <- screen_wilcox(clindf, "timepoint", mbdf)
head(tabdutch)
write.csv2(tabdutch, file = "results/wilcoxon_dutchtime_16s.csv")

path <- 'timepoint_dutch_16s'
dir.create(path, showWarnings = FALSE)
dir.create("timepoint_dutch_16s/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$timepoint)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))


## SAS
otu <- mbc[which(rownames(mbc) %in% sas$sampleID),]
mb1 <- otu[str_detect(rownames(otu), "HELIBA"),]
tk1 <- apply(mb1, 2, function(x) sum(x > 5) > (0.30*length(x)))
mb2 <- otu[str_detect(rownames(otu), "HELIFU"),]
tk2 <- apply(mb2, 2, function(x) sum(x > 5) > (0.30*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% dplyr::select(all_of(names(tk[which(tk == TRUE)])))
clindf <- sas %>% filter(sampleID %in% rownames(mbdf))
clindf <- clindf[match(rownames(mbdf), clindf$sampleID),]

all(clindf$sampleID == rownames(mbdf)) # TRUE
clindf$sampleID; rownames(mbdf)

tabsas <- screen_wilcox(clindf, "timepoint", mbdf)
head(tabsas)
write.csv2(tabsas, file = "results/wilcoxon_sastime_16s.csv")

path <- 'timepoint_sas_16s'
dir.create(path, showWarnings = FALSE)
dir.create("timepoint_sas_16s/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$timepoint)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Baseline - ethnic differences
otu <- mbc[which(rownames(mbc) %in% ba$sampleID),]
mb1 <- otu[rownames(otu) %in% dutchids,]
tk1 <- apply(mb1, 2, function(x) sum(x > 5) > (0.30*length(x)))
mb2 <- otu[rownames(otu) %in% sasids,]
tk2 <- apply(mb2, 2, function(x) sum(x > 5) > (0.30*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% dplyr::select(all_of(names(tk[which(tk == TRUE)])))
clindf <- ba %>% filter(sampleID %in% rownames(mbdf))
clindf <- clindf[match(rownames(mbdf), clindf$sampleID),]

all(clindf$sampleID == rownames(mbdf)) # TRUE
clindf$sampleID; rownames(mbdf)

tabba <- screen_wilcox(clindf, "EthnicityTot", mbdf)
head(tabba)
write.csv2(tabba, file = "results/wilcoxon_ethbase_16s.csv")

path <- 'eth_base_16s'
dir.create(path, showWarnings = FALSE)
dir.create("eth_base_16s/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$EthnicityTot)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Follow-up ethnic differences
otu <- mbc[which(rownames(mbc) %in% fu$sampleID),]
mb1 <- otu[rownames(otu) %in% dutchids,]
tk1 <- apply(mb1, 2, function(x) sum(x > 5) > (0.30*length(x)))
mb2 <- otu[rownames(otu) %in% sasids,]
tk2 <- apply(mb2, 2, function(x) sum(x > 5) > (0.30*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% dplyr::select(all_of(names(tk[which(tk == TRUE)])))
clindf <- fu %>% filter(sampleID %in% rownames(mbdf))
clindf <- clindf[match(rownames(mbdf), clindf$sampleID),]

all(clindf$sampleID == rownames(mbdf)) # TRUE
clindf$sampleID; rownames(mbdf)

tabfu <- screen_wilcox(clindf, "EthnicityTot", mbdf)
head(tabfu)
write.csv2(tabfu, file = "results/wilcoxon_ethfu_16s.csv")
plot_sig(tab = tabfu, clindf = clindf, mbdf = mbdf)

path <- 'eth_fu_16s'
dir.create(path, showWarnings = FALSE)
dir.create("eth_fu_16s/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$EthnicityTot)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))
