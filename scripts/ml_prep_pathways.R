# Prediction models timepoints / ethnicity: XGBoost input

library(dplyr)
library(phyloseq)
library(readr)
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
    tot <- left_join(clindf, mbdf, by = "ID")
    res <- c()
    for(a in colnames(mbdf)[2:ncol(mbdf)]){
        tot$dep <- tot[[a]]
        test <- wilcox.test(dep ~ var, data = tot)
        rowtest <- cbind(a, test$p.value)
        colnames(rowtest) <- c("species", "pvalue")
        res <- rbind(res, rowtest)
    }
    res <- as.data.frame(res)
    res <- res %>% mutate(group = var,
                          pvalue = as.numeric(pvalue),
                          padj = p.adjust(pvalue, method = "BH")) %>% 
        arrange(pvalue)
    return(res)
}


## Open dataframe
sg <- read_csv2("data/shotgun/humann/")
sg$ID <- sg$...1 
sg$...1 <- NULL
dim(sg)
sg <- sg %>% filter(!ID %in% "HELIFU_103370") # this one has no baseline sample in set
dim(sg)
df <- readRDS('data/clinicaldata_long.RDS') %>% filter(sampleID %in% sg$ID)
df <- df %>% mutate(timepoint = case_when(
    timepoint == "baseline" ~ 0,
    timepoint == "follow-up" ~1),
    EthnicityTot = case_when(
        EthnicityTot == "Dutch" ~ 0,
        EthnicityTot == "South-Asian Surinamese" ~ 1
    )
)
dutch <- df %>% filter(Ethnicity == "Dutch") # 237 subjects
dutchids <- dutch$sampleID
sas <- df %>% filter(Ethnicity == "Surinamese") # 238 subjects
sasids <- sas$sampleID
ba <- df %>% filter(str_detect(sampleID, "HELIBA"))
fu <- df %>% filter(str_detect(sampleID, "HELIFU"))

## All
otu <- sg[which(sg$ID %in% df$sampleID),]
mb1 <- otu[str_detect(otu$ID, "HELIBA"),]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
mb2 <- otu[str_detect(otu$ID, "HELIFU"),]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
summary(tk)
mbdf <- otu %>% select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- df %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(mbdf$ID, clindf$sampleID),]
all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID_baseline; mbdf$ID

tab <- screen_wilcox(clindf, "timepoint", mbdf)
head(tab)
tab %>% filter(padj < 0.05)
write.csv2(tab, file = "results/wilcoxon_timepoint_shotgun.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'timepoint_shotgun'
dir.create(path)
dir.create("timepoint_shotgun/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$timepoint)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Dutch
otu <- sg[which(sg$ID %in% dutch$sampleID),]
mb1 <- otu[str_detect(otu$ID, "HELIBA"),]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
mb2 <- otu[str_detect(otu$ID, "HELIFU"),]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% dplyr::select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- dutch %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(mbdf$ID, clindf$sampleID),]

all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID; mbdf$ID

tabdutch <- screen_wilcox(clindf, "timepoint", mbdf)
head(tabdutch)
write.csv2(tabdutch, file = "results/wilcoxon_dutchtime_shotgun.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'timepoint_dutch'
dir.create(path, showWarnings = FALSE)
dir.create("timepoint_dutch/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$timepoint)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))


## SAS
otu <- sg[which(sg$ID %in% sas$sampleID),]
mb1 <- otu[str_detect(otu$ID, "HELIBA"),]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
mb2 <- otu[str_detect(otu$ID, "HELIFU"),]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% dplyr::select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- sas %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(mbdf$ID, clindf$sampleID),]

all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID; mbdf$ID

tabsas <- screen_wilcox(clindf, "timepoint", mbdf)
head(tabsas)
write.csv2(tabsas, file = "results/wilcoxon_sastime_shotgun.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'timepoint_sas'
dir.create(path, showWarnings = FALSE)
dir.create("timepoint_sas/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$timepoint)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Baseline - ethnic differences
otu <- sg[which(sg$ID %in% ba$sampleID),]
mb1 <- otu[otu$ID %in% dutchids,]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
mb2 <- otu[otu$ID %in% sasids,]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% dplyr::select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- ba %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(mbdf$ID, clindf$sampleID),]

all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID; mbdf$ID

tabba <- screen_wilcox(clindf, "EthnicityTot", mbdf)
head(tabba)
write.csv2(tabba, file = "results/wilcoxon_ethbase_shotgun.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'eth_base'
dir.create(path, showWarnings = FALSE)
dir.create("eth_base/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$EthnicityTot)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Follow-up ethnic differences
otu <- sg[which(sg$ID %in% fu$sampleID),]
mb1 <- otu[otu$ID %in% dutchids,]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
mb2 <- otu[otu$ID %in% sasids,]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% dplyr::select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- fu %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(mbdf$ID, clindf$sampleID),]

all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID; mbdf$ID

tabfu <- screen_wilcox(clindf, "EthnicityTot", mbdf)
head(tabfu)
write.csv2(tabsas, file = "results/wilcoxon_ethfu_shotgun.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'eth_fu'
dir.create(path, showWarnings = FALSE)
dir.create("eth_fu/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$EthnicityTot)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))
