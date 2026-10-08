# Prediction models timepoints / ethnicity: XGBoost input

library(dplyr)
library(phyloseq)
library(readr)
rm(list=ls())
dir.create("results/3_species_change/5_mlmodels", showWarnings = FALSE, recursive = TRUE)

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
pw <- read_delim("data/shotgun/humann/merged_tables_renorm_unstratified.tsv", delim = "\t") %>% 
    as.data.frame(.)
rownames(pw) <- pw$`# Pathway`
pw$`# Pathway` <- NULL
colnames(pw) <- str_remove(colnames(pw), "_Abundance-CPM")
pw <- as.data.frame(t(as.matrix(pw)))
pw$ID <- rownames(pw)
pw <- pw %>% filter(!ID %in% c("HELIFU_103370", "HELIBA_103370"))
dim(pw) # 950 samples
df <- readRDS('data/clinicaldata_long.RDS') %>% filter(sampleID %in% pw$ID)
df <- df %>% mutate(timepoint = case_when(
    timepoint == "baseline" ~ 0,
    timepoint == "follow-up" ~ 1),
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

sums <- as.data.frame(colSums(pw[,1:582]))
colnames(sums) <- "prev"
sums <- sums %>% mutate(prev = prev / 950)
hist(sums$prev)
median(sums$prev)
min(sums$prev)
tk <- apply(pw, 2, function(x) sum(x > 15) > (0.25*length(x)))
summary(tk)
# sumsprev <- sums %>% filter(prev > 50) %>% print()

## All
pw2 <- pw[which(pw$ID %in% df$sampleID),]
mb1 <- pw2[str_detect(pw2$ID, "HELIBA"),1:ncol(pw2)-1]
tk1 <- apply(mb1, 2, function(x) sum(x > 200) > (0.3*length(x)))
mb2 <- pw2[str_detect(pw2$ID, "HELIFU"),1:ncol(pw2)-1]
tk2 <- apply(mb2, 2, function(x) sum(x > 200) > (0.3*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
summary(tk)
mbdf <- pw2 %>% select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- df %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(mbdf$ID, clindf$sampleID),]
all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID_baseline; mbdf$ID

tab <- screen_wilcox(clindf, "timepoint", mbdf)
head(tab)
tab %>% filter(padj < 0.05)
write.csv2(tab, file = "results/3_species_change/5_mlmodels/wilcoxon_timepoint_pathways.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'timepoint_pathways'
dir.create(path, showWarnings = FALSE)
dir.create("timepoint_pathways/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$timepoint)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Dutch
pw2 <- pw[which(pw$ID %in% dutch$sampleID),]
mb1 <- pw2[str_detect(pw2$ID, "HELIBA"),1:ncol(pw2)-1]
tk1 <- apply(mb1, 2, function(x) sum(x > 200) > (0.3*length(x)))
mb2 <- pw2[str_detect(pw2$ID, "HELIFU"),1:ncol(pw2)-1]
tk2 <- apply(mb2, 2, function(x) sum(x > 200) > (0.3*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
summary(tk)
mbdf <- pw2 %>% dplyr::select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- dutch %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(mbdf$ID, clindf$sampleID),]

all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID; mbdf$ID

tabdutch <- screen_wilcox(clindf, "timepoint", mbdf)
head(tabdutch)
write.csv2(tabdutch, file = "results/3_species_change/5_mlmodels/wilcoxon_dutchtime_pathways.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'timepoint_dutch_pathways'
dir.create(path, showWarnings = FALSE)
dir.create("timepoint_dutch_pathways/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$timepoint)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))


## SAS
pw2 <- pw[which(pw$ID %in% sas$sampleID),]
mb1 <- pw2[str_detect(pw2$ID, "HELIBA"),1:ncol(pw2)-1]
tk1 <- apply(mb1, 2, function(x) sum(x > 200) > (0.3*length(x)))
mb2 <- pw2[str_detect(pw2$ID, "HELIFU"),1:ncol(pw2)-1]
tk2 <- apply(mb2, 2, function(x) sum(x > 200) > (0.3*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
summary(tk)
mbdf <- pw2 %>% dplyr::select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- sas %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(mbdf$ID, clindf$sampleID),]

all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID; mbdf$ID

tabsas <- screen_wilcox(clindf, "timepoint", mbdf)
head(tabsas)
write.csv2(tabsas, file = "results/3_species_change/5_mlmodels/wilcoxon_sastime_pathways.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'timepoint_sas_pathways'
dir.create(path, showWarnings = FALSE)
dir.create("timepoint_sas_pathways/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$timepoint)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Baseline - ethnic differences
pw2 <- pw[which(pw$ID %in% ba$sampleID),]
mb1 <- pw2[str_detect(pw2$ID, "HELIBA"),1:ncol(pw2)-1]
tk1 <- apply(mb1, 2, function(x) sum(x > 200) > (0.3*length(x)))
mb2 <- pw2[str_detect(pw2$ID, "HELIFU"),1:ncol(pw2)-1]
tk2 <- apply(mb2, 2, function(x) sum(x > 200) > (0.3*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
summary(tk)
mbdf <- pw2 %>% dplyr::select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- ba %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(mbdf$ID, clindf$sampleID),]

all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID; mbdf$ID

tabba <- screen_wilcox(clindf, "EthnicityTot", mbdf)
head(tabba)
write.csv2(tabba, file = "results/3_species_change/5_mlmodels/wilcoxon_ethbase_pathways.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'eth_base_pathways'
dir.create(path, showWarnings = FALSE)
dir.create("eth_base_pathways/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$EthnicityTot)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Follow-up ethnic differences
pw2 <- pw[which(pw$ID %in% fu$sampleID),]
mb1 <- pw2[str_detect(pw2$ID, "HELIBA"),1:ncol(pw2)-1]
tk1 <- apply(mb1, 2, function(x) sum(x > 200) > (0.3*length(x)))
mb2 <- pw2[str_detect(pw2$ID, "HELIFU"),1:ncol(pw2)-1]
tk2 <- apply(mb2, 2, function(x) sum(x > 200) > (0.3*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
summary(tk)
mbdf <- pw2 %>% dplyr::select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- fu %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(mbdf$ID, clindf$sampleID),]

all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID; mbdf$ID

tabfu <- screen_wilcox(clindf, "EthnicityTot", mbdf)
head(tabfu)
write.csv2(tabsas, file = "results/3_species_change/5_mlmodels/wilcoxon_ethfu_pathways.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'eth_fu_pathways'
dir.create(path, showWarnings = FALSE)
dir.create("eth_fu_pathways/input_data", showWarnings = FALSE)
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$EthnicityTot)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))
