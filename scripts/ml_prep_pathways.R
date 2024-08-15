# Prediction models new diagnoses: XGBoost input

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
    clindf$ID <- clindf$sampleID_baseline
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
df <- readRDS('data/clinicaldata_wide.RDS')
dm <- df %>% filter(!is.na(DM_new)) %>% 
    mutate(DM_new = case_when(DM_new == "Yes" ~ 1, DM_new == "No" ~ 0))
ht <- df %>% filter(!is.na(HT_new)) %>% 
    mutate(HT_new = case_when(HT_new == "Yes" ~ 1, HT_new == "No" ~ 0))
metsyn <- df %>% filter(!is.na(MetSyn_new)) %>% 
    mutate(MetSyn_new = case_when(MetSyn_new == "Yes" ~ 1, MetSyn_new == "No" ~ 0))
lld <- df %>% filter(!is.na(LLD_new)) %>% 
    mutate(LLD_new = case_when(LLD_new == "Yes" ~ 1, LLD_new == "No" ~ 0))



## Diabetes
# set.seed(999) # to balance the groups, select a sample of controls (same size as cases)
# dmsub <- c(dm$sampleID_baseline[which(dm$DM_new == 1 & dm$sampleID_baseline %in% sg$ID)],
#            sample(dm$sampleID_baseline[which(dm$DM_new == 0 & dm$sampleID_baseline %in% sg$ID)], size = 33))
# dm <- dm %>% filter(sampleID_baseline %in% dmsub)
otu <- sg[which(sg$ID %in% dm$sampleID_baseline),]
mb1 <- otu[otu$ID %in% dm$sampleID_baseline[which(dm$DM_new == 0)],]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.25) > (0.20*length(x)))
mb2 <- sg[sg$ID %in% dm$sampleID_baseline[which(dm$DM_new == 1)],]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.25) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- dm %>% filter(sampleID_baseline %in% mbdf$ID)
all(clindf$sampleID_baseline == mbdf$ID) # TRUE
clindf$sampleID_baseline; mbdf$ID

tabdm <- screen_wilcox(clindf, "DM_new", mbdf)
head(tabdm)
write.csv2(tabdm, file = "results/wilcoxon_dmnew_shotgun.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'diabetes_shotgun'
dir.create(path)
dir.create("diabetes_shotgun/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$DM_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Hypertension
otu <- sg[which(sg$ID %in% ht$sampleID_baseline),]
mb1 <- otu[otu$ID %in% ht$sampleID_baseline[which(ht$HT_new == 0)],]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.25) > (0.20*length(x)))
mb2 <- sg[sg$ID %in% ht$sampleID_baseline[which(ht$HT_new == 1)],]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.25) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- ht %>% filter(sampleID_baseline %in% mbdf$ID)
all(clindf$sampleID_baseline == mbdf$ID) # TRUE
clindf$sampleID_baseline; mbdf$ID
mbdf <- mbdf[,2:ncol(mbdf)]

tabht <- screen_wilcox(clindf, "HT_new", mbdf)
write.csv2(tabht, file = "results/wilcoxon_htnew_shotgun.csv")

path <- 'hypertension_shotgun'
dir.create(path)
dir.create("hypertension_shotgun/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$HT_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## MetSyn
otu <- sg[which(sg$ID %in% metsyn$sampleID_baseline),]
mb1 <- otu[otu$ID %in% metsyn$sampleID_baseline[which(metsyn$MetSyn_new == 0)],]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.25) > (0.20*length(x)))
mb2 <- sg[sg$ID %in% metsyn$sampleID_baseline[which(metsyn$MetSyn_new == 1)],]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.25) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- metsyn %>% filter(sampleID_baseline %in% mbdf$ID)
all(clindf$sampleID_baseline == mbdf$ID) # TRUE
clindf$sampleID_baseline; mbdf$ID
mbdf <- mbdf[,2:ncol(mbdf)]

tabms <- screen_wilcox(clindf, "MetSyn_new", mbdf)
write.csv2(tabms, file = "results/wilcoxon_msnew_shotgun.csv")

path <- 'metsyn_shotgun'
dir.create(path)
dir.create("metsyn_shotgun/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$MetSyn_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## LLD
otu <- sg[which(sg$ID %in% lld$sampleID_baseline),]
mb1 <- otu[otu$ID %in% lld$sampleID_baseline[which(lld$LLD_new == 0)],]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.25) > (0.20*length(x)))
mb2 <- sg[sg$ID %in% lld$sampleID_baseline[which(lld$LLD_new == 1)],]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.25) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- lld %>% filter(sampleID_baseline %in% mbdf$ID)
all(clindf$sampleID_baseline == mbdf$ID) # TRUE
clindf$sampleID_baseline; mbdf$ID
mbdf <- mbdf[,2:ncol(mbdf)]

tablld <- screen_wilcox(clindf, "LLD_new", mbdf)
write.csv2(tablld, file = "results/wilcoxon_lldnew_shotgun.csv")

path <- 'lld_shotgun'
dir.create(path)
dir.create("lld_shotgun/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$LLD_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

