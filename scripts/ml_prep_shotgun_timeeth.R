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
sg <- read_csv2("data/shotgun/shotgun_abundance.csv")
sg$ID <- sg$...1 
sg$...1 <- NULL
df <- readRDS('data/clinicaldata_long.RDS') %>% filter(sampleID %in% sg$ID)
df <- df %>% mutate(timepoint = case_when(
    timepoint == "baseline" ~ 0,
    timepoint == "follow-up" ~1
))
dutch <- df %>% filter(Ethnicity == "Dutch") # 238 subjects
sas <- df %>% filter(Ethnicity == "Surinamese") # 238 subjects

## All
otu <- sg[which(sg$ID %in% df$sampleID),]
mb1 <- otu[otu$ID %in% df$sampleID,]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
mb2 <- sg[sg$ID %in% df$`sampleID_follow-up`,]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- df %>% filter(sampleID_baseline %in% mbdf$ID)
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

## Dutch
otu <- sg[which(sg$ID %in% dutch$sampleID),]
mb1 <- otu[otu$ID %in% dutch$sampleID,]
tk1 <- apply(mb1[,2:ncol(mb1)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
mb2 <- sg[sg$ID %in% dutch$sampleID,]
tk2 <- apply(mb2[,2:ncol(mb2)], 2, function(x) sum(x > 0.1) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
mbdf <- otu %>% dplyr::select(ID, all_of(names(tk[which(tk == TRUE)])))
clindf <- dutch %>% filter(sampleID %in% mbdf$ID)
clindf <- clindf[match(clindf$sampleID, mbdf$ID),]
all(clindf$sampleID == mbdf$ID) # TRUE
clindf$sampleID; mbdf$ID

tabdm <- screen_wilcox(clindf, "DM_new", mbdf)
head(tabdm)
write.csv2(tabdm, file = "results/wilcoxon_dutch_shotgun.csv")

mbdf <- mbdf[,2:ncol(mbdf)]
path <- 'diabetes_shotgun'
dir.create(path)
dir.create("diabetes_shotgun/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$DM_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))