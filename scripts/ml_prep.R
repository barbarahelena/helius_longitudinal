# Prediction models new diagnoses: XGBoost input

library(dplyr)
library(phyloseq)
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


## Open dataframe
df <- readRDS('data/clinicaldata_delta.RDS')
dm <- df %>% filter(!is.na(DM_new)) %>% 
    mutate(DM_new = case_when(DM_new == "Yes" ~ 1, DM_new == "No" ~ 0))
ht <- df %>% filter(!is.na(HT_new)) %>% 
    mutate(HT_new = case_when(HT_new == "Yes" ~ 1, HT_new == "No" ~ 0))
metsyn <- df %>% filter(!is.na(MetSyn_new)) %>% 
    mutate(MetSyn_new = case_when(MetSyn_new == "Yes" ~ 1, MetSyn_new == "No" ~ 0))
lld <- df %>% filter(!is.na(LLD_new)) %>% 
    mutate(LLD_new = case_when(LLD_new == "Yes" ~ 1, LLD_new == "No" ~ 0))
mb <- readRDS('data/phyloseq_baseline.RDS')

## Diabetes
mb2 <- prune_samples(mb@sam_data$ID %in% dm$ID, mb)
otu <- t(as(mb2@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.25*length(x)))
mbdf <- otu[,tk]
mbdf <- as.data.frame(mbdf)
dim(mbdf)
clindf <- dm[match(rownames(mbdf), dm$sampleID_baseline), ]
all(clindf$sampleID_baseline == rownames(mbdf)) # TRUE
clindf$sampleID_baseline; rownames(mbdf)
path <- 'diabetes'
dir.create(path)
dir.create("diabetes/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$DM_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Hypertension
mb2 <- prune_samples(mb@sam_data$ID %in% ht$ID, mb)
otu <- t(as(mb2@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.25*length(x)))
mbdf <- otu[,tk]
mbdf <- as.data.frame(mbdf)
dim(mbdf)
clindf <- ht[match(rownames(mbdf), ht$sampleID_baseline), ]
all(clindf$sampleID_baseline == rownames(mbdf)) # TRUE
clindf$sampleID_baseline; rownames(mbdf)
path <- 'hypertension'
dir.create(path)
dir.create("hypertension/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$HT_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## MetSyn
mb2 <- prune_samples(mb@sam_data$ID %in% metsyn$ID, mb)
otu <- t(as(mb2@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.25*length(x)))
mbdf <- otu[,tk]
mbdf <- as.data.frame(mbdf)
dim(mbdf)
clindf <- metsyn[match(rownames(mbdf), metsyn$sampleID_baseline), ]
all(clindf$sampleID_baseline == rownames(mbdf)) # TRUE
clindf$sampleID_baseline; rownames(mbdf)
path <- 'metsyn'
dir.create(path)
dir.create("metsyn/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$MetSyn_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## LLD
mb2 <- prune_samples(mb@sam_data$ID %in% lld$ID, mb)
otu <- t(as(mb2@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.25*length(x)))
mbdf <- otu[,tk]
mbdf <- as.data.frame(mbdf)
dim(mbdf)
clindf <- lld[match(rownames(mbdf), lld$sampleID_baseline), ]
all(clindf$sampleID_baseline == rownames(mbdf)) # TRUE
clindf$sampleID_baseline; rownames(mbdf)
path <- 'lld'
dir.create(path)
dir.create("lld/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$LLD_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))


