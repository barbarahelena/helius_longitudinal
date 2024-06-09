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
sbp <- df %>% filter(!is.na(SBP_delta) & HT_BPMed_baseline == "No" & `AntiHT_follow-up` == "No") 
dbp <- df %>% filter(!is.na(DBP_delta) & HT_BPMed_baseline == "No" & `AntiHT_follow-up` == "No")
hba1c <- df %>% filter(!is.na(HbA1c_delta) & DM_baseline == "No" & `GlucLowDrugs_follow-up` == "Nee")
ldl <- df %>% filter(!is.na(LDL_delta) & LLD_baseline == "No" & `LLD_follow-up` == "No")
mb <- readRDS('data/phyloseq_baseline.RDS')
samdata <- inner_join(as(mb@sam_data, 'data.frame'), df, by = "ID")
rownames(samdata) <- samdata$sampleID_baseline
mb <- phyloseq(mb@otu_table, mb@tax_table, mb@refseq, mb@phy_tree, sample_data(samdata))

## Diabetes
mbdm <- prune_samples(!is.na(mb@sam_data$DM_new), mb)
otu <- t(as(mbdm@otu_table, "matrix"))
mb1 <- prune_samples(mbdm@sam_data$DM_new == "No", mbdm)
tk1 <- apply(t(as(mb1@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
mb2 <- prune_samples(mbdm@sam_data$DM_new == "Yes", mbdm)
tk2 <- apply(t(as(mb2@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
tk1[which(tk1 != tk2)]
tk2[which(tk1 != tk2)]
all(tk[which(tk1 != tk2)] == TRUE) # should be TRUE, since either tk1 or tk2 was TRUE
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
mbht <- prune_samples(!is.na(mb@sam_data$HT_new), mb)
otu <- t(as(mbht@otu_table, "matrix"))
mb1 <- prune_samples(mbht@sam_data$HT_new == "No", mbht)
tk1 <- apply(t(as(mb1@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
mb2 <- prune_samples(mbht@sam_data$HT_new == "Yes", mbht)
tk2 <- apply(t(as(mb2@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
tk1[which(tk1 != tk2)]
tk2[which(tk1 != tk2)]
tk[which(tk1 != tk2)] # should all be true (cause one of tk1 or tk2 was)
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
mbms <- prune_samples(!is.na(mb@sam_data$MetSyn_new), mb)
otu <- t(as(mbms@otu_table, "matrix"))
mb1 <- prune_samples(mbms@sam_data$MetSyn_new == "No", mbms)
tk1 <- apply(t(as(mb1@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
mb2 <- prune_samples(mbms@sam_data$MetSyn_new == "Yes", mbms)
tk2 <- apply(t(as(mb2@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
tk1[which(tk1 != tk2)]
tk2[which(tk1 != tk2)]
tk[which(tk1 != tk2)] # should all be true (cause one of tk1 or tk2 was)
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
mbll <- prune_samples(!is.na(mb@sam_data$LLD_new), mb)
otu <- t(as(mbll@otu_table, "matrix"))
mb1 <- prune_samples(mbll@sam_data$LLD_new == "No", mbll)
tk1 <- apply(t(as(mb1@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
mb2 <- prune_samples(mbll@sam_data$LLD_new == "Yes", mbll)
tk2 <- apply(t(as(mb2@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
tk1[which(tk1 != tk2)]
tk2[which(tk1 != tk2)]
tk[which(tk1 != tk2)] # should all be true (cause one of tk1 or tk2 was)
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

## SBP
mb2 <- prune_samples(mb@sam_data$ID %in% sbp$ID, mb)
otu <- t(as(mb2@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.25*length(x)))
mbdf <- otu[,tk]
mbdf <- as.data.frame(mbdf)
dim(mbdf)
clindf <- sbp[match(rownames(mbdf), sbp$sampleID_baseline), ]
all(clindf$sampleID_baseline == rownames(mbdf)) # TRUE
clindf$sampleID_baseline; rownames(mbdf)
path <- 'sbp'
dir.create(path)
dir.create("sbp/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$SBP_delta)
y
write_y(y, name_y = 'y_reg.txt', file.path(path, 'input_data'))

## DBP
mb2 <- prune_samples(mb@sam_data$ID %in% dbp$ID, mb)
otu <- t(as(mb2@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.25*length(x)))
mbdf <- otu[,tk]
mbdf <- as.data.frame(mbdf)
dim(mbdf)
clindf <- dbp[match(rownames(mbdf), dbp$sampleID_baseline), ]
all(clindf$sampleID_baseline == rownames(mbdf)) # TRUE
clindf$sampleID_baseline; rownames(mbdf)
path <- 'dbp'
dir.create(path)
dir.create("dbp/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$DBP_delta)
y
write_y(y, name_y = 'y_reg.txt', file.path(path, 'input_data'))

## HbA1c
mb2 <- prune_samples(mb@sam_data$ID %in% hba1c$ID, mb)
otu <- t(as(mb2@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.25*length(x)))
mbdf <- otu[,tk]
mbdf <- as.data.frame(mbdf)
dim(mbdf)
clindf <- hba1c[match(rownames(mbdf), hba1c$sampleID_baseline), ]
all(clindf$sampleID_baseline == rownames(mbdf)) # TRUE
clindf$sampleID_baseline; rownames(mbdf)
path <- 'hba1c'
dir.create(path)
dir.create("hba1c/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$HbA1c_delta)
y
write_y(y, name_y = 'y_reg.txt', file.path(path, 'input_data'))

## LDL
mb2 <- prune_samples(mb@sam_data$ID %in% ldl$ID, mb)
otu <- t(as(mb2@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.25*length(x)))
mbdf <- otu[,tk]
mbdf <- as.data.frame(mbdf)
dim(mbdf)
clindf <- ldl[match(rownames(mbdf), ldl$sampleID_baseline), ]
all(clindf$sampleID_baseline == rownames(mbdf)) # TRUE
clindf$sampleID_baseline; rownames(mbdf)
path <- 'ldl'
dir.create(path)
dir.create("ldl/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$LDL_delta)
y
write_y(y, name_y = 'y_reg.txt', file.path(path, 'input_data'))
