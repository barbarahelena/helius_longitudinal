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
    clindf$ID <- clindf$sampleID_baseline
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
    clindf$var <- as.factor(clindf[[tab$group[1]]])
    clindf$ID <- clindf$sampleID_baseline
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
df <- readRDS('data/clinicaldata_wide.RDS')
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
mb <- readRDS('data/16s/phyloseq_withclinpaired.RDS')
mb@sam_data <- NULL

## Diabetes
all(str_detect(sample_names(mb), "HELIBA")) # FALSE -> both baseline + follow-up samples
mbdm <- prune_samples(sample_names(mb) %in% dm$sampleID_baseline, mb)
all(str_detect(sample_names(mbdm), "HELIBA")) # only baseline samples left

otu <- t(as(mbdm@otu_table, "matrix")) # get otu table
# select features that make prevalence cut for DM_new == 0
mb1 <- prune_samples(sample_names(mbdm) %in% dm$sampleID_baseline[which(dm$DM_new == 0)], mbdm)
tk1 <- apply(t(as(mb1@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
# select features that make prevalence cut for DM_new == 1
mb2 <- prune_samples(sample_names(mbdm) %in% dm$sampleID_baseline[which(dm$DM_new == 1)], mbdm)
tk2 <- apply(t(as(mb2@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0 # TRUE if either tk1 or tk2 is TRUE
tk1[which(tk1 != tk2)] # which features have different prev between groups
tk2[which(tk1 != tk2)]
all(tk[which(tk1 != tk2)] == TRUE) # should be TRUE, since either tk1 or tk2 was TRUE

mbdf <- otu[,tk] # select all features that made the cut
mbdf <- as.data.frame(mbdf)
dim(mbdf)
clindf <- dm[match(rownames(mbdf), dm$sampleID_baseline), ] # put IDs in order of phyloseq
all(clindf$sampleID_baseline == rownames(mbdf)) # TRUE
clindf$sampleID_baseline; rownames(mbdf) # check if these are the same

tabledm <- screen_wilcox(clindf, "DM_new", mbdf) # wilcoxons (output table)
plot_sig(tabledm, clindf, mbdf) # plot sig features

path <- 'diabetes'
dir.create(path, showWarnings = FALSE); dir.create("diabetes/input_data", showWarnings = FALSE)
head(mbdf) # check if IDs are rownames, data present and numerics, colnames are ASVnrs
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$DM_new)
y # classification model - check if these are 0 and 1
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## Hypertension
mbht <- prune_samples(sample_names(mb) %in% ht$sampleID_baseline, mb)
otu <- t(as(mbht@otu_table, "matrix"))
mb1 <- prune_samples(sample_names(mbht) %in% dm$sampleID_baseline[which(ht$HT_new == 0)], mbht)
tk1 <- apply(t(as(mb1@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
mb2 <- prune_samples(sample_names(mbht) %in% dm$sampleID_baseline[which(ht$HT_new == 1)], mbht)
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

tableht <- screen_wilcox(clindf, "HT_new", mbdf)
plot_sig(tableht, clindf, mbdf) # plot sig features

path <- 'hypertension'
dir.create(path)
dir.create("hypertension/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$HT_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## MetSyn
mbms <- prune_samples(sample_names(mb) %in% metsyn$sampleID_baseline, mb)
otu <- t(as(mbms@otu_table, "matrix"))
mb1 <- prune_samples(sample_names(mbms) %in% dm$sampleID_baseline[which(metsyn$MetSyn_new == 0)], mbms)
tk1 <- apply(t(as(mb1@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
mb2 <- prune_samples(sample_names(mbms) %in% dm$sampleID_baseline[which(metsyn$MetSyn_new == 1)], mbms)
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

tablems <- screen_wilcox(clindf, "MetSyn_new", mbdf)

path <- 'metsyn'
dir.create(path)
dir.create("metsyn/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$MetSyn_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## LLD
mbll <- prune_samples(sample_names(mb) %in% lld$sampleID_baseline, mb)
otu <- t(as(mbll@otu_table, "matrix"))
mb1 <- prune_samples(sample_names(mbll) %in% dm$sampleID_baseline[which(lld$LLD_new == 0)], mbll)
tk1 <- apply(t(as(mb1@otu_table, "matrix")), 2, function(x) sum(x > 5) > (0.20*length(x)))
mb2 <- prune_samples(sample_names(mbll) %in% dm$sampleID_baseline[which(lld$LLD_new == 1)], mbll)
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

tablelld <- screen_wilcox(clindf, "LLD_new", mbdf)

path <- 'lld'
dir.create(path)
dir.create("lld/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(clindf$LLD_new)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))

## SBP
mb2 <- prune_samples(sample_names(mb) %in% sbp$sampleID_baseline, mb)
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
mb2 <- prune_samples(sample_names(mb) %in% dbp$sampleID_baseline, mb)
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
mb2 <- prune_samples(sample_names(mb) %in% hba1c$sampleID_baseline, mb)
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
mb2 <- prune_samples(sample_names(mb) %in% ldl$sampleID_baseline, mb)
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
