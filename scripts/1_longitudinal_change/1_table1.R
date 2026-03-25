# Table 1

# Libraries
library(tableone)
library(dplyr)
library(stringr)

## Data
helius <- readRDS('data/clinicaldata_wide.RDS')
mbids <- read.csv2('data/16s/ids_16s_paired.csv')%>% dplyr::select(ID = x)
mbids$microbiome_16s <- TRUE
mbids2 <- readRDS('data/16s/phyloseq_paired16s.RDS')
mbids2 <- str_c("S", str_extract(str_remove(sample_names(mbids2), "_T1"), "[0-9]+"))
shotids <- read.csv('data/shotgun/helius_ids.csv') %>% dplyr::select(ID = x) %>% filter(!str_detect(ID, "HELIBA"))
shotids$shotgun <- TRUE
shotids$ID <- str_c("S", str_remove_all(shotids$ID, "HELIFU_"))

helius <- left_join(helius, metids, by = "ID") %>% 
    left_join(., mbids, by = 'ID') %>% 
    left_join(., shotids, by = 'ID') %>% 
    mutate(across(c('microbiome_16s', 'shotgun'), ~case_when(
        is.na(.x) ~ FALSE,
        .default = TRUE
    )))
names(helius)

##### Table 1 #####
table1 <- helius %>% 
    dplyr::select(Age_baseline, Sex, Ethnicity, MigrGen, FUtime, BMI_baseline, Smoking_baseline, AlcCons_baseline, 
           DM_baseline, SBP_baseline, DBP_baseline, HT_BPMed_baseline, MetSyn_baseline, 
           PPI_baseline, Metformin_baseline, Statins_baseline, TC_baseline, LDL_baseline, Trig_baseline, 
           HbA1c_baseline, SCORECVDmortNL_baseline,
           metabolomics, microbiome_16s, shotgun) %>% 
    CreateTableOne(data=., test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
table1 <- table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

table2 <- helius %>% 
    dplyr::select(Age_baseline, Sex, Ethnicity, MigrGen, FUtime, BMI_baseline, Smoking_baseline, AlcCons_baseline, 
           DM_baseline, SBP_baseline, DBP_baseline, HT_BPMed_baseline, MetSyn_baseline, 
           PPI_baseline, Metformin_baseline, Statins_baseline, TC_baseline, LDL_baseline, 
           Trig_baseline, HbA1c_baseline, SCORECVDmortNL_baseline,
           metabolomics, microbiome_16s, shotgun) %>% 
    CreateTableOne(data=., strata = 'microbiome_16s', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
table2 <- table2 %>% dplyr::select(`16S` = `TRUE`) %>%  
            mutate(across(everything(.), ~trimws(.x, which = "both")))

table3 <- helius %>% 
    dplyr::select(Age_baseline, Sex, Ethnicity, MigrGen, FUtime, BMI_baseline, Smoking_baseline, AlcCons_baseline, 
           DM_baseline, SBP_baseline, DBP_baseline, HT_BPMed_baseline, MetSyn_baseline, 
           PPI_baseline, Metformin_baseline, Statins_baseline, TC_baseline, 
           LDL_baseline, Trig_baseline, HbA1c_baseline, SCORECVDmortNL_baseline,
           metabolomics, microbiome_16s, shotgun) %>% 
    CreateTableOne(data=., strata = 'shotgun', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
table3 <- table3 %>% dplyr::select(`shotgun` = `TRUE`) %>% 
    mutate(across(everything(.), ~trimws(.x, which = "both")))

table <- cbind(table1,table2,table3)
write.csv2(as.data.frame(table), 'results/tables/table1.csv')

table1eth <- helius %>% 
    dplyr::select(Age_baseline, Sex, EthnicityTot, MigrGen, FUtime, BMI_baseline, Smoking_baseline, AlcCons_baseline, 
           DM_baseline, SBP_baseline, DBP_baseline, HT_BPMed_baseline, MetSyn_baseline, 
           PPI_baseline, Metformin_baseline, Statins_baseline, TC_baseline, LDL_baseline, Trig_baseline, 
           HbA1c_baseline, SCORECVDmortNL_baseline,
           metabolomics, microbiome_16s, shotgun) %>% 
    CreateTableOne(data=., test = FALSE, strata = "EthnicityTot") %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
table1eth <- table1eth %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(as.data.frame(table1eth), 'results/tables/suppltable1.csv')

## Table 2
helius <- readRDS('data/clinicaldata_long.RDS')
table1 <- helius %>% 
    dplyr::select(Age, Sex, Ethnicity, MigrGen, FUtime, BMI, Smoking, AlcCons, 
           DM, SBP, DBP, HT_BPMed, MetSyn, 
           PPI, Metformin, Statins, TC, LDL, Trig, HbA1c, SCORECVDmortNL,
           timepoint) %>% 
    CreateTableOne(data=., strata = 'timepoint', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig")) %>% 
    as.data.frame(.)
table1 <- table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(as.data.frame(table1), 'results/tables/followuptable.csv')

## Table 3
helius <- readRDS('data/clinicaldata_long.RDS')
shotids <- read.csv('data/shotgun/helius_ids.csv') %>% dplyr::select(sampleID = x)
shotids$shotgun <- TRUE

helius <- left_join(helius, shotids, by = "sampleID") %>% 
    mutate(across(c('shotgun'), ~case_when(
        is.na(.x) ~ FALSE,
        .default = TRUE
    )))

table1 <- helius %>% 
    filter(shotgun == TRUE) %>% 
    dplyr::select(Age, Sex, Ethnicity, MigrGen, FUtime, BMI, Smoking, AlcCons, 
           DM, SBP, DBP, HT_BPMed, MetSyn, 
           PPI, Metformin, Statins, TC, LDL, Trig, HbA1c, SCORECVDmortNL,
           timepoint) %>% 
    CreateTableOne(data=., strata = 'timepoint', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig")) %>% 
    as.data.frame(.)
table1 <- table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(as.data.frame(table1), 'results/tables/followuptable_shotgun.csv')

table1 <- helius %>% 
    filter(shotgun == TRUE) %>% droplevels(.) %>% 
    filter(timepoint == "baseline") %>% 
    dplyr::select(Age, Sex, Ethnicity, MigrGen, FUtime, BMI, Smoking, AlcCons, 
           DM, SBP, DBP, HT_BPMed, MetSyn, 
           PPI, Metformin, Statins, TC, LDL, Trig, HbA1c, SCORECVDmortNL,
           timepoint) %>% 
    CreateTableOne(data=., strata = 'Ethnicity', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig")) %>% 
    as.data.frame(.)
table1 <- table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(as.data.frame(table1), 'results/tables/followuptable_shotgun_baseline.csv')

helius <- readRDS('data/clinicaldata_long.RDS')
mbids <- read.csv2('data/16s/ids_16s_paired.csv') %>% dplyr::select(ID = x) %>% 
    mutate(sampleID_baseline = str_c("HELIBA_",str_remove(ID, "S")),
          `sampleID_follow-up` = str_c("HELIFU_",str_remove(ID, "S"))) %>% 
    pivot_longer(., 2:3, names_to = c("var", "timepoint"), values_to = "sampleID", names_sep = "_") 
mbids$microbiome_16s <- TRUE
helius <- left_join(helius, mbids, by = c("sampleID", "ID", "timepoint")) %>% 
    mutate(across(c('microbiome_16s', 'shotgun'), ~case_when(
        is.na(.x) ~ FALSE,
        .default = TRUE
    )))

table1 <- helius %>% 
    filter(microbiome_16s == TRUE) %>% 
    dplyr::select(Age, Sex, Ethnicity, FUtime, BMI, Smoking, AlcCons, 
                  DM, SBP, DBP, HT_BPMed, MetSyn, HT_BPMed,
                  PPI, Metformin, Statins, TC, LDL, Trig, HbA1c, SCORECVDmortNL, timepoint) %>% 
    CreateTableOne(data=., strata = 'timepoint', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig")) %>% 
    as.data.frame(.)
table1 <- table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(as.data.frame(table1), 'results/tables/followuptable_16s.csv')

table1 <- helius %>% 
    filter(shotgun == TRUE) %>% droplevels(.) %>% 
    filter(timepoint == "baseline") %>% 
    dplyr::select(Age, Sex, Ethnicity, MigrGen, FUtime, BMI, Smoking, AlcCons, 
                  DM, SBP, DBP, HT_BPMed, MetSyn, 
                  PPI, Metformin, Statins, TC, LDL, Trig, HbA1c, SCORECVDmortNL,
                  timepoint) %>% 
    CreateTableOne(data=., strata = 'Ethnicity', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig")) %>% 
    as.data.frame(.)
table1 <- table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(as.data.frame(table1), 'results/tables/followuptable_shotgun_baseline.csv')

#### Baseline prediction tables ####
## Diabetes
dfdelta <- readRDS('data/clinicaldata_delta.RDS') %>% 
    left_join(., metids, by = "ID") %>% 
    left_join(., mbids, by = 'ID') %>% 
    left_join(., shotids, by = 'ID') %>% 
    mutate(across(c('microbiome_16s', 'shotgun'), ~case_when(
        is.na(.x) ~ FALSE,
        .default = TRUE
    )))

dmnew <- dfdelta %>% filter(!is.na(DM_new))
tabledm <- dmnew %>% 
    dplyr::select(Age_baseline, Sex, Ethnicity, FUtime, BMI_baseline, BMI_delta, Smoking_baseline, AlcCons_baseline, 
           DM_baseline, SBP_baseline, DBP_baseline, MetSyn_baseline, HT_BPMed_baseline, 
           PPI_baseline, Metformin_baseline, Statins_baseline, 
           TC_baseline, LDL_baseline, Trig_baseline, HbA1c_baseline, SCORECVDmortNL_baseline,
           microbiome_16s, shotgun, DM_new) %>% 
    CreateTableOne(data=., strata = 'DM_new', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
tabledm <- tabledm %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

htnew <- dfdelta %>% filter(!is.na(HT_new))
tableht <- htnew %>% 
    dplyr::select(Age_baseline, Sex, Ethnicity, FUtime, BMI_baseline, BMI_delta, Smoking_baseline, AlcCons_baseline, 
                  DM_baseline, SBP_baseline, DBP_baseline, HT_BPMed_baseline, 
                  MetSyn_baseline, PPI_baseline, Metformin_baseline, Statins_baseline, 
                  TC_baseline, LDL_baseline, Trig_baseline, HbA1c_baseline, SCORECVDmortNL_baseline,
                  microbiome_16s, shotgun, HT_new) %>% 
    CreateTableOne(data=., strata = 'HT_new', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
tableht <- tableht %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

metsynnew <- dfdelta %>% filter(!is.na(MetSyn_new))
tablemetsyn <- metsynnew %>% 
    dplyr::select(Age_baseline, Sex, Ethnicity, FUtime, BMI_baseline, BMI_delta, Smoking_baseline, AlcCons_baseline, 
                  DM_baseline, SBP_baseline, DBP_baseline, MetSyn_baseline, HT_BPMed_baseline, 
                  PPI_baseline, Metformin_baseline, Statins_baseline, 
                  TC_baseline, LDL_baseline, Trig_baseline, HbA1c_baseline, SCORECVDmortNL_baseline,
                  microbiome_16s, shotgun, MetSyn_new) %>% 
    CreateTableOne(data=., strata = 'MetSyn_new', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
tablemetsyn <- tablemetsyn %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

lldnew <- dfdelta %>% filter(!is.na(LLD_new))
tablelld <- metsynnew %>% 
    dplyr::select(Age_baseline, Sex, Ethnicity, FUtime, BMI_baseline, BMI_delta, Smoking_baseline, AlcCons_baseline, 
                  DM_baseline, SBP_baseline, DBP_baseline, MetSyn_baseline, HT_BPMed_baseline, 
                  PPI_baseline, Metformin_baseline, Statins_baseline, 
                  TC_baseline, LDL_baseline, Trig_baseline, HbA1c_baseline, SCORECVDmortNL_baseline,
                  microbiome_16s, shotgun, LLD_new) %>% 
    CreateTableOne(data=., strata = 'LLD_new', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
tablelld <- tablelld %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

tablepred <- cbind(tabledm, tableht, tablemetsyn, tablelld)
write.csv2(as.data.frame(tablepred), 'results/tables/tables_prediction.csv')