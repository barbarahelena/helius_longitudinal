# Table 1

# Libraries
library(tableone)
library(dplyr)
library(stringr)
library(phyloseq)

## Data
helius <- readRDS('data/clinicaldata/clinicaldata_long.RDS')
mbids <- read.csv2('data/16s/ids_16s_paired.csv') |> select(ID = x)
mbids$microbiome_16s <- TRUE
shotids <- read.csv('data/shotgun/shotgunseq_ids.csv') %>% dplyr::select(ID = x)
shotids$shotgun <- TRUE
shotids$ID <- str_c("S", shotids$ID)

helius <- left_join(mbids, helius, by = "ID") |> 
    left_join(shotids, by = 'ID') |> 
    mutate(across(c('microbiome_16s', 'shotgun'), ~case_when(
        is.na(.x) ~ FALSE,
        .default = TRUE
    )))
names(helius)
helius |> filter(microbiome_16s == TRUE) |> nrow()

##### Table 1 #####
table1 <- helius %>% 
    dplyr::select(Age, Sex, EthnicityTot, MigrGen, FUtime, BMI, Smoking, AlcCons, 
           DM, SBP, DBP, HT_BPMed, MetSyn, Dyslipidemia,
           PPI, Metformin, Statins, TC, LDL, Trig, 
           HbA1c, SCORECVDmortNL,
           microbiome_16s, shotgun, timepoint) %>% 
    CreateTableOne(data=., strata = 'timepoint', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig")) %>% 
    as.data.frame(.)
table1 <- table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(as.data.frame(table1), 'results/tables/table1_16s_timepoints.csv')

table2 <- helius %>% filter(shotgun == TRUE) %>%
    dplyr::select(Age, Sex, Ethnicity, MigrGen, FUtime, BMI, Smoking, AlcCons, 
           DM, SBP, DBP, HT_BPMed, MetSyn, 
           PPI, Metformin, Statins, TC, 
           LDL, Trig, HbA1c, SCORECVDmortNL,
           timepoint) %>% 
    CreateTableOne(data=., strata = 'timepoint', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig")) %>% 
    as.data.frame(.)
table2 <- table2 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(as.data.frame(table2), 'results/tables/table2_shotgunset.csv')
