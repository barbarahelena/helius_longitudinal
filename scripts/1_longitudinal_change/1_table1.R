# Table 1

# Libraries
library(tableone)
library(dplyr)
library(stringr)
library(phyloseq)

## Data
helius <- readRDS('data/clinicaldata/clinicaldata_long.RDS')
mbids <- read.csv2('data/16s/ids_16s_paired.csv') |> dplyr::select(ID = x)
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
    dplyr::select(Age, Sex, EthnicityTot, MigrGen, FUtime, BMI, Smoking_current, AlcCons,
           Alcohol, ExerciseNorm, DiscrMean_baseline,
           DM, HT_BPMed, MetSyn, Dyslipidemia,
           SBP, DBP,
           PPI, Metformin, Statins, AntiHT, GlucLowDrugs, PsychoMed, Cortico,
           TC, LDL, Trig, HbA1c,
           TotalCalories, Carbohydrates, Protein, Fiber, FattyAcids, Protein_animal, Sodium_g,
           microbiome_16s, shotgun, timepoint) %>%
    CreateTableOne(data=., strata = 'timepoint', test = FALSE) %>% 
    print(nonnormal=c("Trig")) %>% 
    as.data.frame(.)
table1 <- table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(as.data.frame(table1), 'results/tables/table1_16s_timepoints.csv')

##### Table 1 per EthnicityTot stratum #####
ethnicity_strata <- unique(na.omit(helius$EthnicityTot))

for (eth in ethnicity_strata) {
    eth_label <- gsub("[^A-Za-z0-9]", "_", eth)
    t <- helius %>%
        filter(EthnicityTot == eth) %>%
        dplyr::select(Age, Sex, MigrGen, FUtime, BMI, Smoking_current, AlcCons,
               Alcohol, ExerciseNorm, DiscrMean_baseline,
               DM, HT_BPMed, MetSyn, Dyslipidemia,
               SBP, DBP,
               PPI, Metformin, Statins, AntiHT, GlucLowDrugs, PsychoMed, Cortico,
               TC, LDL, Trig, HbA1c,
               TotalCalories, Carbohydrates, Protein, Fiber, FattyAcids, Protein_animal, Sodium_g,
               timepoint) %>%
        CreateTableOne(data = ., strata = 'timepoint', test = FALSE) %>%
        print(nonnormal = c("Trig")) %>%
        as.data.frame(.)
    t <- t %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
    write.csv2(t, paste0('results/tables/table1_16s_', eth_label, '.csv'))
}

table2 <- helius %>% filter(shotgun == TRUE) %>%
    dplyr::select(Age, Sex, Ethnicity, MigrGen, FUtime, BMI, Smoking_current, AlcCons,
           Alcohol, ExerciseNorm, DiscrMean_baseline,
           DM, HT_BPMed, MetSyn, Dyslipidemia,
           SBP, DBP,
           PPI, Metformin, Statins, AntiHT, GlucLowDrugs, PsychoMed, Cortico,
           TC, LDL, Trig, HbA1c,
           TotalCalories, Carbohydrates, Protein, Fiber, FattyAcids, Protein_animal, Sodium_g,
           timepoint) %>% 
    CreateTableOne(data=., strata = 'timepoint', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig")) %>% 
    as.data.frame(.)
table2 <- table2 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(as.data.frame(table2), 'results/tables/table2_shotgunset.csv')

##### Table 2 per EthnicityTot stratum #####
for (eth in c("Dutch", "South-Asian Surinamese")) {
    eth_label <- gsub("[^A-Za-z0-9]", "_", eth)
    t <- helius %>%
        filter(shotgun == TRUE, EthnicityTot == eth) %>%
        dplyr::select(Age, Sex, MigrGen, FUtime, BMI, Smoking_current, AlcCons,
               Alcohol, ExerciseNorm, DiscrMean_baseline,
               DM, HT_BPMed, MetSyn, Dyslipidemia,
               SBP, DBP,
               PPI, Metformin, Statins, AntiHT, GlucLowDrugs, PsychoMed, Cortico,
               TC, LDL, Trig, HbA1c,
               TotalCalories, Carbohydrates, Protein, Fiber, FattyAcids, Protein_animal, Sodium_g,
               timepoint) %>%
        CreateTableOne(data = ., strata = 'timepoint', test = FALSE) %>%
        print(nonnormal = c("SCORECVDmortNL", "Trig")) %>%
        as.data.frame(.)
    t <- t %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
    write.csv2(t, paste0('results/tables/table2_shotgunset_', eth_label, '.csv'))
}
