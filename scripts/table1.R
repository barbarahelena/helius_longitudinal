# Table 1

# libraries
library(tableone)
library(dplyr)
library(ggplot2)
library(ggsci)

theme_Publication <- function(baselinese_size=14, baselinese_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(baselinese_size=baselinese_size, baselinese_family=baselinese_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.baselineckground = element_rect(colour = NA),
                plot.baselineckground = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.baselineckground=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 


## Data
helius <- readRDS('data/clinicaldata_wide.RDS')
metids <- read.csv2('data/metabolomics/ids_metabolomics.csv', header = TRUE) %>% dplyr::select(ID = x)
metids <- metids %>% filter(!str_detect(ID, "HELICOV"))
metids$metabolomics <- TRUE
metids$ID <- str_c("S", str_remove_all(metids$ID, "HELIFU_")) 
mbids <- read.csv2('data/16s/ids_16s_paired.csv')%>% dplyr::select(ID = x)
mbids$microbiome_16s <- TRUE
shotids <- read.csv('data/shotgun/helius_ids.csv') %>% dplyr::select(ID = x) %>% filter(!str_detect(ID, "HELIbaseline"))
shotids$shotgun <- TRUE
shotids$ID <- str_c("S", str_remove_all(shotids$ID, "HELIFU_"))

helius <- left_join(helius, metids, by = "ID") %>% 
    left_join(., mbids, by = 'ID') %>% 
    left_join(., shotids, by = 'ID') %>% 
    mutate(across(c('metabolomics', 'microbiome_16s', 'shotgun'), ~case_when(
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
    CreateTableOne(data=., strata = 'metabolomics', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
table2 <- table2 %>% dplyr::select(`metabolomics` = `TRUE`) %>% 
    mutate(across(everything(.), ~trimws(.x, which = "both")))

table3 <- helius %>% 
    dplyr::select(Age_baseline, Sex, Ethnicity, MigrGen, FUtime, BMI_baseline, Smoking_baseline, AlcCons_baseline, 
           DM_baseline, SBP_baseline, DBP_baseline, HT_BPMed_baseline, MetSyn_baseline, 
           PPI_baseline, Metformin_baseline, Statins_baseline, TC_baseline, LDL_baseline, 
           Trig_baseline, HbA1c_baseline, SCORECVDmortNL_baseline,
           metabolomics, microbiome_16s, shotgun) %>% 
    CreateTableOne(data=., strata = 'microbiome_16s', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
table3 <- table3 %>% dplyr::select(`16S` = `TRUE`) %>%  
            mutate(across(everything(.), ~trimws(.x, which = "both")))

table4 <- helius %>% 
    dplyr::select(Age_baseline, Sex, Ethnicity, MigrGen, FUtime, BMI_baseline, Smoking_baseline, AlcCons_baseline, 
           DM_baseline, SBP_baseline, DBP_baseline, HT_BPMed_baseline, MetSyn_baseline, 
           PPI_baseline, Metformin_baseline, Statins_baseline, TC_baseline, 
           LDL_baseline, Trig_baseline, HbA1c_baseline, SCORECVDmortNL_baseline,
           metabolomics, microbiome_16s, shotgun) %>% 
    CreateTableOne(data=., strata = 'shotgun', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
table4 <- table4 %>% dplyr::select(`shotgun` = `TRUE`) %>% 
    mutate(across(everything(.), ~trimws(.x, which = "both")))

table <- cbind(table1,table2,table3,table4)

write.csv2(as.data.frame(table), 'results/table1.csv')


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

write.csv2(as.data.frame(table1), 'results/followuptable.csv')

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
write.csv2(as.data.frame(table1), 'results/followuptable_shotgun.csv')

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
write.csv2(as.data.frame(table1), 'results/followuptable_shotgun_baseline.csv')

#### Baseline prediction tables ####
## Diabetes
dfdelta <- readRDS('data/clinicaldata_delta.RDS') %>% 
    left_join(., metids, by = "ID") %>% 
    left_join(., mbids, by = 'ID') %>% 
    left_join(., shotids, by = 'ID') %>% 
    mutate(across(c('metabolomics', 'microbiome_16s', 'shotgun'), ~case_when(
        is.na(.x) ~ FALSE,
        .default = TRUE
    )))

dmnew <- dfdelta %>% filter(!is.na(DM_new))
tabledm <- dmnew %>% 
    dplyr::select(Age_baseline, Sex, Ethnicity, FUtime, BMI_baseline, BMI_delta, Smoking_baseline, AlcCons_baseline, 
           DM_baseline, SBP_baseline, DBP_baseline, MetSyn_baseline, HT_BPMed_baseline, 
           PPI_baseline, Metformin_baseline, Statins_baseline, 
           TC_baseline, LDL_baseline, Trig_baseline, HbA1c_baseline, SCORECVDmortNL_baseline,
           metabolomics, microbiome_16s, shotgun, DM_new) %>% 
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
                  metabolomics, microbiome_16s, shotgun, HT_new) %>% 
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
                  metabolomics, microbiome_16s, shotgun, MetSyn_new) %>% 
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
                  metabolomics, microbiome_16s, shotgun, LLD_new) %>% 
    CreateTableOne(data=., strata = 'LLD_new', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig")) %>% 
    as.data.frame(.)
tablelld <- tablelld %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

tablepred <- cbind(tabledm, tableht, tablemetsyn, tablelld)
write.csv2(as.data.frame(tablepred), 'results/tables_prediction.csv')

#### BP plots ####
pl3 <- ggplot(helius, aes(x=timepoint, y=SBP))+
    geom_violin(aes(fill=timepoint))+
    scale_fill_manual(values = pal_lancet()(2), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'mmHg', title = "Systolic blood pressure")+
    ggpubr::stat_compare_means(method = "t.test", tip.length = 0, label="p.signif", hide.ns = TRUE)
pl3

pl4 <- ggplot(helius, aes(x=timepoint, y=DBP))+
    geom_violin(aes(fill=timepoint))+
    scale_fill_manual(values = pal_lancet()(2), guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'mmHg', title = "Diastolic blood pressure")+
    ggpubr::stat_compare_means(hide.ns = TRUE)
pl4

ggarrange(pl3, pl4, labels = c("A", "B"))
ggsave("results/bloodpressure.pdf", height = 5, width = 8)
