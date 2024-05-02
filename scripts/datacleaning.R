#### Data cleaning

## Libraries
library(tidyverse)
library(dplyr)
library(rio) 
library(haven)
library(phyloseq)
library(forcats)
library(stringr)

## Open HELIUS clinical data
df <- haven::read_sav("data/210517_HELIUS data Ulrika Boulund_2.sav") # more subjects than other set
# df2 <- haven::read_sav("data/210517_HELIUS data Ulrika Boulund.sav")
df2 <- rio::import("data/clinical_data_long_03_07_2023_N3443.csv") %>% 
    mutate(ID = str_c("S",as.character(Heliusnr)))

# Change type of variable 
yesnosmall <- function(x) fct_recode(x, "No"="nee", "Yes"="ja")
yesnocaps <- function(x) fct_recode(x, "No"="Nee", "Yes"="Ja")
zero_one <- function(x) fct_recode(x, "No"="0", "Yes"="1")

## Clean HELIUS dataframe
df_new <- df %>% 
    dplyr::select(# Demographics
                  ID=Heliusnr, Age_BA=H1_lft, Age_FU = H2_lft,
                  Sex=H1_geslacht, EthnicityTot=H1_EtnTotaal, 
                  Ethnicity = H1_etniciteit,
                  MigrGen = H1_MigrGeneratie, FUtime=H2_fu_time,
                  Edu = H1_Opleid, 
                  # Discrimination - not available at H2
                  DiscrSum_BA = H1_Discr_sumscore, DiscrMean_BA = H1_Discr_meanscore,
                  DiscrBin_BA = H1_AnyDiscrim, ResDuration_BA = H1_ResidenceDuration,
                  # Movement, smoking, alcohol - exercise not available at H2
                  ExerciseMinweek_BA = H1_Squash_totmwk,
                  ExerciseNorm_BA = H1_Squash_rlbew, ExerciseScore_BA = H1_Squash_totscor,
                  Smoking_BA = H1_Roken, PackYears_BA = H1_PackYears,
                  SmokingComb = H2_Roken_H1combined, Alcohol_BA = H1_AlcoholJN, Alcohol_FU = H2_AlcoholJN, 
                  AlcCons_FU=H2_AlcoholConsumption, AlcCons_BA = H1_AlcoholConsumption,
                  # CVR scores - not available of H2
                  SCORECVDmort_BA = H1_SCORE_CVDmort, SCORECVDmortNL_BA = H1_SCORE_CVDmort_NL,
                  SCORECVDtotNL_BA = H1_SCORE_CVDtot_NL,
                  Fram_BA = H1_FramScore, FramCVD_BA = H1_Fram_CVD,
                  # Physical exam
                  # Physical_FU=H2_LichamelijkOnderzoekJN,
                  Fatperc_BA = H1_BIA_FatPercent, WHR_BA = H1_LO_WHR, # not available of FU
                  BMI_FU=H2_LO_BMI, SBP_FU=H2_LO_GemBPSysZit, 
                  DBP_FU=H2_LO_GemBPDiaZit, HR_FU=H2_LO_GemBPHRZit,
                  BMI_BA = H1_LO_BMI, SBP_BA = H1_LO_GemBPSysZit,
                  DBP_BA = H1_LO_GemBPDiaZit, 
                  # Diagnoses
                  MetSyn_BA = H1_MetSyn_MetabolicSyndrome, 
                  DM_BA = H1_Diabetes_GlucMed, DM_FU = H2_DM_GlucMed,
                  # Cardiometabolic drugs
                  PPI_BA=H1_ProtPumpInh_Moritz, PPI_FU=H2_ProtPumpInh_Ulrika,
                  Metformin_BA = H1_Metformin_Moritz, Metformin_FU = H2_Metformin_Ulrika,
                  Statins_BA = H1_Statines,
                  # Other drugs
                  PsychoMed_BA = H1_Psychotroop, 
                  AB_FU = H2_Antibiotica, AB_BA = H1_Antibiotica,
                  # Lab
                  Trig_BA = H1_Lab_UitslagTRIG, TC_BA = H1_Lab_UitslagCHOL, HDL_BA=H1_Lab_UitslagHDLS,
                  LDL_BA = H1_Lab_uitslagRLDL, HbA1c_BA = H1_Lab_UitslagIH1C,
                  TC_FU = H2_Lab_UitslagCHOL, LDL_FU=H2_Lab_UitslagRLDL, HDL_FU=H2_Lab_UitslagHDLS, 
                  Trig_FU=H2_Lab_UitslagTRIG, HbA1c_FU=H2_Lab_UitslagIH1C, Glucose_FU = H2_Lab_UitslagGLUC,
                  ASAT_FU = H2_Lab_UitslagROT, ALAT_FU = H2_Lab_UitslagRPT, Trombo_FU = H2_Lab_UitslagHTRO,
                  GGT_FU = H2_Lab_UitslagRGGT,
                  # Fecal sample 
                  SampleAB_BA=H1_Feces_q2, SampleDiarrhoea_BA=H1_Feces_q3,
                  GeneticData = GeneticsGSA_QC
    )
df_new2 <- df_new %>% 
    mutate(across(where(is.character), ~na_if(., c("Missing", "Missing: n.v.t.", "niet ingevuld","nvt", 
                                       "No lab result", "Missing: not applicable",
                                       "Missing: not measured", "missing",
                                       "See comments lab results", "Low (<1 mmol/L)",
                                       "Low (<0,08 mmol/L)", "Low (<0,10 mmol/L)",
                                       "Smoking status unknown", "Number or duration unknown"))),
           ID = str_c("S",as.character(ID)),
           across(c("Sex","Ethnicity", "Smoking_BA", "SmokingComb", "EthnicityTot",
                    "MigrGen", "Edu", "AlcCons_FU", "AlcCons_BA", "PsychoMed_BA",
                    "AB_FU", "AB_BA", "PPI_BA","PPI_FU", "Metformin_BA",
                    "Metformin_FU", "Statins_BA", "MetSyn_BA", "DM_BA", "DM_FU",
                    "Alcohol_BA", "Alcohol_FU", "ExerciseNorm_BA","DiscrBin_BA",
                    "SampleAB_BA", "SampleDiarrhoea_BA"
                    ), as_factor),
           across(c("PsychoMed_BA", "Metformin_BA", "Metformin_FU",
                    "MetSyn_BA","DM_BA", "DM_FU",
                    "DiscrBin_BA", "Alcohol_FU", "Alcohol_BA"), yesnocaps),
           across(c("ExerciseNorm_BA", "Statins_BA", "AB_BA", "AB_FU"), yesnosmall),
           Sex = fct_recode(Sex, "Male" = "1", "Female" = "2"),
           across(c("SampleAB_BA", "SampleDiarrhoea_BA"), zero_one),
           EthnicityTot = forcats::fct_recode(EthnicityTot,
                                              "Other"="Other/unkown",
                                              "Other"="Other/unknown Surinamese"),
    ) %>%
    droplevels(.) %>% 
    mutate(
           across(where(is.numeric), as.numeric), # all other vars to numeric, do this last,
           #ID = str_c("S", ID)
    ) %>%
    # remove unused levels
    droplevels(.)
dim(df_new2)

# Clean phyloseq object
heliusmb <- readRDS("data/16s/phyloseq/rarefied/phyloseq_rarefied.RDS")
all(str_detect(sample_names(heliusmb), "_T1")) # TRUE
any(str_detect(sample_names(heliusmb), "_T2")) # FALSE (hence no duplicated sample IDs)
sample_names(heliusmb) <- str_replace(sample_names(heliusmb), "_T1", "") # remove _T1
sample_names(heliusmb)
all(sample_sums(heliusmb) == 15000) # all rarefied to 15,000 counts
timepoint <- case_when(
    str_detect(sample_names(heliusmb), "HELIBA") ~ "baseline",
    str_detect(sample_names(heliusmb), "HELIFU") ~ "follow-up"
) %>% as.factor(.)
missingdata <- case_when(
    (sample_names(heliusmb) %in% df_new$ID) == TRUE~ "available",
    (sample_names(heliusmb) %in% df_new$ID) == FALSE ~ "missing"
)
df <- data.frame(row.names = sample_names(heliusmb), 
                 timepoint,
                 missingdata,
                 ID = str_c("S",str_remove(str_remove(sample_names(heliusmb), "HELIFU_"), "HELIBA_")))
heliusmb <- phyloseq(heliusmb@otu_table, heliusmb@tax_table, heliusmb@refseq, heliusmb@phy_tree, sample_data(df))

# Check overlap samples
heliusfu <- prune_samples(heliusmb@sam_data$timepoint == "follow-up", heliusmb)
summary(heliusfu@sam_data$ID %in% df_new$ID)
missingfu <- heliusfu@sam_data$ID[which(!heliusfu@sam_data$ID %in% df_new$ID)]

heliusba <- prune_samples(heliusmb@sam_data$timepoint == "baseline", heliusmb)
summary(heliusba@sam_data$ID %in% df_new$ID)
missingba <- heliusba@sam_data$ID[which(!heliusba@sam_data$ID %in% df_new$ID)]
missingba[which(missingba %in% heliusfu@sam_data$ID)] # 1 subject S207389
any(str_detect(missingfu, "S207389")) # TRUE so the subject is already in the follow-up missing list

# Select samples that are in dataset
heliusmb2 <- prune_samples(heliusmb@sam_data$ID %in% df_new2$ID, heliusmb)
heliusmb2
df_new3 <- df_new2 %>% filter(ID %in% heliusmb2@sam_data$ID)
write.csv2(df_new3$ID, 'data/16s/ids_16s_paired.csv')
all(df_new3$ID %in% heliusmb2@sam_data$ID) # TRUE
all(heliusmb2@sam_data$ID %in% df_new3$ID) # TRUE

# Pivot longer clinical data
names(df_new2)
coldouble <- colnames(df_new2)[which(str_detect(colnames(df_new2), "_FU") | str_detect(colnames(df_new2), "_BA"))]
df_new_long <- df_new2 %>% pivot_longer(., cols = coldouble, names_to = c(".value", "timepoint"),
                                        names_sep = "_")
df_new_long <- df_new_long %>% 
    mutate(timepoint = case_when(
        timepoint == "BA" ~ "baseline",
        timepoint == "FU" ~"follow-up"
    ),
            timepoint = as.factor(timepoint),
    sampleID = str_remove(ID, "S"),
    sampleID = case_when(
        timepoint == "baseline" ~ str_c("HELIBA_", sampleID),
        timepoint == "follow-up" ~ str_c("HELIFU_", sampleID)
    )
    )

samdata <- as(heliusmb2@sam_data, "data.frame")
samdata$sampleID <- rownames(samdata)
df_new_long2 <- left_join(samdata, df_new_long, by = c("sampleID", "timepoint", "ID"))
rownames(df_new_long2) <- df_new_long2$sampleID
dfba <- df_new_long2 %>% filter(timepoint == "baseline")
dffu <- df_new_long2 %>% filter(timepoint == "follow-up") %>% 
    filter(ID %in% dfba$ID)
df_new_long2 <- df_new_long2 %>% filter(ID %in% dffu$ID)

heliusmb2 <- phyloseq(heliusmb2@otu_table, heliusmb2@tax_table, heliusmb2@refseq, heliusmb2@phy_tree, 
                      sample_data(df_new_long2))

## Get delta variables
df_wide <- pivot_wider(df_new_long, id_cols = c(1:9), names_from = "timepoint",
                       values_from = c(11:52))
df_wide <- df_wide[,colSums(is.na(df_wide))<nrow(df_wide)]
df_wide_delta <- df_wide %>% 
    mutate(
        Age_delta = `Age_follow-up` - Age_baseline,
        SBP_delta = `SBP_follow-up` - SBP_baseline,
        DBP_delta = `DBP_follow-up` - DBP_baseline,
        BMI_delta = `BMI_follow-up` - BMI_baseline,
        TC_delta = `TC_follow-up` - TC_baseline,
        HDL_delta = `HDL_follow-up` - HDL_baseline,
        HbA1c_delta = `HbA1c_follow-up` - HbA1c_baseline,
        Trig_delta = `Trig_follow-up` - Trig_baseline,
        LDL_delta = `LDL_follow-up` - LDL_baseline
    ) 
df_deltas <- df_wide_delta %>% select(1:9, contains("delta"))


## Save files
saveRDS(df_new2, file = "data/clinicaldata.RDS")
saveRDS(df_wide, file = "data/clinicaldata_wide.RDS")
saveRDS(df_wide_delta, file = "data/clinicaldata_delta.RDS")
saveRDS(df_new_long, file = "data/clinicaldata_long.RDS")
saveRDS(heliusmb2, file = "data/phyloseq_sampledata.RDS")
