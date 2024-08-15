# Data cleaning of clinical data, 16S, shotgun
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(rio) 
library(haven)
library(phyloseq)

#### Clinical data ####
## Open HELIUS clinical data
df <- haven::read_sav("data/210517_HELIUS data Ulrika Boulund_2.sav") # more subjects than other set
df2 <- haven::read_sav("data/240411_HELIUS data Barbara Verhaar.sav")
df3 <- haven::read_sav("data/220712_HELIUS data Barbara Verhaar - Cov1 groep en datum.sav") # covid dates
df4 <- haven::read_sav("data/EGA_standaardvariabelen.sav")
names(df2)[which(!names(df2) %in% names(df))]
df <- df %>% dplyr::select("Heliusnr", !names(df)[which(names(df) %in% names(df2))])
dftot <- full_join(df2, df, by = "Heliusnr") %>% full_join(., df3, by = "Heliusnr") %>% 
    full_join(., df4)
all(df3$Heliusnr %in% df$Heliusnr) # TRUE covid set does not have new subjects
length(df2$Heliusnr[which(!df2$Heliusnr %in% df$Heliusnr)]) # 2 subjects in original set, not in df2
length(df$Heliusnr[which(!df$Heliusnr %in% df2$Heliusnr)]) # 16 subjects not in original set, but in df2

# Change type of variable 
yesnosmall <- function(x) fct_recode(x, "No"="nee", "Yes"="ja")
yesnocaps <- function(x) fct_recode(x, "No"="Nee", "Yes"="Ja")
zero_one <- function(x) fct_recode(x, "No"="0", "Yes"="1")

## Clean HELIUS dataframe
df_new <- dftot %>% 
    dplyr::select(# Demographics
                  ID=Heliusnr, FUtime=H2_fu_time,
                  Age_BA=H1_lft, Age_FU = H2_lft,
                  Sex=H1_geslacht,  Edu = H1_Opleid, EthnicityTot=H1_EtnTotaal, 
                  Ethnicity = H1_etniciteit, MigrGen = H1_MigrGeneratie, 
                  # Discrimination - not available at H2
                  DiscrSum_BA = H1_Discr_sumscore, DiscrMean_BA = H1_Discr_meanscore,
                  DiscrBin_BA = H1_AnyDiscrim, ResDuration_BA = H1_ResidenceDuration,
                  # Movement, smoking, alcohol - exercise not available at H2
                  ExerciseMinweek_BA = H1_Squash_totmwk,
                  ExerciseNorm_BA = H1_Squash_rlbew, ExerciseScore_BA = H1_Squash_totscor,
                  # Intox
                  Smoking_BA = H1_Roken, PackYears_BA = H1_PackYears, Smoking_FU = H2_Roken_H1combined, 
                  Alcohol_BA = H1_AlcoholJN, Alcohol_FU = H2_AlcoholJN, 
                  AlcCons_BA = H1_AlcoholConsumption, AlcCons_FU=H2_AlcoholConsumption, 
                  # CVR scores - not available of H2
                  SCORECVDmort_BA = H1_SCORE_CVDmort, SCORECVDmort_FU = H2_SCORE1_RM_CVDmort,
                  SCORECVDmortNL_BA = H1_SCORE_CVDmort_NL, SCORECVDmortNL_FU = H2_SCORE1_RM_CVDmort_NL,
                  SCORECVDtotNL_BA = H1_SCORE_CVDtot_NL, SCORECVDtotNL_FU = H2_SCORE1_RM_CVDtot_NL,
                  Fram_BA = H1_FramScore, Fram_FU = H2_FramScore,
                  FramCVD_BA = H1_Fram_CVD, FramCVD_FU = H2_Fram_CVD,
                  # Physical exam
                  Physical_BA = H1_LichamelijkOnderzoekJN, Physical_FU=H2_LichamelijkOnderzoekJN,
                  Fatperc_BA = H1_BIA_FatPercent, 
                  WHR_BA = H1_LO_WHR, WHR_FU = H2_LO_WHR,
                  BMI_BA = H1_LO_BMI, BMI_FU=H2_LO_BMI, 
                  SBP_BA = H1_LO_GemBPSysZit, SBP_FU=H2_LO_GemBPSysZit, 
                  DBP_BA = H1_LO_GemBPDiaZit, DBP_FU=H2_LO_GemBPDiaZit, 
                  HR_BA = H2_LO_GemBPHRZit, HR_FU=H2_LO_GemBPHRZit, 
                  # Hypertension
                  HT_Self_BA = H1_HT_Self, HT_Self_FU = H2_HT_Self_H1combined,
                  HT_BP_BA = H1_HT_BP, HT_BP_FU = H2_HT_BP,
                  HT_SelfBP_BA = H1_HT_SelfBP, HT_SelfBP_FU = H2_HT_Self_H1combined,
                  HT_SelfBPMed_BA = H1_HT_SelfBPMed, HT_SelfBPMed_FU = H2_HT_SelfBPMed_H1combined,
                  HT_BPMed_BA = H1_HT_BPMed, HT_BPMed_FU = H2_HT_BPMed,
                  AntiHT_BA = H1_Antihypertensiva, AntiHT_FU = H2_Antihypertensiva,
                  # Cardiovascular complications
                  Claud_BA = H1_CI_Rose_corrected, Inf_BA = H1_possINF_Rose, 
                  AngPec_BA = H1_AP_Rose_corrected, CVD_BA = H1_CVD_Rose_corrected, ## all missing for FU
                  StrokeLoss_BA = H1_UitvalBer, StrokeLoss_FU = H2_UitvalBer_H1combined,
                  Stroke_Self_BA = H1_CVA_Self, # missing FU
                  IschStroke_BA = H1_Infarct, IschStroke_FU = H2_Infarct_H1combined,
                  CardInt_BA = H1_OperDotByp, CardInt_FU = H2_OperDotByp_H1combined,
                  MI_BA = H1_MI_Self, # missing FU
                  # Metabolic diagnoses
                  MetSyn_BA = H1_MetSyn_MetabolicSyndrome, MetSyn_FU = H2_MetSyn_MetabolicSyndrome,
                  MetSyn_Obesity_BA = H1_MetSyn_CentralObesity, MetSyn_Obesity_FU = H2_MetSyn_CentralObesity, 
                  MetSyn_HighTG_BA = H1_MetSyn_HighTriglyceride, MetSyn_HighTG_FU = H2_MetSyn_HighTriglyceride, 
                  MetSyn_LowHDL_BA =H1_MetSyn_LowHDL, MetSyn_LowHDL_FU =H2_MetSyn_LowHDL,
                  MetSyn_HighBP_BA = H1_MetSyn_HighBP, MetSyn_HighBP_FU = H2_MetSyn_HighBP, 
                  MetSyn_HighGluc_BA = H1_MetSyn_HighGluc, MetSyn_HighGluc_FU = H2_MetSyn_HighGluc,
                  DM_BA = H1_Diabetes_GlucMed, DM_FU = H2_DM_GlucMed,
                  Alb_FU = H2_Microalbuminurie, Albstage_FU = H2_ACR_KDIGO,
                  # Cardiometabolic drugs
                  AntiHTTot_BA = H1_AntihypertensivaC02, AntiHTTot_FU = H2_AntihypertensivaC02, 
                  Diuretics_BA = H1_AntihypertensivaC03, Diuretics_FU = H2_AntihypertensivaC03,
                  CalciumAnt_BA = H1_AntihypertensivaC08, CalciumAnt_FU = H2_AntihypertensivaC08, 
                  BetaBlocker_BA = H1_AntihypertensivaC07, BetaBlocker_FU = H2_AntihypertensivaC07,
                  RAASi_BA = H1_AntihypertensivaC09, RAASi_FU = H2_AntihypertensivaC09,
                  LLD_BA = H1_Antilipaemica, LLD_FU = H2_Antilipaemica,
                  PPI_BA=H1_ProtPumpInh_Moritz, PPI_FU=H2_ProtPumpInh_Ulrika,
                  Metformin_BA = H1_Metformin_Moritz, Metformin_FU = H2_Metformin_Ulrika,
                  GlucLowDrugs_BA = H1_Diabetesmiddelen, GlucLowDrugs_FU = H2_Diabetesmiddelen,
                  Statins_BA = H1_Statines, # missing FU
                  # Other drugs
                  PsychoMed_BA = H1_Psychotroop, 
                  AB_BA = H1_Antibiotica, AB_FU = H2_Antibiotica, 
                  Cortico_BA = H1_Corticosteroiden, Cortico_FU = H2_Corticosteroiden,
                  # Lab
                  Trig_BA = H1_Lab_UitslagTRIG, Trig_FU=H2_Lab_UitslagTRIG, 
                  TC_BA = H1_Lab_UitslagCHOL, TC_FU = H2_Lab_UitslagCHOL, 
                  HDL_BA=H1_Lab_UitslagHDLS, HDL_FU=H2_Lab_UitslagHDLS, 
                  LDL_BA = H1_Lab_uitslagRLDL, LDL_FU=H2_Lab_UitslagRLDL, 
                  HbA1c_BA = H1_Lab_UitslagIH1C, HbA1c_FU=H2_Lab_UitslagIH1C, 
                  Glucose_FU = H2_Lab_UitslagGLUC, # missing BA
                  ASAT_FU = H2_Lab_UitslagROT, ALAT_FU = H2_Lab_UitslagRPT, # missing BA
                  GGT_FU = H2_Lab_UitslagRGGT, Trombo_FU = H2_Lab_UitslagHTRO, # missing BA
                  Microalb_BA = H1_Lab_UitslagMIAL, Microalb_FU = H2_Lab_UitslagMIAL,
                  UACR_BA = H1_Lab_uitslagMIKR, UACR_FU = H2_Lab_UitslagMIKR,
                  Kreat_BA = H1_Lab_UitslagKREA_HP, Kreat_FU = H2_Lab_UitslagKREA_HP,
                  UKreat_BA = H1_Lab_UitslagKREA_UP, UKreat_FU = H1_Lab_UitslagKREA_UP,
                  #eGFR_CKDEPI_BA = H1_CKDEPI_eGFR, eGFR_MDRD_BA = H1_MDRD_eGFR, missing BA 
                  eGFR_CKDEPI_FU = H2_CKDEPI_eGFR, eGFR_MDRD_FU = H2_MDRD_eGFR,
                  # Diet at baseline
                  TotalCalories_BA = ENKcal_Sum, Protein_BA = Prot_Sum, 
                  FattyAcids_BA = FattyAcidsTot_Sum, Fiber_BA = Fibre_Sum, 
                  Carbohydrates_BA = Carbo_Sum, Protein_animal_BA = Prot_ani_Sum, 
                  Sodium_g_BA = Natrium_intake_totaal_gram, Sodium_mmol_BA = Natrium_intake_totaal_mmol, 
                  # Fecal sample 
                  SampleAB_BA=H1_Feces_q2, SampleDiarrhoea_BA=H1_Feces_q3,
                  GeneticData = GeneticsGSA_QC
    )

df_new2 <- df_new %>% 
    mutate(across(where(is.character), ~na_if(., c("Missing", "Missing: n.v.t.", "niet ingevuld","nvt", 
                                       "No lab result", "Missing: not applicable", "missing",
                                       "Missing: not measured", "missing",
                                       "See comments lab results", "Low (<1 mmol/L)",
                                       "Low (<0,08 mmol/L)", "Low (<0,10 mmol/L)",
                                       "Smoking status unknown", "Number or duration unknown",
                                       "Rookstatus onbekend", "Rookduur en/of aantal onbekend"))),
           sampleID_BA = str_c("HELIBA_", ID),
           sampleID_FU = str_c("HELIFU_", ID),
           ID = str_c("S",as.character(ID)),
           across(c("FUtime", "Age_BA", "Age_FU", "DiscrSum_BA", "DiscrMean_BA",
                    "ExerciseScore_BA", "ExerciseMinweek_BA", "ResDuration_BA",
                    "SCORECVDmort_BA", "SCORECVDmort_FU",
                    "SCORECVDmortNL_BA", "SCORECVDmortNL_FU", "SCORECVDtotNL_BA", "SCORECVDtotNL_FU",
                    "Fram_BA", "Fram_FU", "FramCVD_BA", "FramCVD_FU", "Fatperc_BA", "WHR_BA",
                    "WHR_FU", "BMI_BA", "BMI_FU", "SBP_BA", "SBP_FU", "DBP_BA", "DBP_FU", 
                    "HR_BA", "HR_FU", "Trig_BA", "Trig_FU", "TC_BA", "TC_FU", "HDL_BA", "HDL_FU",
                    "LDL_BA", "LDL_FU", "HbA1c_BA", "HbA1c_FU", "Glucose_FU", "ASAT_FU", "ALAT_FU",
                    "GGT_FU", "Trombo_FU", "Microalb_BA", "Microalb_FU", "UACR_BA", "UACR_FU",
                    "Kreat_BA", "Kreat_FU", "UKreat_BA", "UKreat_FU", "eGFR_CKDEPI_FU", "eGFR_MDRD_FU",
                    "TotalCalories_BA", "Protein_BA", "FattyAcids_BA", "Fiber_BA", "Carbohydrates_BA",
                    "Protein_animal_BA", "Sodium_g_BA", "Sodium_mmol_BA"), as.numeric), 
           across(where(haven::is.labelled), ~haven::as_factor(.x, levels = "labels")),
           across(c("HT_BP_BA", "HT_BPMed_BA", "HT_SelfBP_BA", "HT_SelfBPMed_BA", "HT_Self_BA",
                    "Alcohol_FU", "LLD_BA", "LLD_FU", "AB_BA", "AB_FU",
                    "RAASi_BA", "RAASi_FU", "BetaBlocker_BA", "BetaBlocker_FU",
                    "CalciumAnt_BA", "CalciumAnt_FU", "Diuretics_BA", "Diuretics_FU",
                    "AntiHT_BA", "AntiHT_FU", "AntiHTTot_BA", "AntiHTTot_FU", "Alb_FU",
                    "DM_BA", "MI_BA", "CardInt_BA", "CardInt_FU",
                    "IschStroke_BA", "IschStroke_FU","StrokeLoss_BA", "StrokeLoss_FU",
                    "Stroke_Self_BA","CVD_BA", "AngPec_BA", "Inf_BA", "Claud_BA"), yesnocaps),
           across(c("Alcohol_BA", "MetSyn_BA", "MetSyn_Obesity_BA",
                    "MetSyn_HighTG_BA", "MetSyn_LowHDL_BA", "MetSyn_HighGluc_BA",
                    "MetSyn_HighGluc_BA", "MetSyn_HighBP_BA"), yesnosmall),
           across(c("ExerciseNorm_BA", "Statins_BA"), ~fct_recode(.x, "No"="no", "Yes"="yes")),
           Sex = fct_recode(Sex, "Male" = "man", "Female" = "vrouw"),
           EthnicityTot = forcats::fct_recode(EthnicityTot, "Dutch" = "NL", "South-Asian Surinamese" = "Hind",
                                              "African Surinamese" = "Creools", "South-Asian Surinamese" = "Javaans",
                                              "Ghanaian" = "Ghanees", "Turkish" = "Turks", "Moroccan" = "Marokkaans",
                                              "Other"="Anders/onbekend",
                                              "Other"="Sur anders/onbekend"),
           across(c("Smoking_BA", "Smoking_FU"), ~fct_recode(.x, "Yes" = "Ja", 
                                                             "Former smoking" = "Nee, maar vroeger wel",
                                                             "Never" = "Nee, ik heb nooit gerookt",
                                                             "Never" = "Nee, nooit gerookt")),
         AlcCons_FU = fct_recode(AlcCons_FU,"low (men 0-4 gl/w, women 0-2 gl/w)" = "Low (men 0-4 gl/w, women 0-2 gl/w)",
                                "moderate (men 5-14 gl/w, women 3-7 gl/w)"="Moderate (men 5-14 gl/w, women 3-7 gl/w)",
                                "high (men >14 gl/w, women >7 gl/wk)" = "High (men >14 gl/w, women >7 gl/wk)"),
         Ethnicity = fct_recode(Ethnicity, "Dutch" = "Nederlands", "Surinamese" = "Surinaams",
                                "Turkish" = "Turks", "Moroccan" = "Marokkaans", "Ghanaian" = "Ghanees",
                                "Unknown" = "Onbekend", "Other" = "Anders"),
         FUtime = as.numeric(dmonths(FUtime), "years"),
         AgeDecade_BA = cut(Age_BA, 
                           breaks = c(18, 40, 50, 71), 
                           labels = paste(c(18, 40, 50), c(39, 49, 70), sep = "-"), 
                           right = FALSE),
         AgeDecade_FU = cut(Age_FU, 
                            breaks = c(18, 40, 50, 78),
                            labels = paste(c(18, 40, 50), c(39, 49, 78), sep = "-"), 
                            right = FALSE)
    ) %>%
    droplevels(.)
dim(df_new2)

# Pivot longer clinical data
names(df_new2)
coldouble <- colnames(df_new2)[which(str_detect(colnames(df_new2), "_FU") | str_detect(colnames(df_new2), "_BA"))]
df_new_long <- df_new2 %>% pivot_longer(., cols = all_of(coldouble), names_to = c(".value", "timepoint"),
                                        names_pattern = "(.*)_([A-Z]+)$")
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

# Get delta variables
df_wide <- pivot_wider(df_new_long, id_cols = c(1:8), names_from = "timepoint",
                       values_from = c(10:ncol(df_new_long)))
df_wide <- df_wide[,colSums(is.na(df_wide))<nrow(df_wide)]

deltafactor <- function(x,y){
    case_when(x == "No" & y == "No" ~ "No",
              x == "No" & y == "Yes" ~ "Yes",
              .default = NA)
}
df_wide_delta <- df_wide %>% 
    mutate(
        Age_delta = `Age_follow-up` - Age_baseline,
        SBP_delta = `SBP_follow-up` - SBP_baseline,
        DBP_delta = `DBP_follow-up` - DBP_baseline,
        BMI_delta = `BMI_follow-up` - BMI_baseline,
        WHR_delta = `WHR_follow-up` - WHR_baseline,
        TC_delta = `TC_follow-up` - TC_baseline,
        HDL_delta = `HDL_follow-up` - HDL_baseline,
        HbA1c_delta = `HbA1c_follow-up` - HbA1c_baseline,
        Trig_delta = `Trig_follow-up` - Trig_baseline,
        LDL_delta = `LDL_follow-up` - LDL_baseline,
        Fram_delta = `Fram_follow-up` - Fram_baseline,
        FramCVD_delta = `FramCVD_follow-up` - FramCVD_baseline,
        SCORECVDmort_delta = `SCORECVDmort_follow-up` - SCORECVDmort_baseline,
        SCORECVDtotNL_delta = `SCORECVDtotNL_follow-up` - SCORECVDtotNL_baseline,
        SCORECVDmortNL_delta = `SCORECVDmortNL_follow-up` - SCORECVDmortNL_baseline
    ) %>% 
    mutate(
        DM_new = deltafactor(DM_baseline, `DM_follow-up`),
        HT_new = deltafactor(HT_BPMed_baseline, `HT_BPMed_follow-up`),
        StrokeLoss_new = deltafactor(StrokeLoss_baseline, `StrokeLoss_follow-up`),
        IschStroke_new = deltafactor(IschStroke_baseline, `IschStroke_follow-up`),
        MetSyn_new = deltafactor(MetSyn_baseline, `MetSyn_follow-up`),
        CardInt_new = deltafactor(CardInt_baseline, `CardInt_follow-up`),
        LLD_new = deltafactor(LLD_baseline, `LLD_follow-up`)
        )

coldouble <- colnames(df_wide_delta)[which(str_detect(colnames(df_wide_delta), "_follow-up") | 
                                               str_detect(colnames(df_wide_delta), "_baseline"))]
coldouble <- coldouble[c(1:2, 5:length(coldouble))]
df_long <- df_wide_delta %>% 
    pivot_longer(., cols = all_of(coldouble), names_to = c(".value", "timepoint"),
                                        names_pattern = "(.*)_([a-z-]+)$")

## Save files
saveRDS(df_wide_delta, file = "data/clinicaldata_wide.RDS")
saveRDS(df_long, file = "data/clinicaldata_long.RDS")


#### 16S ####
# Overlap GWAS and 16s/shotgun
gwas <- rio::import("data/GWAS_ids.txt")
gwas$ID <- str_c("S", gwas$IID)
summary(df_new2$ID %in% gwas$ID)
heliussg <- rio::import("data/shotgun/combined_table.tsv")
sgids <- colnames(heliussg)[which(str_detect(colnames(heliussg), "HELIBA"))]
sgids <- str_c("S", str_remove_all(sgids, "HELIBA_"))
summary(sgids %in% gwas$ID)

# Clean phyloseq object
heliusmb <- readRDS("data/16s/phyloseq/rarefied/phyloseq_rarefied.RDS")
all(str_detect(sample_names(heliusmb), "_T1")) # TRUE
any(str_detect(sample_names(heliusmb), "_T2")) # FALSE (hence no duplicated sample IDs)
sample_names(heliusmb) <- str_replace(sample_names(heliusmb), "_T1", "") # remove _T1
sample_names(heliusmb)
idsfu <- sample_names(heliusmb)[which(str_detect(sample_names(heliusmb), "HELIFU_"))]
all(sample_sums(heliusmb) == 15000) # all rarefied to 15,000 counts
dffu <- df_new2 %>% filter(sampleID_FU %in% sample_names(heliusmb))
table(dffu$AgeDecade_FU, dffu$Ethnicity, dffu$Sex)
table(dffu$EthnicityTot)
df_new3 <- df_new2 %>% filter(sampleID_BA %in% sample_names(heliusmb) | sampleID_FU %in% sample_names(heliusmb))
write.csv2(df_new3$ID, 'data/16s/ids_16s.csv')

# Select samples that are in clinical dataset
heliusmb2 <- prune_samples(sample_names(heliusmb) %in% df_new_long$sampleID, heliusmb)
heliusmb2

# Make baseline mb set for CBS enviroment
cbs_ids <- rio::import("data/240411_HELIUS data Barbara Verhaar_Heliusnrs.csv") %>% 
    dplyr::select(ID = V1) %>% 
    mutate(ID = str_remove(ID, "_T1"))
heliusmb_baseline <- prune_samples(str_detect(sample_names(heliusmb), "HELIBA_"), heliusmb)
heliusmb_baseline <- prune_samples(sample_names(heliusmb) %in% cbs_ids$ID, heliusmb)
heliusmb_baseline <- prune_taxa(taxa_sums(heliusmb_baseline) > 0, heliusmb_baseline)
sample_names(heliusmb_baseline) <- str_replace(sample_names(heliusmb_baseline), "HELIBA_", "") # remove HELIBA_
heliusmb_asv <- as.data.frame(t(as(heliusmb_baseline@otu_table, "matrix")))
write.csv2(heliusmb_asv, "data/CBS/asv_table_heliusba.csv")
heliusmb_tax <- as.data.frame(as(heliusmb_baseline@tax_table, "matrix"))
write.csv2(heliusmb_tax, "data/CBS/tax_table_heliusba.csv")
Biostrings::writeXStringSet(heliusmb_baseline@refseq, "data/CBS/asvs_baseline.fna", append=FALSE,
                            compress=FALSE, compression_level=NA, format="fasta")
ape::write.tree(heliusmb_baseline@phy_tree, "data/CBS/tree_heliusba.tree")

# code to put files back together (test for CBS later)
# asvs <- read.csv2("data/CBS/asv_table_heliusba.csv")
# rownames(asvs) <- asvs$X
# asvs$X <- NULL
# taxs <- read.csv2("data/CBS/tax_table_heliusba.csv")
# rownames(taxs) <- taxs$ASV <- taxs$X
# taxs$X <- NULL
# refseqs <- Biostrings::readDNAStringSet("data/CBS/asvs_baseline.fna")
# trb <- ape::read.tree("data/CBS/tree_heliusba.tree")
# phynew <- phyloseq(otu_table(asvs, taxa_are_rows = FALSE), tax_table(as.matrix(taxs)), refseq(refseqs), phy_tree(trb))

# How many samples have paired data (baseline + follow-up data)
summary(str_detect(sample_names(heliusmb), "HELIBA")) # 2966 follow-up samples total
summary(str_detect(sample_names(heliusmb), "HELIFU")) # 2966 follow-up samples total
summary(str_detect(sample_names(heliusmb2), "HELIFU")) # 1829 follow-up samples with paired clinical data
fusamples <- sample_names(prune_samples(str_detect(sample_names(heliusmb2), "HELIFU"), heliusmb2))
basamples <- sample_names(prune_samples(str_detect(sample_names(heliusmb2), "HELIBA"), heliusmb2))
basamples <- str_remove(basamples, "HELIBA_")
fusamples <- str_remove(fusamples, "HELIFU_")
idspaired <- str_c("S",fusamples[which(fusamples %in% basamples)]) # 1825 with paired baseline-FU data
write.csv2(idspaired, 'data/16s/ids_16s_paired.csv')

## Select paired samples
idspaired2 <- c(str_replace(idspaired, "S", "HELIBA_"), str_replace(idspaired, "S", "HELIFU_"))
heliusmbpaired <- prune_samples(sample_names(heliusmb2) %in% idspaired2, heliusmb2)

## Add sample data to phyloseq
df_all <- df_long %>% filter(sampleID %in% sample_names(heliusmb))
rownames(df_all) <- df_all$sampleID
heliusmb2@sam_data <- sample_data(df_all) # for all subject with paired clin data
heliusmbpaired@sam_data <- sample_data(df_all %>% filter(sampleID %in% sample_names(heliusmbpaired))) # paired 16s

## Save
saveRDS(heliusmb2, file = "data/16s/phyloseq_withclinpaired.RDS")
saveRDS(heliusmbpaired, file = "data/16s/phyloseq_paired16s.RDS")


#### Shotgun ####
abundance <- rio::import('data/shotgun/combined_table_fixedlab.tsv') 
abundance <- abundance[,which(colnames(abundance) != "HELIBA_103370")] # only NAs
abundance2 <- abundance %>% filter(str_detect(clade_name, "s__") & !str_detect(clade_name, "t__"))
colnames(abundance2) <- c(colnames(abundance2)[1], str_replace(colnames(abundance2)[2:ncol(abundance2)], "_T1", ""))

# Short labeling bugs
clade <- abundance2$clade_name
cladesplit <- str_split(clade, "\\|", n = 8, simplify = TRUE)
cladesplit <- as.data.frame(cladesplit[,-8])
colnames(cladesplit) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
cladesplit <- cladesplit %>% mutate(across(everything(.), ~str_remove_all(.x, "[a-z]__")))
cladesplit$rowname <- clade
saveRDS(cladesplit, "data/shotgun/shotgun_taxtable.RDS")
rownames(abundance2) <- cladesplit$Species[match(cladesplit$rowname, abundance2$clade_name)]
abundance2$clade_name <- NULL
abundance2 <- t(as.matrix(abundance2))
write.csv2(abundance2, "data/shotgun/shotgun_abundance.csv", row.names = TRUE)
saveRDS(abundance2, "data/shotgun/shotgun_abundance.RDS")

# Clinical data filtered for available shotgun
rownames(abundance2)
dfshot <- df_new2 %>% 
    filter(sampleID_BA %in% rownames(abundance2) | sampleID_FU %in% rownames(abundance2)) %>% 
    droplevels(.)
table(dfshot$AgeDecade_BA, dfshot$Ethnicity, dfshot$Sex)
