# HELIUS Longitudinal Microbiome Study

Analysis code for the longitudinal follow-up of the HELIUS cohort, examining gut microbiome composition and dynamics across six ethnic groups in Amsterdam.

## Study Overview

The [HELIUS study](https://www.heliusstudy.nl) (Healthy Life in an Urban Setting) is a large multi-ethnic cohort based at Amsterdam UMC. This repository contains analysis scripts for the longitudinal component, linking baseline and follow-up gut microbiome data (16S rRNA and whole-genome shotgun metagenomics) to cardiometabolic outcomes including incident type 2 diabetes and hypertension.

**Ethnic groups**: Dutch, South-Asian Surinamese, African Surinamese, Turkish, Moroccan, Ghanaian
**Microbiome data**: 16S rRNA amplicon sequencing + whole-genome shotgun metagenomics
**Functional annotations**: CAZymes, biosynthetic gene clusters (BGCs), antimicrobial resistance genes (ARGs), metabolic pathways

## Repository Structure

```
scripts/
├── 0_run_workflows/                # Pipeline execution scripts
│   ├── 1_run_16s_pipeline/         #   16S amplicon pipeline
│   ├── 2_run_shotgun_pipelines/    #   Shotgun metagenomics pipelines
│   └── 3_data_cleaning/            #   Master data cleaning and preprocessing
├── 1_longitudinal_change/          # Baseline vs follow-up microbiome changes
├── 2_cmb_microbiome/               # Comorbidity–microbiome associations
├── 3_species_change/               # Species-level longitudinal analysis
│   ├── 1_comparison_16s/           #   16S diversity comparisons
│   ├── 2_species/                  #   LMM species trajectories
│   ├── 3_strain_stability/         #   Strain sharing analysis
│   ├── 4_alistipes_anno/           #   Alistipes functional gene annotation
│   ├── 4_odoribacter_anno/         #   Odoribacter functional gene annotation
│   └── 5_mlmodels/                 #   XGBoost ethnicity prediction models
├── 4_functional_change/            # Functional genomics (CAZymes, BGCs)
└── 5_arg/                          # Antimicrobial resistance genes

data/
├── 16s/                            # 16S phyloseq objects (rarefied, paired, with metadata)
├── shotgun/                        # Shotgun taxonomic profiles, functional annotations
├── clinicaldata/                   # All clinical data: raw SPSS files, processed RDS objects,
│                                   #   alpha diversity, variable dictionaries
├── metabolomics/                   # Metabolite abundance data
├── gwas/                           # GWAS-related data
├── EGA/                            # European Genome Archive submission files
└── otherpeople/                    # Data shared by collaborators

results/
├── tables/                         # Table 1 and supplementary tables
├── 1_longitudinal_change/
├── 2_cmb_microbiome/
├── 3_species_change/
├── 4_functional_change/
└── 5_arg/
```

## Analysis Chapters

### 1. Longitudinal Change
Describes cohort characteristics and overall microbiome shifts between baseline and follow-up. Includes alpha diversity trajectories (Shannon, species richness, Faith's PD), ordination (PCA/NMDS), and multi-panel figure assembly.

### 2. Comorbidity & Microbiome
Associations between microbiome composition and cardiometabolic comorbidities, including incident type 2 diabetes and strain stability analyses.

### 3. Species-level Change
Linear mixed models (`lmer`) for species trajectories over time, stratified by ethnicity. Includes strain sharing analyses, functional annotation of key taxa (*Alistipes*, *Odoribacter*), and XGBoost machine learning models predicting ethnicity from microbiome features.

### 4. Functional Change
Longitudinal LMMs for CAZyme families (carbohydrate-active enzymes) and biosynthetic gene clusters (gutSMASH), with cross-sectional comparisons and beta diversity analyses.

### 5. Antimicrobial Resistance Genes
Descriptive QC, cross-sectional comparisons, and longitudinal LMMs for ARG classes and subclasses.

## Environment Setup

Dependencies are managed with [Pixi](https://pixi.sh). To reproduce the environment:

```bash
# Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# Install all dependencies
pixi install

# Run individual chapters
pixi run lc-all       # Chapter 1: longitudinal change
pixi run cmb-all      # Chapter 2: comorbidity–microbiome
pixi run sc-all       # Chapter 3: species change (includes ML)
pixi run fc-all       # Chapter 4: functional change
pixi run arg-all      # Chapter 5: antimicrobial resistance genes

# Run everything
pixi run all
```

Key R packages: `tidyverse`, `ggpubr`, `phyloseq`, `vegan`, `lme4`, `emmeans`, `ComplexHeatmap`, `tableone`, `mixOmics`, `circlize`
Key Python packages: `xgboost`, `scikit-learn`, `pandas`

## Statistical Methods

- **Alpha diversity**: Wilcoxon signed-rank test (paired baseline vs follow-up); Kruskal-Wallis across ethnicities
- **Beta diversity**: PERMANOVA (vegan); Bray-Curtis dissimilarity
- **Species trajectories**: Linear mixed models with ethnicity × timepoint interaction — `lmer(abundance ~ ethnicity * timepoint + (1|ID))`
- **Post-hoc tests**: `emmeans` with Tukey correction
- **Machine learning**: XGBoost with Wilcoxon pre-screening of features; cross-validation and permutation testing

## Data Availability

Raw sequencing data are deposited in the European Genome Archive (EGA). Clinical data are available on reasonable request via Amsterdam UMC, subject to data access agreements.

## Contact

Barbara Verhaar — Amsterdam UMC
