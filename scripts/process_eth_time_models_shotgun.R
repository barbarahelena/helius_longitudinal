## Process results machine learning with shotgun data

library(tidyverse)
library(ggpubr)
library(ggsci)
library(stringr)

options(scipen=999)

source("scripts/functions_shotgun.R")

#### Composition ####
## Plot assembled figure composition timepoints
path_true <- 'timepoint_shotgun/output_XGB_class_timepoints_2024_08_19__16-28-05'
data_path <- 'timepoint_shotgun/input_data'
labels <- c("Follow-up", "Baseline")

pl2 <- plot_feature_importance_shotgun(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/shotgun_comp_timepoint.pdf",
       width = 14, height = 18)

## Plot assembled figure composition timepoints - Dutch
path_true <- 'timepoint_dutch/output_XGB_class_dutch_timepoints_2024_08_19__18-02-25'
data_path <- 'timepoint_dutch/input_data'
labels <- c("Follow-up", "Baseline")

pl2 <- plot_feature_importance_shotgun(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/shotgun_comp_timepoint_dutch.pdf",
       width = 14, height = 18)

## Plot assembled figure composition timepoints - SAS
path_true <- 'timepoint_sas/output_XGB_class_sas_timepoints_2024_08_19__19-33-47'
data_path <- 'timepoint_sas/input_data'
labels <- c("Follow-up", "Baseline")

pl2 <- plot_feature_importance_shotgun(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/shotgun_comp_timepoint_sas.pdf",
       width = 14, height = 18)

## Plot assembled figure composition ethnicity - baseline
path_true <- 'eth_base/output_XGB_class_eth_baseline_2024_08_19__21-02-46'
data_path <- 'eth_base/input_data'
labels <- c("SAS", "Dutch")

pl2 <- plot_feature_importance_shotgun(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/shotgun_ethnicity_baseline.pdf",
       width = 14, height = 18)

## Plot assembled figure composition ethnicity - baseline
path_true <- 'eth_fu/output_XGB_class_eth_followup_2024_08_20__09-08-32'
data_path <- 'eth_fu/input_data'
labels <- c("SAS", "Dutch")

pl2 <- plot_feature_importance_shotgun(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/shotgun_ethnicity_followup.pdf",
       width = 14, height = 18)


#### Pathways ####
## Plot assembled figure pathway timepoints
path_true <- 'timepoint_pathways/output_XGB_class_timepoints_2024_08_20__13-22-30'
data_path <- 'timepoint_pathways/input_data'
labels <- c("Follow-up", "Baseline")

pl2 <- plot_feature_importance_pathways(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_pathways(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/pathways_timepoint.pdf",
       width = 14, height = 18)

## Plot assembled figure composition timepoints - Dutch
path_true <- 'timepoint_dutch_pathways/output_XGB_class_dutch_timepoints_2024_08_20__14-55-16'
data_path <- 'timepoint_dutch_pathways/input_data'
labels <- c("Follow-up", "Baseline")

pl2 <- plot_feature_importance_pathways(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_pathways(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/pathways_timepoint_dutch.pdf",
       width = 14, height = 18)

## Plot assembled figure composition timepoints - SAS
path_true <- 'timepoint_sas_pathways/output_XGB_class_sas_timepoints_2024_08_20__16-22-16'
data_path <- 'timepoint_sas_pathways/input_data'
labels <- c("Follow-up", "Baseline")

pl2 <- plot_feature_importance_pathways(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_pathways(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/pathways_timepoint_sas.pdf",
       width = 14, height = 18)

## Plot assembled figure pathways ethnicity - baseline
path_true <- 'eth_base_pathways/output_XGB_class_eth_baseline_2024_08_30__21-45-31'
data_path <- 'eth_base_pathways/input_data'
labels <- c("SAS", "Dutch")

pl2 <- plot_feature_importance_pathways(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_pathways(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/pathways_ethnicity_baseline.pdf",
       width = 14, height = 18)

## Plot assembled figure composition ethnicity - baseline
path_true <- 'eth_fu_pathways/output_XGB_class_eth_followup_2024_08_30__23-09-43'
data_path <- 'eth_fu_pathways/input_data'
labels <- c("SAS", "Dutch")

pl2 <- plot_feature_importance_pathways(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_pathways(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/pathways_ethnicity_followup.pdf",
       width = 14, height = 18)

