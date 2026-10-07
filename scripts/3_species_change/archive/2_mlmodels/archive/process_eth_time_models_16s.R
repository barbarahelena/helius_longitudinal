## Process results machine learning with 16s data

library(tidyverse)
library(ggpubr)
library(ggsci)
library(stringr)
library(magick)
library(cowplot)

options(scipen=999)
dir.create("results/3_species_change/5_mlmodels", showWarnings = FALSE, recursive = TRUE)

source("scripts/functions.R")

## Plot assembled figure composition timepoints
path_true <- 'timepoint_16s/output_XGB_class_timepoints_16s_2024_08_19__16-27-25'
data_path <- 'timepoint_16s/input_data'
labels <- c("Follow-up", "Baseline")

pl2 <- plot_feature_importance_class(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/5_mlmodels/16s_comp_timepoint.pdf",
       width = 14, height = 18)

## Plot assembled figure composition timepoints - Dutch
path_true <- 'timepoint_dutch_16s/output_XGB_class_dutch_timepoints_16s_2024_08_19__17-52-06'
data_path <- 'timepoint_dutch_16s/input_data'
labels <- c("Follow-up", "Baseline")

pl2 <- plot_feature_importance_class(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/5_mlmodels/16s_comp_timepoint_dutch.pdf",
       width = 14, height = 18)

## Plot assembled figure composition timepoints - SAS
path_true <- 'timepoint_sas_16s/output_XGB_class_sas_timepoints_16s_2024_08_19__19-14-17'
data_path <- 'timepoint_sas_16s/input_data'
labels <- c("Follow-up", "Baseline")

pl2 <- plot_feature_importance_class(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/5_mlmodels/16s_comp_timepoint_sas.pdf",
       width = 14, height = 18)

## Plot assembled figure composition ethnicity - baseline
path_true <- 'eth_base_16s/output_XGB_class_eth_baseline_16s_2024_08_19__20-32-34'
data_path <- 'eth_base_16s/input_data'
labels <- c("SAS", "Dutch")

pl2 <- plot_feature_importance_class(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/5_mlmodels/16s_ethnicity_baseline.pdf",
       width = 14, height = 18)

## Plot assembled figure composition ethnicity - baseline
path_true <- 'eth_fu_16s/output_XGB_class_eth_followup_16s_2024_08_19__21-53-15'
data_path <- 'eth_fu_16s/input_data'
labels <- c("SAS", "Dutch")

pl2 <- plot_feature_importance_class(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/5_mlmodels/16s_ethnicity_followup.pdf",
       width = 14, height = 18)

