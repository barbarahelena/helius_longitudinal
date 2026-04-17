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
path_true <- 'timepoint_dutch_16s/output_XGB_class_eth_baseline_2024_08_19__21-02-46'
data_path <- 'timepoint_dutch_16s/input_data'
labels <- c("Follow-up", "Baseline")

pl2 <- plot_feature_importance_class(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/5_mlmodels/16s_comp_timepoint.pdf",
       width = 14, height = 18)

## Plot assembled figure composition timepoints - SAS
path_true <- 'timepoint_sas_16s/output_XGB_class_eth_baseline_2024_08_19__21-02-46'
data_path <- 'timepoint_sas_16s/input_data'
labels <- c("Baseline", "Follow-up")

pl2 <- plot_feature_importance_class(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/5_mlmodels/16s_comp_timepoint.pdf",
       width = 14, height = 18)

## Plot assembled figure composition ethnicity - baseline
path_true <- 'eth_base_16s/output_XGB_class_eth_baseline_2024_08_19__21-02-46'
data_path <- 'eth_base_16s/input_data'
labels <- c("Dutch", "SAS")

pl2 <- plot_feature_importance_class(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/5_mlmodels/16s_comp_timepoint.pdf",
       width = 14, height = 18)

## Plot assembled figure composition ethnicity - baseline
path_true <- 'eth_fu_16s/output_XGB_class_eth_baseline_2024_08_19__21-02-46'
data_path <- 'eth_fu_16s/input_data'
labels <- c("Dutch", "SAS")

pl2 <- plot_feature_importance_class(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/5_mlmodels/16s_comp_timepoint.pdf",
       width = 14, height = 18)

