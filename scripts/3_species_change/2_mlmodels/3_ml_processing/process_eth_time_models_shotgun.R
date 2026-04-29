## Process results machine learning with shotgun data

library(tidyverse)
library(ggpubr)
library(ggsci)
library(stringr)
library(magick)
library(cowplot)

# options(scipen=999)
dir.create("results/3_species_change/2_mlmodels", showWarnings = FALSE, recursive = TRUE)

## ---- Helper theme ----
theme_Publication <- function(base_size = 12, base_family = "sans") {
    library(grid)
    library(ggthemes)
    suppressWarnings(
        theme_foundation(base_size = base_size, base_family = base_family) +
            theme(
                plot.title        = element_text(face = "bold", size = rel(1.0), hjust = 0.5,
                                                 family = "Helvetica"),
                text              = element_text(family = "Helvetica"),
                panel.background  = element_rect(colour = NA, fill = NA),
                plot.background   = element_rect(colour = NA, fill = NA),
                panel.border      = element_rect(colour = NA),
                axis.title        = element_text(face = "bold", size = rel(1)),
                axis.title.y      = element_text(angle = 90, vjust = 2),
                axis.title.x      = element_text(vjust = -0.2),
                axis.text         = element_text(),
                axis.line.x       = element_line(colour = "black"),
                axis.ticks.x      = element_line(),
                axis.ticks.y      = element_blank(),
                panel.grid.major  = element_line(colour = "#f0f0f0"),
                panel.grid.minor  = element_blank(),
                legend.key        = element_rect(colour = NA),
                legend.position   = "bottom",
                legend.key.size   = unit(0.2, "cm"),
                legend.spacing    = unit(0, "cm"),
                plot.margin       = unit(c(10, 5, 5, 5), "mm"),
                strip.background  = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
                strip.text        = element_text(face = "bold")
            )
    )
}

## ---- Feature importance bar chart ----
plot_feature_importance_shotgun <- function(path_true, top_n) {
    cols <- list(low = "#ECE7F2", mid = "#0570B0", high = "#034E7B")
    r <- rio::import(file.path(path_true, "feature_importance.txt"))
    r <- r %>% arrange(-RelFeatImp)
    r <- r[1:top_n, ]
    r <- r %>% mutate(FeatName = factor(make.unique(FeatName), levels = rev(make.unique(FeatName))))
    ggplot(data = r, aes(y = RelFeatImp, x = FeatName, fill = RelFeatImp)) +
        theme_Publication() +
        scale_fill_gradient2(low = cols$low, mid = cols$mid, high = cols$high,
                             space = "Lab", midpoint = 50, guide = "none") +
        geom_bar(stat = "identity") +
        coord_flip() +
        ylab("Relative Importance (%)") +
        xlab("") +
        theme(axis.text.x    = element_text(size = 10),
              axis.text.y    = element_text(size = 8),
              legend.key.size = unit(0.5, "cm"),
              legend.position = "right")
}

## ---- Top-N feature boxplots (faceted) ----
plot_features_top_shotgun <- function(input_path, output_path, top_n = 20, nrow = 4, labels) {
    plot_path <- file.path(output_path, "plots")
    dir.create(plot_path, showWarnings = FALSE)
    r <- rio::import(file.path(output_path, "feature_importance.txt"))
    r <- r %>% arrange(-RelFeatImp)
    input_data    <- rio::import(file.path(input_path, "X_data.txt"))
    feature_names <- read.csv(file.path(input_path, "feat_ids.txt"), sep = "\t", header = FALSE)
    names(input_data) <- feature_names$V1
    if (top_n > ncol(input_data)) {
        cat("\n\nRequested no. of features is higher than total number of features in model.\nShowing all features in model.\n\n")
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[!features_tk %in% c("random_variable1", "random_variable2")]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    y <- rio::import(file.path(input_path, "y_binary.txt"))
    dd$y <- factor(ifelse(y$V1 == 1, labels[1], labels[2]), levels = labels)
    df <- dd %>%
        pivot_longer(-y, names_to = "features", values_to = "values") %>%
        mutate(features = fct_inorder(as.factor(features)))
    colorguide <- pal_simpsons()(2)
    comps <- list(c(labels[1], labels[2]))
    ggplot(df, aes(x = y, y = values + 0.001)) +
        geom_boxplot(aes(fill = y)) +
        scale_fill_manual(values = colorguide, guide = "none") +
        theme_Publication() +
        theme(legend.position = "none",
              strip.text = element_text(face = "bold", size = rel(0.8))) +
        labs(x = "Group", y = "Relative abundance (%)") +
        ggpubr::stat_compare_means(comparisons = comps, paired = FALSE, size = rel(3.0)) +
        scale_y_log10() +
        facet_wrap(~features, nrow = nrow, scales = "free")
}

## Helper: pick the most recent XGBoost output folder
latest_output <- function(base_dir, pattern) {
    dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
    dirs <- dirs[grepl(pattern, basename(dirs))]
    dirs[order(basename(dirs), decreasing = TRUE)][1]
}

#### Ethnicity prediction — supplementary figures ####
## Baseline
path_true <- latest_output("results/3_species_change/2_mlmodels/eth_base", "output_XGB")
data_path <- "results/3_species_change/2_mlmodels/eth_base/input_data"
labels <- c("SAS", "Dutch")

pl2 <- plot_feature_importance_shotgun(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"),
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/2_mlmodels/shotgun_ethnicity_baseline.pdf",
       width = 14, height = 18)

## Follow-up
path_true <- latest_output("results/3_species_change/2_mlmodels/eth_fu", "output_XGB")
data_path <- "results/3_species_change/2_mlmodels/eth_fu/input_data"
labels <- c("SAS", "Dutch")

pl2 <- plot_feature_importance_shotgun(path_true, 20)
svg_grob <- ggdraw() + draw_image(image_read_pdf(file.path(path_true, "Plot_AUC.pdf"), density = 200))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(svg_grob, pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"),
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/3_species_change/2_mlmodels/shotgun_ethnicity_followup.pdf",
       width = 14, height = 18)


#### Figure 3A — Ethnicity prediction panel (baseline + follow-up) ####

path_eth_base <- latest_output("results/3_species_change/2_mlmodels/eth_base", "output_XGB")
path_eth_fu   <- latest_output("results/3_species_change/2_mlmodels/eth_fu", "output_XGB")

roc_base <- ggdraw() +
    draw_image(image_read_pdf(file.path(path_eth_base, "Plot_AUC.pdf"), density = 200),
               y = 0, height = 0.93) +
    draw_label("Baseline ethnic differences", x = 0.5, y = 0.97, size = 12, fontface = "bold")

roc_fu <- ggdraw() +
    draw_image(image_read_pdf(file.path(path_eth_fu, "Plot_AUC.pdf"), density = 200),
               y = 0, height = 0.93) +
    draw_label("Follow-up ethnic differences", x = 0.5, y = 0.97, size = 12, fontface = "bold")

pl_fig3_A <- roc_base
pl_fig3_B <- roc_fu
ggsave(ggarrange(pl_fig3_A, pl_fig3_B, nrow = 2, labels = c("A", "B")),
       filename = "results/3_species_change/2_mlmodels/fig3_eth_roc.pdf",
       width = 14, height = 12)

#### Supplemental figure — ML AUROC panels ####
suppl_ml <- ggarrange(
  pl_fig3_A, pl_fig3_B,
  ncol   = 2,
  labels = c("A", "B")
)
ggsave(suppl_ml,
       filename = "results/3_species_change/suppl_figure_ml_auroc.pdf",
       width = 10, height = 5, device = cairo_pdf)
cat("Supplemental ML AUROC figure saved to: results/3_species_change/suppl_figure_ml_auroc.pdf\n")

#### Overlap: top features baseline vs follow-up ####

feat_base <- rio::import(file.path(path_eth_base, "feature_importance.txt")) %>%
    arrange(-RelFeatImp) %>%
    filter(!FeatName %in% c("random_variable1", "random_variable2")) %>%
    slice(1:20) %>%
    pull(FeatName)

feat_fu <- rio::import(file.path(path_eth_fu, "feature_importance.txt")) %>%
    arrange(-RelFeatImp) %>%
    filter(!FeatName %in% c("random_variable1", "random_variable2")) %>%
    slice(1:20) %>%
    pull(FeatName)

overlap_feat   <- intersect(feat_base, feat_fu)
only_base_feat <- setdiff(feat_base, feat_fu)
only_fu_feat   <- setdiff(feat_fu, feat_base)

cat("Total features baseline model:   ", length(feat_base), "\n")
cat("Total features follow-up model:  ", length(feat_fu), "\n")
cat("Overlapping features:            ", length(overlap_feat), "\n")
cat("Baseline only:                   ", length(only_base_feat), "\n")
cat("Follow-up only:                  ", length(only_fu_feat), "\n")
cat("\nOverlapping features (ranked by mean importance across both):\n")

overlap_df <- bind_rows(
    rio::import(file.path(path_eth_base, "feature_importance.txt")) %>% mutate(timepoint = "baseline"),
    rio::import(file.path(path_eth_fu,   "feature_importance.txt")) %>% mutate(timepoint = "follow-up")
) %>%
    filter(FeatName %in% overlap_feat) %>%
    group_by(FeatName) %>%
    mutate(mean_imp = mean(RelFeatImp)) %>%
    ungroup() %>%
    arrange(-mean_imp)

print(overlap_df %>% dplyr::select(FeatName, timepoint, RelFeatImp, mean_imp))

write.csv2(overlap_df, "results/3_species_change/2_mlmodels/ml_overlap_baseline_followup.csv",
           row.names = FALSE)
