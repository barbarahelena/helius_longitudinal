## Alistipes putredinis strain retention (inStrain) — baseline vs follow-up
## Reads the per-participant popANI comparisons produced on Snellius by:
##   scripts/0_run_workflows/2_run_shotgun_pipelines/strain_stability/
##     1_make_instrain_manifest.R      (participant manifest, run locally)
##     2_run_instrain_compare.sh       (SLURM array, run on Snellius)
## (copy instrain_ap/ back from Snellius into data/shotgun/ first) and
## reports strain retention per ethnicity and per clade — the direct test of
## whether the clade shared at baseline and follow-up (necessarily identical,
## since each participant has one MAG; see 3_draw_tree.R / 4_qc_metadata_plots.R)
## reflects genuine strain persistence rather than an artefact of one MAG per
## participant.
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggpubr)
library(ggsci)

# Same Dutch/SAS colours used throughout 4_alistipes_anno/ (utils.R
# jco_palette()) — duplicated rather than sourced so this script stays
# independently reviewable; keep in sync with utils.R by hand.
jco_palette <- function() {
    cols <- pal_jco()(2)
    names(cols) <- c("Dutch", "South-Asian Surinamese")
    cols
}

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
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
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
}

#### Paths ####
in_dir       <- "data/shotgun/instrain_ap/compare"
# instrain_manifest.csv is produced by 1_make_instrain_manifest.R, which
# writes into the shared strain_stability results folder, not here.
manifest_dir <- "results/3_species_change/5_strain_stability"
out_dir      <- "results/3_species_change/4_alistipes_anno"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### Constants ####
# Conventional threshold for "same strain" (Olm et al. 2021, Science).
POPANI_SAME_STRAIN <- 0.99999
# Minimum fraction of the genome compared for the popANI call to be trusted.
MIN_GENOME_COMPARED <- 0.5
# Same eligibility rule as MIN_COMPLETENESS in utils.R (3_draw_tree.R,
# 2_vfdb_comparison.R): the tree-derived clade shown alongside these results
# only exists for bins >= 80% complete, so bins below that are excluded here
# too, for consistency, even though the manifest itself doesn't require it.
# Duplicated rather than sourced from utils.R so this script/PR stays
# independently reviewable — keep this value in sync with utils.R by hand.
MIN_COMPLETENESS <- 80

#### Load ####
files <- list.files(in_dir, pattern = "_genomeWide_compare\\.tsv$", full.names = TRUE)
if (length(files) == 0)
    stop("No inStrain compare output found in ", in_dir,
         " — copy instrain_ap/ back from Snellius first.")
cat("Comparisons found:", length(files), "\n")

cmp <- map_dfr(files, read_tsv, show_col_types = FALSE) %>%
    mutate(subject_id = as.character(subject_id))

# inStrain names this column percent_genome_compared in some versions
if (!"percent_genome_compared" %in% names(cmp) && "percent_compared" %in% names(cmp))
    cmp <- cmp %>% rename(percent_genome_compared = percent_compared)

manifest <- read.csv(file.path(manifest_dir, "instrain_manifest.csv"), colClasses = c(subject_id = "character"))

#### Join clade and ethnicity ####
tip_meta <- readRDS("results/3_species_change/4_alistipes_anno/tip_meta_clades.RDS") %>%
    mutate(subject_id = as.character(subject_id)) %>%
    dplyr::select(subject_id, clade, EthnicityTot)

res <- cmp %>%
    left_join(manifest, by = "subject_id") %>%
    left_join(tip_meta, by = "subject_id") %>%
    mutate(
        enough_compared    = percent_genome_compared >= MIN_GENOME_COMPARED,
        eligible_completeness = completeness >= MIN_COMPLETENESS,
        same_strain        = popANI >= POPANI_SAME_STRAIN
    )

cat("\nComparisons with >=", MIN_GENOME_COMPARED * 100, "% of the genome compared:",
    sum(res$enough_compared), "of", nrow(res), "\n")
cat("Comparisons with MAG completeness >=", MIN_COMPLETENESS, "% (same rule as",
    "3_draw_tree.R / 2_vfdb_comparison.R):", sum(res$eligible_completeness), "of", nrow(res), "\n")

valid <- res %>% filter(enough_compared, eligible_completeness)
cat("Valid comparisons meeting both criteria:", nrow(valid), "\n")

#### Strain retention ####
cat("\npopANI summary (valid comparisons):\n")
print(summary(valid$popANI))
cat("\nSame strain at both timepoints (popANI >=", POPANI_SAME_STRAIN, "):",
    sum(valid$same_strain), "of", nrow(valid),
    sprintf("(%.1f%%)\n", 100 * mean(valid$same_strain)))

cat("\nBy ethnicity:\n")
by_eth <- valid %>%
    filter(!is.na(EthnicityTot)) %>%
    group_by(EthnicityTot) %>%
    summarise(n = n(), n_same = sum(same_strain),
              pct_same = round(100 * mean(same_strain), 1),
              median_popANI = median(popANI), .groups = "drop")
print(by_eth)

if (nrow(by_eth) == 2 && all(by_eth$n > 0)) {
    ft <- fisher.test(matrix(c(by_eth$n_same, by_eth$n - by_eth$n_same), nrow = 2))
    cat("Fisher exact (same strain x ethnicity): OR =", round(ft$estimate, 2),
        " p =", signif(ft$p.value, 3), "\n")
    wt <- wilcox.test(popANI ~ EthnicityTot, data = valid %>% filter(!is.na(EthnicityTot)), exact = FALSE)
    cat("Wilcoxon popANI by ethnicity: p =", signif(wt$p.value, 3), "\n")
}

cat("\nBy clade:\n")
by_clade <- valid %>%
    filter(!is.na(clade)) %>%
    group_by(clade) %>%
    summarise(n = n(), n_same = sum(same_strain),
              pct_same = round(100 * mean(same_strain), 1),
              median_popANI = median(popANI), .groups = "drop")
print(by_clade)

write.csv(res, file.path(out_dir, "instrain_strain_retention.csv"), row.names = FALSE)
write.csv(by_eth, file.path(out_dir, "instrain_retention_by_ethnicity.csv"), row.names = FALSE)

#### Plot ####
# popANI is almost always squeezed into [0.999, 1], so a violin/linear scale
# is uninformative and boxplot+jitter+violin all overplot each other at the
# same handful of positions. Instead plot the genetic distance (1 - popANI)
# on a log10 scale, which spreads out "same strain" vs "diverged" comparisons,
# and drop the violin (no meaningful density shape with this few, spike-like
# values) so the boxplot and jitter don't visually duplicate.
dist_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)

pl_popani <- ggplot(valid %>% filter(!is.na(EthnicityTot)) %>%
                        mutate(popANI_dist = pmax(1 - popANI, 1e-6)),
                    aes(x = EthnicityTot, y = popANI_dist, fill = EthnicityTot)) +
    geom_hline(yintercept = 1 - POPANI_SAME_STRAIN, linetype = "dashed", colour = "firebrick") +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.6, shape = 21, colour = "black") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    scale_y_log10(breaks = dist_breaks,
                  labels = scales::label_number(accuracy = 0.000001)) +
    scale_fill_manual(values = jco_palette(), guide = "none") +
    labs(x = "", y = "Genetic distance (1 - popANI, log scale)",
         title = "A. putredinis strain retention") +
    theme_Publication()

pl_cov <- ggplot(res %>% mutate(popANI_dist = pmax(1 - popANI, 1e-6)),
                 aes(x = percent_genome_compared, y = popANI_dist,
                     colour = enough_compared)) +
    geom_hline(yintercept = 1 - POPANI_SAME_STRAIN, linetype = "dashed", colour = "firebrick") +
    geom_vline(xintercept = MIN_GENOME_COMPARED, linetype = "dashed", colour = "grey50") +
    geom_point(size = 1.8, alpha = 0.8) +
    scale_y_log10(breaks = dist_breaks,
                  labels = scales::label_number(accuracy = 0.000001)) +
    scale_colour_manual(values = c("TRUE" = "#1F78B4", "FALSE" = "grey65"),
                        name = "Included in comparison") +
    labs(x = "Fraction of genome compared", y = "Genetic distance (1 - popANI, log scale)",
         title = "Comparison quality") +
    theme_Publication()

(pl_strain <- ggarrange(pl_popani, pl_cov, ncol = 2, labels = c("A", "B")))
ggsave(file.path(out_dir, "instrain_strain_retention.pdf"), pl_strain, width = 11, height = 5)
cat("\nPlot saved to:", file.path(out_dir, "instrain_strain_retention.pdf"), "\n")
