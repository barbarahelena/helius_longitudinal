## Alistipes putredinis strain retention (inStrain) — baseline vs follow-up
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggpubr)
library(ggsci)

# Color palette
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
# Same QC threshold used for the annotation-track bin set overall
# (filter_samplesheets_by_quality.py: Completeness > 70%, Contamination < 10%)
# rather than the stricter >=80% rule in utils.R (which exists only to match
# external tree membership for 3_draw_tree.R / 2_vfdb_comparison.R). inStrain
# doesn't need tree membership, so it uses the broader QC set, in line with
# the rest of the findings. Participants whose MAG is 70-80% complete will
# therefore have no `clade` (tree tips are still restricted to >=80%; see
# utils.R) but are otherwise included here.
# Duplicated rather than sourced so this script/PR stays independently
# reviewable — keep this value in sync with filter_samplesheets_by_quality.py.
MIN_COMPLETENESS <- 70

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
        eligible_completeness = completeness > MIN_COMPLETENESS,
        same_strain        = popANI >= POPANI_SAME_STRAIN
    )

cat("\nComparisons with >=", MIN_GENOME_COMPARED * 100, "% of the genome compared:",
    sum(res$enough_compared), "of", nrow(res), "\n")
cat("Comparisons with MAG completeness >", MIN_COMPLETENESS, "% (same rule as",
    "filter_samplesheets_by_quality.py):", sum(res$eligible_completeness), "of", nrow(res), "\n")

#### Divergence classification: distinguish likely replacement from in-situ evolution ####
# POPANI_SAME_STRAIN (0.99999) is the conventional short-interval threshold
# (Olm et al. 2021, Science) for calling two samples "the same strain" — it
# assumes a persisting lineage accumulates too few mutations between
# comparisons to cross it. That assumption is much weaker over this cohort's
# ~6-year follow-up: a persisting lineage can plausibly accumulate enough de
# novo mutations in 6 years to drop below 0.99999 without ever being replaced
# by an unrelated strain. popANI alone cannot distinguish "replaced" from
# "evolved in place", but the SCALE of divergence can: an unrelated genomic
# background introduces far more population-level differences than a few
# years of mutation in a single lineage would.
#
# SNPS_PER_MB_DRIFT_MAX is read off this cohort's own distribution of
# population_SNPs among popANI < POPANI_SAME_STRAIN comparisons, which form a
# gradient rather than two clean clusters: ~70% sit within a few dozen
# SNPs/Mb of the threshold (consistent with ordinary within-host evolution
# over 6 years), while a smaller tail shows hundreds to >10,000 SNPs/Mb (far
# more consistent with an unrelated, newly-introduced genome). This is a
# descriptive cutoff read from this dataset's own distribution, not a
# validated threshold from the literature — the resulting counts should be
# treated as approximate, not a precise replacement rate.
SNPS_PER_MB_DRIFT_MAX <- 100

res <- res %>%
    mutate(
        snps_per_mb = population_SNPs / (compared_bases_count / 1e6),
        divergence_class = case_when(
            same_strain                                         ~ "Stable",
            !same_strain & snps_per_mb <  SNPS_PER_MB_DRIFT_MAX  ~ "Modest divergence (consistent with in-situ evolution)",
            !same_strain & snps_per_mb >= SNPS_PER_MB_DRIFT_MAX  ~ "Substantial divergence (consistent with replacement)",
            TRUE ~ NA_character_
        )
    )

valid <- res %>% filter(enough_compared, eligible_completeness)
cat("Valid comparisons meeting both criteria:", nrow(valid), "\n")

cat("\nDivergence classification (valid comparisons):\n")
divergence_summary <- valid %>%
    count(divergence_class) %>%
    mutate(pct = round(100 * n / sum(n), 1))
print(divergence_summary)

cat("\nDivergence classification by ethnicity:\n")
divergence_by_eth <- valid %>%
    filter(!is.na(EthnicityTot)) %>%
    count(EthnicityTot, divergence_class) %>%
    group_by(EthnicityTot) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()
print(divergence_by_eth)

write.csv(divergence_summary, file.path(out_dir, "instrain_divergence_classification.csv"), row.names = FALSE)
write.csv(divergence_by_eth, file.path(out_dir, "instrain_divergence_classification_by_ethnicity.csv"), row.names = FALSE)

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

pl_divergence <- ggplot(
    valid %>% filter(!is.na(divergence_class)) %>%
      mutate(divergence_class = factor(divergence_class,
             levels = c("Stable",
                        "Modest divergence (consistent with in-situ evolution)",
                        "Substantial divergence (consistent with replacement)"))),
    aes(x = divergence_class, fill = divergence_class)
  ) +
    geom_bar() +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.3) +
    scale_fill_manual(values = c("Stable" = "#59A14F",
                                 "Modest divergence (consistent with in-situ evolution)" = "#F28E2B",
                                 "Substantial divergence (consistent with replacement)" = "#E15759"),
                      guide = "none") +
    scale_x_discrete(labels = scales::label_wrap(18)) +
    labs(x = "", y = "n comparisons",
         title = "Divergence classification",
         subtitle = sprintf("Threshold: %d SNPs/Mb among popANI < %g comparisons",
                            SNPS_PER_MB_DRIFT_MAX, POPANI_SAME_STRAIN)) +
    theme_Publication()

(pl_strain <- ggarrange(pl_popani, pl_cov, pl_divergence, ncol = 3, labels = c("A", "B", "C")))
ggsave(file.path(out_dir, "instrain_strain_retention.pdf"), pl_strain, width = 15, height = 5)
cat("\nPlot saved to:", file.path(out_dir, "instrain_strain_retention.pdf"), "\n")
