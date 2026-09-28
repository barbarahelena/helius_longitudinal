## Collect the inStrain baseline vs follow-up comparisons
## Reads the per-participant tables written by 2_run_instrain_compare.sh
## (copy instrain_ap/ back from Snellius into data/shotgun/ first) and
## reports strain retention per ethnicity and per clade.
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggpubr)

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
in_dir  <- "data/shotgun/instrain_ap/compare"
out_dir <- "results/3_species_change/5_strain_stability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### Constants ####
# Conventional threshold for "same strain" (Olm et al. 2021, Science).
POPANI_SAME_STRAIN <- 0.99999
# Minimum fraction of the genome compared for the popANI call to be trusted.
MIN_GENOME_COMPARED <- 0.5

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

manifest <- read.csv(file.path(out_dir, "instrain_manifest.csv"), colClasses = c(subject_id = "character"))

#### Join clade and ethnicity ####
tip_meta <- readRDS("results/3_species_change/4_alistipes_anno/tip_meta_clades.RDS") %>%
    mutate(subject_id = as.character(subject_id)) %>%
    dplyr::select(subject_id, clade, EthnicityTot)

res <- cmp %>%
    left_join(manifest, by = "subject_id") %>%
    left_join(tip_meta, by = "subject_id") %>%
    mutate(
        enough_compared = percent_genome_compared >= MIN_GENOME_COMPARED,
        same_strain     = popANI >= POPANI_SAME_STRAIN
    )

cat("\nComparisons with >=", MIN_GENOME_COMPARED * 100, "% of the genome compared:",
    sum(res$enough_compared), "of", nrow(res), "\n")

valid <- res %>% filter(enough_compared)

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
pl_popani <- ggplot(valid %>% filter(!is.na(EthnicityTot)),
                    aes(x = EthnicityTot, y = popANI, fill = EthnicityTot)) +
    geom_hline(yintercept = POPANI_SAME_STRAIN, linetype = "dashed", colour = "firebrick") +
    geom_violin(alpha = 0.6, colour = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.5) +
    scale_fill_manual(values = c("Dutch" = "#4E79A7",
                                 "South-Asian Surinamese" = "#F28E2B"), guide = "none") +
    labs(x = "", y = "popANI (baseline vs follow-up)",
         title = sprintf("A. putredinis strain retention (n = %d)", nrow(valid)),
         caption = sprintf("Dashed line: popANI = %s, conventional same-strain threshold",
                           POPANI_SAME_STRAIN)) +
    theme_Publication()

pl_cov <- ggplot(res, aes(x = percent_genome_compared, y = popANI,
                          colour = enough_compared)) +
    geom_hline(yintercept = POPANI_SAME_STRAIN, linetype = "dashed", colour = "firebrick") +
    geom_vline(xintercept = MIN_GENOME_COMPARED, linetype = "dashed", colour = "grey50") +
    geom_point(size = 1.8, alpha = 0.8) +
    scale_colour_manual(values = c("TRUE" = "#1F78B4", "FALSE" = "grey65"),
                        name = "Enough genome compared") +
    labs(x = "Fraction of genome compared", y = "popANI",
         title = "Comparison quality") +
    theme_Publication()

(pl_strain <- ggarrange(pl_popani, pl_cov, ncol = 2, labels = c("A", "B")))
ggsave(file.path(out_dir, "instrain_strain_retention.pdf"), pl_strain, width = 11, height = 5)
cat("\nPlot saved to:", file.path(out_dir, "instrain_strain_retention.pdf"), "\n")
