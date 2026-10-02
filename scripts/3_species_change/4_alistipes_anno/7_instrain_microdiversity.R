## Alistipes putredinis within-strain microdiversity — baseline vs follow-up
## Follow-up to 5_instrain_strain_retention.R: that script asks whether the
## population at follow-up is still the SAME strain (popANI). This script
## asks, for strains that ARE retained, whether the strain itself diversified
## differently over time between ethnicities — a distinct question from
## retention, using inStrain profile's nucl_diversity (nucleotide diversity)
## per participant per timepoint rather than the compare-step popANI.
## Reads the per-timepoint genome_info.tsv files produced by the same
## Snellius job as 5_instrain_strain_retention.R:
##   scripts/0_run_workflows/2_run_shotgun_pipelines/strain_stability/
##     2_run_instrain_compare.sh   (writes profile/<subject>_<timepoint>_genome_info.tsv)
## and instrain_strain_retention.csv (same_strain calls) written by
## 5_instrain_strain_retention.R, which must be run first.
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
profile_dir    <- "data/shotgun/instrain_ap/profile"
retention_path <- "results/3_species_change/4_alistipes_anno/instrain_strain_retention.csv"
clin_file      <- "data/clinicaldata/clinicaldata_long.RDS"
out_dir        <- "results/3_species_change/4_alistipes_anno"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### Constants ####
# Same coverage/breadth thresholds used for the compare-step QC in
# 2_run_instrain_compare.sh (MIN_COV/MIN_BREADTH); applied here to each
# timepoint's own profile so a low-coverage sample doesn't inflate/deflate
# its nucl_diversity estimate.
MIN_COV     <- 5
MIN_BREADTH <- 0.5

#### Load per-timepoint profiles ####
files <- list.files(profile_dir, pattern = "_genome_info\\.tsv$", full.names = TRUE)
if (length(files) == 0)
    stop("No inStrain profile output found in ", profile_dir,
         " — copy instrain_ap/ back from Snellius first.")
cat("Profile files found:", length(files), "\n")

prof <- map_dfr(files, read_tsv, show_col_types = FALSE) %>%
    mutate(subject_id = as.character(subject_id))

stopifnot(!any(duplicated(prof[, c("subject_id", "timepoint")])))

wide <- prof %>%
    dplyr::select(subject_id, timepoint, nucl_diversity, coverage, breadth) %>%
    pivot_wider(names_from = timepoint,
                values_from = c(nucl_diversity, coverage, breadth))

cat("Participants with a profile at both timepoints:", sum(complete.cases(wide)), "of", nrow(wide), "\n")

#### QC: trust nucl_diversity only where both timepoints are well covered ####
wide <- wide %>%
    filter(!is.na(nucl_diversity_baseline), !is.na(nucl_diversity_followup)) %>%
    mutate(qc_pass = coverage_baseline >= MIN_COV & coverage_followup >= MIN_COV &
                      breadth_baseline  >= MIN_BREADTH & breadth_followup  >= MIN_BREADTH)
cat("Passing coverage >=", MIN_COV, "x and breadth >=", MIN_BREADTH,
    "at both timepoints:", sum(wide$qc_pass), "of", nrow(wide), "\n")

#### Join ethnicity/clade/covariates and the retention call ####
tip_meta <- readRDS("results/3_species_change/4_alistipes_anno/tip_meta_clades.RDS") %>%
    mutate(subject_id = as.character(subject_id)) %>%
    dplyr::select(subject_id, clade, EthnicityTot, Age, BMI, FUtime)

if (!file.exists(retention_path))
    stop("Run 5_instrain_strain_retention.R first — ", retention_path, " not found.")
retention <- read.csv(retention_path, colClasses = c(subject_id = "character")) %>%
    dplyr::select(subject_id, same_strain, enough_compared, eligible_completeness)

div <- wide %>%
    filter(qc_pass) %>%
    left_join(tip_meta, by = "subject_id") %>%
    left_join(retention, by = "subject_id") %>%
    # Align to the same n=120 population used throughout
    # 5_instrain_strain_retention.R (enough_compared: >=50% of the genome
    # compared between timepoints; eligible_completeness: MAG >70% complete),
    # rather than this script's own, slightly looser profile-based qc_pass
    # (coverage/breadth per sample independently). 3 participants pass
    # qc_pass but fail enough_compared specifically — their baseline and
    # follow-up reads individually cover enough of the genome, but the
    # OVERLAP between the two is <50%, making a baseline-vs-follow-up
    # comparison unreliable even though each timepoint looks fine on its own.
    filter(enough_compared, eligible_completeness) %>%
    mutate(
        delta_nucl_diversity = nucl_diversity_followup - nucl_diversity_baseline,
        log2fc_nucl_diversity = log2(nucl_diversity_followup / nucl_diversity_baseline)
    )

write.csv(div, file.path(out_dir, "instrain_microdiversity.csv"), row.names = FALSE)

#### Overview: retained vs replaced (all QC-passing participants) ####
cat("\nBy retention status (QC-passing participants with a popANI call):\n")
by_retention <- div %>%
    filter(!is.na(same_strain)) %>%
    group_by(same_strain) %>%
    summarise(n = n(),
              median_delta = median(delta_nucl_diversity),
              median_log2fc = median(log2fc_nucl_diversity), .groups = "drop")
print(by_retention)
cat("NOTE: for same_strain == FALSE the follow-up population is a different\n",
    "genotype (replacement), so its nucl_diversity change reflects a new\n",
    "colonist's diversity, not evolution of the baseline strain — the\n",
    "analysis below therefore focuses on same_strain == TRUE.\n", sep = "")

#### Main analysis: retained strains only ####
retained <- div %>% filter(same_strain, !is.na(EthnicityTot))
cat("\nRetained strains with known ethnicity:", nrow(retained), "\n")
print(count(retained, EthnicityTot))

cat("\nPaired baseline vs follow-up nucl_diversity, retained strains:\n")
wt_paired_all <- wilcox.test(retained$nucl_diversity_baseline, retained$nucl_diversity_followup, paired = TRUE)
cat("  All ethnicities: p =", signif(wt_paired_all$p.value, 3), "\n")
for (eth in levels(retained$EthnicityTot)) {
    sub <- retained %>% filter(EthnicityTot == eth)
    if (nrow(sub) >= 3) {
        wt <- wilcox.test(sub$nucl_diversity_baseline, sub$nucl_diversity_followup, paired = TRUE)
        cat("  ", eth, ": n =", nrow(sub), " p =", signif(wt$p.value, 3), "\n")
    }
}

cat("\nChange in nucl_diversity (follow-up - baseline) by ethnicity, retained strains only:\n")
by_eth_delta <- retained %>%
    group_by(EthnicityTot) %>%
    summarise(n = n(), median_delta = median(delta_nucl_diversity),
              median_log2fc = median(log2fc_nucl_diversity), .groups = "drop")
print(by_eth_delta)

if (nrow(by_eth_delta) == 2 && all(by_eth_delta$n >= 3)) {
    wt_delta <- wilcox.test(delta_nucl_diversity ~ EthnicityTot, data = retained, exact = FALSE)
    cat("Wilcoxon, change in nucl_diversity by ethnicity: p =", signif(wt_delta$p.value, 3), "\n")
}

write.csv(by_eth_delta, file.path(out_dir, "instrain_microdiversity_by_ethnicity.csv"), row.names = FALSE)

#### Plot ####
long_retained <- retained %>%
    dplyr::select(subject_id, EthnicityTot, nucl_diversity_baseline, nucl_diversity_followup) %>%
    pivot_longer(starts_with("nucl_diversity"), names_to = "timepoint", values_to = "nucl_diversity") %>%
    mutate(timepoint = factor(if_else(timepoint == "nucl_diversity_baseline", "Baseline", "Follow-up"),
                              levels = c("Baseline", "Follow-up")))

pl_trajectory <- ggplot(long_retained, aes(x = timepoint, y = nucl_diversity)) +
    geom_line(aes(group = subject_id, colour = EthnicityTot), alpha = 0.4) +
    geom_point(aes(fill = EthnicityTot), shape = 21, colour = "black", size = 2, alpha = 0.8) +
    stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.format",
                       comparisons = list(c("Baseline", "Follow-up"))) +
    facet_wrap(~EthnicityTot) +
    scale_y_log10() +
    scale_colour_manual(values = jco_palette(), guide = "none") +
    scale_fill_manual(values = jco_palette(), guide = "none") +
    labs(x = "", y = "Nucleotide diversity (log scale)",
         title = "A. putredinis microdiversity: retained strains") +
    theme_Publication()

pl_delta <- ggplot(retained, aes(x = EthnicityTot, y = delta_nucl_diversity, fill = EthnicityTot)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, shape = 21, colour = "black") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    scale_fill_manual(values = jco_palette(), guide = "none") +
    labs(x = "", y = "Change in nucleotide diversity (follow-up - baseline)",
         title = "B. Change by ethnicity") +
    theme_Publication()

(pl_micro <- ggarrange(pl_trajectory, pl_delta, ncol = 2, widths = c(1.3, 1), labels = c("A", "B")))
ggsave(file.path(out_dir, "instrain_microdiversity.pdf"), pl_micro, width = 11, height = 5)
cat("\nPlot saved to:", file.path(out_dir, "instrain_microdiversity.pdf"), "\n")

#### Cross-sectional: nucl_diversity by ethnicity at each timepoint ####
# Distinct from the paired baseline-vs-follow-up test above (within-person
# change) and from by_eth_delta (change restricted to retained strains) —
# this asks whether Dutch and SAS participants differ in how heterogeneous
# their A. putredinis population is AT a given timepoint, on all QC-passing
# participants (div), not just those with a retained strain.
div_eth <- div %>% filter(!is.na(EthnicityTot))

cat("\nCross-sectional nucl_diversity by ethnicity (all QC-passing participants):\n")
wt_baseline_eth <- wilcox.test(nucl_diversity_baseline ~ EthnicityTot, data = div_eth, exact = FALSE)
wt_followup_eth <- wilcox.test(nucl_diversity_followup ~ EthnicityTot, data = div_eth, exact = FALSE)
cat("  Baseline:  p =", signif(wt_baseline_eth$p.value, 3), "\n")
cat("  Follow-up: p =", signif(wt_followup_eth$p.value, 3), "\n")

cross_sectional_eth <- div_eth %>%
    group_by(EthnicityTot) %>%
    summarise(n = n(),
              median_baseline = median(nucl_diversity_baseline),
              median_followup = median(nucl_diversity_followup), .groups = "drop") %>%
    mutate(wilcox_p_baseline = wt_baseline_eth$p.value,
           wilcox_p_followup = wt_followup_eth$p.value)
print(cross_sectional_eth)
write.csv(cross_sectional_eth, file.path(out_dir, "instrain_microdiversity_crosssectional_ethnicity.csv"),
          row.names = FALSE)

long_div_eth <- div_eth %>%
    dplyr::select(subject_id, EthnicityTot, nucl_diversity_baseline, nucl_diversity_followup) %>%
    pivot_longer(starts_with("nucl_diversity"), names_to = "timepoint", values_to = "nucl_diversity") %>%
    mutate(timepoint = factor(if_else(timepoint == "nucl_diversity_baseline", "Baseline", "Follow-up"),
                              levels = c("Baseline", "Follow-up")))

pl_crosssectional <- ggplot(long_div_eth, aes(x = EthnicityTot, y = nucl_diversity, fill = EthnicityTot)) +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.6, shape = 21, colour = "black") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    facet_wrap(~timepoint) +
    scale_y_log10() +
    scale_fill_manual(values = jco_palette(), guide = "none") +
    labs(x = "", y = "Nucleotide diversity (log scale)",
         title = "A. putredinis microdiversity by ethnicity") +
    theme_Publication()

ggsave(file.path(out_dir, "instrain_microdiversity_crosssectional_ethnicity.pdf"),
       pl_crosssectional, width = 7, height = 5)
cat("Cross-sectional plot saved to:",
    file.path(out_dir, "instrain_microdiversity_crosssectional_ethnicity.pdf"), "\n")

#### Covariate check: does the ethnicity difference in diversity hold once HbA1c is accounted for? ####
# All QC-passing participants (div, not just retained strains) — this is
# about what predicts baseline diversity, not about change over time.
# Ethnicity and HbA1c are both checked because they're entangled here, not
# because HbA1c is assumed to be the explanation: HbA1c differs sharply by
# ethnicity in this subsample (median 36 Dutch vs 41.5 SAS), so a model with
# both terms can't cleanly attribute the diversity association to one or the
# other at this sample size. Age was checked separately and added nothing,
# so it's left out here for a simpler, more interpretable model.
hba1c_lookup <- readRDS(clin_file) %>%
    filter(timepoint == "baseline") %>%
    mutate(subject_id = as.character(str_extract(ID, "[0-9]+$"))) %>%
    distinct(subject_id, HbA1c)

div_hba1c <- div %>% left_join(hba1c_lookup, by = "subject_id")

cat("\nHbA1c by ethnicity (QC-passing participants):\n")
print(div_hba1c %>% filter(!is.na(HbA1c)) %>% group_by(EthnicityTot) %>%
        summarise(n = n(), median_hba1c = median(HbA1c), .groups = "drop"))
cat("Wilcoxon HbA1c ~ ethnicity: p =",
    signif(wilcox.test(HbA1c ~ EthnicityTot, data = div_hba1c, exact = FALSE)$p.value, 3), "\n")

mod_hba1c <- lm(nucl_diversity_baseline ~ EthnicityTot + HbA1c, data = div_hba1c)
cat("\nModel: nucl_diversity_baseline ~ Ethnicity + HbA1c\n")
print(summary(mod_hba1c)$coefficients)
cat("n used (complete cases):", nobs(mod_hba1c), "\n")

write.csv(
  as.data.frame(summary(mod_hba1c)$coefficients) %>% tibble::rownames_to_column("term"),
  file.path(out_dir, "microdiversity_hba1c_model.csv"), row.names = FALSE
)

pl_hba1c <- ggplot(div_hba1c %>% filter(!is.na(HbA1c)),
                    aes(x = HbA1c, y = nucl_diversity_baseline)) +
    geom_point(aes(fill = EthnicityTot), shape = 21, colour = "black", size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", colour = "grey30", se = TRUE) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    scale_y_log10() +
    scale_fill_manual(values = jco_palette(), name = "Ethnicity") +
    labs(x = "HbA1c (mmol/mol)", y = "Baseline nucleotide diversity (log scale)",
         title = "A. putredinis microdiversity vs HbA1c") +
    theme_Publication()

ggsave(file.path(out_dir, "microdiversity_vs_hba1c.pdf"), pl_hba1c, width = 6, height = 5)
cat("HbA1c correlation plot saved to:", file.path(out_dir, "microdiversity_vs_hba1c.pdf"), "\n")
