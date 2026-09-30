## Alistipes putredinis strain retention (inStrain) — between-person background
## Positive control for 5_instrain_strain_retention.R: compares the within-person
## (baseline vs follow-up) popANI distribution against a between-person, same-clade
## background, to test whether "same strain" retention is distinguishable from the
## baseline similarity of unrelated members of the same clade.
## Reads the star-design comparisons produced on Snellius by:
##   scripts/0_run_workflows/2_run_shotgun_pipelines/strain_stability/
##     3_make_between_person_manifest.R  (anchor + other participants per clade, run locally)
##     4_run_instrain_between_person.sh  (SLURM array, run on Snellius)
## (copy instrain_ap_between/ back from Snellius into data/shotgun/ first).
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
in_dir  <- "data/shotgun/instrain_ap_between/compare"
# instrain_between_person_manifest.csv is produced by
# 3_make_between_person_manifest.R, which writes into the shared
# strain_stability results folder, not here.
manifest_dir <- "results/3_species_change/5_strain_stability"
out_dir      <- "results/3_species_change/4_alistipes_anno"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### Constants ####
# Same thresholds as 5_instrain_strain_retention.R — duplicated rather than
# sourced so this script/PR stays independently reviewable; keep in sync by hand.
POPANI_SAME_STRAIN  <- 0.99999  # conventional "same strain" threshold (Olm et al. 2021, Science)
MIN_GENOME_COMPARED <- 0.5      # minimum fraction of the genome compared for the popANI call to be trusted

#### Load between-person comparisons ####
files <- list.files(in_dir, pattern = "_genomeWide_compare\\.tsv$", full.names = TRUE)
if (length(files) == 0)
    stop("No between-person inStrain compare output found in ", in_dir,
         " — copy instrain_ap_between/ back from Snellius first.")
cat("Between-person comparisons found:", length(files), "\n")

btwn <- map_dfr(files, read_tsv, show_col_types = FALSE) %>%
    mutate(anchor_subject_id = as.character(anchor_subject_id),
           other_subject_id  = as.character(other_subject_id))

# inStrain names this column percent_genome_compared in some versions
if (!"percent_genome_compared" %in% names(btwn) && "percent_compared" %in% names(btwn))
    btwn <- btwn %>% rename(percent_genome_compared = percent_compared)

cat("\nCompleted comparisons by clade (of those planned in instrain_between_person_manifest.csv):\n")
print(count(btwn, clade))

manifest_between <- read.csv(file.path(manifest_dir, "instrain_between_person_manifest.csv"))
planned_by_clade <- count(manifest_between, clade, name = "n_planned")
completion <- planned_by_clade %>%
    left_join(count(btwn, clade, name = "n_completed"), by = "clade") %>%
    mutate(n_completed = replace_na(n_completed, 0L))
cat("\nCompletion vs plan:\n")
print(completion)
if (any(completion$n_completed == 0))
    cat("\nNOTE: ", paste(completion$clade[completion$n_completed == 0], collapse = ", "),
        " has/have zero completed between-person comparisons — excluded from the\n",
        "background below, not because of any filter, but because no output exists yet.\n", sep = "")

# All 5 anchors are >90% complete (see instrain_between_person_manifest.csv,
# anchor_completeness) — the MIN_COMPLETENESS>=80 rule from
# 5_instrain_strain_retention.R / utils.R is a no-op here by construction
# (anchors are chosen as the highest-completeness bin per clade), so it is
# not re-applied; this stopifnot makes that explicit.
stopifnot(all(unique(manifest_between[, c("clade", "anchor_completeness")])$anchor_completeness >= 80))

btwn <- btwn %>%
    mutate(
        enough_compared = percent_genome_compared >= MIN_GENOME_COMPARED,
        same_strain     = popANI >= POPANI_SAME_STRAIN
    )
cat("\nBetween-person comparisons with >=", MIN_GENOME_COMPARED * 100,
    "% of the genome compared:", sum(btwn$enough_compared), "of", nrow(btwn), "\n")

valid_between <- btwn %>% filter(enough_compared)

#### Load within-person comparisons (5_instrain_strain_retention.R output) ####
within_path <- file.path(out_dir, "instrain_strain_retention.csv")
if (!file.exists(within_path))
    stop("Run 5_instrain_strain_retention.R first — ", within_path, " not found.")

within_res <- read.csv(within_path, colClasses = c(subject_id = "character")) %>%
    filter(enough_compared, eligible_completeness)
cat("\nWithin-person valid comparisons (from 5_instrain_strain_retention.R):", nrow(within_res), "\n")

#### Combine and compare ####
combined <- bind_rows(
    within_res %>% transmute(clade, popANI, same_strain,
                              group = "Within-person (baseline vs follow-up)"),
    valid_between %>% transmute(clade, popANI, same_strain,
                                 group = "Between-person (same clade, baseline)")
) %>%
    mutate(group = factor(group, levels = c("Within-person (baseline vs follow-up)",
                                             "Between-person (same clade, baseline)")))

cat("\nSame strain (popANI >=", POPANI_SAME_STRAIN, ") by group:\n")
by_group <- combined %>%
    group_by(group) %>%
    summarise(n = n(), n_same = sum(same_strain),
              pct_same = round(100 * mean(same_strain), 1),
              median_popANI = median(popANI), .groups = "drop")
print(by_group)

cat("\nBy clade, within group:\n")
by_group_clade <- combined %>%
    filter(!is.na(clade)) %>%
    group_by(group, clade) %>%
    summarise(n = n(), n_same = sum(same_strain),
              pct_same = round(100 * mean(same_strain), 1),
              median_popANI = median(popANI), .groups = "drop")
print(by_group_clade, n = Inf)

wt <- wilcox.test(popANI ~ group, data = combined, exact = FALSE)
cat("\nWilcoxon popANI, within- vs between-person:\n")
cat("  W =", unname(wt$statistic), " p =", signif(wt$p.value, 3), "\n")

same_strain_tab <- table(combined$group, combined$same_strain)
ft <- fisher.test(same_strain_tab)
cat("\nFisher exact (same strain x group): OR =", round(ft$estimate, 2),
    " p =", signif(ft$p.value, 3), "\n")

write.csv(combined, file.path(out_dir, "instrain_within_vs_between.csv"), row.names = FALSE)
write.csv(by_group, file.path(out_dir, "instrain_within_vs_between_summary.csv"), row.names = FALSE)

#### Plot ####
# Same genetic-distance (1 - popANI), log-scale convention as
# 5_instrain_strain_retention.R, now split by within- vs between-person.
dist_breaks <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)

pl_within_vs_between <- ggplot(combined %>% mutate(popANI_dist = pmax(1 - popANI, 1e-6)),
                               aes(x = group, y = popANI_dist, fill = group)) +
    geom_hline(yintercept = 1 - POPANI_SAME_STRAIN, linetype = "dashed", colour = "firebrick") +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.6, shape = 21, colour = "black") +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       comparisons = list(levels(combined$group))) +
    scale_y_log10(breaks = dist_breaks,
                  labels = scales::label_number(accuracy = 0.000001)) +
    scale_fill_manual(values = c("Within-person (baseline vs follow-up)"  = "#1F78B4",
                                 "Between-person (same clade, baseline)" = "grey60"),
                      guide = "none") +
    scale_x_discrete(labels = scales::label_wrap(15)) +
    labs(x = "", y = "Genetic distance (1 - popANI, log scale)",
         title = "A. putredinis strain identity: within- vs between-person") +
    theme_Publication()

ggsave(file.path(out_dir, "instrain_within_vs_between.pdf"), pl_within_vs_between, width = 6, height = 6)
cat("\nPlot saved to:", file.path(out_dir, "instrain_within_vs_between.pdf"), "\n")
