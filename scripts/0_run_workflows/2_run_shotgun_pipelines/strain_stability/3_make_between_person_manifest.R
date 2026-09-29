## Manifest for the inStrain between-person popANI comparison
## Positive control for the within-person strain-retention result (see
## 5_instrain_strain_retention.R in scripts/3_species_change/4_alistipes_anno/):
## if within-person (baseline vs follow-up) popANI is not distinguishable from
## between-person popANI among unrelated members of the same clade, "20% same
## strain" would not indicate real persistence — it would just be the
## background similarity of any two members of that clade. Comparing
## within-person distances against this between-person background settles that.
##
## Design: one reference ("anchor") MAG per named clade — the best-quality bin
## in that clade — against which every other participant in the same clade is
## compared, using their own baseline reads only. This "star" design (anchor vs
## everyone else) avoids the O(n^2) blow-up of comparing every pair within a
## clade (Clade I alone has 80 bins, i.e. 3160 possible pairs) while still
## giving a same-clade, cross-participant background distribution of a size
## comparable to the within-person result (n = 120).
##
## Writes one row per (clade, other participant). Copy the CSV to Snellius and
## run 4_run_instrain_between_person.sh as a job array over its rows.
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)

#### Output folder ####
out_dir <- "results/3_species_change/5_strain_stability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### Constants ####
NAMED_CLADES        <- paste0("Clade ", as.roman(1:5))
MAX_OTHER_PER_CLADE <- 40   # cap per clade, so a large clade (Clade I, n=80)
                             # doesn't dominate total runtime; random sample
ANCHOR_CONTAMINATION_MAX <- 5  # prefer a clean anchor; falls back if none <5%
SEED <- 20260929

batch_files <- sprintf("data/shotgun/alistipes_annotation/bins_alistipes_batch%d.csv", 1:3)
tree_results_dir <- "results/3_species_change/4_alistipes_anno"

#### Load clade assignment, quality, and batch per bin ####
tip_meta <- readRDS(file.path(tree_results_dir, "tip_meta_clades.RDS")) %>%
    dplyr::select(bin_name, subject_id, clade, EthnicityTot)

quality <- readRDS(file.path(tree_results_dir, "bin_quality_clade.RDS")) %>%
    dplyr::select(bin_name, Completeness, Contamination)

batch_of_bin <- map_dfr(seq_along(batch_files), function(i) {
    read.csv(batch_files[i], check.names = FALSE) %>%
        dplyr::select(bin) %>%
        mutate(bin_name = sub("\\.fa$", "", bin), batch = i) %>%
        dplyr::select(bin_name, batch)
})

bins <- tip_meta %>%
    inner_join(quality, by = "bin_name") %>%
    inner_join(batch_of_bin, by = "bin_name") %>%
    filter(clade %in% NAMED_CLADES) %>%
    mutate(subject_id = as.character(subject_id))

cat("Bins in named clades:", nrow(bins), "\n")
print(count(bins, clade))

#### Pick one anchor per clade: best completeness, preferring < 5% contamination ####
pick_anchor <- function(df) {
    clean <- df %>% filter(Contamination < ANCHOR_CONTAMINATION_MAX)
    pool  <- if (nrow(clean) > 0) clean else df
    pool %>% arrange(desc(Completeness), Contamination) %>% slice(1)
}

anchors <- bins %>%
    group_by(clade) %>%
    group_modify(~ pick_anchor(.x)) %>%
    ungroup() %>%
    transmute(clade, anchor_subject_id = subject_id, anchor_bin_name = bin_name,
              anchor_batch = batch, anchor_completeness = Completeness,
              anchor_contamination = Contamination,
              anchor_sample_baseline = str_c("HELIBA_", subject_id))

cat("\nAnchors (one per clade):\n")
print(anchors %>% dplyr::select(clade, anchor_subject_id, anchor_completeness, anchor_contamination))

#### Sample other participants per clade (excluding the anchor's own participant) ####
set.seed(SEED)
others <- bins %>%
    left_join(anchors %>% dplyr::select(clade, anchor_subject_id), by = "clade") %>%
    filter(subject_id != anchor_subject_id) %>%
    distinct(clade, subject_id) %>%
    group_by(clade) %>%
    group_modify(~ slice_sample(.x, n = min(nrow(.x), MAX_OTHER_PER_CLADE))) %>%
    ungroup() %>%
    transmute(clade, other_subject_id = subject_id,
              other_sample_baseline = str_c("HELIBA_", subject_id))

manifest <- others %>%
    left_join(anchors, by = "clade") %>%
    arrange(clade, other_subject_id) %>%
    dplyr::select(clade, anchor_subject_id, anchor_bin_name, anchor_batch,
                  anchor_sample_baseline, anchor_completeness, anchor_contamination,
                  other_subject_id, other_sample_baseline)

stopifnot(!any(manifest$anchor_subject_id == manifest$other_subject_id))

cat("\nBetween-person comparisons per clade:\n")
print(count(manifest, clade))
cat("\nTotal comparisons:", nrow(manifest), "\n")

write.csv(manifest, file.path(out_dir, "instrain_between_person_manifest.csv"), row.names = FALSE)
cat("\nManifest written to:", file.path(out_dir, "instrain_between_person_manifest.csv"), "\n")
cat("Run the array over 1-", nrow(manifest), ".\n", sep = "")
