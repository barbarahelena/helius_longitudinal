## Manifest for the inStrain baseline vs follow-up strain comparison
## Issue: MAG clade "stability" is not measurable from the co-assembled bins
## (one MAG per participant), so strain retention is tested directly with
## inStrain compare (popANI) between each participant's two samples.
##
## Writes one row per participant that has their own Alistipes putredinis MAG
## detected at both timepoints. Copy the CSV to Snellius and run
## 2_run_instrain_compare.sh as a job array over its rows.
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)

#### Output folder ####
out_dir <- "results/3_species_change/5_strain_stability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### Constants ####
# Minimum re-mapping depth at BOTH timepoints for a participant to be worth
# comparing. inStrain needs ~5x to call popANI reliably; below that the
# comparison is reported but flagged as low coverage at the collection step.
MIN_DEPTH_BOTH <- 5

batch_files <- sprintf("data/shotgun/alistipes_annotation/bins_alistipes_batch%d.csv", 1:3)
trans_file  <- "data/shotgun/alistipes_annotation/bin_translation_table.tsv"

#### Bin -> batch, and depth at the participant's own two samples ####
bins <- map_dfr(seq_along(batch_files), function(i) {
    read.csv(batch_files[i], check.names = FALSE) %>%
        dplyr::select(bin, Completeness, Contamination, starts_with("Depth ")) %>%
        mutate(bin_name = sub("\\.fa$", "", bin), batch = i) %>%
        dplyr::select(bin_name, batch, Completeness, Contamination, starts_with("Depth "))
})

trans <- read.delim(trans_file, stringsAsFactors = FALSE) %>%
    mutate(subject_id = as.character(subject_id)) %>%
    dplyr::select(bin_name, subject_id)

depth_own <- bins %>%
    dplyr::select(bin_name, starts_with("Depth ")) %>%
    pivot_longer(starts_with("Depth "), names_to = "sampleID", values_to = "depth") %>%
    mutate(sampleID = sub("^Depth ", "", sampleID),
           depth    = replace_na(depth, 0)) %>%
    inner_join(trans, by = "bin_name") %>%
    filter(sampleID == str_c("HELIBA_", subject_id) |
           sampleID == str_c("HELIFU_", subject_id)) %>%
    mutate(timepoint = if_else(str_starts(sampleID, "HELIBA"), "depth_baseline", "depth_followup")) %>%
    dplyr::select(bin_name, timepoint, depth) %>%
    pivot_wider(names_from = timepoint, values_from = depth, values_fill = 0)

manifest <- bins %>%
    dplyr::select(bin_name, batch, Completeness, Contamination) %>%
    inner_join(trans, by = "bin_name") %>%
    inner_join(depth_own, by = "bin_name") %>%
    filter(depth_baseline > 0, depth_followup > 0) %>%
    transmute(
        subject_id,
        bin_name,
        batch,
        sample_baseline = str_c("HELIBA_", subject_id),
        sample_followup = str_c("HELIFU_", subject_id),
        depth_baseline  = round(depth_baseline, 2),
        depth_followup  = round(depth_followup, 2),
        completeness    = Completeness,
        contamination   = Contamination,
        adequate_depth  = depth_baseline >= MIN_DEPTH_BOTH & depth_followup >= MIN_DEPTH_BOTH
    ) %>%
    arrange(desc(adequate_depth), subject_id)

stopifnot(!any(duplicated(manifest$subject_id)))

cat("Participants with their MAG detected at both timepoints:", nrow(manifest), "\n")
cat("  of which >=", MIN_DEPTH_BOTH, "x at both:", sum(manifest$adequate_depth), "\n")
cat("Bins per batch:\n"); print(table(manifest$batch))

write.csv(manifest, file.path(out_dir, "instrain_manifest.csv"), row.names = FALSE)
cat("\nManifest written to:", file.path(out_dir, "instrain_manifest.csv"), "\n")
cat("Rows are sorted so adequate-depth participants come first; run the array\n")
cat("over 1-", sum(manifest$adequate_depth), " to cover those only.\n", sep = "")
