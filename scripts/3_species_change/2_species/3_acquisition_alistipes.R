## Alistipes putredinis acquisition — read-level (MetaPhlAn), full cohort
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

suppressMessages(library(tidyverse))
source("scripts/3_species_change/3_alistipes_anno/utils.R")

#### Paths ####
abundance_file <- "data/shotgun/shotgun_abundance.RDS"
clin_file      <- "data/clinicaldata/clinicaldata_long.RDS"
results_dir    <- "results/3_species_change/3_alistipes_anno"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

#### 1. Load MetaPhlAn relative abundance for A. putredinis ####
mb <- readRDS(abundance_file)
stopifnot("Alistipes_putredinis" %in% colnames(mb))

ap_df <- tibble(sampleID = rownames(mb), rel_abund = as.numeric(mb[, "Alistipes_putredinis"]))
cat("Samples with an Alistipes_putredinis value:", nrow(ap_df), "\n")

clin <- readRDS(clin_file) %>%
  filter(EthnicityTot %in% c("Dutch", "South-Asian Surinamese")) %>%
  distinct(sampleID, ID, EthnicityTot) %>%
  mutate(EthnicityTot = droplevels(factor(EthnicityTot,
                                          levels = c("Dutch", "South-Asian Surinamese"))))

ap_clin <- ap_df %>%
  inner_join(clin, by = "sampleID") %>%
  mutate(timepoint = if_else(str_starts(sampleID, "HELIBA"), "baseline", "follow-up"))

cat("Samples with abundance value and Dutch/SAS ethnicity:", nrow(ap_clin), "\n")

#### 2. Restrict to participants profiled at BOTH timepoints ####
# A missing sample (not sequenced / not in the MetaPhlAn table at that
# timepoint) is not the same as a true biological absence — only pairs with
# both samples present can be classified as retained/lost/acquired/never.
wide <- ap_clin %>%
  dplyr::select(ID, EthnicityTot, timepoint, rel_abund) %>%
  pivot_wider(names_from = timepoint, values_from = rel_abund)

paired <- wide %>% filter(!is.na(baseline), !is.na(`follow-up`))
cat("Participants with both timepoints profiled:", nrow(paired), "of", nrow(wide), "\n")

#### 3. Detection threshold sensitivity ####
# MetaPhlAn already applies its own internal detection/calling before
# reporting a nonzero value, so > 0 is the primary, least arbitrary threshold.
# Higher cutoffs are checked purely as a sensitivity range.
DETECTION_THRESHOLDS <- c(0, 1e-5, 1e-4, 1e-3)

classify_at <- function(thresh) {
  paired %>%
    mutate(
      has_bl = baseline    > thresh,
      has_fu = `follow-up` > thresh,
      status = case_when(
        has_bl & has_fu  ~ "retained",
        has_bl & !has_fu ~ "lost",
        !has_bl & has_fu ~ "acquired",
        TRUE             ~ "never"
      ),
      threshold = thresh
    )
}

#### 4. Unconditional comparison (mirrors the MAG-based OR for comparison) ####
cat("\n=== Unconditional: acquired vs rest, full cohort, by ethnicity ===\n")
uncond_results <- map_dfr(DETECTION_THRESHOLDS, function(thr) {
  d   <- classify_at(thr)
  tab <- table(d$EthnicityTot, d$status == "acquired")
  ft  <- fisher.test(tab)
  tibble(threshold = thr, OR = unname(ft$estimate), p = ft$p.value)
})
print(uncond_results)

#### 5. Baseline prevalence by ethnicity (full cohort, not MAG-selected) ####
base_prev <- classify_at(0) %>%
  group_by(EthnicityTot) %>%
  summarise(n = n(), pct_baseline_positive = round(100 * mean(has_bl), 1), .groups = "drop")
cat("\nBaseline colonisation prevalence (full cohort, both-timepoint pairs):\n")
print(base_prev)
prev_fisher <- fisher.test(table(classify_at(0)$EthnicityTot, classify_at(0)$has_bl))
cat("Fisher exact, baseline prevalence x ethnicity: p =", signif(prev_fisher$p.value, 3), "\n")

#### 6. Conditioned on at-risk (baseline-negative) ####
cat("\n=== Conditioned on baseline-negative (at-risk population) ===\n")
atrisk_results <- map_dfr(DETECTION_THRESHOLDS, function(thr) {
  d <- classify_at(thr) %>% filter(!has_bl)
  ft <- fisher.test(table(d$EthnicityTot, d$has_fu))
  rates <- d %>% group_by(EthnicityTot) %>%
    summarise(n = n(), pct_acquired = round(100 * mean(has_fu), 1), .groups = "drop")
  tibble(
    threshold   = thr,
    n_dutch     = rates$n[rates$EthnicityTot == "Dutch"],
    n_sas       = rates$n[rates$EthnicityTot == "South-Asian Surinamese"],
    pct_dutch   = rates$pct_acquired[rates$EthnicityTot == "Dutch"],
    pct_sas     = rates$pct_acquired[rates$EthnicityTot == "South-Asian Surinamese"],
    OR          = unname(ft$estimate),
    p           = ft$p.value
  )
})
print(atrisk_results)

write.csv(uncond_results, file.path(results_dir, "acquisition_readlevel_unconditional.csv"), row.names = FALSE)
write.csv(atrisk_results, file.path(results_dir, "acquisition_readlevel_atrisk.csv"), row.names = FALSE)
write.csv(base_prev, file.path(results_dir, "baseline_prevalence_fullcohort.csv"), row.names = FALSE)

#### 7. Plot: acquisition rate among at-risk participants, by ethnicity ####
plot_df <- classify_at(0) %>%
  filter(!has_bl) %>%
  mutate(EthnicityTot = factor(EthnicityTot, levels = c("Dutch", "South-Asian Surinamese")))

# Normal-approximation binomial CI (no extra package dependency); fine for a
# supplementary sanity plot, not for a number quoted in the text — Dutch n=22
# is too small for the normal approximation to be taken literally there.
binom_ci <- function(x) {
  p  <- mean(x); n <- length(x)
  se <- sqrt(p * (1 - p) / n)
  data.frame(y = p, ymin = max(0, p - 1.96 * se), ymax = min(1, p + 1.96 * se))
}

p_atrisk <- ggplot(plot_df, aes(x = EthnicityTot, y = as.numeric(has_fu), fill = EthnicityTot)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = binom_ci, geom = "errorbar", width = 0.2) +
  scale_fill_manual(values = jco_palette(), guide = "none") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    title    = "Acquisition rate among at-risk (baseline-negative) participants",
    subtitle = sprintf("Fisher exact p = %.3f (vs. unconditional p = %.1e on the same full cohort)",
                       atrisk_results$p[atrisk_results$threshold == 0],
                       uncond_results$p[uncond_results$threshold == 0]),
    x = "", y = "% acquiring A. putredinis by follow-up"
  ) +
  theme_Publication()

ggsave(file.path(results_dir, "acquisition_atrisk_readlevel.pdf"), p_atrisk, width = 6, height = 5)

cat("\nDone. Results written to:", results_dir, "\n")
