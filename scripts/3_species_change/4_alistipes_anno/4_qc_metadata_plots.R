## QC and metadata plots — Alistipes putredinis bins
## Reads pre-computed RDS files saved by 3_draw_tree.R
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

source("scripts/3_species_change/4_alistipes_anno/utils.R")

library(ggsci)
library(ggalluvial)   # install.packages("ggalluvial") if needed

#### Paths ####
results_dir <- "results/3_species_change/4_alistipes_anno"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

batch_files <- c(
  "data/shotgun/alistipes_annotation/bins_alistipes_batch1.csv",
  "data/shotgun/alistipes_annotation/bins_alistipes_batch2.csv",
  "data/shotgun/alistipes_annotation/bins_alistipes_batch3.csv"
)

#### Constants ####
PLOT_CLADES <- c("Clade I", "Clade II", "Clade III", "Clade IV")

# Hand-picked high-contrast palette — must match 3_draw_tree.R
clade_pal <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
               "#76B7B2", "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC")

#### Load pre-computed data ####
tip_meta_clades  <- readRDS(file.path(results_dir, "tip_meta_clades.RDS"))
bin_quality_clade <- readRDS(file.path(results_dir, "bin_quality_clade.RDS"))

# Reconstruct clade_colors from the clades present in the data
all_clades   <- sort(unique(tip_meta_clades$clade))
clade_colors <- setNames(clade_pal[seq_along(all_clades)], all_clades)

# Restrict analyses to clades with at least MIN_CLADE_N bins total
MIN_CLADE_N  <- 20
clade_levels <- tip_meta_clades %>%
  count(clade) %>%
  filter(n >= MIN_CLADE_N) %>%
  arrange(clade) %>%
  pull(clade)
cat("Clades included (n >=", MIN_CLADE_N, "):", paste(clade_levels, collapse = ", "), "\n")

clade_fill_cols <- clade_colors[clade_levels]
clade_pairs_all <- combn(clade_levels, 2, simplify = FALSE)

# Helper: pairwise Wilcoxon (BH) list for stat_compare_means
make_clade_pairs <- function(kw_res, var_col, df) {
  if (kw_res$p.value >= 0.05) return(list())
  map_dfr(clade_pairs_all, function(pair) {
    d <- df %>% filter(clade %in% pair)
    p <- wilcox.test(d[[var_col]] ~ d$clade, exact = FALSE)$p.value
    tibble(group1 = pair[1], group2 = pair[2], p_raw = p)
  }) %>%
    mutate(p_adj = p.adjust(p_raw, "BH")) %>%
    filter(p_adj < 0.05) %>%
    { map2(.$group1, .$group2, c) }
}

#### 1. Clinical variables per clade (Age, BMI, Sex) — baseline bins only ####
clin_clade <- tip_meta_clades %>%
  filter(timepoint == "baseline", clade %in% clade_levels) %>%  mutate(clade = factor(clade, levels = clade_levels))

aov_age <- aov(Age ~ clade, data = clin_clade)
aov_bmi <- aov(BMI ~ clade, data = clin_clade)
p_age_anova <- summary(aov_age)[[1]][["Pr(>F)"]][1]
p_bmi_anova <- summary(aov_bmi)[[1]][["Pr(>F)"]][1]
cat("ANOVA Age ~ clade: p =", signif(p_age_anova, 3), "\n")
cat("ANOVA BMI ~ clade: p =", signif(p_bmi_anova, 3), "\n")

# Tukey post-hoc — extract significant pairs for annotation
tukey_age_pairs <- if (p_age_anova < 0.05) {
  tk <- as.data.frame(TukeyHSD(aov_age)$clade)
  tk$pair <- rownames(tk)
  tk %>% filter(`p adj` < 0.05) %>%
    rowwise() %>%
    mutate(g = list(strsplit(pair, "-")[[1]])) %>%
    pull(g)
} else list()

tukey_bmi_pairs <- if (p_bmi_anova < 0.05) {
  tk <- as.data.frame(TukeyHSD(aov_bmi)$clade)
  tk$pair <- rownames(tk)
  tk %>% filter(`p adj` < 0.05) %>%
    rowwise() %>%
    mutate(g = list(strsplit(pair, "-")[[1]])) %>%
    pull(g)
} else list()

p_age_clade <- ggplot(clin_clade, aes(x = clade, y = Age, fill = clade)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6, alpha = 0.7) +
  scale_fill_manual(values = clade_fill_cols, guide = "none") +
  labs(title = "Age per clade",
       subtitle = sprintf("ANOVA p = %s", signif(p_age_anova, 3)),
       x = "", y = "Age (years)") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
if (length(tukey_age_pairs) > 0)
  p_age_clade <- p_age_clade +
    stat_compare_means(comparisons = tukey_age_pairs, method = "t.test",
                       label = "p.format", tip.length = 0)

p_bmi_clade <- ggplot(clin_clade, aes(x = clade, y = BMI, fill = clade)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6, alpha = 0.7) +
  scale_fill_manual(values = clade_fill_cols, guide = "none") +
  labs(title = "BMI per clade",
       subtitle = sprintf("ANOVA p = %s", signif(p_bmi_anova, 3)),
       x = "", y = "BMI (kg/m²)") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
if (length(tukey_bmi_pairs) > 0)
  p_bmi_clade <- p_bmi_clade +
    stat_compare_means(comparisons = tukey_bmi_pairs, method = "t.test",
                       label = "p.format", tip.length = 0)

sex_tab   <- table(clade = clin_clade$clade, sex = clin_clade$Sex)
chisq_sex <- chisq.test(sex_tab)
cat("Chi-square Sex ~ clade: p =", signif(chisq_sex$p.value, 3), "\n")

sex_long <- clin_clade %>%
  filter(!is.na(Sex)) %>%
  count(clade, Sex) %>%
  group_by(clade) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_sex_clade <- ggplot(sex_long, aes(x = clade, y = prop, fill = Sex)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Male" = "royalblue", "Female" = "firebrick1")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title    = "Sex distribution per clade",
       subtitle = paste0("Chi-square p = ", signif(chisq_sex$p.value, 3)),
       x = "", y = "Proportion of bins") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_clin_clade <- ggarrange(p_age_clade, p_bmi_clade, p_sex_clade,
                           ncol = 3, nrow = 1, labels = c("A", "B", "C"))

ggsave(
  file.path(results_dir, "clin_vars_per_clade.pdf"),
  p_clin_clade,
  width  = 12,
  height = 5
)
cat("Clinical variable plots saved to:", file.path(results_dir, "clin_vars_per_clade.pdf"), "\n")

#### 2. Bin completeness / contamination per ethnicity ####
# All clades included here — this is a QC check for assembly bias, not a
# biological analysis, so restricting to PLOT_CLADES would be misleading.
quality_eth <- bin_quality_clade %>%
  filter(!is.na(EthnicityTot)) %>%
  mutate(EthnicityTot = factor(EthnicityTot,
                               levels = c("Dutch", "South-Asian Surinamese")))

cat("\nBins with quality + ethnicity:", nrow(quality_eth), "\n")
cat("Completeness range:", round(range(quality_eth$Completeness, na.rm = TRUE), 1), "\n")

wx_comp <- wilcox.test(Completeness ~ EthnicityTot, data = quality_eth, exact = FALSE)
wx_cont <- wilcox.test(Contamination ~ EthnicityTot, data = quality_eth, exact = FALSE)
cat("Wilcoxon completeness p =", signif(wx_comp$p.value, 3), "\n")
cat("Wilcoxon contamination p =", signif(wx_cont$p.value, 3), "\n")

p_comp <- ggplot(quality_eth, aes(x = EthnicityTot, y = Completeness, fill = EthnicityTot)) +
  geom_boxplot(outlier.size = 0.8, width = 0.5) +
  scale_fill_manual(values = jco_palette(), guide = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x.npc = "center") +
  labs(x = "", y = "Completeness (%)") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_cont <- ggplot(quality_eth, aes(x = EthnicityTot, y = Contamination, fill = EthnicityTot)) +
  geom_boxplot(outlier.size = 0.8, width = 0.5) +
  scale_fill_manual(values = jco_palette(), guide = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x.npc = "center") +
  labs(x = "", y = "Contamination (%)") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_quality <- ggarrange(p_comp, p_cont, ncol = 2, labels = c("A", "B"))

ggsave(
  file.path(results_dir, "bin_quality_per_ethnicity.pdf"),
  p_quality,
  width  = 7,
  height = 5
)
cat("Bin quality per ethnicity saved to:",
    file.path(results_dir, "bin_quality_per_ethnicity.pdf"), "\n")

#### 3. Bin completeness / contamination per clade ####
quality_clade <- bin_quality_clade %>%
  filter(clade %in% clade_levels) %>%
  mutate(clade = factor(clade, levels = clade_levels))

kw_comp_clade <- kruskal.test(Completeness ~ clade, data = quality_clade)
kw_cont_clade <- kruskal.test(Contamination ~ clade, data = quality_clade)
cat("\nKruskal-Wallis Completeness ~ clade: p =", signif(kw_comp_clade$p.value, 3), "\n")
cat("Kruskal-Wallis Contamination ~ clade: p =", signif(kw_cont_clade$p.value, 3), "\n")

comp_pairs_clade <- make_clade_pairs(kw_comp_clade, "Completeness", quality_clade)
cont_pairs_clade <- make_clade_pairs(kw_cont_clade, "Contamination", quality_clade)

p_comp_clade <- ggplot(quality_clade, aes(x = clade, y = Completeness, fill = clade)) +
  geom_boxplot(outlier.size = 0.8, width = 0.5, alpha = 0.7) +
  scale_fill_manual(values = clade_colors, guide = "none") +
  labs(title = "Completeness per clade", x = "", y = "Completeness (%)") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
if (length(comp_pairs_clade) > 0)
  p_comp_clade <- p_comp_clade +
    stat_compare_means(comparisons = comp_pairs_clade, method = "wilcox.test",
                       label = "p.format", tip.length = 0)

p_cont_clade <- ggplot(quality_clade, aes(x = clade, y = Contamination, fill = clade)) +
  geom_boxplot(outlier.size = 0.8, width = 0.5, alpha = 0.7) +
  scale_fill_manual(values = clade_colors, guide = "none") +
  labs(title = "Contamination per clade", x = "", y = "Contamination (%)") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
if (length(cont_pairs_clade) > 0)
  p_cont_clade <- p_cont_clade +
    stat_compare_means(comparisons = cont_pairs_clade, method = "wilcox.test",
                       label = "p.format", tip.length = 0)

p_quality_clade <- ggarrange(p_comp_clade, p_cont_clade,
                              ncol = 2, labels = c("A", "B"))

ggsave(
  file.path(results_dir, "bin_quality_per_clade.pdf"),
  p_quality_clade,
  width  = 8,
  height = 5
)
cat("Bin quality per clade saved to:",
    file.path(results_dir, "bin_quality_per_clade.pdf"), "\n")

#### 4. Bin abundance per ethnicity × timepoint ####
abund_long <- tip_meta_clades %>%
  filter(!is.na(EthnicityTot), clade %in% clade_levels) %>%
  dplyr::select(bin_name, EthnicityTot, depth_baseline, depth_followup) %>%
  pivot_longer(cols = c(depth_baseline, depth_followup),
               names_to  = "timepoint",
               values_to = "depth") %>%
  mutate(
    timepoint    = recode(timepoint,
                          depth_baseline = "Baseline",
                          depth_followup = "Follow-up"),
    timepoint    = factor(timepoint, levels = c("Baseline", "Follow-up")),
    EthnicityTot = factor(EthnicityTot,
                          levels = c("Dutch", "South-Asian Surinamese"))
  )

wx_abund <- abund_long %>%
  group_by(timepoint) %>%
  summarise(
    p = wilcox.test(depth ~ EthnicityTot, exact = FALSE)$p.value,
    .groups = "drop"
  )
cat("\nWilcoxon abundance Dutch vs SAS:\n")
print(wx_abund)

p_abund <- ggplot(abund_long,
                  aes(x = EthnicityTot, y = depth, fill = EthnicityTot)) +
  geom_boxplot(outlier.size = 0.6, width = 0.55) +
  scale_fill_manual(values = jco_palette(), guide = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x.npc = "center") +
  facet_wrap(~ timepoint) +
  labs(x = "", y = "Sequencing depth") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(
  file.path(results_dir, "bin_abundance_per_ethnicity.pdf"),
  p_abund,
  width  = 7,
  height = 5
)
cat("Bin abundance plot saved to:",
    file.path(results_dir, "bin_abundance_per_ethnicity.pdf"), "\n")

#### 5. Clade composition by ethnicity × timepoint ####
# MAGs were co-assembled from pooled baseline + follow-up reads per participant.
# Strain switches would produce two distinct bins with divergent depth profiles
# (one baseline-dominant, one follow-up-dominant), potentially in different clades.
# The dominant clade at each timepoint is therefore the clade of the highest-depth
# bin at that timepoint — a valid proxy even from co-assembled data.
#
# -- 5a. Stacked bar: clade composition per ethnicity × timepoint --

# Ethnicity lookup per subject — derived from bins that have clinical data
# (follow-up bins often have EthnicityTot = NA because clin sampleIDs are HELIBA_*)
eth_lookup <- tip_meta_clades %>%
  filter(!is.na(EthnicityTot)) %>%
  mutate(subject_id = str_extract(sampleID, "\\d+$")) %>%
  distinct(subject_id, EthnicityTot)

subj_depth_all <- tip_meta_clades %>%
  filter(clade %in% clade_levels) %>%
  mutate(subject_id = str_extract(sampleID, "\\d+$")) %>%
  dplyr::select(-EthnicityTot) %>%
  left_join(eth_lookup, by = "subject_id")

# Dominant clade per subject at baseline = clade of the bin with max depth_baseline
clade_at_bl <- subj_depth_all %>%
  filter(depth_baseline > 0) %>%
  group_by(subject_id, EthnicityTot) %>%
  slice_max(depth_baseline, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(subject_id, EthnicityTot, timepoint = "baseline", clade)

# Dominant clade per subject at follow-up = clade of the bin with max depth_followup
clade_at_fu <- subj_depth_all %>%
  filter(depth_followup > 0) %>%
  group_by(subject_id, EthnicityTot) %>%
  slice_max(depth_followup, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(subject_id, EthnicityTot, timepoint = "follow-up", clade)

eth_tp_clade <- bind_rows(clade_at_bl, clade_at_fu) %>%
  filter(!is.na(EthnicityTot)) %>%
  mutate(
    EthnicityTot = factor(EthnicityTot, levels = c("Dutch", "South-Asian Surinamese")),
    timepoint    = factor(timepoint,    levels = c("baseline", "follow-up")),
    clade        = factor(clade,        levels = clade_levels)
  ) %>%
  count(EthnicityTot, timepoint, clade) %>%
  group_by(EthnicityTot, timepoint) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Chi-square per timepoint
for (tp in c("baseline", "follow-up")) {
  mat <- eth_tp_clade %>%
    filter(timepoint == tp) %>%
    dplyr::select(EthnicityTot, clade, n) %>%
    pivot_wider(names_from = EthnicityTot, values_from = n, values_fill = 0L) %>%
    tibble::column_to_rownames("clade") %>%
    as.matrix()
  chi <- chisq.test(mat)
  cat(sprintf("\nChi-square clade × ethnicity at %s: chi2 = %.2f, df = %d, p = %.4f\n",
              tp, chi$statistic, chi$parameter, chi$p.value))
}

p_eth_tp <- ggplot(eth_tp_clade,
                   aes(x = timepoint, y = prop, fill = clade)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = clade_fill_cols, name = "Clade") +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~ EthnicityTot, ncol = 2) +
  labs(
    title = "Clade composition per ethnicity and timepoint",
    x     = "",
    y     = "Proportion of bins"
  ) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(
  file.path(results_dir, "clade_by_ethnicity_timepoint.pdf"),
  p_eth_tp,
  width  = 8,
  height = 5
)
cat("Clade × ethnicity × timepoint plot saved to:",
    file.path(results_dir, "clade_by_ethnicity_timepoint.pdf"), "\n")

# -- 5b. Within-participant clade stability --
# MAGs were co-assembled from pooled baseline + follow-up reads per participant.
# This actually *enables* strain-switch detection: if a participant harboured a
# different Alistipes strain at follow-up, the assembler would produce two distinct
# bins — one recruiting mainly baseline reads (high depth_baseline, low
# depth_followup) and one recruiting mainly follow-up reads (the reverse).
# Such a pair would appear as a clade transition below.
#
# Conversely, a single bin with depth > 0 at both timepoints means the same
# strain was dominant at both — genuine longitudinal stability.
#
# Strategy: per participant, find the highest-depth bin at each timepoint
# (using depth_baseline / depth_followup). If those two bins are in different
# clades, that is a real clade transition.

# Diagnostic: check how many bins have each depth > 0
cat("\nBins with depth_baseline > 0:", sum(tip_meta_clades$depth_baseline > 0, na.rm = TRUE), "\n")
cat("Bins with depth_followup  > 0:", sum(tip_meta_clades$depth_followup  > 0, na.rm = TRUE), "\n")
cat("Bins with both depths     > 0:", sum(tip_meta_clades$depth_baseline > 0 &
                                           tip_meta_clades$depth_followup  > 0, na.rm = TRUE), "\n")

# Per-subject clade sets derived from depth columns — use ALL clades (no MIN_CLADE_N
# filter) so that participants whose dominant bin is in a small clade are not excluded
subj_depth_allclades <- tip_meta_clades %>%
  mutate(subject_id = str_extract(sampleID, "\\d+$")) %>%
  dplyr::select(-EthnicityTot) %>%
  left_join(eth_lookup, by = "subject_id")

# Clades present at baseline (depth_baseline > 0) per subject
clades_bl <- subj_depth_allclades %>%
  filter(depth_baseline > 0) %>%
  group_by(subject_id, EthnicityTot) %>%
  slice_max(depth_baseline, n = 1, with_ties = FALSE) %>%
  summarise(clades_bl = list(unique(clade)), .groups = "drop")

# Clades present at follow-up (depth_followup > 0) per subject
clades_fu <- subj_depth_allclades %>%
  filter(depth_followup > 0) %>%
  group_by(subject_id) %>%
  slice_max(depth_followup, n = 1, with_ties = FALSE) %>%
  summarise(clades_fu = list(unique(clade)), .groups = "drop")

paired_both <- clades_bl %>%
  inner_join(clades_fu, by = "subject_id") %>%
  mutate(
    # Use dominant clade (first in list, which = highest-depth bin picked earlier)
    clade_baseline = map_chr(clades_bl, 1),
    clade_followup = map_chr(clades_fu, 1),
    same_clade     = clade_baseline == clade_followup,
    EthnicityTot   = factor(EthnicityTot, levels = c("Dutch", "South-Asian Surinamese"))
  )

cat("\nParticipants with Alistipes detected at both timepoints:", nrow(paired_both), "\n")
cat("Clade stability overall:",
    sum(paired_both$same_clade), "of", nrow(paired_both),
    sprintf("(%.1f%%) — NOTE: excludes participants with baseline-only or follow-up-only detection\n",
            100 * mean(paired_both$same_clade)))

# Full detection picture: retained, lost, acquired
# n_baseline = detected at baseline (baseline-only + both)
# n_lost     = baseline-only (had Alistipes, lost it)
# n_retained = both timepoints (kept Alistipes)
# n_acquired = follow-up-only (gained Alistipes)
all_subj <- bind_rows(
  clades_bl %>% transmute(subject_id, EthnicityTot, has_bl = TRUE),
  clades_fu %>% transmute(subject_id, has_fu = TRUE)
) %>%
  group_by(subject_id) %>%
  summarise(
    has_bl       = any(!is.na(has_bl) & has_bl),
    has_fu       = any(!is.na(has_fu) & has_fu),
    EthnicityTot = first(na.omit(EthnicityTot)),
    .groups      = "drop"
  ) %>%
  mutate(
    detection = case_when(
      has_bl & has_fu  ~ "retained",
      has_bl & !has_fu ~ "lost",
      !has_bl & has_fu ~ "acquired"
    ),
    EthnicityTot = factor(EthnicityTot, levels = c("Dutch", "South-Asian Surinamese"))
  )

cat("\nAlistipes detection outcome overall:\n")
print(table(all_subj$detection))
cat("\nAlistipes detection outcome by ethnicity:\n")
print(table(all_subj$EthnicityTot, all_subj$detection))

cat("\nClade stability by ethnicity (detected at both timepoints only):\n")
print(
  paired_both %>%
    group_by(EthnicityTot) %>%
    summarise(n          = n(),
              n_stable   = sum(same_clade),
              pct_stable = round(100 * mean(same_clade), 1),
              .groups    = "drop")
)

# McNemar / exact test: is stability different between ethnicities?
stab_tab <- table(ethnicity = paired_both$EthnicityTot,
                  stable    = paired_both$same_clade)
cat("\nStability table (ethnicity × same_clade):\n")
print(stab_tab)
if (all(dim(stab_tab) == c(2, 2))) {
  ft <- fisher.test(stab_tab)
  cat("Fisher exact test p =", signif(ft$p.value, 3), "\n")
}

# Logistic regression on same_clade is not informative when nearly all
# participants are stable — instead, model detection pattern (baseline only,
# follow-up only, or both) as the outcome, which captures who gains/loses
# Alistipes over time and whether this differs by ethnicity and clade.

detect_pat <- subj_depth_all %>%
  group_by(subject_id, EthnicityTot, clade) %>%
  summarise(
    has_bl = any(depth_baseline > 0),
    has_fu = any(depth_followup  > 0),
    .groups = "drop"
  ) %>%
  filter(has_bl | has_fu) %>%
  mutate(
    pattern = case_when(
      has_bl & has_fu  ~ "both",
      has_bl & !has_fu ~ "baseline only",
      !has_bl & has_fu ~ "follow-up only"
    ),
    pattern      = factor(pattern,
                          levels = c("baseline only", "both", "follow-up only")),
    EthnicityTot = factor(EthnicityTot,
                          levels = c("Dutch", "South-Asian Surinamese")),
    clade        = factor(clade, levels = clade_levels)
  )

cat("\nDetection pattern overall:\n")
print(table(detect_pat$pattern))
cat("\nDetection pattern by ethnicity:\n")
print(table(detect_pat$EthnicityTot, detect_pat$pattern))

# Chi-square / Fisher: does detection pattern differ by ethnicity?
# Use Fisher's exact (simulated p-value) because several cells have n < 5
det_tab <- table(ethnicity = detect_pat$EthnicityTot,
                 pattern   = detect_pat$pattern)
ft_det <- fisher.test(det_tab, simulate.p.value = TRUE, B = 10000)
cat(sprintf("\nFisher exact detection pattern × ethnicity: p = %.4f\n",
            ft_det$p.value))

# More targeted: is "follow-up only" (newly acquired) enriched in SAS?
# Binary logistic regression: fu_only ~ EthnicityTot
# (each participant contributes one row, outcome = gained Alistipes at FU)
detect_subj <- detect_pat %>%
  filter(!is.na(EthnicityTot)) %>%
  # one row per participant: take their dominant clade pattern
  group_by(subject_id, EthnicityTot) %>%
  # if any clade is "follow-up only" and none are "both", classify as gained
  summarise(
    fu_only  = any(pattern == "follow-up only") & !any(pattern == "both"),
    bl_only  = any(pattern == "baseline only")  & !any(pattern == "both"),
    .groups  = "drop"
  ) %>%
  left_join(
    tip_meta_clades %>%
      filter(timepoint == "baseline") %>%
      mutate(subject_id = str_extract(sampleID, "\\d+$")) %>%
      distinct(subject_id, Age, BMI, Sex),
    by = "subject_id"
  )

# Simple 2x2 Fisher exact: is "follow-up only" enriched in SAS vs Dutch?
fu_tab <- table(ethnicity = detect_subj$EthnicityTot,
                fu_only   = detect_subj$fu_only)
cat("\nFollow-up only detection by ethnicity:\n")
print(fu_tab)
ft_fu <- fisher.test(fu_tab)
cat(sprintf("Fisher exact test (fu_only × ethnicity): OR = %.2f, p = %.4f\n",
            ft_fu$estimate, ft_fu$p.value))

# Which clade(s) are being newly acquired at follow-up?
cat("\nClade of newly acquired (follow-up only) bins, by ethnicity:\n")
fu_clade_tab <- detect_pat %>%
  filter(pattern == "follow-up only", !is.na(EthnicityTot)) %>%
  count(EthnicityTot, clade) %>%
  arrange(EthnicityTot, desc(n))
print(fu_clade_tab)

p_fu_clade <- ggplot(fu_clade_tab, aes(x = clade, y = n, fill = EthnicityTot)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = jco_palette(), name = "Ethnicity") +
  labs(
    title = "Clade of newly acquired Alistipes at follow-up",
    x     = "",
    y     = "Number of participants"
  ) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(results_dir, "fu_only_clade.pdf"),
  p_fu_clade,
  width  = 7,
  height = 5
)
cat("Follow-up acquisition by clade plot saved to:",
    file.path(results_dir, "fu_only_clade.pdf"), "\n")

# Fisher exact: does the acquired clade differ between Dutch and SAS?
fu_clade_mat <- fu_clade_tab %>%
  pivot_wider(names_from = EthnicityTot, values_from = n, values_fill = 0L) %>%
  tibble::column_to_rownames("clade") %>%
  as.matrix()
if (ncol(fu_clade_mat) == 2 && nrow(fu_clade_mat) > 1) {
  ft_fu_clade <- fisher.test(fu_clade_mat, simulate.p.value = TRUE, B = 10000)
  cat(sprintf("Fisher exact acquired clade × ethnicity: p = %.4f\n", ft_fu_clade$p.value))
}

# Fisher exact: does detection pattern differ by clade?
det_clade_tab <- table(clade   = detect_pat$clade,
                       pattern = detect_pat$pattern)
ft_det_clade <- fisher.test(det_clade_tab, simulate.p.value = TRUE, B = 10000)
cat(sprintf("Fisher exact detection pattern × clade: p = %.4f\n",
            ft_det_clade$p.value))

# Plot: detection pattern per clade, faceted by ethnicity
detect_sum <- detect_pat %>%
  filter(!is.na(EthnicityTot)) %>%
  count(EthnicityTot, clade, pattern) %>%
  group_by(EthnicityTot, clade) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_detect <- ggplot(detect_sum, aes(x = clade, y = prop, fill = pattern)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c("baseline only"   = "#4E79A7",
               "both"            = "#59A14F",
               "follow-up only"  = "#F28E2B"),
    name   = "Detection"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~ EthnicityTot, ncol = 2) +
  labs(
    title = "Alistipes detection pattern per clade and ethnicity",
    x     = "",
    y     = "Proportion of participants"
  ) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(results_dir, "detection_pattern_per_clade.pdf"),
  p_detect,
  width  = 10,
  height = 5
)
cat("Detection pattern plot saved to:",
    file.path(results_dir, "detection_pattern_per_clade.pdf"), "\n")

# Sankey-style alluvial plot of baseline → follow-up clade per participant
alluvial_df <- paired_both %>%
  count(EthnicityTot, clade_baseline, clade_followup) %>%
  mutate(
    clade_baseline = factor(clade_baseline, levels = clade_levels),
    clade_followup = factor(clade_followup, levels = clade_levels)
  )

p_alluvial <- ggplot(alluvial_df,
                     aes(axis1 = clade_baseline, axis2 = clade_followup,
                         y = n)) +
  ggalluvial::geom_alluvium(aes(fill = clade_baseline), width = 1/4, alpha = 0.7) +
  ggalluvial::geom_stratum(width = 1/4, fill = "grey90", colour = "grey40") +
  geom_text(stat = ggalluvial::StatStratum, aes(label = after_stat(stratum)),
            size = 3.5) +
  scale_fill_manual(values = clade_fill_cols, name = "Baseline clade") +
  scale_x_discrete(limits = c("Baseline", "Follow-up"), expand = c(0.05, 0.05)) +
  facet_wrap(~ EthnicityTot, ncol = 2) +
  labs(
    title = "Within-participant clade transitions (baseline → follow-up)",
    y     = "Number of participants"
  ) +
  theme_Publication()

ggsave(
  file.path(results_dir, "clade_stability_alluvial.pdf"),
  p_alluvial,
  width  = 10,
  height = 6
)
cat("Clade stability alluvial plot saved to:",
    file.path(results_dir, "clade_stability_alluvial.pdf"), "\n")

# Heatmap of transition counts (baseline clade × follow-up clade), per ethnicity
transition_heat <- paired_both %>%
  count(EthnicityTot, clade_baseline, clade_followup) %>%
  mutate(
    clade_baseline = factor(clade_baseline, levels = clade_levels),
    clade_followup = factor(clade_followup, levels = rev(clade_levels))
  )

p_trans_heat <- ggplot(transition_heat,
                       aes(x = clade_baseline, y = clade_followup, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = n), size = 4) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "n participants") +
  facet_wrap(~ EthnicityTot, ncol = 2) +
  labs(
    title = "Clade transition matrix (baseline → follow-up)",
    x     = "Baseline clade",
    y     = "Follow-up clade"
  ) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(results_dir, "clade_transition_heatmap.pdf"),
  p_trans_heat,
  width  = 10,
  height = 5
)
cat("Clade transition heatmap saved to:",
    file.path(results_dir, "clade_transition_heatmap.pdf"), "\n")

# Summary table
write.csv(paired_both %>% dplyr::select(-clades_bl, -clades_fu),
          file.path(results_dir, "clade_stability_per_participant.csv"),
          row.names = FALSE)
cat("Per-participant clade stability table saved to:", results_dir, "\n")

#### 6. Dominant timepoint per clade ####
# Each bin is assigned a dominant timepoint = the sample with the highest depth.
# This shows whether bins in each clade tend to be baseline- or follow-up-dominant.
# Also stratified by ethnicity to test whether the clade × timepoint pattern
# differs between Dutch and South-Asian Surinamese participants.

tp_clade <- tip_meta_clades %>%
  filter(clade %in% clade_levels) %>%
  mutate(
    clade    = factor(clade,    levels = clade_levels),
    timepoint = factor(timepoint, levels = c("baseline", "follow-up"))
  ) %>%
  count(clade, timepoint) %>%
  group_by(clade) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Chi-square: is timepoint distribution independent of clade?
tp_mat <- tp_clade %>%
  dplyr::select(clade, timepoint, n) %>%
  pivot_wider(names_from = timepoint, values_from = n, values_fill = 0L) %>%
  tibble::column_to_rownames("clade") %>%
  as.matrix()
chi_tp <- chisq.test(tp_mat)
cat(sprintf("\nChi-square dominant timepoint × clade: chi2 = %.2f, df = %d, p = %.4f\n",
            chi_tp$statistic, chi_tp$parameter, chi_tp$p.value))

# Stratified by ethnicity: does the timepoint × clade pattern hold within each group?
tp_clade_eth <- tip_meta_clades %>%
  filter(clade %in% clade_levels, !is.na(EthnicityTot)) %>%
  mutate(
    clade     = factor(clade,     levels = clade_levels),
    timepoint = factor(timepoint, levels = c("baseline", "follow-up"))
  )

for (eth in c("Dutch", "South-Asian Surinamese")) {
  mat_eth <- tp_clade_eth %>%
    filter(EthnicityTot == eth) %>%
    count(clade, timepoint) %>%
    pivot_wider(names_from = timepoint, values_from = n, values_fill = 0L) %>%
    tibble::column_to_rownames("clade") %>%
    as.matrix()
  chi_eth <- chisq.test(mat_eth)
  cat(sprintf("  %s — chi2 = %.2f, df = %d, p = %.4f\n",
              eth, chi_eth$statistic, chi_eth$parameter, chi_eth$p.value))
}

# Plot stratified by ethnicity
tp_clade_eth_sum <- tp_clade_eth %>%
  count(EthnicityTot, clade, timepoint) %>%
  group_by(EthnicityTot, clade) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(EthnicityTot = factor(EthnicityTot,
                               levels = c("Dutch", "South-Asian Surinamese")))

p_tp_clade_eth <- ggplot(tp_clade_eth_sum,
                          aes(x = clade, y = prop, fill = timepoint)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("baseline"  = "#4E79A7",
                               "follow-up" = "#F28E2B"),
                    name   = "Dominant timepoint") +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~ EthnicityTot, ncol = 2) +
  labs(
    title = "Dominant timepoint per clade, stratified by ethnicity",
    x     = "",
    y     = "Proportion of bins"
  ) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(results_dir, "dominant_timepoint_per_clade_by_ethnicity.pdf"),
  p_tp_clade_eth,
  width  = 10,
  height = 5
)
cat("Stratified dominant timepoint plot saved to:",
    file.path(results_dir, "dominant_timepoint_per_clade_by_ethnicity.pdf"), "\n")

p_tp_clade <- ggplot(tp_clade, aes(x = clade, y = prop, fill = timepoint)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("baseline"   = "#4E79A7",
                               "follow-up"  = "#F28E2B"),
                    name   = "Dominant timepoint") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = "Dominant timepoint per clade",
    subtitle = sprintf("Chi-square p = %.4f", chi_tp$p.value),
    x        = "",
    y        = "Proportion of bins"
  ) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(results_dir, "dominant_timepoint_per_clade.pdf"),
  p_tp_clade,
  width  = 7,
  height = 5
)
cat("Dominant timepoint per clade plot saved to:",
    file.path(results_dir, "dominant_timepoint_per_clade.pdf"), "\n")

#### 7. Clade abundance (depth) at baseline vs follow-up — LMM ####
# Per participant × clade: model depth ~ time_elapsed * EthnicityTot + (1|subject_id)
# time_elapsed = 0 at baseline, = FUtime (years) at follow-up.
# This naturally handles variable follow-up intervals: the slope represents change
# per year, and the interaction tests whether that rate differs between SAS and Dutch.

library(lme4)
library(lmerTest)   # provides p-values via Satterthwaite df

# FUtime per subject — take from any row (not just baseline) to maximise coverage
futime_lookup <- tip_meta_clades %>%
  mutate(subject_id = str_extract(sampleID, "\\d+$")) %>%
  filter(!is.na(FUtime)) %>%
  distinct(subject_id, FUtime)

# Aggregate to one depth value per subject × clade × timepoint (sum across bins),
# then reshape to long format; set time_elapsed = 0 (baseline) or FUtime (follow-up)
depth_long <- subj_depth_all %>%
  filter(!is.na(EthnicityTot)) %>%
  dplyr::select(subject_id, EthnicityTot, clade, depth_baseline, depth_followup) %>%
  group_by(subject_id, EthnicityTot, clade) %>%
  summarise(depth_baseline = sum(depth_baseline, na.rm = TRUE),
            depth_followup = sum(depth_followup, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_longer(cols = c(depth_baseline, depth_followup),
               names_to  = "timepoint",
               values_to = "depth") %>%
  mutate(
    timepoint    = recode(timepoint,
                          depth_baseline = "Baseline",
                          depth_followup = "Follow-up"),
    timepoint    = factor(timepoint, levels = c("Baseline", "Follow-up")),
    EthnicityTot = factor(EthnicityTot, levels = c("Dutch", "South-Asian Surinamese")),
    clade        = factor(clade, levels = clade_levels)
  ) %>%
  left_join(futime_lookup, by = "subject_id") %>%
  mutate(time_elapsed = if_else(timepoint == "Baseline", 0, FUtime))

cat("FUtime available for",
    sum(!is.na(depth_long$FUtime[depth_long$timepoint == "Baseline"])),
    "of", sum(depth_long$timepoint == "Baseline"), "baseline rows\n")

# LMM per clade: depth ~ time_elapsed * EthnicityTot + (1|subject_id)
cat("\n=== LMM: depth ~ time_elapsed * EthnicityTot + (1|subject_id) per clade ===\n")

lmm_results <- map_dfr(clade_levels, function(cl) {
  # Keep participants with depth > 0 at either timepoint for this clade —
  # includes both stable carriers, those who lost it, and those who acquired it
  carriers <- depth_long %>%
    filter(clade == cl) %>%
    group_by(subject_id) %>%
    filter(any(depth > 0)) %>%
    pull(subject_id) %>%
    unique()
  df_cl <- depth_long %>%
    filter(clade == cl, subject_id %in% carriers, !is.na(time_elapsed))
  cat(sprintf("Clade %s: %d participants with clade detected at either timepoint\n", cl, length(carriers)))
  if (nrow(df_cl) < 10) return(NULL)
  m <- tryCatch(
    lmer(depth ~ time_elapsed * EthnicityTot + (1 | subject_id),
         data = df_cl, REML = FALSE),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  coef_tbl <- as.data.frame(coef(summary(m)))
  coef_tbl$term  <- rownames(coef_tbl)
  coef_tbl$clade <- cl
  coef_tbl
}) %>%
  rename(estimate = Estimate, se = `Std. Error`, df = df, t = `t value`, p = `Pr(>|t|)`) %>%
  dplyr::select(clade, term, estimate, se, t, df, p) %>%
  group_by(term) %>%
  mutate(p_adj = p.adjust(p, "BH")) %>%
  ungroup()

cat("\nAll LMM coefficients (BH-adjusted p across clades per term):\n")
print(lmm_results %>% arrange(term, clade), n = Inf)

# Highlight the key interaction term
cat("\n--- Interaction term: time_elapsed:EthnicityTotSouth-Asian Surinamese ---\n")
print(
  lmm_results %>%
    filter(grepl("time_elapsed.*Ethnicity|Ethnicity.*time_elapsed", term)) %>%
    dplyr::select(clade, estimate, se, t, p, p_adj) %>%
    arrange(p)
)

# Plot: depth by timepoint × ethnicity per clade (boxplot + paired lines)
p_depth_clade <- ggplot(depth_long,
                        aes(x = timepoint, y = depth, fill = EthnicityTot)) +
  geom_boxplot(outlier.size = 0.5, width = 0.5, alpha = 0.7,
               position = position_dodge(0.6)) +
  scale_fill_manual(values = jco_palette(), name = "Ethnicity") +
  facet_wrap(~ clade, nrow = 1) +
  labs(
    title = "Alistipes putredinis clade abundance: baseline vs follow-up",
    x     = "",
    y     = "Sequencing depth"
  ) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text  = element_text(size = 9))

ggsave(
  file.path(results_dir, "clade_abundance_baseline_vs_followup.pdf"),
  p_depth_clade,
  width  = 14,
  height = 5
)
cat("Clade abundance plot saved to:",
    file.path(results_dir, "clade_abundance_baseline_vs_followup.pdf"), "\n")
