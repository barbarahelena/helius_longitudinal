## Cayman Longitudinal Analysis with Linear Mixed Models
library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggrepel)

# Theme
theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family) +
     theme(plot.title = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA, fill = NA),
           plot.background = element_rect(colour = NA, fill = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold", size = rel(0.8)),
           axis.title.y = element_text(angle=90, vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(size = rel(0.7)),
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing  = unit(0, "cm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold"),
           plot.caption = element_text(size = rel(0.5), face = "italic"),
           plot.subtitle = element_text(size=8, hjust = 0.5, face = "italic")))
}

# Data import
df_raw <- rio::import("data/shotgun/cayman_results/families_cpm_table.tsv") |> dplyr::select(-HELIBA_103370, -HELIFU_103370)
head(df_raw)[1:5,1:5]
clinical <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
dir.create("results/4_functional_change/cayman/longitudinal", showWarnings = FALSE, recursive = TRUE)

stats <- rio::import("data/shotgun/cayman_results/sample_statistics.tsv") |> 
  mutate(timepoint = case_when(str_detect(sample, "HELIBA") ~ "Baseline",
                                str_detect(sample, "HELIFU") ~ "Follow-up"),
          timepoint = as.factor(timepoint),
        log_richness= log10(richness)) |> 
  dplyr::rename(sampleID = sample) |> dplyr::select(-timepoint)
stats <- left_join(stats, clinical, by = "sampleID")
head(stats)[1:5,1:10]


# Prevalence filtering (>5%)
rownames(df_raw) <- df_raw$family
family_prevalence <- rowSums(df_raw > 0) / ncol(df_raw) * 100
prevalent_families <- names(family_prevalence)[family_prevalence > 5]

# DATA PREPARATION ----
df_raw <- df_raw %>% filter(family %in% prevalent_families)
rownames(df_raw) <- df_raw$family
df_raw$family <- NULL
df <- as.data.frame(t(as.matrix(df_raw)))
names(df_raw)
dim(df)
head(df)[1:5,1:5]
df$sampleID <- rownames(df)
fam <- ncol(df)

dftot <- left_join(df, clinical) |> droplevels()
summary(as.factor(dftot$timepoint))
summary(dftot$Fiber)

write.csv2(data.frame(family = prevalent_families), "results/4_functional_change/cayman/longitudinal/prevalent_families_list.csv", row.names = FALSE)

# LMM: Richness ~ Ethnicity * timepoint (+ log_depth) ----
model1 <- lmer(log_richness ~ EthnicityTot * timepoint + (1 | ID), data = stats)
res <- summary(model1)
ci  <- confint(model1, method = "Wald")

# Extract all Ethnicity:timepoint interaction terms
term_idx <- grep("^EthnicityTot.*:timepointfollow-up$", rownames(res$coefficients))
ci_idx   <- grep("^EthnicityTot.*:timepointfollow-up$", rownames(ci))

statres <- data.frame(
  term     = rownames(res$coefficients)[term_idx],
  estimate = res$coefficients[term_idx, 1],
  conflow  = if (length(ci_idx) > 0) ci[ci_idx, 1] else NA,
  confhigh = if (length(ci_idx) > 0) ci[ci_idx, 2] else NA,
  pval     = res$coefficients[term_idx, 5],
  stringsAsFactors = FALSE
) %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))

write.csv2(statres, "results/4_functional_change/cayman/longitudinal/lmm_richness_ethnicity_timepoint.csv", row.names = FALSE)

# Plot: spaghetti + group means + interaction p (pick one term) ----
# Choose which interaction to annotate:
# - If you want a specific ethnicity vs reference, set it here:
target_term <- statres$term[1]  # or e.g. "EthnicityTotSouth-Asian Surinamese:timepointfollow-up"

p_to_show <- statres %>%
  filter(term == target_term) %>%
  mutate(sigq = paste0("p=", formatC(pval, format = "e", digits = 2))) %>%
  slice(1)

df_means <- stats %>%
  group_by(EthnicityTot, timepoint) %>%
  summarise(mean = mean(log_richness, na.rm = TRUE),
            sd   = sd(log_richness, na.rm = TRUE),
            n    = sum(!is.na(log_richness)),
            .groups = "drop")

ymax  <- max(stats$log_richness, na.rm = TRUE)
ystat <- ymax + 0.08 * diff(range(stats$log_richness, na.rm = TRUE))

# stat_pvalue_manual wants these columns:
ann <- data.frame(
  group1 = "baseline",
  group2 = "follow-up",
  y.position = ystat,
  sigq = p_to_show$sigq
)

(gg_richness_lmm <- ggplot() +
  geom_line(data = stats, aes(group = ID, stats, x = timepoint, y = log_richness, color = EthnicityTot), alpha = 0.08, linewidth = 0.4) +
  geom_point(data = stats, aes(group = ID, stats, x = timepoint, y = log_richness, color = EthnicityTot), alpha = 0.08, size = 0.7) +
  geom_line(data = df_means, aes(x = timepoint, y = mean, group = EthnicityTot, color = EthnicityTot), linewidth = 0.9) +
  geom_point(data = df_means, aes(x = timepoint, y = mean, color = EthnicityTot), size = 1.6) +
  geom_errorbar(
    data = df_means,
    aes(x = timepoint, ymin = mean - (sd/sqrt(n)), ymax = mean + (sd/sqrt(n)), color = EthnicityTot),
    width = 0.12, linewidth = 0.4
  ) +
  stat_pvalue_manual(ann, label = "sigq", tip.length = 0, bracket.shorten = 0.15, size = 3.5) +
  scale_color_jco() +
  theme_Publication() +
  labs(
    x = "Timepoint",
    y = "log10(Richness)",
    title = "Richness over time by ethnicity",
    color = ""
  ))
ggsave("results/4_functional_change/cayman/longitudinal/richness_lmm_spaghetti_means.pdf",
       gg_richness_lmm, width = 6.5, height = 4.5)

# LMM: %CAZy reads ~ Ethnicity * timepoint (unadjusted + adjusted) ----
# Create baseline covariates for stats
stats_baseline_cov <- stats %>%
  filter(timepoint == "baseline") %>%
  dplyr::select(ID, Age_baseline = Age, BMI_baseline = BMI,
                Smoking_baseline = Smoking, PPI_baseline = PPI) %>%
  distinct(ID, .keep_all = TRUE)

stats_adj <- stats %>%
  left_join(stats_baseline_cov, by = "ID") %>%
  filter(!is.na(Age_baseline), !is.na(Sex), !is.na(BMI_baseline),
         !is.na(Smoking_baseline), !is.na(PPI_baseline))

# Unadjusted
model_cazy_unadj <- lmer(pct_cazy_reads ~ EthnicityTot * timepoint + (1 | ID), data = stats)
res_cazy_unadj <- summary(model_cazy_unadj)
ci_cazy_unadj <- confint(model_cazy_unadj, method = "Wald")

term_idx <- grep("^EthnicityTot.*:timepointfollow-up$", rownames(res_cazy_unadj$coefficients))
ci_idx <- grep("^EthnicityTot.*:timepointfollow-up$", rownames(ci_cazy_unadj))

statres_cazy_unadj <- data.frame(
  term     = rownames(res_cazy_unadj$coefficients)[term_idx],
  estimate = res_cazy_unadj$coefficients[term_idx, 1],
  conflow  = if (length(ci_idx) > 0) ci_cazy_unadj[ci_idx, 1] else NA,
  confhigh = if (length(ci_idx) > 0) ci_cazy_unadj[ci_idx, 2] else NA,
  pval     = res_cazy_unadj$coefficients[term_idx, 5],
  stringsAsFactors = FALSE
) %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv2(statres_cazy_unadj, "results/4_functional_change/cayman/longitudinal/lmm_pctcazy_ethnicity_timepoint_unadjusted.csv", row.names = FALSE)

# Adjusted
model_cazy_adj <- lmer(pct_cazy_reads ~ EthnicityTot * timepoint + Age_baseline +
                         Sex + BMI_baseline + Smoking_baseline + PPI_baseline +
                         (1 | ID), data = stats_adj)
res_cazy_adj <- summary(model_cazy_adj)
ci_cazy_adj <- confint(model_cazy_adj, method = "Wald")

term_idx <- grep("^EthnicityTot.*:timepointfollow-up$", rownames(res_cazy_adj$coefficients))
ci_idx <- grep("^EthnicityTot.*:timepointfollow-up$", rownames(ci_cazy_adj))

statres_cazy_adj <- data.frame(
  term     = rownames(res_cazy_adj$coefficients)[term_idx],
  estimate = res_cazy_adj$coefficients[term_idx, 1],
  conflow  = if (length(ci_idx) > 0) ci_cazy_adj[ci_idx, 1] else NA,
  confhigh = if (length(ci_idx) > 0) ci_cazy_adj[ci_idx, 2] else NA,
  pval     = res_cazy_adj$coefficients[term_idx, 5],
  stringsAsFactors = FALSE
) %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv2(statres_cazy_adj, "results/4_functional_change/cayman/longitudinal/lmm_pctcazy_ethnicity_timepoint_adjusted.csv", row.names = FALSE)

# GENE FAMILY-LEVEL LMM ----
# Prepare abundance matrix: df (samples x gene families), rownames(df) are sample IDs
# Join with clinical by sampleID (rownames)
gene_families <- setdiff(colnames(df), c("sampleID"))

# LMMs for each gene family
statres <- data.frame()
for (gf in gene_families) {
  dftot$mb <- log10(dftot[[gf]] + 1)
    model1 <- lmer(mb ~ EthnicityTot * timepoint + (1|ID), data = dftot)
    res <- summary(model1)
    confint_model1 <- confint(model1, method = "Wald")
    interaction_row <- grep("EthnicityTotSouth-Asian Surinamese:timepointfollow-up", rownames(res$coefficients))
    ci_row <- grep("EthnicityTotSouth-Asian Surinamese:timepointfollow-up", rownames(confint_model1))
    statres <- rbind(statres, data.frame(
        family = gf,
        estimate = res$coefficients[interaction_row, 1],
        conflow = ifelse(length(ci_row) > 0, confint_model1[ci_row, 1], NA),
        confhigh = ifelse(length(ci_row) > 0, confint_model1[ci_row, 2], NA),
        pval = res$coefficients[interaction_row, 5]
      ))
}
statres <- as.data.frame(statres)
statres <- statres %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv2(statres, "results/4_functional_change/cayman/longitudinal/lmm_ethnicity_timepoint_results.csv", row.names = FALSE)

# LONGITUDINAL PLOTS ----
statres_sig <- statres %>% filter(pval < 0.05) %>% arrange(pval)

plist <- list()
for(i in 1:min(nrow(statres_sig), 20)){
  nm <- statres_sig$family[i]
  dftot$mb <- log10(dftot[[nm]] + 1)

  df_means <- dftot %>%
    group_by(EthnicityTot, timepoint) %>%
    summarise(mean = mean(mb, na.rm = TRUE), sd = sd(mb, na.rm = TRUE),
              n = n(), .groups = "drop")

  res_lmm <- statres_sig %>%
    filter(family == nm) %>%
    mutate(group1 = "baseline", group2 = "follow-up",
            sigq = paste0("p=", formatC(pval, format = "e", digits = 2))) %>%
    dplyr::select(-family)

  mbmax <- ifelse(max(dftot$mb, na.rm = TRUE) < 0,
                  max(dftot$mb, na.rm = TRUE) * 0.8,
                  max(dftot$mb, na.rm = TRUE) * 1.2)
  mbstat <- ifelse(max(dftot$mb, na.rm = TRUE) < 0,
                    max(dftot$mb, na.rm = TRUE) * 0.9,
                    max(dftot$mb, na.rm = TRUE) * 1.1)
  mbmin <- min(dftot$mb, na.rm = TRUE)

  pl2 <- ggplot() +
    geom_line(data = dftot, aes(x = timepoint, y = mb, color = EthnicityTot, group = ID),
              alpha = 0.15, linewidth = 0.5) +
    geom_point(data = dftot, aes(x = timepoint, y = mb, color = EthnicityTot),
                alpha = 0.15, size = 0.8) +
    geom_line(data = df_means, aes(x = timepoint, y = mean, color = EthnicityTot, group = EthnicityTot),
              alpha = 1, linewidth = 0.8) +
    geom_point(data = df_means, aes(x = timepoint, y = mean, color = EthnicityTot, group = EthnicityTot),
                alpha = 1, size = 1.3) +
    geom_errorbar(data = df_means,
                  aes(ymin = mean - (sd/sqrt(n)), ymax = mean + (sd/sqrt(n)),
                      x = timepoint, color = EthnicityTot), width = 0.1) +
    stat_pvalue_manual(res_lmm, y.position = mbstat, label = "sigq",
                        tip.length = 0, bracket.shorten = 0.1, size = 4) +
    scale_color_jco() +
    coord_cartesian(ylim = c(mbmin, mbmax)) +
    theme_Publication() +
    labs(x = "Timepoint", y = "log10(RPKM + 1)",
          title = nm,
          color = "")

  plist[[i]] <- pl2
}

n_plots <- length(plist)
n_cols <- 3
n_rows <- ceiling(n_plots / n_cols)

plots <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
                    labels = LETTERS[1:n_plots], nrow = n_rows, ncol = n_cols)

ggsave("results/4_functional_change/cayman/longitudinal/significant_cayman_timepoint_ethnicity.pdf", plots,
        width = 12, height = 4 * n_rows)

# ADJUSTED LMMs ----
# Model: CAZy ~ EthnicityTot * timepoint + Age_baseline + Sex + BMI + Smoking + PPI + (1|ID)
library(emmeans)

# Create baseline covariates (fixed per ID)
baseline_covariates <- dftot %>%
  filter(timepoint == "baseline") %>%
  dplyr::select(ID, Age_baseline = Age, BMI_baseline = BMI,
                Smoking_baseline = Smoking, PPI_baseline = PPI, AB_baseline = AB) %>%
  distinct(ID, .keep_all = TRUE)

# Merge stats for total_reads
dftot_adj <- dftot %>%
  left_join(stats %>% dplyr::select(sampleID, total_reads), by = "sampleID") %>%
  left_join(baseline_covariates, by = "ID") %>%
  filter(!is.na(Age_baseline), !is.na(Sex), !is.na(BMI_baseline),
         !is.na(Smoking_baseline), !is.na(PPI_baseline))

# Adjusted LMMs for each gene family
statres_adj <- data.frame()
for (gf in gene_families) {
  print(gf)
  dftot_adj$mb <- log10(dftot_adj[[gf]] + 1)
  tryCatch({
    model_adj <- lmer(mb ~ EthnicityTot * timepoint + Age_baseline +
                        Sex + BMI_baseline + Smoking_baseline + PPI_baseline +
                        (1|ID), data = dftot_adj)
    res <- summary(model_adj)
    ci <- confint(model_adj, method = "Wald")
    interaction_row <- grep("EthnicityTotSouth-Asian Surinamese:timepointfollow-up", rownames(res$coefficients))
    ci_row <- grep("EthnicityTotSouth-Asian Surinamese:timepointfollow-up", rownames(ci))
    statres_adj <- rbind(statres_adj, data.frame(
      family = gf,
      estimate = res$coefficients[interaction_row, 1],
      conflow = ifelse(length(ci_row) > 0, ci[ci_row, 1], NA),
      confhigh = ifelse(length(ci_row) > 0, ci[ci_row, 2], NA),
      pval = res$coefficients[interaction_row, 5]
    ))
  }, error = function(e) NULL)
}
statres_adj <- as.data.frame(statres_adj) %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv2(statres_adj, "results/4_functional_change/cayman/longitudinal/lmm_ethnicity_timepoint_adjusted.csv", row.names = FALSE)

# FOREST PLOT: Ethnicity x Timepoint interaction estimates ----
# Unadjusted
statres_forest <- statres %>%
  mutate(sig = ifelse(padj < 0.05, "FDR < 0.05", ifelse(pval < 0.05, "p < 0.05", "NS")),
         family = fct_reorder(family, estimate))

(p_forest_unadj <- ggplot(statres_forest %>% filter(pval < 0.05),
                          aes(x = estimate, y = family, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(xmin = conflow, xmax = confhigh), size = 0.4, linewidth = 0.5) +
  scale_color_manual(values = c("FDR < 0.05" = "#E64B35", "p < 0.05" = "#4DBBD5", "NS" = "gray70")) +
  theme_Publication() +
  labs(x = "Estimate (Ethnicity x Timepoint interaction)", y = "",
       title = "Unadjusted LMM: Ethnicity x Timepoint interaction",
       color = ""))
ggsave("results/4_functional_change/cayman/longitudinal/forest_ethnicity_timepoint_unadjusted.pdf",
       p_forest_unadj, width = 8, height = max(4, sum(statres$pval < 0.05) * 0.25 + 2))

# Adjusted
statres_adj_forest <- statres_adj %>%
  mutate(sig = ifelse(padj < 0.05, "FDR < 0.05", ifelse(pval < 0.05, "p < 0.05", "NS")),
         family = fct_reorder(family, estimate))

(p_forest_adj <- ggplot(statres_adj_forest %>% filter(pval < 0.05),
                        aes(x = estimate, y = family, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(xmin = conflow, xmax = confhigh), size = 0.4, linewidth = 0.5) +
  scale_color_manual(values = c("FDR < 0.05" = "#E64B35", "p < 0.05" = "#4DBBD5", "NS" = "gray70")) +
  theme_Publication() +
  labs(x = "Estimate (Ethnicity x Timepoint interaction)", y = "",
       title = "Adjusted LMM: Ethnicity x Timepoint interaction",
       subtitle = "Adjusted for baseline age, sex, BMI, smoking, PPI",
       color = ""))
ggsave("results/4_functional_change/cayman/longitudinal/forest_ethnicity_timepoint_adjusted.pdf",
       p_forest_adj, width = 8, height = max(4, sum(statres_adj$pval < 0.05) * 0.25 + 1))

# Forest plots: FDR significant only ----
if (sum(statres$padj < 0.05) > 0) {
  statres_fdr <- statres %>% filter(padj < 0.05) %>% mutate(family = fct_reorder(family, estimate))
  (p_forest_unadj_fdr <- ggplot(statres_fdr, aes(x = estimate, y = family)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(xmin = conflow, xmax = confhigh), size = 0.4, linewidth = 0.5, color = "#E64B35") +
    theme_Publication() +
    labs(x = "Estimate (Ethnicity x Timepoint interaction)", y = "",
         title = "Unadjusted LMM: FDR significant interactions"))
  ggsave("results/4_functional_change/cayman/longitudinal/forest_ethnicity_timepoint_unadjusted_fdr.pdf",
         p_forest_unadj_fdr, width = 8, height = max(4, nrow(statres_fdr) * 0.25 + 1))
}

if (sum(statres_adj$padj < 0.05) > 0) {
  statres_adj_fdr <- statres_adj %>% filter(padj < 0.05) %>% mutate(family = fct_reorder(family, estimate))
  (p_forest_adj_fdr <- ggplot(statres_adj_fdr, aes(x = estimate, y = family)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(xmin = conflow, xmax = confhigh), size = 0.4, linewidth = 0.5, color = "#E64B35") +
    theme_Publication() +
    labs(x = "Estimate (Ethnicity x Timepoint interaction)", y = "",
         title = "Differential changes in ethnic groups over time",
         subtitle = "Adjusted for baseline age, sex, BMI, smoking, PPI"))
  ggsave("results/4_functional_change/cayman/longitudinal/forest_ethnicity_timepoint_adjusted_fdr.pdf",
         p_forest_adj_fdr, width = 8, height = max(4, nrow(statres_adj_fdr) * 0.25 + 1))
}

# DIET-ADJUSTED LMMs ----
# Model: CAZy ~ EthnicityTot * timepoint + Age_baseline + Sex + BMI + Smoking + PPI + Fiber + (1|ID)

diet_baseline <- dftot %>%
  filter(timepoint == "baseline") %>%
  dplyr::select(ID, Fiber_BL = Fiber, Carbohydrates_BL = Carbohydrates) %>%
  distinct(ID, .keep_all = TRUE)

dftot_diet <- dftot_adj %>%
  left_join(diet_baseline, by = "ID") %>%
  filter(!is.na(Fiber_BL), !is.na(Carbohydrates_BL))

statres_diet <- data.frame(family   = character(),
                           estimate = numeric(),
                           conflow  = numeric(),
                           confhigh = numeric(),
                           pval     = numeric())
for (gf in gene_families) {
  dftot_diet$mb <- log10(dftot_diet[[gf]] + 1)
  tryCatch({
    model_diet <- lmer(mb ~ EthnicityTot * timepoint + Age_baseline +
                         Sex + BMI_baseline + Smoking_baseline + PPI_baseline +
                         Fiber_BL + Carbohydrates_BL + (1|ID), data = dftot_diet)
    res <- summary(model_diet)
    ci  <- confint(model_diet, method = "Wald")
    interaction_row <- grep("EthnicityTotSouth-Asian Surinamese:timepointfollow-up", rownames(res$coefficients))
    ci_row          <- grep("EthnicityTotSouth-Asian Surinamese:timepointfollow-up", rownames(ci))
    if (length(interaction_row) == 0) return(NULL)
    statres_diet <- rbind(statres_diet, data.frame(
      family   = gf,
      estimate = res$coefficients[interaction_row, 1],
      conflow  = ifelse(length(ci_row) > 0, ci[ci_row, 1], NA_real_),
      confhigh = ifelse(length(ci_row) > 0, ci[ci_row, 2], NA_real_),
      pval     = res$coefficients[interaction_row, 5]
    ))
  }, error = function(e) NULL)
}

statres_diet <- as.data.frame(statres_diet) %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))

write.csv2(statres_diet,
           "results/4_functional_change/cayman/longitudinal/lmm_ethnicity_timepoint_dietary.csv",
           row.names = FALSE)

# Forest plot: dietary model, FDR significant only
if (sum(statres_diet$padj < 0.05) > 0) {
  statres_diet_fdr <- statres_diet %>%
    filter(padj < 0.05) %>%
    mutate(family    = fct_reorder(family, estimate),
           direction = ifelse(estimate > 0, "SAS", "Dutch"))

  (p_forest_diet <- ggplot(statres_diet_fdr,
                           aes(x = estimate, y = family, colour = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(xmin = conflow, xmax = confhigh), size = 0.4, linewidth = 0.5) +
    scale_colour_manual(values = c("Dutch" = "#2166AC", "SAS" = "#E6B800"),
                        name   = "",
                        labels = c("Dutch" = "More increase in Dutch",
                                   "SAS"   = "More increase in SAS")) +
    theme_Publication() +
    labs(x        = "Interaction effect (\u00b1 95% CI)",
         y        = "",
         title    = "Differential CAZyme dynamics by ethnicity",
         subtitle = "Adjusted for baseline age, sex, BMI, smoking, PPI, fiber, carbohydrates (FDR < 0.05)"))

  ggsave("results/4_functional_change/cayman/longitudinal/forest_ethnicity_timepoint_dietary.pdf",
         p_forest_diet,
         width  = 8,
         height = max(4, nrow(statres_diet_fdr) * 0.25 + 1))
}

# PREDICTED MARGINAL MEANS PLOTS (emmeans) ----
# Much cleaner than spaghetti plots: shows model-estimated means with CIs
statres_adj_sig <- statres_adj %>% filter(pval < 0.05) %>% arrange(pval)

plist_emm <- list()
for (i in seq_len(min(nrow(statres_adj_sig), 20))) {
  nm <- statres_adj_sig$family[i]
  dftot_adj$mb <- log10(dftot_adj[[nm]] + 1)

  tryCatch({
    model_adj <- lmer(mb ~ EthnicityTot * timepoint + Age_baseline +
                        Sex + BMI_baseline + Smoking_baseline + PPI_baseline +
                        (1|ID), data = dftot_adj)

    # Get predicted marginal means per ethnicity x timepoint
    emm <- emmeans(model_adj, ~ EthnicityTot | timepoint)
    emm_df <- as.data.frame(emm)

    pval_label <- formatC(statres_adj_sig$pval[i], format = "e", digits = 2)

    pl <- ggplot(emm_df, aes(x = timepoint, y = emmean, color = EthnicityTot, group = EthnicityTot)) +
      geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),
                      position = position_dodge(width = 0.3), size = 0.5, linewidth = 0.6) +
      geom_line(position = position_dodge(width = 0.3), linewidth = 0.7) +
      scale_color_jco() +
      theme_Publication() +
      labs(x = "Timepoint", y = "Predicted log10(RPKM + 1)",
           title = nm,
           subtitle = paste0("Ethnicity x Timepoint p=", pval_label),
           color = "")

    plist_emm[[i]] <- pl
  }, error = function(e) NULL)
}

n_plots <- length(plist_emm)
if (n_plots > 0) {
  n_cols <- 3
  n_rows <- ceiling(n_plots / n_cols)
  plots_emm <- ggarrange(plotlist = plist_emm, common.legend = TRUE, legend = "bottom",
                          labels = LETTERS[1:n_plots], nrow = n_rows, ncol = n_cols)
  ggsave("results/4_functional_change/cayman/longitudinal/significant_cayman_adjusted_emmeans.pdf", plots_emm,
         width = 12, height = 4 * n_rows)
}
