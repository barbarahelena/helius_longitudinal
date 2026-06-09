## Cayman Descriptive Analysis & QC Plots
library(tidyverse)
library(ggsci)
library(ggpubr)
library(ggrepel)

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
           strip.text = element_text(face="bold")))
}

df_raw <- rio::import("data/shotgun/cayman_results/families_cpm_table.tsv") |> 
  dplyr::select(-HELIBA_103370, -HELIFU_103370)
head(df_raw)[1:5,1:5]
rownames(df_raw) <- df_raw$family
df_raw$family <- NULL
df <- as.data.frame(t(as.matrix(df_raw)))
names(df_raw)
dim(df)
head(df)[1:5,1:5]
df$sampleID <- rownames(df)
fam <- ncol(df)

clinical <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
stats <- rio::import("data/shotgun/cayman_results/sample_statistics.tsv") |>
  rename(sampleID = sample) |>
  dplyr::select(sampleID, total_reads)
dftot <- left_join(df, clinical) |>
  left_join(stats, by = "sampleID") |>
  droplevels()
summary(as.factor(dftot$timepoint))

# Create cross-sectional output directory
dir.create("results/4_functional_change/cayman", showWarnings = FALSE, recursive = TRUE)
cross_dir <- "results/4_functional_change/cayman/crossectional"
dir.create(cross_dir, showWarnings = FALSE, recursive = TRUE)

# Add timepoint variable to dftot
if (!"timepoint" %in% names(dftot)) {
  dftot <- dftot |> 
    mutate(timepoint = case_when(
      str_detect(sampleID, "HELIBA") ~ "Baseline",
      str_detect(sampleID, "HELIFU") ~ "Follow-up",
      TRUE ~ NA_character_
    ))
}

# Test difference between ethnic groups for each gene family at Baseline and Follow-up
# Assumes ethnicity column is present in dftot (adjust name if needed)
gene_families <- setdiff(names(df), c("sampleID"))
results_list <- list()
for (tp in c("baseline", "follow-up")) {
  df_tp <- dftot %>% filter(timepoint == tp)
  res <- purrr::map_dfr(gene_families, function(gf) {
    # Only test if there is variation
    if (length(unique(df_tp[[gf]])) > 1 && length(unique(df_tp$EthnicityTot)) > 1) {
      kw <- wilcox.test(df_tp[[gf]] ~ df_tp$EthnicityTot)
      tibble(
        gene_family = gf,
        timepoint = tp,
        p_value = kw$p.value
      )
    } else {
      tibble(
        gene_family = gf,
        timepoint = tp,
        p_value = NA_real_
      )
    }
  })
  results_list[[tp]] <- res
}
all_results <- bind_rows(results_list)
all_results <- all_results %>% mutate(p_adj = p.adjust(p_value, method = "fdr")) |> arrange(p_adj)
print(all_results |> filter(p_adj < 0.05)) |> arrange(p_adj)
write.csv2(all_results, file.path(cross_dir, "kruskal_ethnicity_by_genefamily.csv"), row.names = FALSE)

# --- CAYMAN GENE FAMILY ABUNDANCE DIFFERENCES BY ETHNICITY ---
# Prepare wide abundance table (log10 transformed)
gene_families <- setdiff(names(df), c("sampleID"))
df_wide <- df %>%
  dplyr::select(sampleID, all_of(gene_families))
df_tot <- df_wide %>%
  left_join(clinical, by = "sampleID") %>%
  filter(!is.na(EthnicityTot)) %>%
  droplevels()

# Baseline analysis
if (any(df_tot$timepoint == "baseline")) {
  df_baseline <- df_tot %>% filter(timepoint == "baseline")
  statres_baseline <- data.frame()
  for (gf in gene_families) {
    df_baseline$mb <- log10(df_baseline[[gf]] + 1)
      model_baseline <- lm(mb ~ EthnicityTot, data = df_baseline)
      res <- summary(model_baseline)
      if (nrow(res$coefficients) >= 2) {
        confint_baseline <- confint(model_baseline)
        statres_baseline <- rbind(statres_baseline, data.frame(
          gene_family = gf,
          pval = res$coefficients[2, 4],
          estimate = res$coefficients[2, 1],
          se = res$coefficients[2, 2],
          conflow = confint_baseline[2, 1],
          confhigh = confint_baseline[2, 2]
        ))
      }
  }
  statres_baseline <- statres_baseline %>%
    arrange(pval) %>%
    mutate(padj = p.adjust(pval, method = "fdr"))
  write.csv2(statres_baseline, file.path(cross_dir, "baseline_ethnicity_abundance.csv"), row.names = FALSE)
  # Volcano plot
  baseline_sig <- statres_baseline %>%
    mutate(sig_level = case_when(
      padj < 0.05 & abs(estimate) > 0.5 ~ "Significant & Large Effect",
      padj < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"))
  p_baseline_volcano <- ggplot(baseline_sig, aes(x = estimate, y = -log10(pval), color = sig_level)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_text_repel(data = baseline_sig %>% filter(padj < 0.05),
                    aes(label = gene_family), size = 3, max.overlaps = 20,
                    box.padding = 0.5, point.padding = 0.3, color = "black") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Significant & Large Effect" = "red",
                                   "Significant" = "orange",
                                   "Not Significant" = "gray70")) +
    theme_Publication() +
    labs(x = "log10(Abundance) Difference (South-Asian Surinamese - Dutch)",
         y = "-log10(p-value)",
         title = "Baseline Cayman Abundance Differences Between Ethnicities", color = "")
  ggsave(file.path(cross_dir, "baseline_abundance_volcano.pdf"), p_baseline_volcano, width = 7, height = 7, device = cairo_pdf)
}
statres_baseline |> filter(estimate > 0 & padj < 0.05) |> nrow()
statres_baseline |> filter(estimate < 0 & padj < 0.05) |> nrow()

# Follow-up analysis
if (any(df_tot$timepoint == "follow-up")) {
  df_followup <- df_tot %>% filter(timepoint == "follow-up")
  statres_followup <- data.frame()
  for (gf in gene_families) {
    df_followup$mb <- log10(df_followup[[gf]] + 1)
    tryCatch({
      model_followup <- lm(mb ~ EthnicityTot, data = df_followup)
      res <- summary(model_followup)
      if (nrow(res$coefficients) >= 2) {
        confint_followup <- confint(model_followup)
        statres_followup <- rbind(statres_followup, data.frame(
          gene_family = gf,
          pval = res$coefficients[2, 4],
          estimate = res$coefficients[2, 1],
          se = res$coefficients[2, 2],
          conflow = confint_followup[2, 1],
          confhigh = confint_followup[2, 2]
        ))
      }
    }, error = function(e) NULL)
  }
  statres_followup <- statres_followup %>%
    arrange(pval) %>%
    mutate(padj = p.adjust(pval, method = "fdr"))
  write.csv2(statres_followup, file.path(cross_dir, "followup_ethnicity_abundance.csv"), row.names = FALSE)
  # Volcano plot
  followup_sig <- statres_followup %>%
    mutate(sig_level = case_when(
      padj < 0.05 & abs(estimate) > 0.5 ~ "Significant & Large Effect",
      padj < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"))
  p_followup_volcano <- ggplot(followup_sig, aes(x = estimate, y = -log10(pval), color = sig_level)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_text_repel(data = followup_sig %>% filter(padj < 0.05),
                    aes(label = gene_family), size = 3, max.overlaps = 20,
                    box.padding = 0.5, point.padding = 0.3, color = "black") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Significant & Large Effect" = "red",
                                   "Significant" = "orange",
                                   "Not Significant" = "gray70")) +
    theme_Publication() +
    labs(x = "log10(Abundance) Difference (South-Asian Surinamese - Dutch)",
         y = "-log10(p-value)",
         title = "Follow-up Cayman Abundance Differences Between Ethnicities", color = "")
  ggsave(file.path(cross_dir, "followup_abundance_volcano.pdf"), p_followup_volcano, width = 7, height = 7, device = cairo_pdf)
}
statres_followup |> filter(estimate > 0 & padj < 0.05) |> nrow()
statres_followup |> filter(estimate < 0 & padj < 0.05) |> nrow()

# ---------------------------------------------------------------------------
# Save unified cross-sectional results table (for LMM overlap analysis)
# ---------------------------------------------------------------------------
cs_combined <- statres_baseline |>
    dplyr::select(gene_family, padj_baseline = padj,
                  cs_baseline_direction = estimate) |>
    mutate(cs_baseline_direction = ifelse(cs_baseline_direction > 0, "SAS higher", "Dutch higher"),
           sig_baseline = padj_baseline < 0.05) |>
    full_join(
        statres_followup |>
            dplyr::select(gene_family, padj_followup = padj,
                          cs_followup_direction = estimate) |>
            mutate(cs_followup_direction = ifelse(cs_followup_direction > 0, "SAS higher", "Dutch higher"),
                   sig_followup = padj_followup < 0.05),
        by = "gene_family"
    )

cat("\nOf", length(gene_families), "gene families,",
    sum(cs_combined$sig_baseline | cs_combined$sig_followup, na.rm = TRUE),
    "show a significant ethnicity difference at baseline or follow-up (FDR < 0.05)\n")
cat("  Baseline: ", sum(cs_combined$sig_baseline, na.rm = TRUE), "significant\n")
cat("  Follow-up:", sum(cs_combined$sig_followup, na.rm = TRUE), "significant\n")
cat("  Overlap:  ", sum(cs_combined$sig_baseline & cs_combined$sig_followup, na.rm = TRUE),
    "significant at both\n")

write.csv2(cs_combined,
           file.path(cross_dir, "crosssectional_results_combined.csv"),
           row.names = FALSE)

# --- ADJUSTED MODELS: CAZy ~ EthnicityTot + Age + Sex + BMI + DM + PPI + total_reads ---
# Merge total_reads into df_tot
df_tot <- df_tot %>%
  left_join(stats, by = "sampleID")

# Baseline adjusted
if (any(df_tot$timepoint == "baseline")) {
  df_baseline <- df_tot %>% filter(timepoint == "baseline") %>%
    filter(!is.na(Age), !is.na(Sex), !is.na(BMI), !is.na(DM), !is.na(PPI), !is.na(total_reads))
  statres_baseline_adj <- data.frame()
  for (gf in gene_families) {
    df_baseline$mb <- log10(df_baseline[[gf]] + 1)
    tryCatch({
      model <- lm(mb ~ EthnicityTot + Age + Sex + BMI + DM + PPI + total_reads, data = df_baseline)
      res <- summary(model)
      if (nrow(res$coefficients) >= 2) {
        ci <- confint(model)
        statres_baseline_adj <- rbind(statres_baseline_adj, data.frame(
          gene_family = gf,
          pval = res$coefficients[2, 4],
          estimate = res$coefficients[2, 1],
          se = res$coefficients[2, 2],
          conflow = ci[2, 1],
          confhigh = ci[2, 2]
        ))
      }
    }, error = function(e) NULL)
  }
  statres_baseline_adj <- statres_baseline_adj %>%
    arrange(pval) %>%
    mutate(padj = p.adjust(pval, method = "fdr"))
  write.csv2(statres_baseline_adj, file.path(cross_dir, "baseline_ethnicity_abundance_adjusted.csv"), row.names = FALSE)

  baseline_adj_sig <- statres_baseline_adj %>%
    mutate(sig_level = case_when(
      padj < 0.05 & abs(estimate) > 0.5 ~ "Significant & Large Effect",
      padj < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"))
  p_baseline_adj_volcano <- ggplot(baseline_adj_sig, aes(x = estimate, y = -log10(pval), color = sig_level)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_text_repel(data = baseline_adj_sig %>% filter(padj < 0.05),
                    aes(label = gene_family), size = 3, max.overlaps = 20,
                    box.padding = 0.5, point.padding = 0.3, color = "black") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Significant & Large Effect" = "red",
                                   "Significant" = "orange",
                                   "Not Significant" = "gray70")) +
    theme_Publication() +
    labs(x = "log10(Abundance) Difference (South-Asian Surinamese - Dutch)",
         y = "-log10(p-value)",
         title = "Baseline (adjusted) CAZy Differences Between Ethnicities",
         subtitle = "Adjusted for age, sex, BMI, diabetes, PPI, total reads",
         color = "")
  ggsave(file.path(cross_dir, "baseline_abundance_volcano_adjusted.pdf"), p_baseline_adj_volcano, width = 7, height = 7)
}

# Follow-up adjusted
if (any(df_tot$timepoint == "follow-up")) {
  df_followup <- df_tot %>% filter(timepoint == "follow-up") %>%
    filter(!is.na(Age), !is.na(Sex), !is.na(BMI), !is.na(DM), !is.na(PPI), !is.na(total_reads))
  statres_followup_adj <- data.frame()
  for (gf in gene_families) {
    df_followup$mb <- log10(df_followup[[gf]] + 1)
    tryCatch({
      model <- lm(mb ~ EthnicityTot + Age + Sex + BMI + DM + PPI + total_reads, data = df_followup)
      res <- summary(model)
      if (nrow(res$coefficients) >= 2) {
        ci <- confint(model)
        statres_followup_adj <- rbind(statres_followup_adj, data.frame(
          gene_family = gf,
          pval = res$coefficients[2, 4],
          estimate = res$coefficients[2, 1],
          se = res$coefficients[2, 2],
          conflow = ci[2, 1],
          confhigh = ci[2, 2]
        ))
      }
    }, error = function(e) NULL)
  }
  statres_followup_adj <- statres_followup_adj %>%
    arrange(pval) %>%
    mutate(padj = p.adjust(pval, method = "fdr"))
  write.csv2(statres_followup_adj, file.path(cross_dir, "followup_ethnicity_abundance_adjusted.csv"), row.names = FALSE)

  followup_adj_sig <- statres_followup_adj %>%
    mutate(sig_level = case_when(
      padj < 0.05 & abs(estimate) > 0.5 ~ "Significant & Large Effect",
      padj < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"))
  p_followup_adj_volcano <- ggplot(followup_adj_sig, aes(x = estimate, y = -log10(pval), color = sig_level)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_text_repel(data = followup_adj_sig %>% filter(padj < 0.05),
                    aes(label = gene_family), size = 3, max.overlaps = 20,
                    box.padding = 0.5, point.padding = 0.3, color = "black") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Significant & Large Effect" = "red",
                                   "Significant" = "orange",
                                   "Not Significant" = "gray70")) +
    theme_Publication() +
    labs(x = "log10(Abundance) Difference (South-Asian Surinamese - Dutch)",
         y = "-log10(p-value)",
         title = "Follow-up (adjusted) CAZy Differences Between Ethnicities",
         subtitle = "Adjusted for age, sex, BMI, diabetes, PPI, total reads",
         color = "")
  ggsave(file.path(cross_dir, "followup_abundance_volcano_adjusted.pdf"), p_followup_adj_volcano, width = 7, height = 7)
}
