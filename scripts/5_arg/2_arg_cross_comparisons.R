## ARG Baseline Comparisons Between Ethnicities
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
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")))
}

df_raw <- rio::import("data/shotgun/arg/all_samples.merged_arg_counts.tsv")
clinical <- readRDS("data/clinicaldata_long.RDS")
dir.create("results/5_arg", showWarnings = FALSE, recursive = TRUE)
dir.create("results/5_arg/crosssectional", showWarnings = FALSE, recursive = TRUE)

# Get prevalent genes (>5%)
gene_prevalence <- df_raw %>%
  filter(Prevalence == 1) |> 
  group_by(Gene_Symbol, Class, Subclass) %>%
  summarise(n_samples = n_distinct(Sample),
            prevalence_pct = (n_distinct(Sample) / n_distinct(df_raw$Sample)) * 100,
            .groups = "drop")

prevalent_genes <- gene_prevalence %>%
  filter(prevalence_pct > 5) %>%
  pull(Gene_Symbol)

# GENE-LEVEL PREVALENCE ----
# Get all baseline samples first
baseline_samples <- clinical %>% 
  filter(sampleID %in% df_raw$Sample) |> 
  filter(timepoint == "baseline") %>%
  pull(sampleID) %>%
  unique()

prevalence_data <- df_raw %>%
  filter(Gene_Symbol %in% prevalent_genes) %>%
  mutate(sampleID = Sample,
         present = ifelse(Prevalence > 0 & Mapped_Reads > 0, 1, 0)) %>%
  dplyr::select(sampleID, Gene_Symbol, present) %>%
  # Complete all combinations of samples and genes
  complete(sampleID = baseline_samples,
           Gene_Symbol = prevalent_genes,
           fill = list(present = 0))

prevalence_clin <- prevalence_data %>%
  left_join(clinical, by = "sampleID") %>%
  filter(timepoint == "baseline", !is.na(EthnicityTot)) %>%
  droplevels()

# Test each gene
prevalence_results <- data.frame()
for(gene in prevalent_genes){
  gene_data <- prevalence_clin %>% filter(Gene_Symbol == gene)
  if(nrow(gene_data) == 0) next

  prev_by_eth <- gene_data %>%
    group_by(EthnicityTot) %>%
    summarise(n_total = n(), n_present = sum(present),
              prevalence_pct = (n_present / n_total) * 100, .groups = "drop_last")

  if(sum(gene_data$present) == 0 | sum(gene_data$present) == nrow(gene_data)) next

  contingency <- table(gene_data$EthnicityTot, gene_data$present)

  test_result <- chisq.test(contingency)
  pval <- test_result$p.value
  or <- (contingency[1,2] * contingency[2,1]) / (contingency[1,1] * contingency[2,2])

    gene_info <- gene_prevalence %>%
      filter(Gene_Symbol == gene) %>%
      dplyr::select(Class, Subclass, prevalence_pct)

    prevalence_results <- rbind(prevalence_results, data.frame(
      gene = gene,
      class = ifelse(nrow(gene_info) > 0, gene_info$Class, NA),
      subclass = ifelse(nrow(gene_info) > 0, gene_info$Subclass, NA),
      overall_prevalence = ifelse(nrow(gene_info) > 0, gene_info$prevalence_pct, NA),
      prev_group1 = prev_by_eth$prevalence_pct[1],
      prev_group2 = prev_by_eth$prevalence_pct[2],
      prev_diff = prev_by_eth$prevalence_pct[2] - prev_by_eth$prevalence_pct[1],
      odds_ratio = or,
      pval = pval
    ))
}

prevalence_results <- prevalence_results %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"),
         direction = case_when(
           prev_diff > 0 ~ paste0("Higher in ", levels(prevalence_clin$EthnicityTot)[2]),
           prev_diff < 0 ~ paste0("Higher in ", levels(prevalence_clin$EthnicityTot)[1]),
           TRUE ~ "Equal"))

write.csv2(prevalence_results, "results/5_arg/crosssectional/prevalence_ethnicity_comparison.csv", row.names = FALSE)

# Volcano plot
prevalence_results <- prevalence_results %>%
  mutate(sig_level = case_when(
    padj < 0.05 & abs(prev_diff) > 10 ~ "Significant & Large Effect",
    padj < 0.05 ~ "Significant",
    TRUE ~ "Not Significant"))

(p_prev_volcano <- ggplot(prevalence_results, aes(x = prev_diff, y = -log10(pval), color = sig_level)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(data = prevalence_results %>% filter(padj < 0.05),
                  aes(label = gene), size = 3, max.overlaps = 20,
                  box.padding = 0.5, point.padding = 0.3, color = "black") +
  geom_vline(xintercept = c(-10, 10), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Significant & Large Effect" = "red",
                                 "Significant" = "orange",
                                 "Not Significant" = "gray70")) +
  theme_Publication() +
  labs(x = paste0("Prevalence Difference (% in ", levels(prevalence_clin$EthnicityTot)[2],
                  " - % in ", levels(prevalence_clin$EthnicityTot)[1], ")"),
       y = "-log10(p-value)",
       title = "Differential Gene Prevalence Between Ethnicities (baseline)",
       color = ""))
ggsave("results/5_arg/crosssectional/prevalence_volcano_plot.pdf", p_prev_volcano, width = 7, height = 7)

# Top differential genes
prev_sig <- prevalence_results |> filter(pval < 0.05) |> arrange(-abs(prev_diff)) |> 
  mutate(gene = fct_reorder(gene, prev_diff))

prev_sig_long <- prev_sig |> 
  dplyr::select(gene, class, prev_group1, prev_group2) |> 
  pivot_longer(cols = c(prev_group1, prev_group2),
                names_to = "ethnicity", values_to = "prevalence") |> 
  mutate(ethnicity = ifelse(ethnicity == "prev_group1",
                            levels(prevalence_clin$EthnicityTot)[1],
                            levels(prevalence_clin$EthnicityTot)[2]))

(p_prev_top <- ggplot(prev_sig_long, aes(x = gene, y = prevalence, fill = ethnicity)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_jco() +
  theme_Publication() +
  labs(x = "", y = "Prevalence (%)",
        title = "Genes with Differential Prevalence", fill = ""))
ggsave("results/5_arg/crosssectional/prevalence_top_differences.pdf", p_prev_top, width = 10, height = 8)

# FOLLOW-UP GENE-LEVEL PREVALENCE ----
followup_samples_prev <- clinical %>%
  filter(sampleID %in% df_raw$Sample) |>
  filter(timepoint == "follow-up") %>%
  pull(sampleID) %>%
  unique()

prevalence_data_followup <- df_raw %>%
  filter(Gene_Symbol %in% prevalent_genes) %>%
  mutate(sampleID = Sample,
         present = ifelse(Prevalence > 0 & Mapped_Reads > 0, 1, 0)) %>%
  dplyr::select(sampleID, Gene_Symbol, present) %>%
  complete(sampleID = followup_samples_prev,
           Gene_Symbol = prevalent_genes,
           fill = list(present = 0))

prevalence_clin_followup <- prevalence_data_followup %>%
  left_join(clinical, by = "sampleID") %>%
  filter(timepoint == "follow-up", !is.na(EthnicityTot)) %>%
  droplevels()

prevalence_results_followup <- data.frame()
for(gene in prevalent_genes){
  gene_data <- prevalence_clin_followup %>% filter(Gene_Symbol == gene)
  if(nrow(gene_data) == 0) next

  prev_by_eth <- gene_data %>%
    group_by(EthnicityTot) %>%
    summarise(n_total = n(), n_present = sum(present),
              prevalence_pct = (n_present / n_total) * 100, .groups = "drop_last")

  if(sum(gene_data$present) == 0 | sum(gene_data$present) == nrow(gene_data)) next

  contingency <- table(gene_data$EthnicityTot, gene_data$present)

  test_result <- chisq.test(contingency)
  pval <- test_result$p.value
  or <- (contingency[1,2] * contingency[2,1]) / (contingency[1,1] * contingency[2,2])

  gene_info <- gene_prevalence %>%
    filter(Gene_Symbol == gene) %>%
    dplyr::select(Class, Subclass, prevalence_pct)

  prevalence_results_followup <- rbind(prevalence_results_followup, data.frame(
    gene = gene,
    class = ifelse(nrow(gene_info) > 0, gene_info$Class, NA),
    subclass = ifelse(nrow(gene_info) > 0, gene_info$Subclass, NA),
    overall_prevalence = ifelse(nrow(gene_info) > 0, gene_info$prevalence_pct, NA),
    prev_group1 = prev_by_eth$prevalence_pct[1],
    prev_group2 = prev_by_eth$prevalence_pct[2],
    prev_diff = prev_by_eth$prevalence_pct[2] - prev_by_eth$prevalence_pct[1],
    odds_ratio = or,
    pval = pval
  ))
}

prevalence_results_followup <- prevalence_results_followup %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"),
         direction = case_when(
           prev_diff > 0 ~ paste0("Higher in ", levels(prevalence_clin_followup$EthnicityTot)[2]),
           prev_diff < 0 ~ paste0("Higher in ", levels(prevalence_clin_followup$EthnicityTot)[1]),
           TRUE ~ "Equal"))

write.csv2(prevalence_results_followup, "results/5_arg/crosssectional/prevalence_ethnicity_comparison_followup.csv", row.names = FALSE)

prevalence_results_followup <- prevalence_results_followup %>%
  mutate(sig_level = case_when(
    padj < 0.05 & abs(prev_diff) > 10 ~ "Significant & Large Effect",
    padj < 0.05 ~ "Significant",
    TRUE ~ "Not Significant"))

(p_prev_volcano_followup <- ggplot(prevalence_results_followup, aes(x = prev_diff, y = -log10(pval), color = sig_level)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(data = prevalence_results_followup %>% filter(padj < 0.05),
                  aes(label = gene), size = 3, max.overlaps = 20,
                  box.padding = 0.5, point.padding = 0.3, color = "black") +
  geom_vline(xintercept = c(-10, 10), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Significant & Large Effect" = "red",
                                 "Significant" = "orange",
                                 "Not Significant" = "gray70")) +
  theme_Publication() +
  labs(x = paste0("Prevalence Difference (% in ", levels(prevalence_clin_followup$EthnicityTot)[2],
                  " - % in ", levels(prevalence_clin_followup$EthnicityTot)[1], ")"),
       y = "-log10(p-value)",
       title = "Differential Gene Prevalence Between Ethnicities (follow-up)",
       color = ""))
ggsave("results/5_arg/crosssectional/prevalence_volcano_plot_followup.pdf", p_prev_volcano_followup, width = 7, height = 7)

# Heatmap of top genes
top_genes <- prev_sig$gene[1:min(20, nrow(prev_sig))]
prev_matrix <- prevalence_clin %>%
  filter(Gene_Symbol %in% top_genes) %>%
  dplyr::select(sampleID, Gene_Symbol, present, EthnicityTot) %>%
  arrange(EthnicityTot, sampleID) %>%
  mutate(sample_eth = paste0(sampleID, " (", EthnicityTot, ")"))

gene_order <- prev_matrix %>%
  group_by(Gene_Symbol) %>%
  summarise(prev = mean(present), .groups = "drop_last") %>%
  arrange(-prev)

prev_matrix <- prev_matrix %>%
  mutate(Gene_Symbol = factor(Gene_Symbol, levels = gene_order$Gene_Symbol))

(p_prev_heatmap <- ggplot(prev_matrix, aes(x = sample_eth, y = Gene_Symbol, fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "darkblue"), labels = c("Absent", "Present")) +
  facet_grid(~EthnicityTot, scales = "free_x", space = "free_x") +
  theme_Publication() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing = unit(0.5, "lines")) +
  labs(x = "Samples", y = "Gene",
        title = "Presence/Absence Pattern of Top Differential Genes", fill = ""))
ggsave("results/5_arg/crosssectional/prevalence_heatmap.pdf", p_prev_heatmap, width = 12, height = 8)

# CLASS-LEVEL PREVALENCE ----
baseline_samples <- df_raw %>%
  mutate(sampleID = Sample) %>%
  inner_join(clinical %>% dplyr::select(sampleID, timepoint, EthnicityTot), by = "sampleID") %>%
  filter(timepoint == "baseline", !is.na(EthnicityTot)) %>%
  pull(sampleID) %>%
  unique()

class_prevalence_data <- df_raw %>%
  filter(Gene_Symbol %in% prevalent_genes) %>%
  mutate(sampleID = Sample) %>%
  filter(sampleID %in% baseline_samples) %>%
  mutate(present = ifelse(Prevalence > 0 & Mapped_Reads > 0, 1, 0)) %>%
  filter(present == 1) %>%
  group_by(sampleID, Class) %>%
  summarise(class_present = 1, .groups = "drop") %>%
  complete(sampleID = baseline_samples, Class, fill = list(class_present = 0))

class_prev_clin <- class_prevalence_data %>%
  left_join(clinical, by = "sampleID") %>%
  filter(!is.na(EthnicityTot)) %>%
  droplevels()

all_classes <- unique(class_prev_clin$Class)

# Test each class
class_prev_results <- data.frame()
for(arg_class in all_classes){
  class_data <- class_prev_clin %>% filter(Class == arg_class)
  prev_by_eth <- class_data %>% group_by(EthnicityTot) %>%
    summarise(n_total = n(), n_present = sum(class_present),
              prevalence_pct = (n_present / n_total) * 100, .groups = "drop_last")
  if(sum(class_data$class_present) == 0 | sum(class_data$class_present) == nrow(class_data)) next

  contingency <- table(class_data$EthnicityTot, class_data$class_present)

  tryCatch({
    test_result <- chisq.test(contingency)
    pval <- test_result$p.value
    or <- (contingency[1,2] * contingency[2,1]) / (contingency[1,1] * contingency[2,2])

    class_prev_results <- rbind(class_prev_results, data.frame(
      class = arg_class,
      prev_group1 = prev_by_eth$prevalence_pct[1],
      prev_group2 = prev_by_eth$prevalence_pct[2],
      prev_diff = prev_by_eth$prevalence_pct[2] - prev_by_eth$prevalence_pct[1],
      odds_ratio = or,
      pval = pval
    ))
  }, error = function(e) NULL)
}

class_prev_results <- class_prev_results %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"),
         direction = case_when(
           prev_diff > 0 ~ paste0("Higher in ", levels(class_prev_clin$EthnicityTot)[2]),
           prev_diff < 0 ~ paste0("Higher in ", levels(class_prev_clin$EthnicityTot)[1]),
           TRUE ~ "Equal"))
write.csv2(class_prev_results, "results/5_arg/crosssectional/class_prevalence_ethnicity_comparison.csv", row.names = FALSE)

# Class prevalence plots
class_prev_results <- class_prev_results %>%
  mutate(class = fct_reorder(class, prev_diff),
         sig_color = case_when(
           padj < 0.05 ~ "Significant (FDR < 0.05)",
           pval < 0.05 ~ "Nominally Significant (p < 0.05)",
           TRUE ~ "Not Significant"))

p_class_prev <- ggplot(class_prev_results, aes(x = class, y = prev_diff, fill = sig_color)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("Significant (FDR < 0.05)" = "red",
                                "Nominally Significant (p < 0.05)" = "orange",
                                "Not Significant" = "gray70")) +
  theme_Publication() +
  labs(x = "ARG Class",
       y = paste0("Prevalence Difference (% in ", levels(class_prev_clin$EthnicityTot)[2],
                  " - % in ", levels(class_prev_clin$EthnicityTot)[1], ")"),
       title = "ARG Class Prevalence Differences Between Ethnicities", fill = "")
ggsave("results/5_arg/crosssectional/class_prevalence_differences.pdf", p_class_prev, width = 10, height = 8)

class_prev_long <- class_prev_results %>%
  dplyr::select(class, prev_group1, prev_group2) %>%
  pivot_longer(cols = c(prev_group1, prev_group2),
               names_to = "ethnicity", values_to = "prevalence") %>%
  mutate(ethnicity = ifelse(ethnicity == "prev_group1",
                            levels(class_prev_clin$EthnicityTot)[1],
                            levels(class_prev_clin$EthnicityTot)[2]),
         class = fct_reorder(class, prevalence, .fun = mean))

p_class_grouped <- ggplot(class_prev_long, aes(x = class, y = prevalence, fill = ethnicity)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_jco() +
  theme_Publication() +
  labs(x = "ARG Class", y = "Prevalence (% of samples with at least one gene)",
       title = "ARG Class Prevalence by Ethnicity", fill = "")
ggsave("results/5_arg/crosssectional/class_prevalence_grouped.pdf", p_class_grouped, width = 10, height = 8)

# Baseline heatmap
class_matrix <- class_prev_clin %>%
  dplyr::select(sampleID, Class, class_present, EthnicityTot) %>%
  arrange(EthnicityTot, sampleID)

class_order <- class_matrix %>%
  group_by(Class) %>%
  summarise(prev = mean(class_present), .groups = "drop_last") %>%
  arrange(-prev)

class_matrix <- class_matrix %>%
  mutate(Class = factor(Class, levels = class_order$Class))

p_class_heatmap <- ggplot(class_matrix, aes(x = sampleID, y = Class, fill = factor(class_present))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "darkred"),
                    labels = c("Absent", "Present")) +
  facet_grid(~EthnicityTot, scales = "free_x", space = "free_x") +
  theme_Publication() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing = unit(0.5, "lines")) +
  labs(x = "Samples (grouped by ethnicity)", y = "ARG Class",
       title = "ARG Class Presence/Absence Pattern by Ethnicity (Baseline)", fill = "")
ggsave("results/5_arg/crosssectional/class_prevalence_heatmap_baseline.pdf", p_class_heatmap, width = 14, height = 8)

# Follow-up heatmap
followup_samples <- df_raw %>%
  mutate(sampleID = Sample) %>%
  inner_join(clinical %>% dplyr::select(sampleID, timepoint, EthnicityTot), by = "sampleID") %>%
  filter(timepoint == "follow-up", !is.na(EthnicityTot)) %>%
  pull(sampleID) %>%
  unique()

class_prevalence_followup <- df_raw %>%
  filter(Gene_Symbol %in% prevalent_genes) %>%
  mutate(sampleID = Sample) %>%
  filter(sampleID %in% followup_samples) %>%
  mutate(present = ifelse(Prevalence > 0 & Mapped_Reads > 0, 1, 0)) %>%
  filter(present == 1) %>%
  group_by(sampleID, Class) %>%
  summarise(class_present = 1, .groups = "drop") %>%
  complete(sampleID = followup_samples, Class, fill = list(class_present = 0))

class_prev_followup <- class_prevalence_followup %>%
  left_join(clinical, by = "sampleID") %>%
  filter(!is.na(EthnicityTot)) %>%
  droplevels()

class_matrix_followup <- class_prev_followup %>%
  dplyr::select(sampleID, Class, class_present, EthnicityTot) %>%
  arrange(EthnicityTot, sampleID) %>%
  mutate(Class = factor(Class, levels = class_order$Class))

p_class_heatmap_followup <- ggplot(class_matrix_followup, aes(x = sampleID, y = Class, fill = factor(class_present))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "darkred"),
                    labels = c("Absent", "Present")) +
  facet_grid(~EthnicityTot, scales = "free_x", space = "free_x") +
  theme_Publication() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing = unit(0.5, "lines")) +
  labs(x = "Samples (grouped by ethnicity)", y = "ARG Class",
       title = "ARG Class Presence/Absence Pattern by Ethnicity (Follow-up)", fill = "")
ggsave("results/5_arg/crosssectional/class_prevalence_heatmap_followup.pdf", p_class_heatmap_followup, width = 14, height = 8)

# BASELINE ABUNDANCE DIFFERENCES ----
df_wide <- df_raw %>%
  filter(Gene_Symbol %in% prevalent_genes) %>%
  dplyr::select(Sample, Gene_Symbol, RPKM) %>%
  group_by(Sample, Gene_Symbol) %>%
  summarise(RPKM = mean(RPKM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Gene_Symbol, values_from = RPKM,
              values_fill = list(RPKM = 0)) %>%
  mutate(sampleID = Sample) %>%
  dplyr::select(sampleID, everything(), -Sample)

df_tot <- df_wide %>%
  left_join(clinical, by = "sampleID") %>%
  filter(!is.na(EthnicityTot)) %>%
  droplevels()

df_baseline <- df_tot %>% filter(timepoint == "baseline")
gene_cols <- colnames(df_baseline)[2:(length(prevalent_genes) + 1)]

statres_baseline <- data.frame()
for(a in 1:length(gene_cols)){
  mbname <- gene_cols[a]
  df_baseline$mb <- log10(df_baseline[[mbname]] + 1)

  tryCatch({
    model_baseline <- lm(mb ~ EthnicityTot, data = df_baseline)
    res <- summary(model_baseline)

    if(nrow(res$coefficients) >= 2){
      confint_baseline <- confint(model_baseline)
      gene_info <- gene_prevalence %>%
        filter(Gene_Symbol == mbname) %>%
        dplyr::select(Class, Subclass, prevalence_pct)

      statres_baseline <- rbind(statres_baseline, data.frame(
        mbname = mbname,
        class = ifelse(nrow(gene_info) > 0, gene_info$Class, NA),
        subclass = ifelse(nrow(gene_info) > 0, gene_info$Subclass, NA),
        prevalence_pct = ifelse(nrow(gene_info) > 0, gene_info$prevalence_pct, NA),
        pval = res$coefficients[2, 4],
        estimate = res$coefficients[2, 1],
        se = res$coefficients[2, 2],
        conflow = confint_baseline[2, 1],
        confhigh = confint_baseline[2, 2]
      ))
    }
  }, error = function(e) NULL)
}

statres_baseline <- statres_baseline %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))

write.csv2(statres_baseline, "results/5_arg/crosssectional/baseline_ethnicity_abundance.csv", row.names = FALSE)

# Volcano plot for abundance
baseline_sig <- statres_baseline %>%
  mutate(sig_level = case_when(
    padj < 0.05 & abs(estimate) > 0.5 ~ "Significant & Large Effect",
    padj < 0.05 ~ "Significant",
    TRUE ~ "Not Significant"))

p_baseline_volcano <- ggplot(baseline_sig, aes(x = estimate, y = -log10(pval), color = sig_level)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(data = baseline_sig %>% filter(padj < 0.05),
                  aes(label = mbname), size = 3, max.overlaps = 20,
                  box.padding = 0.5, point.padding = 0.3, color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Significant & Large Effect" = "red",
                                 "Significant" = "orange",
                                 "Not Significant" = "gray70")) +
  theme_Publication() +
  labs(x = paste0("log10(RPKM) Difference (", levels(df_baseline$EthnicityTot)[2],
                  " - ", levels(df_baseline$EthnicityTot)[1], ")"),
       y = "-log10(p-value)",
       title = "Baseline ARG Abundance Differences Between Ethnicities", color = "")
ggsave("results/5_arg/crosssectional/baseline_abundance_volcano.pdf", p_baseline_volcano, width = 7, height = 7)

# FOLLOW-UP ABUNDANCE DIFFERENCES ----
df_followup <- df_tot %>% filter(timepoint == "follow-up")

statres_followup <- data.frame()
for(a in 1:length(gene_cols)){
  mbname <- gene_cols[a]
  df_followup$mb <- log10(df_followup[[mbname]] + 1)

  tryCatch({
    model_followup <- lm(mb ~ EthnicityTot, data = df_followup)
    res <- summary(model_followup)

    if(nrow(res$coefficients) >= 2){
      confint_followup <- confint(model_followup)
      gene_info <- gene_prevalence %>%
        filter(Gene_Symbol == mbname) %>%
        dplyr::select(Class, Subclass, prevalence_pct)

      statres_followup <- rbind(statres_followup, data.frame(
        mbname = mbname,
        class = ifelse(nrow(gene_info) > 0, gene_info$Class, NA),
        subclass = ifelse(nrow(gene_info) > 0, gene_info$Subclass, NA),
        prevalence_pct = ifelse(nrow(gene_info) > 0, gene_info$prevalence_pct, NA),
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

write.csv2(statres_followup, "results/5_arg/crosssectional/followup_ethnicity_abundance.csv", row.names = FALSE)

# Volcano plot for follow-up abundance
followup_sig <- statres_followup %>%
  mutate(sig_level = case_when(
    padj < 0.05 & abs(estimate) > 0.5 ~ "Significant & Large Effect",
    padj < 0.05 ~ "Significant",
    TRUE ~ "Not Significant"))

p_followup_volcano <- ggplot(followup_sig, aes(x = estimate, y = -log10(pval), color = sig_level)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(data = followup_sig %>% filter(padj < 0.05),
                  aes(label = mbname), size = 3, max.overlaps = 20,
                  box.padding = 0.5, point.padding = 0.3, color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Significant & Large Effect" = "red",
                                 "Significant" = "orange",
                                 "Not Significant" = "gray70")) +
  theme_Publication() +
  labs(x = paste0("log10(RPKM) Difference (", levels(df_followup$EthnicityTot)[2],
                  " - ", levels(df_followup$EthnicityTot)[1], ")"),
       y = "-log10(p-value)",
       title = "Follow-up ARG Abundance Differences Between Ethnicities", color = "")
ggsave("results/5_arg/crosssectional/followup_abundance_volcano.pdf", p_followup_volcano, width = 7, height = 7)
