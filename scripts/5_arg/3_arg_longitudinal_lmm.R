## ARG Longitudinal Analysis with Linear Mixed Models
library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(lmerTest)

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
           strip.text = element_text(face="bold"),
           plot.caption = element_text(size = rel(0.5), face = "italic"),
           plot.subtitle = element_text(size=8, hjust = 0.5, face = "italic")))
}

df_raw <- rio::import("data/shotgun/arg/all_samples.merged_arg_counts.tsv")
clinical <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
dir.create("results/5_arg/longitudinal", showWarnings = FALSE, recursive = TRUE)

# Get prevalent genes (>5%)
gene_prevalence <- df_raw %>%
  filter(Prevalence > 0) %>%
  mutate(Class = str_to_title(Class),
         Subclass = str_to_title(Subclass)) %>%
  group_by(Gene_Symbol, Class, Subclass) %>%
  summarise(n_samples = n_distinct(Sample),
            prevalence_pct = (n_distinct(Sample) / n_distinct(df_raw$Sample)) * 100,
            .groups = "drop")

prevalent_genes <- gene_prevalence %>%
  filter(prevalence_pct > 0.5) %>%
  pull(Gene_Symbol)

# DATA PREPARATION ----
df_wide <- df_raw %>%
  filter(Gene_Symbol %in% prevalent_genes) %>%
  dplyr::select(Sample, Gene_Symbol, RPKM) %>%
  group_by(Sample, Gene_Symbol) %>%
  summarise(RPKM = mean(RPKM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Gene_Symbol, values_from = RPKM,
              values_fill = list(RPKM = 0)) %>%
  mutate(sampleID = Sample) %>%
  dplyr::select(sampleID, everything(), -Sample)

# Add sequencing depth
sample_depth <- df_raw %>%
  group_by(Sample) %>%
  summarise(sequencing_depth = mean(Total_Reads, na.rm = TRUE), .groups = "drop") %>%
  rename(sampleID = Sample)

df_tot <- left_join(df_wide, clinical, by = "sampleID") %>%
  left_join(sample_depth, by = "sampleID") %>%
  filter(!is.na(EthnicityTot)) %>%
  mutate(timepoint = factor(timepoint, levels = c("baseline", "follow-up")),
         EthnicityTot = factor(EthnicityTot),
         ID = factor(ID),
         log_depth = log10(sequencing_depth)) %>%
  droplevels()

saveRDS(df_tot, "data/arg_prepared_for_lmm.RDS")
write.csv2(data.frame(gene = prevalent_genes), "results/5_arg/longitudinal/prevalent_genes_list.csv", row.names = FALSE)

# TOTAL ARG BURDEN ----
# Using RPM (Reads Per Million) instead of RPKM to avoid gene length bias
total_arg_burden <- df_raw %>%
  group_by(Sample) %>%
  summarise(total_arg_reads = sum(Mapped_Reads, na.rm = TRUE),
            n_arg_genes = sum(Prevalence > 0),
            sequencing_depth = mean(Total_Reads, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(total_arg_rpm = (total_arg_reads / sequencing_depth) * 1e6) %>%
  rename(sampleID = Sample)

arg_burden_clin <- left_join(total_arg_burden, clinical, by = "sampleID") %>%
  filter(!is.na(EthnicityTot)) %>%
  mutate(timepoint = factor(timepoint, levels = c("baseline", "follow-up")),
         EthnicityTot = factor(EthnicityTot),
         ID = factor(ID),
         log_rpm = log10(total_arg_rpm + 1)) %>%
  droplevels()

# Tests (no need for depth adjustment - RPM already normalizes for depth)
arg_baseline <- arg_burden_clin %>% filter(timepoint == "baseline")
model_eth <- lm(log_rpm ~ EthnicityTot, data = arg_baseline)
model_time <- lmer(log_rpm ~ timepoint + (1|ID), data = arg_burden_clin)
model_int <- lmer(log_rpm ~ EthnicityTot * timepoint + (1|ID), data = arg_burden_clin)

res_eth <- summary(model_eth)
res_time <- summary(model_time)
res_int <- summary(model_int)

arg_burden_results <- data.frame(
  test = c("Ethnicity (baseline)", "Timepoint (overall)", "Ethnicity*Timepoint"),
  estimate = c(res_eth$coefficients[2, 1], res_time$coefficients[2, 1], res_int$coefficients[4, 1]),
  se = c(res_eth$coefficients[2, 2], res_time$coefficients[2, 2], res_int$coefficients[4, 2]),
  pval = c(res_eth$coefficients[2, 4], res_time$coefficients[2, 5], res_int$coefficients[4, 5]))
write.csv2(arg_burden_results, "results/5_arg/longitudinal/total_arg_burden_results.csv", row.names = FALSE)

# Plots
p_burden_eth <- ggplot(arg_baseline, aes(x = EthnicityTot, y = total_arg_rpm, fill = EthnicityTot)) +
  geom_violin(colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA) +
  geom_jitter(alpha = 0.2, width = 0.2, size = 1) +
  stat_compare_means(method = "wilcox.test") +
  scale_fill_jco() +
  scale_y_log10() +
  theme_Publication() +
  labs(x = "", y = "Total ARG Burden (RPM, log scale)",
       title = "Total ARG Burden by Ethnicity at Baseline") +
  theme(legend.position = "none")
ggsave("results/5_arg/longitudinal/total_burden_ethnicity_baseline.pdf", p_burden_eth, width = 6, height = 6)

arg_change <- arg_burden_clin %>%
  dplyr::select(ID, timepoint, total_arg_rpm, EthnicityTot) %>%
  pivot_wider(names_from = timepoint, values_from = total_arg_rpm) %>%
  mutate(change = `follow-up` - baseline)

p_burden_time <- ggplot(arg_burden_clin, aes(x = timepoint, y = total_arg_rpm, fill = timepoint)) +
  geom_violin(colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA) +
  geom_line(aes(group = ID), alpha = 0.1, color = "gray50") +
  scale_fill_simpsons() +
  scale_y_log10() +
  annotate("text", x = 1.5, y = max(arg_burden_clin$total_arg_rpm, na.rm = TRUE),
           label = paste0("Paired t-test p = ", formatC(res_time$coefficients[2, 5], format = "e", digits = 2)),
           size = 4) +
  theme_Publication() +
  labs(x = "", y = "Total ARG Burden (RPM, log scale)",
       title = "Total ARG Burden Change Over Time") +
  theme(legend.position = "none")
ggsave("results/5_arg/longitudinal/total_burden_timepoint.pdf", p_burden_time, width = 6, height = 6)

means_int <- arg_burden_clin %>%
  group_by(EthnicityTot, timepoint) %>%
  summarise(mean_rpm = mean(total_arg_rpm), sd_rpm = sd(total_arg_rpm),
            n = n(), se = sd_rpm / sqrt(n), .groups = "drop")

p_burden_int <- ggplot() +
  geom_line(data = arg_burden_clin,
            aes(x = timepoint, y = total_arg_rpm, color = EthnicityTot, group = ID),
            alpha = 0.15) +
  geom_line(data = means_int,
            aes(x = timepoint, y = mean_rpm, color = EthnicityTot, group = EthnicityTot),
            linewidth = 1.5) +
  geom_point(data = means_int,
             aes(x = timepoint, y = mean_rpm, color = EthnicityTot), size = 3) +
  geom_errorbar(data = means_int,
                aes(x = timepoint, ymin = mean_rpm - se, ymax = mean_rpm + se,
                    color = EthnicityTot), width = 0.1, linewidth = 1) +
  scale_color_jco() +
  scale_y_log10() +
  theme_Publication() +
  labs(x = "Timepoint", y = "Total ARG Burden (RPM, log scale)",
       title = "Total ARG Burden: Ethnicity x Timepoint", color = "",
       caption = paste0("Interaction p = ", formatC(res_int$coefficients[4, 5], format = "e", digits = 2)))
ggsave("results/5_arg/longitudinal/total_burden_interaction.pdf", p_burden_int, width = 8, height = 6)

p_change <- ggplot(arg_change, aes(x = EthnicityTot, y = change, fill = EthnicityTot)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_violin(colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA) +
  geom_jitter(alpha = 0.2, width = 0.2, size = 1) +
  stat_compare_means(method = "t.test") +
  scale_fill_jco() +
  theme_Publication() +
  labs(x = "", y = "Change in Total ARG Burden\n(Follow-up - Baseline RPM)",
       title = "Change in Total ARG Burden by Ethnicity") +
  theme(legend.position = "none")
ggsave("results/5_arg/longitudinal/total_burden_change_by_ethnicity.pdf", p_change, width = 6, height = 6)

# ARG DIVERSITY ----
# Richness = number of distinct ARG genes detected per sample
# Shannon  = evenness-weighted diversity from RPKM proportions
arg_div_raw <- df_raw %>%
  filter(Prevalence > 0, Mapped_Reads > 0) %>%
  group_by(Sample, Gene_Symbol) %>%
  summarise(RPKM = sum(RPKM, na.rm = TRUE), .groups = "drop")

arg_diversity <- arg_div_raw %>%
  group_by(Sample) %>%
  summarise(
    richness = n_distinct(Gene_Symbol),
    shannon  = {
      p <- RPKM / sum(RPKM)
      p <- p[p > 0]
      -sum(p * log(p))
    },
    .groups = "drop"
  ) %>%
  rename(sampleID = Sample)

arg_div_clin <- arg_diversity %>%
  left_join(clinical, by = "sampleID") %>%
  filter(!is.na(EthnicityTot)) %>%
  mutate(timepoint   = factor(timepoint, levels = c("baseline", "follow-up")),
         EthnicityTot = factor(EthnicityTot),
         ID           = factor(ID)) %>%
  droplevels()

# LMMs: ethnicity × timepoint interaction
model_rich_int <- lmer(richness ~ EthnicityTot * timepoint + (1|ID), data = arg_div_clin)
model_shan_int <- lmer(shannon  ~ EthnicityTot * timepoint + (1|ID), data = arg_div_clin)

res_rich <- summary(model_rich_int)
res_shan <- summary(model_shan_int)

div_results <- data.frame(
  metric   = c("Richness", "Shannon"),
  estimate = c(res_rich$coefficients[4, 1], res_shan$coefficients[4, 1]),
  se       = c(res_rich$coefficients[4, 2], res_shan$coefficients[4, 2]),
  pval     = c(res_rich$coefficients[4, 5], res_shan$coefficients[4, 5])
)
write.csv2(div_results, "results/5_arg/longitudinal/arg_diversity_lmm_results.csv", row.names = FALSE)

# Exploratory plots
arg_div_plot <- arg_div_clin %>%
  mutate(tp_label = factor(recode(as.character(timepoint),
                                  "baseline"  = "Baseline",
                                  "follow-up" = "Follow-up"),
                           levels = c("Baseline", "Follow-up")))

p_richness <- ggplot(arg_div_plot, aes(x = EthnicityTot, y = richness, fill = EthnicityTot)) +
  geom_violin(alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.22, fill = "white", outlier.shape = NA, colour = "gray30") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5, label.y.npc = 0.97, size = 2.8, hjust = 0.5) +
  facet_wrap(~tp_label) +
  scale_fill_jco() +
  scale_x_discrete(labels = function(x) gsub("South-Asian Surinamese", "South-Asian\nSurinamese", x)) +
  theme_Publication() +
  labs(x = "", y = "ARG Richness (n genes)", title = "ARG Richness by Ethnicity") +
  theme(legend.position = "none")
ggsave("results/5_arg/longitudinal/arg_richness_ethnicity.pdf", p_richness, width = 8, height = 5)

p_shannon <- ggplot(arg_div_plot, aes(x = EthnicityTot, y = shannon, fill = EthnicityTot)) +
  geom_violin(alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.22, fill = "white", outlier.shape = NA, colour = "gray30") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5, label.y.npc = 0.97, size = 2.8, hjust = 0.5) +
  facet_wrap(~tp_label) +
  scale_fill_jco() +
  scale_x_discrete(labels = function(x) gsub("South-Asian Surinamese", "South-Asian\nSurinamese", x)) +
  theme_Publication() +
  labs(x = "", y = "Shannon Diversity", title = "ARG Shannon diversity by ethnicity") +
  theme(legend.position = "none")
ggsave("results/5_arg/longitudinal/arg_shannon_ethnicity.pdf", p_shannon, width = 8, height = 5)

# GENE-LEVEL LMM ----
gene_cols <- colnames(df_tot)[2:(length(prevalent_genes) + 1)]
statres <- data.frame()

for(a in 1:length(gene_cols)){
  mbname <- gene_cols[a]
  df_tot$mb <- log10(df_tot[[mbname]] + 1)
  model1 <- lmer(mb ~ EthnicityTot * timepoint + log_depth + (1|ID), data = df_tot)
  res <- summary(model1)
  confint_model1 <- confint(model1, method = "Wald")
  interaction_row <- grep("EthnicityTot.*:.*timepoint", rownames(res$coefficients))

  if(length(interaction_row) > 0){
    ci_row <- grep("EthnicityTot.*:.*timepoint", rownames(confint_model1))
    gene_info <- gene_prevalence %>%
      filter(Gene_Symbol == mbname) %>%
      dplyr::select(Class, Subclass, prevalence_pct)

    statres <- rbind(statres, data.frame(
      mbname = mbname,
      class = ifelse(nrow(gene_info) > 0, gene_info$Class, NA),
      subclass = ifelse(nrow(gene_info) > 0, gene_info$Subclass, NA),
      prevalence_pct = ifelse(nrow(gene_info) > 0, gene_info$prevalence_pct, NA),
      pval = res$coefficients[interaction_row, 5],
      estimate = res$coefficients[interaction_row, 1],
      conflow = ifelse(length(ci_row) > 0, confint_model1[ci_row, 1], NA),
      confhigh = ifelse(length(ci_row) > 0, confint_model1[ci_row, 2], NA)
    ))
  }
}

statres <- statres %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv2(statres, "results/5_arg/longitudinal/lmm_ethnicity_timepoint_results.csv", row.names = FALSE)

statres_interaction <- statres   # keep before timepoint-only LMM overwrites statres

# LONGITUDINAL PLOTS ----
statres_sig <- statres %>%
  filter(pval < 0.05) %>%
  arrange(pval)

if(nrow(statres_sig) > 0){
  plist <- list()

  for(i in 1:min(nrow(statres_sig), 20)){
    nm <- statres_sig$mbname[i]
    df_tot$mb <- log10(df_tot[[nm]] + 1)

    df_means <- df_tot %>%
      group_by(EthnicityTot, timepoint) %>%
      summarise(mean = mean(mb, na.rm = TRUE), sd = sd(mb, na.rm = TRUE),
                n = n(), .groups = "drop")

    res_lmm <- statres_sig %>%
      filter(mbname == nm) %>%
      mutate(group1 = "baseline", group2 = "follow-up",
             sigq = paste0("p=", formatC(pval, format = "e", digits = 2))) %>%
      dplyr::select(-mbname)

    mbmax <- ifelse(max(df_tot$mb, na.rm = TRUE) < 0,
                    max(df_tot$mb, na.rm = TRUE) * 0.8,
                    max(df_tot$mb, na.rm = TRUE) * 1.2)
    mbstat <- ifelse(max(df_tot$mb, na.rm = TRUE) < 0,
                     max(df_tot$mb, na.rm = TRUE) * 0.9,
                     max(df_tot$mb, na.rm = TRUE) * 1.1)
    mbmin <- min(df_tot$mb, na.rm = TRUE)

    pl2 <- ggplot() +
      geom_line(data = df_tot, aes(x = timepoint, y = mb, color = EthnicityTot, group = ID),
                alpha = 0.15, linewidth = 0.5) +
      geom_point(data = df_tot, aes(x = timepoint, y = mb, color = EthnicityTot),
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
           title = paste0(nm, "\n", statres_sig$subclass[i]),
           color = "")

    plist[[i]] <- pl2
  }

  n_plots <- length(plist)
  n_cols <- 3
  n_rows <- ceiling(n_plots / n_cols)

  plots <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
                     labels = LETTERS[1:n_plots], nrow = n_rows, ncol = n_cols)

  ggsave("results/5_arg/longitudinal/significant_arg_timepoint_ethnicity.pdf", plots,
         width = 12, height = 4 * n_rows)
}


# GENE-LEVEL timepoint ----
gene_cols <- colnames(df_tot)[2:(length(prevalent_genes) + 1)]
statres <- data.frame()
for(a in 1:length(gene_cols)){
  mbname <- gene_cols[a]
  df_tot$mb <- log10(df_tot[[mbname]] + 1)
  model1 <- lmer(mb ~ timepoint + log_depth + (1|ID), data = df_tot)
  res <- summary(model1)
  confint_model1 <- confint(model1, method = "Wald")
  pval_row <- grep("timepoint", rownames(res$coefficients))

  if(length(pval_row) > 0){
    ci_row <- grep("timepoint", rownames(confint_model1))
    gene_info <- gene_prevalence %>%
      filter(Gene_Symbol == mbname) %>%
      dplyr::select(Class, Subclass, prevalence_pct)

    statres <- rbind(statres, data.frame(
      mbname = mbname,
      class = ifelse(nrow(gene_info) > 0, gene_info$Class, NA),
      subclass = ifelse(nrow(gene_info) > 0, gene_info$Subclass, NA),
      prevalence_pct = ifelse(nrow(gene_info) > 0, gene_info$prevalence_pct, NA),
      pval = res$coefficients[pval_row, 5],
      estimate = res$coefficients[pval_row, 1],
      conflow = ifelse(length(ci_row) > 0, confint_model1[ci_row, 1], NA),
      confhigh = ifelse(length(ci_row) > 0, confint_model1[ci_row, 2], NA)
    ))
  }
}

statres <- statres %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv2(statres, "results/5_arg/longitudinal/lmm_timepoint_results.csv", row.names = FALSE)

# LONGITUDINAL PLOTS ----
statres_sig <- statres |> filter(pval < 0.05) |> arrange(pval)

plist <- list()
for(i in 1:min(nrow(statres_sig), 20)){
  nm <- statres_sig$mbname[i]
  df_tot$mb <- log10(df_tot[[nm]] + 1)

  df_means <- df_tot %>%
    group_by(timepoint) %>%
    summarise(mean = mean(mb, na.rm = TRUE), sd = sd(mb, na.rm = TRUE),
              n = n(), .groups = "drop")

  res_lmm <- statres_sig %>%
    filter(mbname == nm) %>%
    mutate(group1 = "baseline", group2 = "follow-up",
            sigq = paste0("p=", formatC(pval, format = "e", digits = 2))) %>%
    dplyr::select(-mbname)

  mbmax <- ifelse(max(df_tot$mb, na.rm = TRUE) < 0,
                  max(df_tot$mb, na.rm = TRUE) * 0.8,
                  max(df_tot$mb, na.rm = TRUE) * 1.2)
  mbstat <- ifelse(max(df_tot$mb, na.rm = TRUE) < 0,
                    max(df_tot$mb, na.rm = TRUE) * 0.9,
                    max(df_tot$mb, na.rm = TRUE) * 1.1)
  mbmin <- min(df_tot$mb, na.rm = TRUE)

  pl2 <- ggplot() +
    geom_line(data = df_tot, aes(x = timepoint, y = mb, group = ID),
              alpha = 0.15, linewidth = 0.5, color = pal_nejm()(6)[6]) +
    geom_point(data = df_tot, aes(x = timepoint, y = mb),
                alpha = 0.15, size = 0.8, color = pal_nejm()(6)[6]) +
    geom_line(data = df_means, aes(x = timepoint, y = mean, group = 1),
              alpha = 1, linewidth = 1.0, color = pal_nejm()(3)[3]) +
    geom_point(data = df_means, aes(x = timepoint, y = mean),
                alpha = 1, size = 1.3, color = pal_nejm()(3)[3]) +
    geom_errorbar(data = df_means,
                  aes(ymin = mean - (sd/sqrt(n)), ymax = mean + (sd/sqrt(n)),
                      x = timepoint), width = 0.1, color = pal_nejm()(3)[3]) +
    stat_pvalue_manual(res_lmm, y.position = mbstat, label = "sigq",
                        tip.length = 0, bracket.shorten = 0.1, size = 4) +
    coord_cartesian(ylim = c(mbmin, mbmax)) +
    theme_Publication() +
    labs(x = "Timepoint", y = "log10(RPKM + 1)",
          title = paste0(nm), subtitle = paste0(statres_sig$subclass[i]),
          color = "")

  plist[[i]] <- pl2
}

n_plots <- length(plist)
n_cols <- 3
n_rows <- ceiling(n_plots / n_cols)

plots <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
                    labels = LETTERS[1:n_plots], nrow = n_rows, ncol = n_cols)

ggsave("results/5_arg/longitudinal/significant_arg_timepoint.pdf", plots,
        width = 12, height = 4 * n_rows)

## ── Figure 5 panels: shared aesthetics ───────────────────────────────────────
BASE_SIZE  <- 11
tp_colors  <- pal_simpsons()(8)[7:8]
names(tp_colors) <- c("baseline", "follow-up")
ETH_DUTCH  <- levels(arg_burden_clin$EthnicityTot)[1]
ETH_SAS    <- levels(arg_burden_clin$EthnicityTot)[2]
eth_colors <- c("#2166AC", "#E6B800")
names(eth_colors) <- c(ETH_DUTCH, ETH_SAS)
tp_labels  <- c("baseline" = "Baseline", "follow-up" = "Follow-up")

## ── Panel A: Total ARG Burden Over Time ──────────────────────────────────────
pl_A <- ggplot(arg_burden_clin,
               aes(x = timepoint, y = log_rpm, fill = timepoint,
                   alpha = timepoint)) +
  geom_violin(colour = NA) +
  geom_boxplot(width = 0.22, fill = "white", outlier.shape = NA, colour = "gray30",
               alpha = 1) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.8, hjust = 0.5,
           label = "p = 5.4e-15", size = 3.2) +
  scale_fill_manual(values = tp_colors, guide = "none") +
  scale_alpha_manual(values = c("baseline" = 0.60, "follow-up" = 0.90),
                     guide = "none") +
  scale_x_discrete(labels = tp_labels) +
  theme_Publication(base_size = BASE_SIZE) +
  labs(x = "", y = "Total ARG Burden (log\u2081\u2080 RPM)",
       title = "Total ARG burden over time")

## ── Panel B: ARG Burden by Ethnicity × Timepoint ─────────────────────────────
arg_burden_B <- arg_burden_clin %>%
  mutate(tp_label = factor(recode(as.character(timepoint),
                                  "baseline"  = "Baseline",
                                  "follow-up" = "Follow-up"),
                           levels = c("Baseline", "Follow-up")))

pl_B <- ggplot(arg_burden_B,
               aes(x = EthnicityTot, y = log_rpm, fill = EthnicityTot)) +
  geom_violin(alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.22, fill = "white", outlier.shape = NA, colour = "gray30") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5, label.y.npc = 0.97,
                     size = 2.8, hjust = 0.5) +
  facet_wrap(~tp_label) +
  scale_fill_manual(values = eth_colors, guide = "none") +
  scale_x_discrete(labels = function(x)
    gsub("South-Asian Surinamese", "South-Asian\nSurinamese", x)) +
  theme_Publication(base_size = BASE_SIZE) +
  labs(x = "", y = "Total ARG Burden (log\u2081\u2080 RPM)",
       title = "ARG burden by ethnicity")

## ── Panel B_rich / B_shan: ARG Diversity by Ethnicity × Timepoint ────────────
arg_div_fig <- arg_div_clin %>%
  mutate(tp_label = factor(recode(as.character(timepoint),
                                  "baseline"  = "Baseline",
                                  "follow-up" = "Follow-up"),
                           levels = c("Baseline", "Follow-up")))

pl_B_rich <- ggplot(arg_div_fig, aes(x = EthnicityTot, y = richness, fill = EthnicityTot)) +
  geom_violin(alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.22, fill = "white", outlier.shape = NA, colour = "gray30") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5, label.y.npc = 0.97, size = 2.8, hjust = 0.5) +
  facet_wrap(~tp_label) +
  scale_fill_manual(values = eth_colors, guide = "none") +
  scale_x_discrete(labels = function(x) gsub("South-Asian Surinamese", "South-Asian\nSurinamese", x)) +
  theme_Publication(base_size = BASE_SIZE) +
  labs(x = "", y = "ARG Richness (n genes)", title = "ARG Richness by Ethnicity")

pl_B_shan <- ggplot(arg_div_fig, aes(x = EthnicityTot, y = shannon, fill = EthnicityTot)) +
  geom_violin(alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.22, fill = "white", outlier.shape = NA, colour = "gray30") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.5, label.y.npc = 0.97, size = 2.8, hjust = 0.5) +
  facet_wrap(~tp_label) +
  scale_fill_manual(values = eth_colors, guide = "none") +
  scale_x_discrete(labels = function(x) gsub("South-Asian Surinamese", "South-Asian\nSurinamese", x)) +
  theme_Publication(base_size = BASE_SIZE) +
  labs(x = "", y = "Shannon Diversity", title = "ARG Shannon diversity by ethnicity")

## ── Panels F–: Key gene box + violin (FDR < 0.05 interaction) ────────────────
key_genes <- statres_interaction %>%
  filter(as.numeric(padj) < 0.05) %>%
  arrange(as.numeric(pval)) %>%
  pull(mbname)

gene_panels <- lapply(key_genes, function(nm) {
  gr      <- statres_interaction %>% filter(mbname == nm)
  sub_cls <- gr$subclass[1]
  int_p   <- as.numeric(gr$pval[1])
  p_str   <- if (int_p < 0.001) sprintf("interaction p = %.2e", int_p) else
                                 sprintf("interaction p = %.3f", int_p)

  plot_dat <- df_tot %>%
    mutate(mb        = log10(.data[[nm]] + 1),
           timepoint = factor(timepoint, levels = c("baseline", "follow-up"),
                              labels = c("Baseline", "Follow-up")))

  ggplot(plot_dat, aes(x = EthnicityTot, y = mb, fill = EthnicityTot)) +
    geom_violin(alpha = 0.75, colour = NA) +
    geom_boxplot(width = 0.20, fill = "white", outlier.shape = NA, colour = "gray30") +
    facet_wrap(~timepoint) +
    scale_fill_manual(values = eth_colors, guide = "none") +
    scale_x_discrete(labels = function(x)
      gsub("South-Asian Surinamese", "South-Asian\nSurinamese", x)) +
    theme_Publication(base_size = BASE_SIZE) +
    labs(x = "", y = "log\u2081\u2080(RPKM + 1)",
         title    = paste0(nm, "\n", sub_cls),
         subtitle = p_str)
})

