## Figure 5 — Antimicrobial resistance gene (ARG) dynamics
##
## Panels:
##   A — Total ARG burden over time (paired violin + spaghetti, lines coloured by ethnicity)
##   B — Total ARG burden by ethnicity × timepoint (2×2 violin; Wilcoxon p per timepoint)
##   C — Prevalence volcano: baseline gene prevalence differences Dutch vs SAS
##   D — Dumbbell: top FDR-significant differential genes (% prevalence Dutch vs SAS)
##   E — Slope chart: mean ± 95% CI for key declining/increasing genes, by ethnicity

library(tidyverse)
library(ggpubr)
library(ggrepel)

## ── Source subscripts ─────────────────────────────────────────────────────────
source("scripts/5_arg/1_arg_descriptive_qc.R")
source("scripts/5_arg/2_arg_cross_comparisons.R")
source("scripts/5_arg/3_arg_longitudinal_lmm.R")
## After sourcing:
##   arg_burden_clin   — sample-level burden (total_arg_rpm, log_rpm, EthnicityTot, ID, timepoint)
##   prevalence_results — gene-level ethnicity prevalence comparison (script 2)
##   prevalence_clin   — baseline prevalence + clinical data (script 2)
##   statres           — timepoint-only LMM results with padj (script 3)
##   df_tot            — gene-level wide data (script 3)

## ── Shared aesthetics ─────────────────────────────────────────────────────────
BASE_SIZE <- 11

# Ethnicity colours: Dutch = blue, South-Asian Surinamese = gold
# Dutch is assumed to be levels(EthnicityTot)[1] (alphabetical default in HELIUS).
ETH_DUTCH  <- levels(arg_burden_clin$EthnicityTot)[1]
ETH_SAS    <- levels(arg_burden_clin$EthnicityTot)[2]
eth_colors <- c("#2166AC", "#E6B800")
names(eth_colors) <- c(ETH_DUTCH, ETH_SAS)

tp_labels <- c("baseline" = "Baseline", "follow-up" = "Follow-up")

## ── Panel A: Total ARG Burden Over Time ───────────────────────────────────────
pl_A <- ggplot(arg_burden_clin,
               aes(x = timepoint, y = log_rpm, fill = timepoint)) +
  geom_violin(alpha = 0.50, colour = NA) +
  geom_boxplot(width = 0.22, fill = "white", outlier.shape = NA, colour = "gray30") +
  geom_line(aes(group = ID, colour = EthnicityTot),
            alpha = 0.12, linewidth = 0.30) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.8, hjust = 0.5,
           label = "p = 5.4\u00d710\u207b\u00b9\u2075", size = 3.2) +
  scale_fill_manual(
    values = c("baseline" = "#D9E8F5", "follow-up" = "#A8C8E8"),
    guide  = "none") +
  scale_colour_manual(values = eth_colors, name = "") +
  scale_x_discrete(labels = tp_labels) +
  theme_Publication(base_size = BASE_SIZE) +
  labs(x = "", y = "Total ARG Burden (log\u2081\u2080 RPM)",
       title = "Total ARG Burden Over Time") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))

## ── Panel B: Ethnicity × Timepoint 2×2 violins ───────────────────────────────
arg_burden_B <- arg_burden_clin %>%
  mutate(tp_label = factor(recode(as.character(timepoint),
                                  "baseline"  = "Baseline",
                                  "follow-up" = "Follow-up"),
                           levels = c("Baseline", "Follow-up")))

pl_B <- ggplot(arg_burden_B,
               aes(x = EthnicityTot, y = log_rpm, fill = EthnicityTot)) +
  geom_violin(alpha = 0.60, colour = NA) +
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
       title = "ARG Burden by Ethnicity")

## ── Panel C: Prevalence Volcano (Dutch as reference) ─────────────────────────
# prev_diff in script 2 = group2 - group1.
# If Dutch = group1, Dutch - SAS = -prev_diff  (sign_mult = -1).
dutch_idx <- which(levels(prevalence_clin$EthnicityTot) == ETH_DUTCH)
sign_mult  <- if (dutch_idx == 1L) -1 else 1

n_tested <- nrow(prevalence_results)
n_sig    <- sum(prevalence_results$padj < 0.05, na.rm = TRUE)

volcano_dat <- prevalence_results %>%
  mutate(prev_diff_dutch = prev_diff * sign_mult)

pl_C <- ggplot(volcano_dat,
               aes(x = prev_diff_dutch, y = -log10(pval), colour = sig_level)) +
  geom_point(alpha = 0.55, size = 1.8) +
  geom_text_repel(
    data            = filter(volcano_dat, sig_level == "Significant & Large Effect"),
    aes(label = gene),
    size            = 3.2, colour = "black",
    max.overlaps    = 30, box.padding = 0.5, point.padding = 0.3,
    min.segment.length = 0.2) +
  geom_text_repel(
    data            = filter(volcano_dat, sig_level == "Significant"),
    aes(label = gene),
    size            = 2.4, colour = "gray30",
    max.overlaps    = 20, box.padding = 0.4, point.padding = 0.3,
    min.segment.length = 0.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "gray50") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
           label = paste0("N tested = ", n_tested, "\nFDR sig = ", n_sig),
           size = 2.8, colour = "gray40") +
  scale_colour_manual(
    values = c("Significant & Large Effect" = "red",
               "Significant"               = "orange",
               "Not Significant"           = "gray70"),
    name = "") +
  theme_Publication(base_size = BASE_SIZE) +
  labs(x = "Prevalence Difference (% Dutch \u2212 % South-Asian Surinamese)",
       y = "\u2212log\u2081\u2080(p-value)",
       title = "ARG Prevalence Differences (Baseline)")

## ── Panel D: Dumbbell — top differential genes, Dutch vs SAS prevalence ───────
top_diff <- prevalence_results %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(prev_diff))) %>%
  head(20)

# Assign Dutch / SAS labels to group1 / group2
if (dutch_idx == 1L) {
  top_diff <- top_diff %>% rename(prev_dutch = prev_group1, prev_sas = prev_group2)
} else {
  top_diff <- top_diff %>% rename(prev_dutch = prev_group2, prev_sas = prev_group1)
}

top_diff <- top_diff %>%
  mutate(gene = fct_reorder(gene, prev_dutch))

top_diff_long <- top_diff %>%
  pivot_longer(c(prev_dutch, prev_sas),
               names_to = "eth_key", values_to = "prevalence") %>%
  mutate(ethnicity = if_else(eth_key == "prev_dutch", ETH_DUTCH, ETH_SAS))

pl_D <- ggplot(top_diff_long, aes(y = gene, x = prevalence, colour = ethnicity)) +
  geom_line(aes(group = gene), colour = "gray75", linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "gray55") +
  scale_colour_manual(values = eth_colors, name = "") +
  scale_x_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%"),
                     breaks = seq(0, 100, 25)) +
  theme_Publication(base_size = BASE_SIZE) +
  labs(x = "Prevalence (% of samples)", y = "",
       title = "Top Differential ARG Genes (Baseline)")

## ── Panel E: Slope chart — key genes mean ± 95% CI by ethnicity ───────────────
fdr_sig <- statres %>%
  filter(padj < 0.05) %>%
  arrange(padj)

top_declines <- head(fdr_sig$mbname, 3)

gene_G <- if ("catB" %in% fdr_sig$mbname) {
  "catB"
} else {
  fdr_sig %>% filter(estimate > 0) %>% slice_min(padj, n = 1) %>% pull(mbname)
}
if (length(gene_G) == 0 || is.na(gene_G[1])) gene_G <- fdr_sig$mbname[4]

key_genes <- c(top_declines[1:3], gene_G[1])

slope_data <- lapply(key_genes, function(nm) {
  gr      <- statres %>% filter(mbname == nm)
  sub_cls <- if (!is.na(gr$subclass[1])) gr$subclass[1] else gr$class[1]
  fdr_p   <- gr$padj[1]
  fdr_str <- if (!is.na(fdr_p) && fdr_p < 0.001) {
    sprintf("%.2e", fdr_p)
  } else {
    sprintf("%.3f", fdr_p)
  }
  gene_lab <- paste0(nm, "\n", sub_cls, "\nFDR p = ", fdr_str)

  df_tot %>%
    mutate(mb = log10(.data[[nm]] + 1)) %>%
    group_by(EthnicityTot, timepoint) %>%
    summarise(mean  = mean(mb, na.rm = TRUE),
              se    = sd(mb, na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    mutate(ci_lo    = mean - 1.96 * se,
           ci_hi    = mean + 1.96 * se,
           gene     = nm,
           gene_lab = gene_lab)
}) %>%
  bind_rows() %>%
  mutate(gene_lab  = factor(gene_lab, levels = unique(gene_lab)),
         timepoint = factor(timepoint, levels = c("baseline", "follow-up")))

pl_E <- ggplot(slope_data,
               aes(x = timepoint, y = mean, colour = EthnicityTot,
                   group = EthnicityTot)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.10, linewidth = 0.7) +
  facet_wrap(~gene_lab, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = eth_colors, name = "") +
  scale_x_discrete(labels = tp_labels) +
  theme_Publication(base_size = BASE_SIZE) +
  labs(x = "", y = "log\u2081\u2080(RPKM + 1)",
       title = "Key ARG Gene Trajectories (mean \u00b1 95% CI)")

## ── Assemble Figure 5 ────────────────────────────────────────────────────────
top_row <- ggarrange(
  pl_A, pl_B, pl_C,
  ncol          = 3,
  widths        = c(1, 1.2, 1.2),
  labels        = c("A", "B", "C"),
  font.label    = list(size = 14, face = "bold"),
  common.legend = FALSE
)

bot_row <- ggarrange(
  pl_D, pl_E,
  ncol          = 2,
  widths        = c(1, 2),
  labels        = c("D", "E"),
  font.label    = list(size = 14, face = "bold"),
  common.legend = FALSE
)

fig5 <- ggarrange(
  top_row, bot_row,
  nrow    = 2,
  heights = c(1, 1.2)
)

dir.create("results/5_arg", showWarnings = FALSE, recursive = TRUE)
ggsave(fig5, filename = "results/5_arg/figure5.pdf",
       width = 18, height = 14, device = cairo_pdf)
