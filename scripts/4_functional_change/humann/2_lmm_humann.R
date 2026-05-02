# HUMAnN pathways — longitudinal LMM (ethnicity x timepoint interaction)

library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggrepel)
library(aplot)

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
              axis.title.y = element_text(angle=90, vjust=2),
              axis.title.x = element_text(vjust=-0.2),
              axis.text = element_text(size = rel(0.7)),
              axis.line = element_line(colour="black"),
              axis.ticks = element_line(),
              panel.grid.major = element_line(colour="#f0f0f0"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.key.size = unit(0.2, "cm"),
              legend.spacing  = unit(0, "cm"),
              strip.background = element_rect(colour="#f0f0f0", fill="#f0f0f0"),
              strip.text = element_text(face="bold"),
              plot.caption  = element_text(size = rel(0.5), face = "italic"),
              plot.subtitle = element_text(size=8, hjust=0.5, face="italic")))
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
raw <- read_tsv("data/shotgun/humann/merged_tables_renorm_unstratified.tsv",
                show_col_types = FALSE)

df <- raw |>
    pivot_longer(-`# Pathway`, names_to = "sampleID", values_to = "cpm") |>
    pivot_wider(names_from = `# Pathway`, values_from = cpm) |>
    mutate(sampleID = str_remove(sampleID, "_Abundance-CPM")) |>
    filter(!sampleID %in% c("HELIBA_103370", "HELIFU_103370"))

pathway_cols <- setdiff(names(df), "sampleID")
df_mat       <- as.matrix(df[, pathway_cols])
df_rel_mat   <- df_mat / rowSums(df_mat, na.rm = TRUE)
df_rel       <- as.data.frame(df_rel_mat)
df_rel$sampleID <- df$sampleID

# Filter: relative abundance >= 0.005 in >= 10% of samples
prev_threshold  <- 0.25
abund_threshold <- 0.0025
keep_pw      <- colMeans(df_rel[, pathway_cols] >= abund_threshold, na.rm = TRUE) >= prev_threshold
pathway_cols <- pathway_cols[keep_pw]
df_rel       <- df_rel[, c("sampleID", pathway_cols)]

pseudocount <- min(df_rel[, pathway_cols][df_rel[, pathway_cols] > 0], na.rm = TRUE) * 100 / 2

clinical <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
df_clin  <- df_rel |> left_join(clinical, by = "sampleID") |> droplevels()
df_clin  <- df_clin |> filter(!is.na(EthnicityTot))
table(df_clin$EthnicityTot)

dir.create("results/4_functional_change/humann/longitudinal", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# LMM for each pathway — ethnicity x timepoint interaction
# ---------------------------------------------------------------------------
statres <- data.frame()
for (pw in pathway_cols) {
    df_clin$mb <- log10(df_clin[[pw]] * 100 + pseudocount)
    model1 <- lmer(mb ~ EthnicityTot * timepoint + FUtime + (1|ID), data = df_clin)
    res    <- summary(model1)
    confint_model1 <- confint(model1, method = "Wald")
    interaction_row <- grep("EthnicityTotSouth-Asian Surinamese:timepoint", rownames(res$coefficients))
    ci_row          <- grep("EthnicityTotSouth-Asian Surinamese:timepoint", rownames(confint_model1))
    statres <- rbind(statres, data.frame(
        pathway  = pw,
        estimate = res$coefficients[interaction_row, 1],
        conflow  = ifelse(length(ci_row) > 0, confint_model1[ci_row, 1], NA),
        confhigh = ifelse(length(ci_row) > 0, confint_model1[ci_row, 2], NA),
        pval     = res$coefficients[interaction_row, 5]
    ))
}

statres <- statres |>
    arrange(pval) |>
    mutate(padj = p.adjust(pval, method = "fdr"))

write.csv2(statres,
           "results/4_functional_change/humann/longitudinal/lmm_ethnicity_timepoint_results.csv",
           row.names = FALSE)

# ---------------------------------------------------------------------------
# Longitudinal line plots — top 20 significant pathways
# ---------------------------------------------------------------------------
statres_sig <- statres |> filter(padj < 0.05) |> arrange(pval)

plist <- list()
for (i in seq_len(min(nrow(statres_sig), 20))) {
    nm <- statres_sig$pathway[i]
    df_clin$mb <- log10(df_clin[[nm]] * 100 + pseudocount)

    df_means <- df_clin |>
        group_by(EthnicityTot, timepoint) |>
        summarise(mean = mean(mb, na.rm = TRUE), sd = sd(mb, na.rm = TRUE),
                  n = n(), .groups = "drop")

    res_lmm <- statres_sig |>
        filter(pathway == nm) |>
        mutate(group1 = "baseline", group2 = "follow-up",
               sigq = paste0("p=", formatC(pval, format = "e", digits = 2))) |>
        dplyr::select(-pathway)

    mbmax  <- ifelse(max(df_clin$mb, na.rm = TRUE) < 0,
                     max(df_clin$mb, na.rm = TRUE) * 0.8,
                     max(df_clin$mb, na.rm = TRUE) * 1.2)
    mbstat <- ifelse(max(df_clin$mb, na.rm = TRUE) < 0,
                     max(df_clin$mb, na.rm = TRUE) * 0.9,
                     max(df_clin$mb, na.rm = TRUE) * 1.1)
    mbmin  <- min(df_clin$mb, na.rm = TRUE)

    nm_label <- str_remove(nm, "^[A-Z0-9_-]+: ")

    pl <- ggplot() +
        geom_line(data = df_clin, aes(x = timepoint, y = mb, color = EthnicityTot, group = ID),
                  alpha = 0.07, linewidth = 0.5) +
        geom_point(data = df_clin, aes(x = timepoint, y = mb, color = EthnicityTot),
                   alpha = 0.07, size = 0.8) +
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
        labs(x = "Timepoint", y = "log10(abundance %)", title = nm_label, color = "")

    plist[[i]] <- pl
}

n_plots <- length(plist)
if (n_plots > 0) {
    n_cols <- 3
    n_rows <- ceiling(n_plots / n_cols)
    plots  <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
                        labels = LETTERS[1:n_plots], nrow = n_rows, ncol = n_cols)
    ggsave("results/4_functional_change/humann/longitudinal/significant_humann_timepoint_ethnicity.pdf",
           plots, width = 12, height = 4 * n_rows, device = cairo_pdf)
} else {
    message("No FDR-significant pathways — longitudinal plot skipped.")
}

# ---------------------------------------------------------------------------
# Forest plot (FDR < 0.05)
# ---------------------------------------------------------------------------
humann_fdr <- statres |>
    filter(padj < 0.05) |>
    mutate(label     = str_remove(pathway, "^[A-Z0-9_-]+: "),
           label     = str_trunc(label, 55),
           label     = fct_reorder(label, estimate),
           direction = ifelse(estimate > 0, "SAS", "Dutch"))

if (nrow(humann_fdr) >= 2) {
    pl_forest <- ggplot(humann_fdr, aes(x = estimate, y = label, colour = direction)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        geom_pointrange(aes(xmin = conflow, xmax = confhigh), size = 0.45, linewidth = 0.55) +
        scale_colour_manual(values = c("Dutch" = "#2166AC", "SAS" = "#E6B800"),
                            name = "",
                            labels = c("Dutch" = "More increase in Dutch",
                                       "SAS"   = "More increase in SAS")) +
        theme_Publication() +
        labs(x = "Interaction effect (± 95% CI)", y = "",
             title    = "HUMAnN pathway shifts",
             subtitle = "LMM ethnicity × timepoint, FDR < 0.05")
} else {
    pl_forest <- ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No FDR-significant\nHUMAnN pathways",
                 hjust = 0.5, vjust = 0.5, size = 4)
}

# ---------------------------------------------------------------------------
# Top-3 FDR-significant box-violin plots
# ---------------------------------------------------------------------------
make_boxviolin <- function(df, pathway_name, pval, y_label) {
    df$mb <- log10(df[[pathway_name]] * 100 + pseudocount)
    pval_label <- formatC(pval, format = "e", digits = 2)
    ggplot(df, aes(x = EthnicityTot, y = mb, fill = EthnicityTot)) +
        geom_violin(colour = NA, aes(alpha = timepoint)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        facet_wrap(~timepoint) +
        scale_fill_jco(guide = "none") +
        scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
        theme_Publication() +
        labs(x = "", y = y_label,
             subtitle = paste0("Ethnicity × Timepoint p=", pval_label))
}

humann_top3 <- statres |> filter(padj < 0.05) |> arrange(padj) |> slice_head(n = 3)

if (nrow(humann_top3) > 0) {
    pl_top <- lapply(seq_len(nrow(humann_top3)), function(i) {
        nm       <- humann_top3$pathway[i]
        nm_clean <- str_trunc(str_remove(nm, "^[A-Z0-9_-]+: "), 40)
        make_boxviolin(df_clin, nm, humann_top3$pval[i],
                       y_label = "log10(abundance % + pseudocount)") +
            labs(title = nm_clean) +
            theme(plot.title = element_text(face = "bold", size = rel(0.75), hjust = 0.5))
    })
    pl_top_panel <- ggarrange(plotlist = pl_top, ncol = length(pl_top),
                              labels = LETTERS[seq_len(nrow(humann_top3))])
    ggsave("results/4_functional_change/humann/longitudinal/top3_humann_boxviolin.pdf",
           pl_top_panel, width = 6 * nrow(humann_top3), height = 5, device = cairo_pdf)
} else {
    message("No FDR-significant pathways — top-3 box-violin plot skipped.")
}

ggsave("results/4_functional_change/humann/longitudinal/forest_humann.pdf",
       pl_forest, width = 10, height = max(4, nrow(humann_fdr) * 0.3 + 2), device = cairo_pdf)

# ---------------------------------------------------------------------------
# Overlap: cross-sectional significant vs LMM interaction significant
# ---------------------------------------------------------------------------
cs_results <- read.csv2("results/4_functional_change/humann/crosssectional_wilcox_results.csv",
                        stringsAsFactors = FALSE) |>
    mutate(across(starts_with("padj_"), as.numeric))

# Cross-sectional direction: which group has higher median log-abundance per timepoint
cs_direction <- df_clin |>
    dplyr::select(sampleID, EthnicityTot, timepoint, all_of(pathway_cols)) |>
    pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "relab") |>
    mutate(log_relab = log10(relab * 100 + pseudocount)) |>
    group_by(pathway, timepoint, EthnicityTot) |>
    summarise(median_log = median(log_relab, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = EthnicityTot, values_from = median_log) |>
    mutate(cs_direction = ifelse(`South-Asian Surinamese` > Dutch, "SAS higher", "Dutch higher")) |>
    dplyr::select(pathway, timepoint, cs_direction)

cs_direction_baseline <- cs_direction |>
    filter(timepoint == "baseline") |>
    dplyr::select(pathway, cs_baseline_direction = cs_direction)

cs_direction_followup <- cs_direction |>
    filter(timepoint == "follow-up") |>
    dplyr::select(pathway, cs_followup_direction = cs_direction)

# LMM interaction direction
lmm_annotated <- statres |>
    mutate(
        sig_lmm_interaction = padj < 0.05,
        lmm_direction = ifelse(estimate > 0, "SAS increases more", "Dutch increases more")
    )

# Combine all
overlap <- cs_results |>
    left_join(lmm_annotated |> dplyr::select(pathway, sig_lmm_interaction, lmm_direction, estimate, padj),
              by = "pathway") |>
    left_join(cs_direction_baseline, by = "pathway") |>
    left_join(cs_direction_followup, by = "pathway") |>
    mutate(
        pattern = case_when(
            # Both CS timepoints significant
            sig_lmm_interaction & sig_baseline & sig_followup &
                cs_baseline_direction == "SAS higher"  & lmm_direction == "SAS increases more" ~ "divergence (SAS)",
            sig_lmm_interaction & sig_baseline & sig_followup &
                cs_baseline_direction == "Dutch higher" & lmm_direction == "Dutch increases more" ~ "divergence (Dutch)",
            sig_lmm_interaction & sig_baseline & sig_followup &
                cs_baseline_direction == "SAS higher"  & lmm_direction == "Dutch increases more" ~ "convergence (SAS shrinks)",
            sig_lmm_interaction & sig_baseline & sig_followup &
                cs_baseline_direction == "Dutch higher" & lmm_direction == "SAS increases more" ~ "convergence (Dutch shrinks)",
            # Baseline only
            sig_lmm_interaction & sig_baseline & !sig_followup &
                cs_baseline_direction == "SAS higher"  & lmm_direction == "Dutch increases more" ~ "convergence (SAS shrinks)",
            sig_lmm_interaction & sig_baseline & !sig_followup &
                cs_baseline_direction == "Dutch higher" & lmm_direction == "SAS increases more" ~ "convergence (Dutch shrinks)",
            sig_lmm_interaction & sig_baseline & !sig_followup ~ "divergence (lost significance)",
            # Follow-up only
            sig_lmm_interaction & !sig_baseline & sig_followup &
                cs_followup_direction == "SAS higher"  & lmm_direction == "SAS increases more" ~ "emergence (SAS)",
            sig_lmm_interaction & !sig_baseline & sig_followup &
                cs_followup_direction == "Dutch higher" & lmm_direction == "Dutch increases more" ~ "emergence (Dutch)",
            sig_lmm_interaction & !sig_baseline & sig_followup ~ "emergence (other)",
            # LMM sig, neither CS timepoint sig
            sig_lmm_interaction & !sig_baseline & !sig_followup ~ "LMM only",
            TRUE ~ NA_character_
        )
    )

cat("\n--- Cross-sectional vs LMM interaction overlap ---\n")
lmm_sig_n <- sum(overlap$sig_lmm_interaction, na.rm = TRUE)
cat("LMM interaction-significant pathways (FDR < 0.05):", lmm_sig_n, "\n")
cat("Of those, also sig at baseline in cross-sectional: ",
    sum(overlap$sig_lmm_interaction & overlap$sig_baseline,  na.rm = TRUE), "\n")
cat("Of those, also sig at follow-up in cross-sectional: ",
    sum(overlap$sig_lmm_interaction & overlap$sig_followup,  na.rm = TRUE), "\n")
cat("Of those, sig at both timepoints in cross-sectional:",
    sum(overlap$sig_lmm_interaction & overlap$sig_baseline & overlap$sig_followup, na.rm = TRUE), "\n")
cat("\nPattern classification (LMM-significant pathways):\n")
print(table(overlap$pattern[overlap$sig_lmm_interaction], useNA = "ifany"))

write.csv2(overlap,
           "results/4_functional_change/humann/longitudinal/crosssectional_vs_lmm_overlap.csv",
           row.names = FALSE)

# ---------------------------------------------------------------------------
# Forest plot + cross-sectional heatmap (aplot)
# ---------------------------------------------------------------------------
if (nrow(humann_fdr) >= 2) {
    pw_vec <- as.character(humann_fdr$pathway)

    # log10 fold change: log10(SAS median) − log10(Dutch median) = log10(SAS/Dutch)
    heatmap_diff <- df_clin |>
        dplyr::select(sampleID, EthnicityTot, timepoint, all_of(pw_vec)) |>
        pivot_longer(cols = all_of(pw_vec), names_to = "pathway", values_to = "relab") |>
        mutate(log_relab = log10(relab * 100 + pseudocount)) |>
        group_by(pathway, timepoint, EthnicityTot) |>
        summarise(median_log = median(log_relab, na.rm = TRUE), .groups = "drop") |>
        pivot_wider(names_from = EthnicityTot, values_from = median_log) |>
        mutate(log2fc = `South-Asian Surinamese` - Dutch)

    cs_sig_mask <- overlap |>
        filter(pathway %in% pw_vec) |>
        dplyr::select(pathway, sig_baseline, sig_followup) |>
        mutate(across(c(sig_baseline, sig_followup), as.logical)) |>
        left_join(
            read.csv2("results/4_functional_change/humann/crosssectional_wilcox_results.csv",
                      stringsAsFactors = FALSE) |>
                dplyr::select(pathway, padj_baseline, padj_follow.up),
            by = "pathway"
        )

    heatmap_data <- heatmap_diff |>
        left_join(cs_sig_mask, by = "pathway") |>
        mutate(
            diff_display = case_when(
                timepoint == "baseline"  & sig_baseline  ~ log2fc,
                timepoint == "follow-up" & sig_followup  ~ log2fc,
                TRUE ~ NA_real_
            ),
            padj_cs = case_when(
                timepoint == "baseline"  ~ padj_baseline,
                timepoint == "follow-up" ~ padj_follow.up
            ),
            star = case_when(
                !is.na(diff_display) & padj_cs < 0.001 ~ "***",
                !is.na(diff_display) & padj_cs < 0.01  ~ "**",
                !is.na(diff_display) & padj_cs < 0.05  ~ "*",
                TRUE ~ ""
            )
        ) |>
        left_join(humann_fdr |> dplyr::select(pathway, label), by = "pathway") |>
        mutate(timepoint = factor(timepoint, levels = c("baseline", "follow-up")))

    abs_lim <- max(abs(heatmap_data$diff_display), na.rm = TRUE)

    pl_heatmap <- ggplot(heatmap_data,
                         aes(x = timepoint, y = label, fill = diff_display)) +
        geom_tile(color = "white", linewidth = 0.4) +
        geom_text(aes(label = star), color = "black", size = 2.5, vjust = 0.75) +
        scale_fill_gradient2(
            low      = "#2166AC",
            mid      = "white",
            high     = "#E6B800",
            na.value = "grey93",
            limits   = c(-abs_lim, abs_lim),
            name     = "log₁₀ FC\n(SAS / Dutch)",
            guide    = guide_colorbar(barheight = unit(6, "cm"),
                                      barwidth  = unit(0.5, "cm"))
        ) +
        theme_Publication() +
        theme(
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y  = element_blank(),
            legend.position = "right"
        ) +
        labs(x = "", y = "",
             caption = "* q<0.05  ** q<0.01  *** q<0.001")

    pl_combined <- pl_forest |> aplot::insert_right(pl_heatmap, width = 0.25)

    combined_height <- max(4, nrow(humann_fdr) * 0.3 + 2)
    cairo_pdf("results/4_functional_change/humann/longitudinal/forest_heatmap_humann.pdf",
              width = 14, height = combined_height)
    print(pl_combined)
    dev.off()
}
