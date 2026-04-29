# HUMAnN pathways — cross-sectional ethnicity comparison

# Libraries
library(tidyverse)
library(ggpubr)
library(ggsci)

# Theme
theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(),
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
                strip.text = element_text(face="bold")
        ))
}

# Data — pathways in rows, samples in columns; transpose to samples x pathways
raw <- read_tsv("data/shotgun/humann/merged_tables_renorm_unstratified.tsv",
                show_col_types = FALSE)

# Pivot to samples x pathways
df <- raw |>
    pivot_longer(-`# Pathway`, names_to = "sampleID", values_to = "cpm") |>
    pivot_wider(names_from = `# Pathway`, values_from = cpm)

# Strip "_Abundance-CPM" suffix to match clinical sampleIDs
df <- df |>
    mutate(sampleID = str_remove(sampleID, "_Abundance-CPM"))

# Remove known outlier samples
df <- df |> filter(!sampleID %in% c("HELIBA_103370", "HELIFU_103370"))

# Convert to relative proportions within each sample
pathway_cols <- setdiff(names(df), "sampleID")
df_mat <- as.matrix(df[, pathway_cols])
row_sums <- rowSums(df_mat, na.rm = TRUE)
df_rel_mat <- df_mat / row_sums
df_rel <- as.data.frame(df_rel_mat)
df_rel$sampleID <- df$sampleID

# Filter: keep pathways with relative abundance >= 0.005 in >= 10% of samples
prev_threshold  <- 0.25
abund_threshold <- 0.0025
keep_pw <- colMeans(df_rel[, pathway_cols] >= abund_threshold, na.rm = TRUE) >= prev_threshold
pathway_cols <- pathway_cols[keep_pw]
df_rel <- df_rel[, c("sampleID", pathway_cols)]

pseudocount <- min(df_rel[, pathway_cols][df_rel[, pathway_cols] > 0], na.rm = TRUE) * 100 / 2

# Join clinical metadata
clinical <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
df_clin  <- df_rel |> left_join(clinical, by = "sampleID") |> droplevels()
df_clin  <- df_clin |> filter(!is.na(EthnicityTot))
table(df_clin$EthnicityTot)

# ---------------------------------------------------------------------------
# Test all filtered pathways by Wilcoxon, then keep significant ones for plots
# ---------------------------------------------------------------------------
wilcox_res <- df_clin |>
    dplyr::select(sampleID, EthnicityTot, timepoint, all_of(pathway_cols)) |>
    pivot_longer(cols = all_of(pathway_cols), names_to = "pathway", values_to = "abundance") |>
    group_by(timepoint, pathway) |>
    summarise(
        pvalue = {
            g <- split(abundance, EthnicityTot)
            wilcox.test(g[[1]], g[[2]])$p.value
        },
        .groups = "drop"
    ) |>
    group_by(timepoint) |>
    mutate(padj = p.adjust(pvalue, method = "BH")) |>
    ungroup()

sig_pathways <- wilcox_res |>
    filter(padj < 0.05) |>
    distinct(pathway) |>
    pull(pathway)

sig_baseline <- wilcox_res |> filter(timepoint == "baseline",  padj < 0.05) |> pull(pathway)
sig_followup <- wilcox_res |> filter(timepoint == "follow-up", padj < 0.05) |> pull(pathway)
cat("\nOf", length(pathway_cols), "filtered pathways,", length(sig_pathways),
    "show a significant ethnicity difference at baseline or follow-up (padj < 0.05)\n")
cat("  Baseline: ", length(sig_baseline), "significant\n")
cat("  Follow-up:", length(sig_followup), "significant\n")
cat("  Overlap:  ", length(intersect(sig_baseline, sig_followup)), "significant at both\n")

# Save cross-sectional results table
wilcox_wide <- wilcox_res |>
    dplyr::select(pathway, timepoint, padj) |>
    pivot_wider(names_from = timepoint, values_from = padj,
                names_prefix = "padj_") |>
    mutate(
        sig_baseline = `padj_baseline`  < 0.05,
        sig_followup = `padj_follow-up` < 0.05
    )
dir.create("results/4_functional_change/humann", showWarnings = FALSE, recursive = TRUE)
write.csv2(wilcox_wide,
           "results/4_functional_change/humann/crosssectional_wilcox_results.csv",
           row.names = FALSE)

plot_pathways <- sig_pathways

# ---------------------------------------------------------------------------
# Violin plots for significant pathways
# ---------------------------------------------------------------------------
df_pw_long <- df_clin |>
    dplyr::select(sampleID, EthnicityTot, timepoint, all_of(plot_pathways)) |>
    pivot_longer(cols = all_of(plot_pathways), names_to = "pathway", values_to = "relative_abundance") |>
    mutate(
        pathway_label = str_remove(pathway, "^[A-Z0-9_-]+: "),
        log_relab     = log10(relative_abundance * 100 + pseudocount)
    )

plot_list <- lapply(plot_pathways, \(pw) {
    sub <- df_pw_long |> filter(pathway == pw)
    ggplot(sub, aes(x = EthnicityTot, y = log_relab, fill = EthnicityTot)) +
        geom_violin(colour = NA) +
        geom_boxplot(fill = "white", width = 0.25) +
        geom_point(alpha = 0.2, size = 0.5) +
        scale_fill_jco(guide = "none") +
        stat_compare_means(size = 2.5) +
        facet_wrap(~ timepoint) +
        labs(title = unique(sub$pathway_label), x = "", y = "Relative abundance (log10)") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

dir.create("results/4_functional_change/humann", showWarnings = FALSE, recursive = TRUE)

n_pages <- ceiling(length(plot_list) / 15)
for (i in seq_len(n_pages)) {
    idx <- ((i - 1) * 15 + 1) : min(i * 15, length(plot_list))
    ggarrange(plotlist = plot_list[idx], nrow = 5, ncol = 3, labels = LETTERS[idx])
    ggsave(sprintf("results/4_functional_change/humann/humann_pw_page%d.pdf", i),
           width = 15, height = 30, device = cairo_pdf)
}
