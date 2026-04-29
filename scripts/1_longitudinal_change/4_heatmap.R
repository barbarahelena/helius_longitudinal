# Heatmap: 16S core OTUs (>95% prevalence in any ethnic group) across ethnic groups
# Replicating HELIUS Fig. 3 style — baseline, follow-up, and delta versions
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

library(tidyverse)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(ggpubr)
library(Cairo)
library(grid)

theme_Publication <- function(base_size = 14, base_family = "sans") {
    library(grid)
    library(ggthemes)
    suppressWarnings(theme_foundation(base_size = base_size, base_family = base_family) +
        theme(
            plot.title   = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
            text         = element_text(),
            panel.background  = element_rect(colour = NA, fill = NA),
            plot.background   = element_rect(colour = NA, fill = NA),
            panel.border      = element_rect(colour = NA),
            axis.title        = element_text(face = "bold", size = rel(0.8)),
            axis.title.y      = element_text(angle = 90, vjust = 2),
            axis.title.x      = element_text(vjust = -0.2),
            axis.text         = element_text(size = rel(0.7)),
            axis.text.x       = element_text(angle = 0),
            axis.line         = element_line(colour = "black"),
            axis.ticks        = element_line(),
            panel.grid.major  = element_line(colour = "#f0f0f0"),
            panel.grid.minor  = element_blank(),
            legend.key        = element_rect(colour = NA),
            legend.position   = "bottom",
            legend.key.size   = unit(0.2, "cm"),
            legend.spacing    = unit(0, "cm"),
            plot.margin       = unit(c(10, 5, 5, 5), "mm"),
            strip.background  = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text        = element_text(face = "bold"),
            plot.caption      = element_text(size = rel(0.5), face = "italic")
        ))
}

resultsfolder <- "results/1_longitudinal_change/heatmap"
dir.create(resultsfolder, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
ps  <- readRDS("data/16s/phyloseq_withclinpaired.RDS")
sdf <- as.data.frame(sample_data(ps))  # rownames = sa1, sa2, ...; sampleID col = HELIBA_xxx

otu_counts <- as(otu_table(ps), "matrix")        # taxa x samples
otu_rel    <- sweep(otu_counts, 2, colSums(otu_counts), "/")  # relative abundance

# Taxonomy labels
tt <- as.data.frame(tax_table(ps)) |>
    rownames_to_column("ASV") |>
    mutate(label = case_when(
        !is.na(Genus) & !is.na(Species) ~ paste(Genus, Species),
        !is.na(Genus)                   ~ Genus,
        !is.na(Family)                  ~ paste0(Family, " spp."),
        TRUE                            ~ ASV
    ))

# ---- Build tidy data frame (samples x ASVs + metadata) ----
otu_df <- as.data.frame(t(otu_rel))
otu_df$sampleID <- rownames(otu_df)
otu_df <- left_join(otu_df,
                    sdf |> dplyr::select(sampleID, EthnicityTot, timepoint),
                    by = "sampleID")

eth_levels <- c("Dutch", "Ghanaian", "Moroccan", "African Surinamese",
                "South-Asian Surinamese", "Turkish")
eth_abbr   <- c("D", "G", "M", "AS", "SAS", "T")
all_groups <- c("All", eth_levels)
all_abbr   <- c("All", eth_abbr)

# ---- Identify core ASVs: >95% prevalent in any ethnic group at baseline ----
asv_cols <- setdiff(names(otu_df), c("sampleID", "EthnicityTot", "timepoint"))

otu_base <- otu_df |> filter(timepoint == "baseline", EthnicityTot %in% eth_levels)

max_prev <- setNames(rep(0, length(asv_cols)), asv_cols)
for (eth in eth_levels) {
    sub_vals  <- otu_base |> filter(EthnicityTot == eth) |> dplyr::select(all_of(asv_cols))
    prev      <- colMeans(sub_vals > 0)
    max_prev  <- pmax(max_prev, prev)
}
core_asvs_base <- names(max_prev[max_prev >= 0.95])
cat("Core ASVs (>95% in any ethnic group at baseline):", length(core_asvs_base), "\n")

# ---- Also identify core ASVs at follow-up ----
otu_fu_eth <- otu_df |> filter(timepoint == "follow-up", EthnicityTot %in% eth_levels)

max_prev_fu <- setNames(rep(0, length(asv_cols)), asv_cols)
for (eth in eth_levels) {
    sub_vals      <- otu_fu_eth |> filter(EthnicityTot == eth) |> dplyr::select(all_of(asv_cols))
    prev          <- colMeans(sub_vals > 0)
    max_prev_fu   <- pmax(max_prev_fu, prev)
}
core_asvs_fu <- names(max_prev_fu[max_prev_fu >= 0.95])
cat("Core ASVs (>95% in any ethnic group at follow-up):", length(core_asvs_fu), "\n")

# Union of both cores as the row set for the two-panel figure
core_asvs_union <- union(core_asvs_base, core_asvs_fu)
cat("Union core ASVs (rows in two-panel figure):", length(core_asvs_union), "\n")

# Baseline-only core kept for single-timepoint and delta figures
core_asvs <- core_asvs_base

row_labels_vec <- tt |>
    filter(ASV %in% core_asvs_union) |>
    dplyr::select(ASV, label) |>
    deframe()

# ---- Helper: compute median and prevalence matrices per group ----
compute_group_matrices <- function(otu_tp_df, asvs) {
    med_mat  <- matrix(NA_real_, nrow = length(asvs), ncol = length(all_groups),
                       dimnames = list(asvs, all_abbr))
    prev_mat <- matrix(NA_real_, nrow = length(asvs), ncol = length(all_groups),
                       dimnames = list(asvs, all_abbr))

    for (j in seq_along(all_groups)) {
        sub <- if (all_groups[j] == "All") {
            otu_tp_df |> filter(EthnicityTot %in% eth_levels)
        } else {
            otu_tp_df |> filter(EthnicityTot == all_groups[j])
        }
        vals <- sub |> dplyr::select(all_of(asvs))
        med_mat[asvs, j]  <- sapply(vals, median)
        prev_mat[asvs, j] <- sapply(vals, function(x) mean(x > 0))
    }
    list(median = med_mat, prevalence = prev_mat)
}

# Single-panel baseline heatmap uses baseline core only
base_mats <- compute_group_matrices(otu_df |> filter(timepoint == "baseline"), core_asvs)

# Two-panel figure uses the union core, each timepoint masked by its own prevalence
base_mats_union <- compute_group_matrices(otu_df |> filter(timepoint == "baseline"), core_asvs_union)
fu_mats_union   <- compute_group_matrices(otu_df |> filter(timepoint == "follow-up"),  core_asvs_union)

# ---- Color scale for abundance (log10; dark red = high, green = low) ----
col_fun_abund <- colorRamp2(
    c(-3.1, -2.5, -2.0, -1.5, -1.0, -0.7),
    c("#1A9850", "#66BD63", "#FEE08B", "#FDAE61", "#F46D43", "#A50026")
)

lgd_abund <- Legend(
    col_fun   = col_fun_abund,
    title     = "Median rel. ab. (log10)",
    at        = c(-3.1, -2.5, -2.0, -1.5, -1.0, -0.7),
    labels    = c("0.0008", "0.003", "0.01", "0.03", "0.1", "0.2")
)

# ---- Helper: mask cells below prevalence threshold and log10-transform ----
mask_and_log <- function(med_mat, prev_mat, threshold = 0.95) {
    log_mat <- log10(pmax(med_mat, 1e-6))
    log_mat[prev_mat < threshold] <- NA
    log_mat
}

# ---- Helper: build ComplexHeatmap ----
make_core_heatmap <- function(med_mat, prev_mat, col_fun, row_labs,
                              title = "", show_legend = FALSE,
                              name = "Median rel.ab. (log10)") {
    plot_mat <- mask_and_log(med_mat, prev_mat)

    # For clustering: replace NAs with the minimum observed value so hclust works
    cluster_mat <- plot_mat
    cluster_mat[is.na(cluster_mat)] <- min(plot_mat, na.rm = TRUE)

    Heatmap(
        plot_mat,
        name              = name,
        col               = col_fun,
        rect_gp           = gpar(col = "white", lwd = 1.5),
        na_col            = "white",
        cluster_rows      = TRUE,
        cluster_columns   = FALSE,
        clustering_distance_rows = function(x) dist(cluster_mat),
        show_row_names    = TRUE,
        show_column_names = TRUE,
        row_names_side    = "left",
        show_row_dend     = FALSE,
        row_names_gp      = gpar(fontsize = 9),
        column_names_gp   = gpar(fontsize = 11, fontface = "bold"),
        row_labels        = row_labs[rownames(plot_mat)],
        row_names_rot     = 0,
        column_names_rot  = 0,
        column_title      = title,
        column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
        show_heatmap_legend = show_legend,
        width             = unit(7, "cm")
    )
}

# ---- Figure 1: Baseline heatmap ----
ht_base <- make_core_heatmap(base_mats$median, base_mats$prevalence,
                              col_fun_abund, row_labels_vec,
                              title = "Baseline")

CairoPDF(file.path(resultsfolder, "heatmap_16s_core_baseline.pdf"), width = 9, height = 11)
draw(ht_base,
     annotation_legend_list = list(lgd_abund),
     padding = unit(c(5, 5, 5, 60), "mm"))
dev.off()

# ---- Figure 2: Baseline + Follow-up side-by-side (union core, each panel masked independently) ----
ht_base_union <- make_core_heatmap(base_mats_union$median, base_mats_union$prevalence,
                                    col_fun_abund, row_labels_vec,
                                    title = "Baseline")
ht_fu <- make_core_heatmap(fu_mats_union$median, fu_mats_union$prevalence,
                            col_fun_abund, row_labels_vec,
                            title = "Follow-up")

CairoPDF(file.path(resultsfolder, "heatmap_16s_core_baseline_vs_followup.pdf"), width = 15, height = 11)
draw(ht_base_union + ht_fu,
     annotation_legend_list = list(lgd_abund),
     padding = unit(c(5, 5, 5, 60), "mm"))
dev.off()