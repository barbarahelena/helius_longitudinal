# Gutsmash 

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
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
} 

# Data
df <- rio::import("data/shotgun/gutsmash_results/population_pathways.tsv")
dim(df)
head(df)[1:5,1:5]
df <- df |> filter(! sample %in% c("HELIBA_103370", "HELIFU_103370"))
head(df)
rownames(df) <- df$sample
df$sample <- NULL
df <- as.matrix(df)
df_rel <- df / rowSums(df)
df_rel <- as.data.frame(df_rel)

# Filter: keep pathways with abundance >= 0.01 in at least 15% of subjects
prev_threshold <- 0.15
abund_threshold <- 0.01
keep_pw <- colMeans(df_rel >= abund_threshold) >= prev_threshold
df_rel <- df_rel[, keep_pw]

pathway_cols <- colnames(df_rel)
df_rel$sampleID <- rownames(df_rel)

clinical <- readRDS("data/clinicaldata_long.RDS")
df_clin  <- df_rel |> left_join(clinical, by = "sampleID") |> droplevels()
df_clin <- df_clin |> filter(!is.na(EthnicityTot))
table(df_clin$EthnicityTot) # sample per ethnicity

# ---------------------------------------------------------------------------
# Violin plots: top 30 pathways by ethnicity (ggarrange)
# ---------------------------------------------------------------------------
# Rank pathways by mean relative abundance
pw_prev_wide <- tibble(
    pathway = pathway_cols,
    mean_abund = colMeans(df_clin[, pathway_cols], na.rm = TRUE)
) |> arrange(desc(mean_abund))

top30 <- head(pw_prev_wide$pathway, 60)

wilcox_res <- df_clin |>
    select(sampleID, EthnicityTot, timepoint, all_of(top30)) |>
    pivot_longer(cols = all_of(top30), names_to = "pathway", values_to = "abundance") |>
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

cat("\nOf the top 60 pathways,", length(sig_pathways),
    "show a significant ethnicity difference at baseline or follow-up (padj < 0.05)\n")

top30 <- sig_pathways

# Pivot to long for the significant pathways
df_pw_long <- df_clin |>
    select(sampleID, EthnicityTot, timepoint, all_of(top30)) |>
    pivot_longer(cols = all_of(top30), names_to = "pathway", values_to = "relative_abundance") |>
    mutate(pathway_label = str_remove_all(pathway, "[\\[\\]']"),
                log_relab = log10(relative_abundance + 0.01))

# One violin per pathway
plot_list <- list()
plot_list <- lapply(top30, \(pw) {
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

dir.create("results/4_functional_change/gutsmash", showWarnings = FALSE, recursive = TRUE)
# 15 per page, dynamic number of pages
n_pages <- ceiling(length(plot_list) / 15)
for (i in seq_len(n_pages)) {
    idx <- ((i - 1) * 15 + 1) : min(i * 15, length(plot_list))
    ggarrange(plotlist = plot_list[idx], nrow = 5, ncol = 3, labels = LETTERS[idx])
    ggsave(sprintf("results/4_functional_change/gutsmash/gutsmashpw_page%d.pdf", i), width = 15, height = 30)
}
