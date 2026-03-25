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
pathway_cols <- colnames(df_rel)
df_rel$sampleID <- rownames(df_rel)

clinical <- readRDS("data/clinicaldata_long.RDS")
df_clin  <- df_rel |> left_join(clinical, by = "sampleID") |> droplevels()
df_clin <- df_clin |> filter(!is.na(EthnicityTot))
table(df_clin$EthnicityTot) # sample per ethnicity

# ---------------------------------------------------------------------------
# Per-ethnicity pathway summary
# ---------------------------------------------------------------------------
summarise_pw <- function(sub, cols) {
    mat <- sub[, cols]
    data.frame(
        pathway      = cols,
        n_samples    = nrow(mat),
        n_nonzero    = colSums(mat > 0),
        prevalence   = colMeans(mat > 0) * 100,
        mean_abund   = colMeans(mat),
        median_abund = apply(mat, 2, median),
        stringsAsFactors = FALSE
    )
}

pw_summary <- split(df_clin, df_clin$EthnicityTot) |>
    lapply(\(sub) {
        res <- summarise_pw(sub, pathway_cols)
        res$EthnicityTot <- as.character(sub$EthnicityTot[1])
        res
    }) |>
    bind_rows()

# ---------------------------------------------------------------------------
# Plot 1: prevalence histogram, faceted by ethnicity
# ---------------------------------------------------------------------------
(p_prev <- ggplot(pw_summary, aes(x = prevalence)) +
    geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
    geom_vline(xintercept = c(25, 50, 75), linetype = "dashed", colour = "grey40") +
    facet_wrap(~ EthnicityTot, scales = "free_y") +
    labs(x = "Prevalence (% of samples with value > 0)",
         y = "Number of pathways",
         title = "Pathway prevalence by ethnicity") +
    theme_Publication())

# ---------------------------------------------------------------------------
# Plot 2: prevalence vs mean abundance, faceted by ethnicity
# ---------------------------------------------------------------------------
(p_scatter <- ggplot(pw_summary, aes(x = prevalence, y = mean_abund, colour = median_abund)) +
    geom_point(size = 1, alpha = 0.7) +
    scale_color_viridis_c(name = "Median\nabundance") +
    geom_vline(xintercept = c(25, 50), linetype = "dashed", colour = "grey40") +
    facet_wrap(~ EthnicityTot, scales = "free_y") +
    labs(x = "Prevalence (%)",
         y = "Mean abundance",
         title = "Prevalence vs abundance by ethnicity") +
    theme_Publication())

# ---------------------------------------------------------------------------
# Plot 3: threshold grid per ethnicity
# ---------------------------------------------------------------------------
thresh_grid <- expand.grid(
    min_prevalence = c(10, 15, 20, 25),
    min_abundance  = c(0, 0.01, 0.05, 0.1, 0.5, 1)
)

thresholds <- split(df_clin, df_clin$EthnicityTot) |>
    lapply(\(sub) {
        mat  <- sub[, pathway_cols]
        res  <- thresh_grid
        res$n_passing <- mapply(function(prev, abund) {
            sum(colMeans(mat > abund) * 100 > prev)
        }, res$min_prevalence, res$min_abundance)
        res$EthnicityTot <- as.character(sub$EthnicityTot[1])
        res
    }) |>
    bind_rows()

(p_thresh <- ggplot(thresholds, aes(
        x = factor(min_abundance), y = n_passing,
        colour = factor(min_prevalence), group = factor(min_prevalence))) +
    geom_line() + geom_point(size = 2.5) +
    scale_color_atlassian() +
    facet_wrap(~ EthnicityTot, scales = "free_y") +
    labs(x = "Minimum abundance threshold",
         y = "Number of pathways retained",
         colour = "Min prevalence (%)",
         title = "Pathways retained at various cutoffs, by ethnicity") +
    theme_Publication())

# ---------------------------------------------------------------------------
# Wide comparison table: prevalence & mean abundance side-by-side
# ---------------------------------------------------------------------------
pw_prev_wide <- pw_summary |>
    select(pathway, EthnicityTot, prevalence, mean_abund) |>
    pivot_wider(names_from = EthnicityTot,
                values_from = c(prevalence, mean_abund),
                names_sep = "_")

# minimum prevalence across all ethnic groups
prev_cols <- grep("^prevalence_", colnames(pw_prev_wide), value = TRUE)
pw_prev_wide$min_prevalence <- apply(pw_prev_wide[, prev_cols], 1, min)
pw_prev_wide <- pw_prev_wide |> arrange(desc(min_prevalence))

cat("\n--- Top 30 pathways by minimum prevalence across all ethnicities ---\n")
print(head(pw_prev_wide, 60))
