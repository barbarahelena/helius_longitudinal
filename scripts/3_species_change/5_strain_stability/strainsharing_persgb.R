## Strain sharing — per-SGB logistic regression
## For each SGB (n > 50), logistic regression of binary strain sharing (0/1)
## against covariate blocks. Outputs ORs, forest plots, and a heatmap.
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

#### Libraries ####
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(broom)

#### Theme ####
theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0),
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
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
}

#### Output folder ####
resultsfolder <- "results/3_species_change/5_strain_stability/covariates"
dir.create(resultsfolder, showWarnings = FALSE, recursive = TRUE)

#### Data ####
df <- rio::import("data/shotgun/strainsharing_merged.csv")
colnames(df) <- str_remove(colnames(df), "sharing_")

dfsame <- df %>% filter(str_remove(sampleid_1, "HELIBA_") == str_remove(sampleid_2, "HELIFU_"))

clin <- readRDS("data/clinicaldata/clinicaldata_long.RDS")

#### Identify SGBs and build long dataset ####
cladesplit <- readRDS("data/shotgun/shotgun_taxtable_sgb.RDS")

dfsame_long_all <- dfsame %>%
    filter(str_detect(sampleid_1, "HELIBA_")) %>%
    dplyr::rename(sampleID = sampleid_1) %>%
    dplyr::select(sampleID, starts_with("t__")) %>%
    pivot_longer(cols = starts_with("t__"), names_to = "SGB_col", values_to = "shared") %>%
    mutate(SGB = str_remove(SGB_col, "t__"),
           shared = as.integer(shared)) %>%
    dplyr::select(-SGB_col)

sgb_n <- dfsame_long_all %>%
    group_by(SGB) %>%
    summarise(n = sum(!is.na(shared)), .groups = "drop") %>%
    filter(n > 237)

sgbs_use <- sgb_n$SGB
message(length(sgbs_use), " SGBs retained (n > 237)")

cladesplit_match <- cladesplit %>%
    mutate(SGB_num = str_extract(SGB, "[0-9]+")) %>%
    mutate(SGB_key = sgbs_use[match(SGB_num, str_extract(sgbs_use, "[0-9]+"))]) %>%
    filter(!is.na(SGB_key)) %>%
    dplyr::select(SGB = SGB_key, Species) %>%
    distinct(SGB, .keep_all = TRUE)

dfsame_long <- dfsame_long_all %>%
    filter(SGB %in% sgbs_use) %>%
    left_join(clin, by = "sampleID") %>%
    left_join(cladesplit_match, by = "SGB") %>%
    mutate(
        EthnicityTot = relevel(factor(EthnicityTot), ref = "Dutch"),
        Age_z        = as.numeric(scale(Age)),
        FUtime_z     = as.numeric(scale(FUtime))
    )

message("dfsame_long rows: ", nrow(dfsame_long),
        " | SGBs: ", n_distinct(dfsame_long$SGB))

#### Per-SGB logistic regression ####
# Outcome:   shared (0/1) — was this strain retained between timepoints?
# Predictor: EthnicityTot (South-Asian Surinamese vs Dutch as reference)
# Adjustment: Age (z-scored), FUtime (z-scored)
# Model:     logistic regression, one model per SGB

message("Running per-SGB logistic regressions: shared ~ EthnicityTot + Age_z + FUtime_z ...")
persgb_results <- purrr::map_dfr(sgbs_use, function(sgb) {
    df_sgb <- dfsame_long %>%
        filter(SGB == sgb, !is.na(shared), !is.na(EthnicityTot),
               !is.na(Age_z), !is.na(FUtime_z))

    if (nrow(df_sgb) < 30 || length(unique(df_sgb$shared)) < 2) {
        return(tibble(SGB = sgb, n = nrow(df_sgb)))
    }

    tryCatch({
        mod <- glm(shared ~ EthnicityTot + Age_z + FUtime_z, data = df_sgb, family = binomial())
        tidy(mod, conf.int = TRUE, exponentiate = TRUE) %>%
            filter(str_detect(term, "EthnicityTot")) %>%
            mutate(term = str_remove(term, "EthnicityTot"),
                   SGB  = sgb,
                   n    = nrow(df_sgb))
    }, error = function(e) {
        tibble(SGB = sgb, n = nrow(df_sgb))
    })
})

persgb_results <- persgb_results %>%
    filter(!is.na(p.value)) %>%
    mutate(qval = p.adjust(p.value, method = "BH"),
           sig  = case_when(qval < 0.001 ~ "***", qval < 0.01 ~ "**",
                            qval <= 0.05 ~ "*", TRUE ~ "")) %>%
    left_join(cladesplit_match, by = "SGB")

write.csv2(persgb_results,
           file.path(resultsfolder, "persgb_ethnicity.csv"),
           row.names = FALSE)

message("SGBs significant (FDR < 0.05): ", sum(persgb_results$qval <= 0.05, na.rm = TRUE),
        " out of ", nrow(persgb_results))

#### Sanity check: chi-square test per SGB (ethnicity vs shared) ####
chisq_results <- purrr::map_dfr(sgbs_use, function(sgb) {
    df_sgb <- dfsame_long %>%
        filter(SGB == sgb, !is.na(shared),
               EthnicityTot %in% c("Dutch", "South-Asian Surinamese"))

    if (nrow(df_sgb) < 10) return(NULL)

    tbl <- table(df_sgb$EthnicityTot, df_sgb$shared)
    if (any(dim(tbl) < 2)) return(NULL)

    tryCatch({
        ct <- chisq.test(tbl)
        tibble(
            SGB     = sgb,
            n       = nrow(df_sgb),
            p.value = ct$p.value
        )
    }, error = function(e) NULL)
}) %>%
    left_join(cladesplit_match, by = "SGB") %>%
    mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        sig   = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
                          p.adj <= 0.05 ~ "*", TRUE ~ "")
    )

write.csv2(chisq_results,
           file.path(resultsfolder, "persgb_chisq_ethnicity.csv"),
           row.names = FALSE)

message("Chi-square significant (FDR < 0.05): ", sum(chisq_results$p.adj <= 0.05, na.rm = TRUE),
        " out of ", nrow(chisq_results))

# Compare overlap with logistic regression results
overlap <- inner_join(
    persgb_results %>% filter(qval <= 0.05) %>% dplyr::select(SGB, Species),
    chisq_results  %>% filter(p.adj <= 0.05) %>% dplyr::select(SGB),
    by = "SGB"
)
message("SGBs significant in both logistic regression and chi-square: ", nrow(overlap))

#### Plots ####

## Forest plot: OR for South-Asian Surinamese vs Dutch, per SGB
## Species ordered by OR; significant hits highlighted
df_plot <- persgb_results %>%
    filter(!is.na(estimate)) %>%
    mutate(Species = str_replace_all(Species, "_", " ")) %>%
    add_count(Species, name = "n_species") %>%
    mutate(Species = case_when(
        n_species > 1 ~ paste0(Species, " (SGB", SGB, ")"),
        TRUE          ~ Species
    )) %>%
    dplyr::select(-n_species) %>%
    mutate(Species = fct_reorder(Species, estimate))

pl_persgb <- ggplot(df_plot, aes(x = estimate, y = Species, color = qval <= 0.05)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                  linewidth = 0.5, width = 0.3) +
    geom_point(size = 2) +
    scale_color_manual(
        values = c("TRUE" = "#E18727FF", "FALSE" = "grey60"),
        labels = c("TRUE" = "FDR < 0.05", "FALSE" = "FDR >= 0.05"),
        name = NULL) +
    labs(x = "OR South-Asian Surinamese vs Dutch (95% CI)",
         y = NULL,
         title = "Ethnicity and per-SGB strain sharing",
         caption = "Logistic regression adjusted for Age and follow-up time. Outcome: strain retained (1) vs lost (0).") +
    theme_Publication() +
    theme(legend.position = "bottom")

ggsave(file.path(resultsfolder, "persgb_ethnicity_forest.pdf"),
       plot = pl_persgb,
       width = 7, height = max(5, nrow(df_plot) * 0.25),
       device = cairo_pdf)

message("Done. Results saved to: ", resultsfolder)
