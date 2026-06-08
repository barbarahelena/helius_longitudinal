## Strain sharing — ethnic differences: logistic regression-based SGB selection + attenuation
## 1. Logistic regression of shared ~ EthnicityTot for all SGBs (n > 50)
## 2. Select SGBs with FDR < 0.05 → dumbbell plot
## 3. Attenuation analysis: ethnicity OR before/after covariate adjustment
## 4. Per-covariate effects within ethnically-differential SGBs
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
resultsfolder <- "results/3_species_change/5_strain_stability/ethnic_covariates"
dir.create(resultsfolder, showWarnings = FALSE, recursive = TRUE)

#### Data ####
df <- rio::import("data/shotgun/strainsharing_merged.csv")
colnames(df) <- str_remove(colnames(df), "sharing_")

dfsame <- df %>%
    filter(str_remove(sampleid_1, "HELIBA_") == str_remove(sampleid_2, "HELIFU_"))

clin    <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
cladesplit <- readRDS("data/shotgun/shotgun_taxtable_sgb.RDS")

#### Build long dataset — all SGBs with n > 50 ####
dfsame_long_all <- dfsame %>%
    filter(str_detect(sampleid_1, "HELIBA_")) %>%
    dplyr::rename(sampleID = sampleid_1) %>%
    dplyr::select(sampleID, starts_with("t__")) %>%
    pivot_longer(cols = starts_with("t__"), names_to = "SGB_col", values_to = "shared") %>%
    mutate(SGB = str_remove(SGB_col, "t__")) %>%
    dplyr::select(-SGB_col)

sgb_n <- dfsame_long_all %>%
    group_by(SGB) %>%
    summarise(n = sum(!is.na(shared)), .groups = "drop") %>%
    filter(n > 50)
sgbs_use <- sgb_n$SGB
message(length(sgbs_use), " SGBs with n > 50")

# Match taxonomy
cladesplit_match <- cladesplit %>%
    mutate(SGB_num = str_extract(SGB, "[0-9]+")) %>%
    mutate(SGB_key = sgbs_use[match(SGB_num, str_extract(sgbs_use, "[0-9]+"))]) %>%
    filter(!is.na(SGB_key)) %>%
    dplyr::select(SGB = SGB_key, Species) %>%
    mutate(Species = str_replace_all(Species, "_", " ")) %>%
    distinct(SGB, .keep_all = TRUE)

dfsame_long <- dfsame_long_all %>%
    filter(SGB %in% sgbs_use) %>%
    left_join(clin, by = "sampleID") %>%
    left_join(cladesplit_match, by = "SGB") %>%
    mutate(
        EthnicityTot = relevel(factor(EthnicityTot), ref = "Dutch"),
        Alcohol      = relevel(factor(Alcohol),      ref = "No"),
        PPI          = relevel(factor(PPI),          ref = "No"),
        FUtime_z     = as.numeric(scale(FUtime))
    )

message("dfsame_long rows: ", nrow(dfsame_long),
        " | SGBs present: ", n_distinct(dfsame_long$SGB))

#### 1. Logistic regression: shared ~ EthnicityTot per SGB ####
eth_results <- purrr::map_dfr(sgbs_use, function(sgb) {
    df_sgb <- dfsame_long %>%
        filter(SGB == sgb, !is.na(shared), !is.na(EthnicityTot))
    if (nrow(df_sgb) < 20 || length(unique(df_sgb$shared)) < 2) return(NULL)
    tryCatch({
        m <- glm(shared ~ EthnicityTot, data = df_sgb, family = binomial())
        tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
            filter(str_detect(term, "EthnicityTot")) %>%
            slice(1) %>%
            mutate(SGB = sgb, n = nrow(df_sgb))
    }, error = function(e) NULL)
}) %>%
    left_join(cladesplit_match, by = "SGB") %>%
    mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        sig   = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR >= 0.05")
    )

write.csv2(eth_results,
           file.path(resultsfolder, "ethnicity_OR_allsgbs.csv"),
           row.names = FALSE)

message(sum(eth_results$p.adj < 0.05, na.rm = TRUE),
        " SGBs with significant ethnicity effect (FDR < 0.05)")

# Select ethnically-differential SGBs
sgbs_eth_sig <- eth_results %>%
    filter(p.adj < 0.05) %>%
    arrange(estimate) %>%
    mutate(Species = factor(Species, levels = unique(Species)))

#### 2. Dumbbell plot: per-ethnicity strain sharing % for significant SGBs ####
# Calculate per-ethnicity sharing proportions
dfsame_long_sig <- dfsame_long %>%
    filter(SGB %in% sgbs_eth_sig$SGB)

eth_sharing_perc <- dfsame_long_sig %>%
    filter(!is.na(shared), EthnicityTot %in% c("Dutch", "South-Asian Surinamese")) %>%
    group_by(SGB, Species, EthnicityTot) %>%
    summarise(
        sharing_perc = mean(shared, na.rm = TRUE) * 100,
        n            = sum(!is.na(shared)),
        .groups = "drop"
    ) %>%
    mutate(Species = factor(Species, levels = levels(sgbs_eth_sig$Species)))

leg <- c("Dutch" = pal_simpsons()(2)[2], "South-Asian Surinamese" = pal_simpsons()(1))

pl_dumbbell <- ggplot(eth_sharing_perc, aes(x = Species)) +
    geom_segment(data = eth_sharing_perc %>%
                     pivot_wider(id_cols = c(SGB, Species),
                                 names_from = EthnicityTot,
                                 values_from = sharing_perc),
                 aes(y = Dutch, yend = `South-Asian Surinamese`),
                 color = "darkgrey") +
    geom_point(aes(y = sharing_perc, color = EthnicityTot), size = 3) +
    scale_color_manual(values = leg, name = NULL) +
    labs(y = "% subjects with stable strain",
         x = NULL,
         title = "Ethnically-differential strain stability\n(FDR < 0.05)") +
    theme_Publication() +
    coord_flip()

ggsave(file.path(resultsfolder, "dumbbell_ethnic_sgbs.pdf"),
       plot = pl_dumbbell,
       width = 6, height = max(4, nrow(sgbs_eth_sig) * 0.4),
       device = cairo_pdf)

#### 3. Attenuation analysis: ethnicity OR before/after covariate adjustment ####
## Model 1: shared ~ EthnicityTot (unadjusted)
## Model 2: shared ~ EthnicityTot + FUtime_z + Alcohol + PPI (adjusted)

models_def <- list(
    "Unadjusted"                        = "shared ~ EthnicityTot",
    "Adjusted (FUtime + Alcohol + PPI)" = "shared ~ EthnicityTot + FUtime_z + Alcohol + PPI"
)
adj_vars <- c("FUtime_z", "Alcohol", "PPI")

attenuation_results <- purrr::map_dfr(sgbs_eth_sig$SGB, function(sgb) {
    purrr::map_dfr(names(models_def), function(mod_name) {
        df_sgb <- dfsame_long %>%
            filter(SGB == sgb, !is.na(shared), !is.na(EthnicityTot))
        if (mod_name != "Unadjusted") {
            df_sgb <- df_sgb %>%
                filter(complete.cases(dplyr::select(., all_of(adj_vars))))
        }
        if (nrow(df_sgb) < 20 || length(unique(df_sgb$shared)) < 2) return(NULL)
        tryCatch({
            m <- glm(as.formula(models_def[[mod_name]]), data = df_sgb, family = binomial())
            tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
                filter(str_detect(term, "EthnicityTot")) %>%
                slice(1) %>%
                mutate(SGB = sgb, model = mod_name, n = nrow(df_sgb))
        }, error = function(e) NULL)
    })
}) %>%
    left_join(cladesplit_match, by = "SGB") %>%
    mutate(
        p.adj   = p.adjust(p.value, method = "BH"),
        sig     = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR >= 0.05"),
        model   = factor(model, levels = names(models_def)),
        Species = factor(Species, levels = levels(sgbs_eth_sig$Species))
    )

write.csv2(attenuation_results,
           file.path(resultsfolder, "attenuation_ethnicity_OR.csv"),
           row.names = FALSE)

n_sig_unadj <- attenuation_results %>%
    filter(model == "Unadjusted", p.adj < 0.05) %>%
    nrow()
n_sig_adj <- attenuation_results %>%
    filter(model == "Adjusted (FUtime + Alcohol + PPI)", p.adj < 0.05) %>%
    nrow()
message("SGBs significant after attenuation: unadjusted = ", n_sig_unadj,
        ", adjusted = ", n_sig_adj)

#### 4. Per-covariate effects within ethnically-differential SGBs ####
## Outcome: binary shared (0/1); adjusted for EthnicityTot
covars <- data.frame(
    var   = c("FUtime_z", "Alcohol", "PPI"),
    label = c("Follow-up time (per SD)", "Alcohol use", "PPI use"),
    stringsAsFactors = FALSE
)

percovar_results <- purrr::map_dfr(seq_len(nrow(covars)), function(i) {
    var_name  <- covars$var[i]
    var_label <- covars$label[i]
    purrr::map_dfr(sgbs_eth_sig$SGB, function(sgb) {
        df_sgb <- dfsame_long %>%
            filter(SGB == sgb, !is.na(shared), !is.na(.data[[var_name]]),
                   !is.na(EthnicityTot)) %>%
            mutate(predictor = .data[[var_name]])
        if (nrow(df_sgb) < 20 || length(unique(df_sgb$shared)) < 2) return(NULL)
        tryCatch({
            m <- glm(shared ~ predictor + EthnicityTot, data = df_sgb, family = binomial())
            tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
                filter(term == "predictor" | str_detect(term, "^predictor")) %>%
                slice(1) %>%
                mutate(SGB = sgb, label = var_label, n = nrow(df_sgb))
        }, error = function(e) NULL)
    })
}) %>%
    left_join(cladesplit_match, by = "SGB") %>%
    mutate(
        p.adj   = p.adjust(p.value, method = "BH"),
        sig     = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR >= 0.05"),
        label   = factor(label, levels = covars$label),
        Species = factor(Species, levels = levels(sgbs_eth_sig$Species))
    )

write.csv2(percovar_results,
           file.path(resultsfolder, "percovar_ethnic_sgbs.csv"),
           row.names = FALSE)

#### Plots ####
n_sig    <- nrow(sgbs_eth_sig)
pl_h     <- max(4, n_sig * 0.4)

## Plot: Attenuation forest — ethnicity OR before/after adjustment
pl_atten <- ggplot(attenuation_results,
                   aes(x = estimate, y = Species, color = model)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                  position = position_dodge(width = 0.5),
                  linewidth = 0.6, width = 0.3) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    scale_color_manual(
        values = c("Unadjusted"                        = "grey50",
                   "Adjusted (FUtime + Alcohol + PPI)" = "#E18727FF"),
        name = NULL) +
    labs(x = "OR South-Asian Surinamese vs Dutch (95% CI)",
         y = NULL,
         title = "Ethnic differences in strain stability:\nattenuation by covariates") +
    theme_Publication() +
    theme(legend.position = "bottom")

ggsave(file.path(resultsfolder, "attenuation_ethnicity_forest.pdf"),
       plot = pl_atten, width = 7, height = pl_h, device = cairo_pdf)

## Plot: Per-covariate forest plots (adjusted for ethnicity), one panel per covariate
pl_list_cov <- lapply(levels(percovar_results$label), function(lbl) {
    df_plot <- percovar_results %>% filter(label == lbl, !is.na(estimate))
    ggplot(df_plot, aes(x = estimate, y = Species, color = sig)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
        geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                      linewidth = 0.5, width = 0.3) +
        geom_point(size = 2.5) +
        scale_color_manual(
            values = c("FDR < 0.05" = "#E18727FF", "FDR >= 0.05" = "grey55"),
            name = NULL) +
        labs(x = "OR (95% CI)", y = NULL, title = lbl) +
        theme_Publication() +
        theme(legend.position = "bottom",
              axis.text.y  = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y  = element_blank())
})

pl_covars_combined <- ggarrange(
    pl_atten,
    ggarrange(plotlist = pl_list_cov, ncol = length(pl_list_cov),
              common.legend = TRUE, legend = "bottom"),
    ncol = 2, widths = c(1.5, length(pl_list_cov))
)

ggsave(file.path(resultsfolder, "attenuation_percovar_combined.pdf"),
       plot = pl_covars_combined,
       width = 5 + 4 * length(pl_list_cov), height = pl_h,
       device = cairo_pdf)

message("Done. Results saved to: ", resultsfolder)
