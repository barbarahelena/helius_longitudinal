## Strain sharing — overall covariate analysis
## Linear regression of sharing_perc against covariate blocks (Risk factors, Disease, Medication, Microbiota)
## Three-panel plot: forest plot + per-ethnicity companion + N bar
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

transfnum <- function(var) {
    var1 <- as.numeric(gsub(",", ".", gsub("\\.", "", var)))
    return(var1)
}

#### Output folder ####
resultsfolder <- "results/3_species_change/5_strain_stability/covariates"
dir.create(resultsfolder, showWarnings = FALSE, recursive = TRUE)

#### Data ####
df <- rio::import("data/shotgun/strainsharing_merged.csv")
thres <- rio::import("data/shotgun/thresholds_merged.csv") %>%
    mutate(
        across(c("n_markers", "n_samples", "aln_length", "avg_gap_prop",
                 "threshold_value", "max_youden", "false_positive_rate", "false_negative_rate"),
               transfnum))
thres$n_markers <- NULL
thres$n_markers <- thres$n_samples
thres$n_samples <- NULL

colnames(df) <- str_remove(colnames(df), "sharing_")
sharing_sum <- apply(df[,3:ncol(df)], 1, function(x) sum(x, na.rm = TRUE))
df$sharing_sum <- sharing_sum
strain_total <- apply(df[,3:ncol(df)], 1, function(x) sum(!is.na(x), na.rm = TRUE))
df$strain_total <- strain_total
df$sharing_perc <- (sharing_sum / strain_total) * 100

dfsame <- df %>% filter(str_remove(sampleid_1, "HELIBA_") == str_remove(sampleid_2, "HELIFU_"))
df <- NULL

#### Merge with clinical data ####
clin <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
dfsh <- dfsame %>%
    filter(str_detect(sampleid_1, "HELIBA_")) %>%
    dplyr::select(sampleID = sampleid_1, sharing_perc)
dftot <- inner_join(clin, dfsh, by = "sampleID")

bray_alphadiv <- readRDS("data/shotgun/alphabetadiversity_shotgun.RDS") %>%
    dplyr::select(sampleID, shannon, shannon_delta, richness, richness_delta, distance)
dftot <- dftot %>% left_join(., bray_alphadiv, by = "sampleID")

helius_wide <- readRDS("data/clinicaldata/clinicaldata_wide.RDS")
nutrient_adj <- helius_wide %>%
    dplyr::select(ID, TotalCalories_baseline,
                  Protein_baseline_adj, FattyAcids_baseline_adj,
                  Carbohydrates_baseline_adj, Fiber_baseline_adj, Sodium_g_baseline_adj)
dftot <- dftot %>% left_join(nutrient_adj, by = "ID")

#### 1. Overall strain sharing — covariate analysis ####
# Linear regression of sharing_perc per covariate block.
# Continuous predictors are z-scored so effect sizes are comparable.
# Adjustment: FUtime + EthnicityTot (dropped when FUtime is the test variable).

## Variable metadata — defines covariate blocks, labels, and types
var_meta_ss <- data.frame(
    var   = c("Age", "Sex", "BMI", "FUtime", "Alcohol", "Smoking_current",
              "DM", "HT_BPMed", "MetSyn", "Dyslipidemia",
              "Statins", "Metformin", "AntiHT", "PPI",
              "shannon", "richness", "distance",
              "TotalCalories_baseline",
              "Protein_baseline_adj", "FattyAcids_baseline_adj", "Carbohydrates_baseline_adj",
              "Fiber_baseline_adj", "Sodium_g_baseline_adj"),
    label = c("Age (per SD)", "Sex (female)", "BMI (per SD)", "Follow-up time (per SD)",
              "Alcohol use", "Current smoking",
              "Diabetes", "Hypertension", "MetSyn", "Dyslipidemia",
              "Statins", "Metformin", "Antihypertensives", "PPI",
              "Shannon index (per SD)", "Richness (per SD)", "Bray-Curtis dissimilarity (per SD)",
              "Total calories (per SD)",
              "Protein (per SD)", "Fatty acids (per SD)", "Carbohydrates (per SD)",
              "Fiber (per SD)", "Sodium (per SD)"),
    group = c(rep("Risk factors", 6),
              rep("Disease", 4),
              rep("Medication", 4),
              rep("Microbiota", 3),
              rep("Diet", 6)),
    type  = c("continuous", "binary", "continuous", "continuous", "binary", "binary",
              "binary", "binary", "binary", "binary",
              "binary", "binary", "binary", "binary",
              "continuous", "continuous", "continuous",
              rep("continuous", 6)),
    stringsAsFactors = FALSE
)

## Analysis dataset: z-score continuous predictors, relevel binary to "No" / "Male" / "Dutch"
binary_vars    <- var_meta_ss$var[var_meta_ss$type == "binary"]
binary_vars_no <- binary_vars[binary_vars != "Sex"]

dftot_scaled <- dftot %>%
    mutate(
        across(all_of(var_meta_ss$var[var_meta_ss$type == "continuous"]), ~ as.numeric(scale(.))),
        across(all_of(intersect(binary_vars_no, names(dftot))),           ~ relevel(factor(.), ref = "No")),
        Sex          = relevel(factor(Sex), ref = "Male"),
        EthnicityTot = relevel(factor(EthnicityTot), ref = "Dutch")
    )

## Helper: linear regression, adjusted for FUtime + EthnicityTot
extract_lm_ss <- function(data, predictor_var, label, group) {
    adj <- c("FUtime")
    adj <- adj[adj != predictor_var]
    rhs <- paste(c("predictor", adj), collapse = " + ")
    fml <- as.formula(paste("sharing_perc ~", rhs))

    df_use <- data %>%
        filter(!is.na(.data[[predictor_var]]), !is.na(sharing_perc)) %>%
        mutate(predictor = .data[[predictor_var]])
    if (length(adj) > 0) {
        df_use <- df_use %>% filter(complete.cases(dplyr::select(., all_of(adj))))
    }
    if (nrow(df_use) < 20) return(NULL)

    tryCatch({
        m  <- lm(fml, data = df_use)
        cf <- coef(summary(m))
        ci <- confint(m)
        data.frame(
            label = label, group = group,
            estimate  = cf[2, "Estimate"],
            conf.low  = ci[2, 1],
            conf.high = ci[2, 2],
            p.value   = cf[2, "Pr(>|t|)"],
            n         = nrow(df_use),
            stringsAsFactors = FALSE
        )
    }, error = function(e) NULL)
}

## Helper: per-ethnicity linear regression (unadjusted, stratified)
extract_lm_ss_eth <- function(data, predictor_var, label, group) {
    eth_list <- c("Dutch", "South-Asian Surinamese")
    purrr::map_dfr(eth_list, function(eth) {
        sub <- data %>%
            filter(EthnicityTot == eth, !is.na(.data[[predictor_var]]), !is.na(sharing_perc)) %>%
            mutate(predictor = .data[[predictor_var]])
        if (nrow(sub) < 10) return(NULL)
        if (is.factor(sub$predictor) && sum(sub$predictor != levels(sub$predictor)[1]) < 5) return(NULL)
        tryCatch({
            m  <- lm(sharing_perc ~ predictor, data = sub)
            cf <- coef(summary(m))
            ci <- confint(m)
            data.frame(
                label = label, group = group, ethnicity = eth,
                estimate  = cf[2, "Estimate"],
                conf.low  = ci[2, 1],
                conf.high = ci[2, 2],
                p.value   = cf[2, "Pr(>|t|)"],
                stringsAsFactors = FALSE
            )
        }, error = function(e) NULL)
    })
}

## Run main models
ss_effects <- purrr::map_dfr(seq_len(nrow(var_meta_ss)), function(i) {
    extract_lm_ss(dftot_scaled, var_meta_ss$var[i], var_meta_ss$label[i], var_meta_ss$group[i])
}) %>%
    mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        sig   = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR \u2265 0.05"),
        label = factor(label, levels = rev(var_meta_ss$label)),
        group = factor(group, levels = c("Risk factors", "Disease", "Medication", "Microbiota", "Diet"))
    )

write.csv2(ss_effects, file.path(resultsfolder, "overall_lm_covariates.csv"), row.names = FALSE)

## Run per-ethnicity models
ss_eth_effects <- purrr::map_dfr(seq_len(nrow(var_meta_ss)), function(i) {
    extract_lm_ss_eth(dftot_scaled, var_meta_ss$var[i], var_meta_ss$label[i], var_meta_ss$group[i])
}) %>%
    mutate(
        p.adj     = p.adjust(p.value, method = "BH"),
        sig       = ifelse(p.adj < 0.05, "FDR < 0.05", "FDR \u2265 0.05"),
        label     = factor(label, levels = levels(ss_effects$label)),
        group     = factor(group, levels = levels(ss_effects$group)),
        ethnicity = factor(ethnicity, levels = c("Dutch", "South-Asian Surinamese"))
    )

write.csv2(ss_eth_effects, file.path(resultsfolder, "overall_lm_covariates_ethnicity.csv"), row.names = FALSE)

## N bar data
n_total_ss <- nrow(dftot_scaled)
bar_data_ss <- purrr::map_dfr(seq_len(nrow(var_meta_ss)), function(i) {
    col  <- var_meta_ss$var[i]
    x    <- dftot[[col]]   # use original (unscaled) for counts
    if (col == "Sex") {
        n_female <- sum(as.character(x) == "Female", na.rm = TRUE)
        n_male   <- sum(as.character(x) == "Male",   na.rm = TRUE)
        n_tot    <- n_female + n_male
        data.frame(label = var_meta_ss$label[i], group = var_meta_ss$group[i],
                   category = c("Female", "Male"),
                   pct      = c(n_female / n_tot * 100, n_male / n_tot * 100),
                   n_label  = c(n_female, NA))
    } else if (var_meta_ss$type[i] == "binary") {
        n_yes <- sum(as.character(x) == "Yes", na.rm = TRUE)
        n_no  <- sum(as.character(x) == "No",  na.rm = TRUE)
        n_tot <- n_yes + n_no
        data.frame(label = var_meta_ss$label[i], group = var_meta_ss$group[i],
                   category = c("Yes", "No"),
                   pct      = c(n_yes / n_tot * 100, n_no / n_tot * 100),
                   n_label  = c(n_yes, NA))
    } else {
        n_valid <- sum(!is.na(x))
        data.frame(label = var_meta_ss$label[i], group = var_meta_ss$group[i],
                   category = "Non-missing",
                   pct      = n_valid / n_total_ss * 100,
                   n_label  = n_valid)
    }
}) %>%
    mutate(label    = factor(label,    levels = levels(ss_effects$label)),
           group    = factor(group,    levels = levels(ss_effects$group)),
           category = factor(category, levels = c("No", "Male", "Yes", "Female", "Non-missing")))

#### Plots — overall analysis ####
eth_colors_ss <- c("Dutch" = "#709AE1FF", "South-Asian Surinamese" = "#FED439FF")

## Main forest plot (ORs, faceted by covariate block)
pl_ss_main <- ggplot(ss_effects, aes(x = estimate, y = label, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                  orientation = "y", linewidth = 0.8, width = 0.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("FDR < 0.05" = "#F05C3BFF", "FDR \u2265 0.05" = "grey55"),
                       name = NULL) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "\u03b2 coefficient (95% CI)", y = NULL,
         title = "" ) +
    theme_Publication() +
    theme(legend.position = "bottom",
          strip.text.y = element_blank(),
          strip.background.y = element_blank())

## Per-ethnicity companion (dots per ethnicity, same y-axis)
pl_ss_eth <- ggplot(ss_eth_effects, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(fill = ethnicity), shape = 21, color = "black",
               size = 2.5, stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = eth_colors_ss, name = NULL) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "\u03b2", y = NULL, title = "Baseline predictors of strain stability") +
    theme_Publication() +
    theme(legend.position = "bottom",
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          strip.text.y = element_blank(),
          strip.background.y = element_blank())

## N bar (% Yes / % non-missing per predictor)
pl_ss_bar <- ggplot(bar_data_ss, aes(x = pct, y = label, fill = category)) +
    geom_col(width = 0.6) +
    geom_text(data = bar_data_ss %>% filter(!is.na(n_label)),
              aes(label = paste0("n=", n_label), x = 101),
              hjust = 0, size = 2.8, color = "grey30") +
    scale_fill_manual(values = c("Yes"    = "#6F99ADFF", "No"   = "#d0d8e4",
                                 "Female" = "#6F99ADFF", "Male" = "#d0d8e4",
                                 "Non-missing" = "#6F99ADFF"),
                      guide = "none") +
    scale_x_continuous(limits = c(0, 140), breaks = c(0, 50, 100),
                       expand = expansion(mult = c(0, 0))) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "% / % non-missing", y = NULL, title = " ") +
    theme_Publication() +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          legend.position = "bottom")

(pl_ss_combined <- pl_ss_main + pl_ss_eth + pl_ss_bar + plot_layout(widths = c(0.4,0.4,0.3)))
ggsave(file.path(resultsfolder, "overall_strainsharing_covariates.pdf"),
       plot = pl_ss_combined, width = 12, height = 10, device = cairo_pdf)

message("Done. Results saved to: ", resultsfolder)
