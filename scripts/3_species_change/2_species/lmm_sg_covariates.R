# LMMs - Covariate attenuation analysis
# For each species with significant Ethnicity×timepoint interaction,
# progressively add covariate blocks and track how the interaction estimate changes.
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)

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

# ── Data ───────────────────────────────────────────────────────────────────────
# Use pcdiet version for DietPC1/DietPC2 availability
df <- readRDS("data/clinicaldata_long_pcdiet.RDS")
mb <- readRDS("data/shotgun/shotgun_abundance.RDS")

otu <- mb[which(rownames(mb) %in% df$sampleID), ]
mb1 <- otu[str_detect(rownames(otu), "HELIBA"), ]
tk1 <- apply(mb1[, 2:ncol(mb1)], 2, function(x) sum(x > 0.1) > (0.20 * length(x)))
mb2 <- otu[str_detect(rownames(otu), "HELIFU"), ]
tk2 <- apply(mb2[, 2:ncol(mb2)], 2, function(x) sum(x > 0.1) > (0.20 * length(x)))
tk <- Reduce(`+`, list(tk1, tk2)) > 0
mb <- mb[, tk == TRUE]
mb <- apply(mb, 2, function(x) log10(x + 0.01))
mb <- as.data.frame(mb)
mb$sampleID <- rownames(mb)

# Statins is baseline-only — carry forward per person so it remains usable
# as a time-invariant covariate without dropping follow-up rows
df <- df %>%
    group_by(ID) %>%
    mutate(Statins = if_else(is.na(Statins),
                             Statins[timepoint == "baseline"][1],
                             Statins)) %>%
    ungroup()

df_tot <- left_join(mb, df, by = "sampleID")

# ── Significant species from main analysis ─────────────────────────────────────
statres <- read.csv2("results/3_species_change/2_species/lmer/lmm_results.csv")
sig_species <- statres %>%
    mutate(qval = as.numeric(as.character(qval))) %>%
    filter(qval <= 0.05) %>%
    pull(mbname) %>%
    unique()

# ── Cumulative covariate model specifications ──────────────────────────────────
# Each model adds one block on top of the previous, so attenuation is directly
# attributable to the newly added block.
# Note: Diet PCs (~50% missingness) will reduce n substantially in the last model.
model_specs <- list(
    "Base"            = "Ethnicity * timepoint",
    "+Age, Sex"       = "Ethnicity * timepoint + Age + Sex",
    "+BMI"            = "Ethnicity * timepoint + Age + Sex + BMI",
    "+Discrimination" = "Ethnicity * timepoint + Age + Sex + BMI + DiscrMean_baseline",
    "+Alcohol"        = "Ethnicity * timepoint + Age + Sex + BMI + DiscrMean_baseline + AlcCons",
    "+Medications"    = "Ethnicity * timepoint + Age + Sex + BMI + DiscrMean_baseline + AlcCons + Statins + Metformin + AntiHT + PPI",
    "+Diet PCs"       = "Ethnicity * timepoint + Age + Sex + BMI + DiscrMean_baseline + AlcCons + Statins + Metformin + AntiHT + PPI + DietPC1 + DietPC2"
)

# ── Run all models for each significant species ────────────────────────────────
results <- map_dfr(sig_species, function(sp) {
    df_tot$microbe <- df_tot[[sp]]

    map_dfr(names(model_specs), function(model_name) {
        formula_str <- paste0("microbe ~ ", model_specs[[model_name]], " + (1|ID)")
        tryCatch({
            model    <- lmer(as.formula(formula_str), data = df_tot)
            coefs    <- summary(model)$coefficients
            inter_idx <- grep("Ethnicity.*timepoint", rownames(coefs))
            if (length(inter_idx) == 0) stop("interaction term not found")

            est <- coefs[inter_idx, 1]
            se  <- coefs[inter_idx, 2]
            pv  <- coefs[inter_idx, 5]

            tibble(
                species  = sp,
                model    = model_name,
                estimate = est,
                conflow  = est - 1.96 * se,
                confhigh = est + 1.96 * se,
                pval     = pv,
                n_obs    = nobs(model)
            )
        }, error = function(e) {
            tibble(species = sp, model = model_name,
                   estimate = NA_real_, conflow = NA_real_, confhigh = NA_real_,
                   pval = NA_real_, n_obs = NA_integer_)
        })
    })
})

# ── Tidy for plotting ──────────────────────────────────────────────────────────
model_order <- names(model_specs)

results_plot <- results %>%
    mutate(
        model = factor(model, levels = rev(model_order)),
        sig   = case_when(
            pval < 0.001 ~ "***",
            pval < 0.01  ~ "**",
            pval <= 0.05 ~ "*",
            TRUE         ~ ""
        ),
        species_clean = str_replace_all(species, "_", " ")
    )

# ── Forest plots: one panel per species, rows = model specifications ───────────
jco_blue <- pal_jco()(1)

pl_list <- lapply(unique(results_plot$species_clean), function(sp_clean) {
    sp_data  <- results_plot %>% filter(species_clean == sp_clean)
    x_range  <- range(c(sp_data$conflow, sp_data$confhigh), na.rm = TRUE)
    label_x  <- x_range[2] + diff(x_range) * 0.08

    ggplot(sp_data, aes(x = estimate, y = model)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
        geom_errorbar(aes(xmin = conflow, xmax = confhigh),
                      height = 0.2, linewidth = 0.5, color = jco_blue) +
        geom_point(size = 2.5, color = jco_blue) +
        geom_text(aes(x = label_x, label = sig),
                  hjust = 0, size = 4, color = "black") +
        theme_Publication() +
        labs(
            x     = "Interaction estimate (± 95% CI)",
            y     = NULL,
            title = sp_clean
        )
})

n_sp   <- length(pl_list)
n_cols <- min(3, n_sp)
n_rows <- ceiling(n_sp / n_cols)

dir.create("results/3_species_change/2_species/lmer", recursive = TRUE, showWarnings = FALSE)

(plots <- ggarrange(plotlist = pl_list, ncol = n_cols, nrow = n_rows))
ggsave(plots,
       filename = "results/3_species_change/2_species/lmer/lmm_covariate_attenuation.pdf",
       width = 5 * n_cols, height = 4 * n_rows)

write.csv2(results_plot %>% dplyr::select(-species_clean),
           "results/3_species_change/2_species/lmer/lmm_covariate_attenuation.csv",
           row.names = FALSE)
