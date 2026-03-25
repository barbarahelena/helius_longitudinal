## Gutsmash Longitudinal Analysis with Linear Mixed Models
library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggrepel)

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
           axis.title.y = element_text(angle=90, vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(size = rel(0.7)),
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing  = unit(0, "cm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold"),
           plot.caption = element_text(size = rel(0.5), face = "italic"),
           plot.subtitle = element_text(size=8, hjust = 0.5, face = "italic")))
}

# Data import
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
dim(df_rel)
pathway_cols <- colnames(df_rel) # for later use
df_rel$sampleID <- rownames(df_rel)

clinical <- readRDS("data/clinicaldata_long.RDS")
dir.create("results/4_functional_change/gutsmash/longitudinal", showWarnings = FALSE, recursive = TRUE)
df_clin  <- df_rel |> left_join(clinical, by = "sampleID") |> droplevels()
df_clin <- df_clin |> filter(!is.na(EthnicityTot))
table(df_clin$EthnicityTot) # sample per ethnicity

# GENE FAMILY-LEVEL LMM ----
# Prepare abundance matrix: df (samples x gene families), rownames(df) are sample IDs
# Join with clinical by sampleID (rownames)
# LMMs for each gene family
statres <- data.frame()
for (gf in pathway_cols) {
    df_clin$mb <- log10(df_clin[[gf]] + 1)
    model1 <- lmer(mb ~ EthnicityTot * timepoint + (1|ID), data = df_clin)
    res <- summary(model1)
    confint_model1 <- confint(model1, method = "Wald")
    interaction_row <- grep("EthnicityTotSouth-Asian Surinamese:timepointfollow-up", rownames(res$coefficients))
    ci_row <- grep("EthnicityTotSouth-Asian Surinamese:timepointfollow-up", rownames(confint_model1))
    statres <- rbind(statres, data.frame(
        family = gf,
        estimate = res$coefficients[interaction_row, 1],
        conflow = ifelse(length(ci_row) > 0, confint_model1[ci_row, 1], NA),
        confhigh = ifelse(length(ci_row) > 0, confint_model1[ci_row, 2], NA),
        pval = res$coefficients[interaction_row, 5]
      ))
}
statres <- as.data.frame(statres)
statres <- statres %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv2(statres, "results/4_functional_change/gutsmash/longitudinal/lmm_ethnicity_timepoint_results.csv", row.names = FALSE)

# LONGITUDINAL PLOTS ----
statres_sig <- statres %>% filter(padj < 0.05) %>% arrange(pval)

plist <- list()
for(i in 1:min(nrow(statres_sig), 20)){
  nm <- statres_sig$family[i]
  df_clin$mb <- log10(df_clin[[nm]] + 1)

  df_means <- df_clin %>%
    group_by(EthnicityTot, timepoint) %>%
    summarise(mean = mean(mb, na.rm = TRUE), sd = sd(mb, na.rm = TRUE),
              n = n(), .groups = "drop")

  res_lmm <- statres_sig %>%
    filter(family == nm) %>%
    mutate(group1 = "baseline", group2 = "follow-up",
            sigq = paste0("p=", formatC(pval, format = "e", digits = 2))) %>%
    dplyr::select(-family)

  mbmax <- ifelse(max(df_clin$mb, na.rm = TRUE) < 0,
                  max(df_clin$mb, na.rm = TRUE) * 0.8,
                  max(df_clin$mb, na.rm = TRUE) * 1.2)
  mbstat <- ifelse(max(df_clin$mb, na.rm = TRUE) < 0,
                    max(df_clin$mb, na.rm = TRUE) * 0.9,
                    max(df_clin$mb, na.rm = TRUE) * 1.1)
  mbmin <- min(df_clin$mb, na.rm = TRUE)

  pl2 <- ggplot() +
    geom_line(data = df_clin, aes(x = timepoint, y = mb, color = EthnicityTot, group = ID),
              alpha = 0.15, linewidth = 0.5) +
    geom_point(data = df_clin, aes(x = timepoint, y = mb, color = EthnicityTot),
                alpha = 0.15, size = 0.8) +
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
    labs(x = "Timepoint", y = "log10(abundance % + 0.01)",
          title = paste0(nm, "\n", statres_sig$subclass[i]),
          color = "")

  plist[[i]] <- pl2
}

n_plots <- length(plist)
n_cols <- 3
n_rows <- ceiling(n_plots / n_cols)

plots <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
                    labels = LETTERS[1:n_plots], nrow = n_rows, ncol = n_cols)
plots
ggsave("results/4_functional_change/gutsmash/longitudinal/significant_gutsmash_timepoint_ethnicity.pdf", plots,
        width = 12, height = 4 * n_rows)
