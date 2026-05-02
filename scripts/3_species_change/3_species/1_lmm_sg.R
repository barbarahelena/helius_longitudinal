# LMMs
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(afex)
library(aplot)

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
df <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
mb <- readRDS("data/shotgun/shotgun_abundance.RDS")
otu <- mb[which(rownames(mb) %in% df$sampleID),]
dutch_ids <- df %>% filter(Ethnicity == "Dutch") %>% pull(sampleID)
sas_ids   <- df %>% filter(Ethnicity == "Surinamese") %>% pull(sampleID)
mb1 <- otu[rownames(otu) %in% dutch_ids,]
tk1 <- apply(mb1, 2, function(x) sum(x > 0.1) > (0.20*length(x)))
mb2 <- otu[rownames(otu) %in% sas_ids,]
tk2 <- apply(mb2, 2, function(x) sum(x > 0.1) > (0.20*length(x)))
tk <- Reduce(`+`,list(tk1,tk2)) > 0
summary(tk)
mb <- mb[,tk==TRUE]
dim(mb)
mean(mb, na.rm = TRUE)[1:5]
mb <- apply(mb, 2, function(x) log10(x + 0.01))
mean(mb, na.rm = TRUE)[1:5]
mb <- as.data.frame(mb)
mb$sampleID <- rownames(mb)

# Metadata
df_tot <- left_join(mb, df, by = c("sampleID"))

statres <- c()
for(i in c(1:(ncol(mb)-1))) {
    df_tot$microbe <- df_tot[,i]
    mbname <- colnames(df_tot)[i]
    model1 <- lmer(microbe ~ Ethnicity*timepoint + FUtime + (1|ID), data = df_tot)
    res <- summary(model1)
    confint_model1 <- confint(model1)
    estimate <- as.numeric(format(round(res$coefficients[5,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model1[7,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model1[7,2], 3), nsmall = 3))
    pval <- format(round(res$coefficients[5,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    sig <- case_when(
        pval < 0.0001 ~ paste0("****"),
        pval < 0.001 ~paste0("***"),
        pval < 0.01 ~paste0("**"),
        pval <= 0.05 ~paste0("*"),
        pval > 0.05 ~paste0("")
    )
    statres_line <- cbind(mbname, group1 = "baseline", group2 = "follow-up", pval, 
                          sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
}

statres <- as.data.frame(statres)
statres <- statres %>% arrange(pval, group1) %>% 
    mutate(qval = p.adjust(pval, method = "fdr")) |> 
    mutate(sigq = case_when(
        qval < 0.0001 ~ paste0("****"),
        qval < 0.001 ~paste0("***"),
        qval < 0.01 ~paste0("**"),
        qval <= 0.05 ~paste0("*"),
        qval > 0.05 ~paste0("")
    ))
head(statres, n = 20)
maxsig <- statres %>% filter(qval <= 0.05) %>% filter(!duplicated(mbname))

plist <- list()
for(i in 1:nrow(maxsig)){
    nm <- maxsig$mbname[i]
    pval <- maxsig$pval[i]
    df_tot$mb <- df_tot[,maxsig$mbname[i]]
    df_means <- df_tot %>% group_by(Ethnicity, timepoint) %>% 
        summarise(mean = mean(mb), sd = sd(mb), n = length(mb), .groups = "drop_last")
    res_lmm <- statres %>% filter(mbname == nm) %>% dplyr::select(-mbname) %>% filter(sig != "")
    if(max(df_tot$mb) < 0) mbmax <- max(df_tot$mb*0.8) else mbmax <- max(df_tot$mb*1.2)
    if(max(df_tot$mb) < 0) mbstat <- max(df_tot$mb*0.9) else mbstat <- max(df_tot$mb*0.7)
    mbmin <- min(df_tot$mb)
    pl2 <- ggplot() +
        geom_line(data = df_tot, aes(x = timepoint, y = mb,
                  color = Ethnicity, group = ID), alpha = 0.05, linewidth = 0.5) +
        geom_point(data = df_tot, aes(x = timepoint, y = mb,
                  color = Ethnicity, group = Ethnicity), alpha = 0.05, size = 0.8) +
        geom_line(data = df_means, aes(x = timepoint, y = mean, 
                  color = Ethnicity, group = Ethnicity), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = timepoint, y = mean, 
                  color = Ethnicity, group = Ethnicity), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = mean - (sd/sqrt(n)),
                          ymax = mean + (sd/sqrt(n)),
                          x = timepoint,
                          color = Ethnicity), width=0.1) +
        stat_pvalue_manual(res_lmm, y.position = mbstat, label = "{sigq}", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5) +
        scale_color_jco() + 
        coord_cartesian(ylim = c(mbmin,mbmax)) +
        theme_Publication() +
        labs(x = "Timepoint", y = "log10(abundance+0.01)", title = nm, color = "")
        plist[[i]] <- pl2
}

dir.create("results/3_species_change/3_species/lmer", recursive = TRUE, showWarnings = FALSE)

(plots <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:8],
          nrow = 3, ncol = 3))
ggsave(plots, filename = "results/3_species_change/3_species/lmer/lmer_plots.pdf", width = 12, height = 13)
write.csv2(statres, "results/3_species_change/3_species/lmer/lmm_results.csv")

#### Figure 3B — Forest plot: species with significant ethnicity × timepoint interaction ####

lmm_sig <- statres %>%
    filter(sigq != "") %>%
    mutate(
        estimate = as.numeric(estimate),
        conflow  = as.numeric(conflow),
        confhigh = as.numeric(confhigh),
        mbname   = str_replace_all(mbname, "_", " "),
        mbname   = factor(mbname, levels = mbname[order(as.numeric(estimate))]),
        direction = ifelse(estimate > 0, "SAS more increase", "Dutch more increase")
    )

dir_colors <- c("Dutch more increase" = pal_jco()(2)[1], "SAS more increase" = pal_jco()(2)[2])

(pl_fig3_C <- ggplot(lmm_sig, aes(x = estimate, y = mbname, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conflow, xmax = confhigh), orientation = "y", linewidth = 0.5, width = 0.2) +
    geom_point(size = 3.5) +
    scale_color_manual(values = dir_colors, name = NULL) +
    labs(x = "Interaction effect (\u00b1 95% CI)",
         y = NULL,
         title = "Species changing differently\nby ethnicity over time") +
    theme_Publication() +
    theme(legend.position = "bottom", plot.title = element_text(size = rel(1.2))))

ggsave(pl_fig3_C, filename = "results/3_species_change/3_species/lmer/lmm_species_forest.pdf",
       width = 10, height = 10)

#### Baseline differential abundance between ethnicities ####

df_base <- df_tot %>% filter(str_detect(sampleID, "HELIBA"))

statres_base <- c()
for(i in c(1:(ncol(mb)-1))) {
    df_base$microbe <- df_base[, i]
    mbname <- colnames(df_tot)[i]
    model_base <- lm(microbe ~ Ethnicity, data = df_base)
    res <- summary(model_base)
    confint_base <- confint(model_base)
    estimate  <- as.numeric(format(round(res$coefficients[2, 1], 3), nsmall = 3))
    conflow   <- as.numeric(format(round(confint_base[2, 1], 3), nsmall = 3))
    confhigh  <- as.numeric(format(round(confint_base[2, 2], 3), nsmall = 3))
    pval      <- as.numeric(format(round(res$coefficients[2, 4], 3), nsmall = 3))
    statres_line <- cbind(mbname, pval, estimate, conflow, confhigh)
    statres_base <- rbind(statres_base, statres_line)
}

statres_base <- as.data.frame(statres_base) %>%
    mutate(across(c(pval, estimate, conflow, confhigh), as.numeric)) %>%
    arrange(pval) %>%
    mutate(
        qval = p.adjust(pval, method = "fdr"),
        sigq = case_when(
            qval < 0.0001 ~ "****",
            qval < 0.001  ~ "***",
            qval < 0.01   ~ "**",
            qval <= 0.05  ~ "*",
            qval > 0.05   ~ ""
        )
    )

write.csv2(statres_base, "results/3_species_change/3_species/lmer/lm_baseline_ethnicity_results.csv")

base_sig <- statres_base %>%
    filter(sigq != "") %>%
    mutate(
        mbname    = str_replace_all(mbname, "_", " "),
        mbname    = factor(mbname, levels = mbname[order(estimate)]),
        direction = ifelse(estimate > 0, "Higher in SAS", "Higher in Dutch")
    )

dir_colors_base <- c("Higher in Dutch" = pal_jco()(2)[1], "Higher in SAS" = pal_jco()(2)[2])

(pl_baseline_eth <- ggplot(base_sig, aes(x = estimate, y = mbname, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conflow, xmax = confhigh), orientation = "y", linewidth = 0.5, width = 0.2) +
    geom_point(size = 3.5) +
    scale_color_manual(values = dir_colors_base, name = NULL) +
    labs(x = "Ethnicity effect (± 95% CI)",
         y = NULL,
         title = "Species differing at baseline\nbetween ethnicities") +
    theme_Publication() +
    theme(legend.position = "bottom"))

ggsave(pl_baseline_eth, filename = "results/3_species_change/3_species/lmer/lm_baseline_ethnicity_forest.pdf",
       width = 11, height = 13)

#### Follow-up differential abundance between ethnicities ####

df_fu <- df_tot %>% filter(str_detect(sampleID, "HELIFU"))

statres_fu <- c()
for(i in c(1:(ncol(mb)-1))) {
    df_fu$microbe <- df_fu[, i]
    mbname <- colnames(df_tot)[i]
    model_fu <- lm(microbe ~ Ethnicity + FUtime, data = df_fu)
    res <- summary(model_fu)
    confint_fu <- confint(model_fu)
    estimate  <- as.numeric(format(round(res$coefficients[2, 1], 3), nsmall = 3))
    conflow   <- as.numeric(format(round(confint_fu[2, 1], 3), nsmall = 3))
    confhigh  <- as.numeric(format(round(confint_fu[2, 2], 3), nsmall = 3))
    pval      <- as.numeric(format(round(res$coefficients[2, 4], 3), nsmall = 3))
    statres_line <- cbind(mbname, pval, estimate, conflow, confhigh)
    statres_fu <- rbind(statres_fu, statres_line)
}

statres_fu <- as.data.frame(statres_fu) %>%
    mutate(across(c(pval, estimate, conflow, confhigh), as.numeric)) %>%
    arrange(pval) %>%
    mutate(
        qval = p.adjust(pval, method = "fdr"),
        sigq = case_when(
            qval < 0.0001 ~ "****",
            qval < 0.001  ~ "***",
            qval < 0.01   ~ "**",
            qval <= 0.05  ~ "*",
            qval > 0.05   ~ ""
        )
    )

write.csv2(statres_fu, "results/3_species_change/3_species/lmer/lm_followup_ethnicity_results.csv")

fu_sig <- statres_fu %>%
    filter(sigq != "") %>%
    mutate(
        mbname    = str_replace_all(mbname, "_", " "),
        mbname    = factor(mbname, levels = mbname[order(estimate)]),
        direction = ifelse(estimate > 0, "Higher in SAS", "Higher in Dutch")
    )
nrow(fu_sig)

(pl_fu_eth <- ggplot(fu_sig, aes(x = estimate, y = mbname, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(xmin = conflow, xmax = confhigh), orientation = "y", linewidth = 0.5, width = 0.2) +
    geom_point(size = 3.5) +
    scale_color_manual(values = dir_colors_base, name = NULL) +
    labs(x = "Ethnicity effect (± 95% CI)",
         y = NULL,
         title = "Species differing at follow-up\nbetween ethnicities") +
    theme_Publication() +
    theme(legend.position = "bottom"))

ggsave(pl_fu_eth, filename = "results/3_species_change/3_species/lmer/lm_followup_ethnicity_forest.pdf",
       width = 11, height = 13)

#### Heatmap: cross-sectional ethnic differences for LMM-significant species ####

lmm_sig_names_orig <- statres %>% filter(sigq != "") %>% pull(mbname)

heatmap_data <- bind_rows(
    statres_base %>%
        filter(mbname %in% lmm_sig_names_orig) %>%
        dplyr::select(mbname, estimate, qval) %>%
        mutate(timepoint = "baseline"),
    statres_fu %>%
        filter(mbname %in% lmm_sig_names_orig) %>%
        dplyr::select(mbname, estimate, qval) %>%
        mutate(timepoint = "follow-up")
) %>%
    mutate(
        diff_display = ifelse(qval < 0.05, estimate, NA_real_),
        star = case_when(
            qval < 0.001 ~ "***",
            qval < 0.01  ~ "**",
            qval < 0.05  ~ "*",
            TRUE ~ ""
        ),
        mbname_clean = str_replace_all(mbname, "_", " "),
        mbname_clean = factor(mbname_clean, levels = levels(lmm_sig$mbname)),
        timepoint    = factor(timepoint, levels = c("baseline", "follow-up"))
    )

abs_lim <- max(abs(heatmap_data$diff_display), na.rm = TRUE)
if (is.na(abs_lim) || abs_lim == 0) abs_lim <- 1

pl_heatmap_species <- ggplot(heatmap_data,
                              aes(x = timepoint, y = mbname_clean, fill = diff_display)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = star), color = "black", size = 2.5, vjust = 0.75) +
    scale_fill_gradient2(
        low      = "#2166AC",
        mid      = "white",
        high     = "#E6B800",
        na.value = "grey93",
        limits   = c(-abs_lim, abs_lim),
        name     = "Effect\n(SAS vs Dutch)",
        guide    = guide_colorbar(barheight = unit(6, "cm"), barwidth = unit(0.5, "cm"))
    ) +
    theme_Publication() +
    theme(
        axis.text.x  = element_text(angle = 45, hjust = 1),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position = "right"
    ) +
    labs(x = "", y = "", caption = "* q<0.05  ** q<0.01  *** q<0.001")

pl_forest_heatmap <- pl_fig3_C |> aplot::insert_right(pl_heatmap_species, width = 0.25)

combined_height <- max(4, nrow(lmm_sig) * 0.3 + 2)
cairo_pdf("results/3_species_change/3_species/lmer/forest_heatmap_species.pdf",
          width = 9, height = combined_height)
print(pl_forest_heatmap)
dev.off()

#### Overlap: baseline vs follow-up significant species ####

base_sig_names <- statres_base %>% filter(sigq != "") %>% pull(mbname)
fu_sig_names   <- statres_fu   %>% filter(sigq != "") %>% pull(mbname)

overlap_names  <- intersect(base_sig_names, fu_sig_names)
only_base      <- setdiff(base_sig_names, fu_sig_names)
only_fu        <- setdiff(fu_sig_names, base_sig_names)

cat("Significant at baseline only:    ", length(only_base), "\n")
cat("Significant at follow-up only:   ", length(only_fu), "\n")
cat("Significant at both timepoints:  ", length(overlap_names), "\n")
cat("\nOverlapping species:\n")
cat(paste0("  ", overlap_names), sep = "\n")

overlap_df <- bind_rows(
    statres_base %>% filter(mbname %in% overlap_names) %>% mutate(timepoint = "baseline"),
    statres_fu   %>% filter(mbname %in% overlap_names) %>% mutate(timepoint = "follow-up")
) %>%
    mutate(
        mbname    = str_replace_all(mbname, "_", " "),
        direction = ifelse(estimate > 0, "Higher in SAS", "Higher in Dutch")
    )
nrow(overlap_df)

write.csv2(overlap_df, "results/3_species_change/3_species/lmer/lm_overlap_baseline_followup.csv")

#### LMM interaction species × baseline/follow-up cross-sectional models ####

lmm_sig_names <- statres %>% filter(sigq != "") %>% pull(mbname)

in_base <- intersect(lmm_sig_names, base_sig_names)
in_fu   <- intersect(lmm_sig_names, fu_sig_names)
in_both <- intersect(in_base, in_fu)
in_none <- setdiff(lmm_sig_names, union(base_sig_names, fu_sig_names))

cat("\nOf the", length(lmm_sig_names), "LMM interaction-significant species:\n")
cat("  Also significant at baseline:         ", length(in_base),  "-", paste(in_base,  collapse = ", "), "\n")
cat("  Also significant at follow-up:        ", length(in_fu),    "-", paste(in_fu,    collapse = ", "), "\n")
cat("  Significant at both timepoints:       ", length(in_both),  "-", paste(in_both,  collapse = ", "), "\n")
cat("  Not significant at either timepoint:  ", length(in_none),  "-", paste(in_none,  collapse = ", "), "\n")
