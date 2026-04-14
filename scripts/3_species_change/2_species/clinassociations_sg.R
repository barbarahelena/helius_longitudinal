  # Clinical associations with specific species
  # Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

  # Libraries
  library(tidyverse)
  library(Cairo)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggsci)
  library(ggpubr)

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
                  # legend.direction = "horizontal",
                  legend.key.size= unit(0.2, "cm"),
                  legend.spacing  = unit(0, "cm"),
                  # legend.title = element_text(face="italic"),
                  plot.margin=unit(c(10,5,5,5),"mm"),
                  strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                  strip.text = element_text(face="bold"),
                  plot.caption = element_text(size = rel(0.5), face = "italic")
          ))
      
  } 

  # Data
  df <- readRDS("data/clinicaldata/clinicaldata_wide.RDS")
  mb <- readRDS("data/shotgun/shotgun_abundance.RDS")
  sp <- rio::import("results/3_species_change/2_species/lmer/lmm_results.csv")
  sp <- sp |> filter(sigq != "")
  sp
  mb <- mb[,sp$mbname]
  dim(mb)
  head(mb)
  df <- df |> filter(sampleID_baseline %in% rownames(mb)| `sampleID_follow-up`  %in%  rownames(mb))

  mb <- as.data.frame(mb)
  mb$sampleID <- rownames(mb)
  head(mb)
  for (a in 1:(ncol(mb)-2)){
    microbe <- names(mb)[a]
    mb1 <- mb |> dplyr::select(sampleID, a) |>
                mutate(timepoint = case_when(str_detect(sampleID, "HELIBA") ~ "baseline",
                                            str_detect(sampleID, "HELIFU") ~ "follow-up"),
                        ID = str_c("S", str_remove(str_remove(sampleID, "HELIBA_"), "HELIFU_"))) |> 
                pivot_wider(id_cols = ID, names_from = c("timepoint"), values_from = microbe) |> 
                filter(!is.na(baseline) & !is.na(`follow-up`)) |> 
                filter(baseline != 0) |> 
                mutate(change_cat = case_when(baseline < `follow-up` ~ "increase", 
                                          `follow-up` < baseline ~ "decrease", 
                                          baseline == `follow-up` ~ "equal"), 
                        change_cat = as.factor(change_cat),
                      change = `follow-up` - baseline)
    summary(mb1$change)
    tot <- full_join(mb1, df, by = "ID")
    head(tot)
    names(tot)

    ## Change in these microbes and clinical changes
    pl <- ggplot(data = tot, aes(x = change, y = SBP_delta, color = Ethnicity)) +
            scale_color_jco() +
            geom_point(alpha = 0.5) +
            stat_cor(method = "spearman", color = "darkgrey") +
            geom_smooth(method = "lm", formula = y ~ x, color = "darkgrey") +
            theme_Publication() +
            labs(x = "Delta microbe", y = "Delta SBP", title = microbe)
    print(pl)

    pl <- ggplot(data = tot, aes(x = change, y = BMI_delta, color = Ethnicity)) +
            scale_color_jco() +
            geom_point(alpha = 0.5) +
            stat_cor(method = "spearman", color = "darkgrey") +
            geom_smooth(method = "lm", formula = y ~ x, color = "darkgrey") +
            theme_Publication() +
            labs(x = "Delta microbe", y = "Delta BMI", title = microbe)
    print(pl)

    pl <- ggplot(data = tot |> filter(!is.na(change_cat)), aes(x = change_cat, y = BMI_delta)) +
            scale_fill_jco() +
            geom_violin(colour = NA, aes(fill = change_cat), alpha = 0.5) +
            geom_boxplot(fill = "white") +
            stat_compare_means() +
            geom_smooth(method = "lm", formula = y ~ x, color = "darkgrey") +
            theme_Publication() +
            labs(x = "Change category", y = "Delta BMI", title = microbe)
    print(pl)

    pl <- ggplot(data = tot, aes(x = change, y = Trig_delta, color = Ethnicity)) +
            scale_color_jco() +
            geom_point(alpha = 0.5) +
            geom_smooth(method = "lm", formula = y ~ x, color = "darkgrey") +
            stat_cor(method = "spearman", color = "darkgrey") +
            theme_Publication() +
            labs(x = "Delta microbe", y = "Delta Triglycerides", title = microbe)
    print(pl)

    pl <- ggplot(data = tot, aes(x = change, y = HbA1c_delta, color = Ethnicity)) +
            scale_color_jco() +
            geom_point(alpha = 0.5) +
            stat_cor(method = "spearman", color = "darkgrey") +
            geom_smooth(method = "lm", formula = y ~ x, color = "darkgrey") +
            theme_Publication() +
            labs(x = "Delta microbe", y = "Delta HbA1c", title = microbe)
    print(pl)

    ## Baseline
    pl <- ggplot(data = tot, aes(x = baseline, y = SBP_delta, color = Ethnicity)) +
            scale_color_jco() +
            geom_point(alpha = 0.5) +
            stat_cor(method = "spearman", color = "darkgrey") +
            geom_smooth(method = "lm", formula = y ~ x, color = "darkgrey") +
            theme_Publication() +
            labs(x = "Baseline microbe", y = "Delta SBP", title = microbe)
    print(pl)

    pl <- ggplot(data = tot, aes(x = baseline, y = BMI_delta, color = Ethnicity)) +
            scale_color_jco() +
            geom_point(alpha = 0.5) +
            stat_cor(method = "spearman", color = "darkgrey") +
            geom_smooth(method = "lm", formula = y ~ x, color = "darkgrey") +
            theme_Publication() +
            labs(x = "Baseline microbe", y = "Delta BMI", title = microbe)
    print(pl)

    pl <- ggplot(data = tot, aes(x = baseline, y = Trig_delta, color = Ethnicity)) +
            scale_color_jco() +
            geom_point(alpha = 0.5) +
            geom_smooth(method = "lm", formula = y ~ x, color = "darkgrey") +
            stat_cor(method = "spearman", color = "darkgrey") +
            theme_Publication() +
            labs(x = "Baseline microbe", y = "Delta Triglycerides", title = microbe)
    print(pl)

    pl <- ggplot(data = tot, aes(x = baseline, y = HbA1c_delta, color = Ethnicity)) +
            scale_color_jco() +
            geom_point(alpha = 0.5) +
            stat_cor(method = "spearman", color = "darkgrey") +
            geom_smooth(method = "lm", formula = y ~ x, color = "darkgrey") +
            theme_Publication() +
            labs(x = "Baseline microbe", y = "Delta HbA1c", title = microbe)
    print(pl)

  }

  # heatmap of all correlations
  cor_pval <- function(x, y) {
    idx <- complete.cases(x, y)
    if (sum(idx) == 0) return(c("estimate.rho" = NA, p.value = NA))
    res <- suppressWarnings(cor.test(x[idx], y[idx], method = "spearman"))
    return(c(estimate = res$estimate, p.value = res$p.value))
  }

  # Initialize matrices
  cor_mat <- matrix(nrow = ncol(mb)-2, ncol = 4)
  pval_mat <- matrix(nrow = ncol(mb)-2, ncol = 4)
  rownames(cor_mat) <- rownames(pval_mat) <- names(mb)[1:(ncol(mb)-2)]
  colnames(cor_mat) <- colnames(pval_mat) <- c("baseline_SBP", "baseline_BMI", "baseline_Trig", "baseline_HbA1c")
  mb <- mb[rownames(mb) != "HELIFU_103370",]
  for (a in 1:(ncol(mb)-2)){
    microbe <- names(mb)[a]
    
    mb1 <- mb |>
      dplyr::select(sampleID, a) |>
      mutate(
        timepoint = case_when(
          str_detect(sampleID, "HELIBA") ~ "baseline",
          str_detect(sampleID, "HELIFU") ~ "follow-up"
        ),
        ID = str_c("S", str_remove(str_remove(sampleID, "HELIBA_"), "HELIFU_"))
      ) |>
      pivot_wider(id_cols = ID, names_from = timepoint, values_from = microbe) |>
      filter(!is.na(baseline) | !is.na(`follow-up`)) |>
      mutate(change = `follow-up` - baseline)
    
    tot <- inner_join(mb1, df, by = "ID")
    if(any(is.na(tot[[2]]))) print("There are NAs!")

    vars <- list(
      baseline_SBP = cor_pval(tot$baseline, tot$SBP_delta),
      baseline_BMI = cor_pval(tot$baseline, tot$BMI_delta),
      baseline_Trig = cor_pval(tot$baseline, tot$Trig_delta),
      baseline_HbA1c = cor_pval(tot$baseline, tot$HbA1c_delta)
    )
    
    cor_mat[microbe, ] <- sapply(vars, function(v) v[["estimate.rho"]])
    pval_mat[microbe, ] <- sapply(vars, function(v) v[["p.value"]])
  }

  col_fun <- circlize::colorRamp2(
    c(-0.3, 0, 0.3),
    c(pal_nejm()(6)[6], "white", pal_nejm()(3)[3])
  )

  colnames(cor_mat) <- c("ΔSBP", "ΔBMI", "ΔTriglycerides", "ΔHbA1c")

  # Heatmap
  heatmap_bugs <- Heatmap(
    as.matrix(cor_mat),
    name = "Spearman\nCorrelation",
    col = col_fun,
    rect_gp = gpar(col = "white", lwd = 2),
    na_col = "grey95",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    show_row_dend = FALSE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 12),
    row_names_rot = 0,
    column_names_rot = 45,
    show_heatmap_legend = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      pval <- pval_mat[i, j]
      if (!is.na(pval)) {
        sig <- ""
        if (pval < 0.001) sig <- "***"
        else if (pval < 0.01) sig <- "**"
        else if (pval < 0.05) sig <- "*"
        if (sig != "") {
          grid.text(sig, x, y, gp = gpar(fontsize = 16), vjust = 0.75)
        }
      }
    }
  )

  # Build legends and stack vertically
  lgd_cor <- Legend(
    col_fun = col_fun,
    title = "Spearman\nCorrelation",
    at = c(-0.3, -0.15, 0, 0.15, 0.3),
    labels = c("-0.30", "-0.15", "0", "0.15", "0.30")
  )

  lgd_sig <- Legend(
    pch = c("*", "**", "***"), type = "points",
    labels = c("p < 0.05", "p < 0.01", "p < 0.001"),
    legend_gp = gpar(fontsize = 10)
  )

  lgd_packed <- packLegend(lgd_cor, lgd_sig, direction = "vertical", gap = unit(4, "mm"))

  heatmap_gg <- as_ggplot(grid.grabExpr(
    draw(
      heatmap_bugs,
      annotation_legend_list = list(lgd_packed),
      padding = unit(c(5, 25, 5, 5), "mm") # bottom, left, top, right
    )
  ))
  pl_fig3_D <- heatmap_gg

  #### Heatmap stratified by ethnicity (Dutch vs South-Asian Surinamese) ####

  ethnicities <- c("Dutch", "South-Asian Surinamese")

  cor_list  <- setNames(lapply(ethnicities, function(e) {
    m <- matrix(nrow = ncol(mb)-2, ncol = 4)
    rownames(m) <- names(mb)[1:(ncol(mb)-2)]
    colnames(m) <- c("baseline_SBP", "baseline_BMI", "baseline_Trig", "baseline_HbA1c")
    m
  }), ethnicities)

  pval_list <- cor_list  # same structure

  for (a in 1:(ncol(mb)-2)) {
    microbe <- names(mb)[a]

    mb1 <- mb |>
      dplyr::select(sampleID, a) |>
      mutate(
        timepoint = case_when(
          str_detect(sampleID, "HELIBA") ~ "baseline",
          str_detect(sampleID, "HELIFU") ~ "follow-up"
        ),
        ID = str_c("S", str_remove(str_remove(sampleID, "HELIBA_"), "HELIFU_"))
      ) |>
      pivot_wider(id_cols = ID, names_from = timepoint, values_from = microbe) |>
      filter(!is.na(baseline) | !is.na(`follow-up`)) |>
      mutate(change = `follow-up` - baseline)

    tot <- inner_join(mb1, df, by = "ID")

    for (eth in ethnicities) {
      tot_eth <- tot |> filter(EthnicityTot == eth)

      vars <- list(
        baseline_SBP   = cor_pval(tot_eth$baseline, tot_eth$SBP_delta),
        baseline_BMI   = cor_pval(tot_eth$baseline, tot_eth$BMI_delta),
        baseline_Trig  = cor_pval(tot_eth$baseline, tot_eth$Trig_delta),
        baseline_HbA1c = cor_pval(tot_eth$baseline, tot_eth$HbA1c_delta)
      )

      cor_list[[eth]][microbe, ]  <- sapply(vars, function(v) v[["estimate.rho"]])
      pval_list[[eth]][microbe, ] <- sapply(vars, function(v) v[["p.value"]])
    }
  }

  # Rename columns for display
  for (eth in ethnicities) {
    colnames(cor_list[[eth]])  <- c("ΔSBP", "ΔBMI", "ΔTriglycerides", "ΔHbA1c")
    colnames(pval_list[[eth]]) <- c("ΔSBP", "ΔBMI", "ΔTriglycerides", "ΔHbA1c")
  }

  make_strat_heatmap <- function(eth, show_legend = FALSE) {
    cm <- cor_list[[eth]]
    pm <- pval_list[[eth]]
    Heatmap(
      as.matrix(cm),
      name = paste0("Spearman\nCorrelation (", eth, ")"),
      col = col_fun,
      rect_gp = gpar(col = "white", lwd = 2),
      na_col = "grey95",
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_side = "left",
      show_row_dend = FALSE,
      row_names_gp = gpar(fontsize = 10),
      column_names_gp = gpar(fontsize = 12),
      row_names_rot = 0,
      column_names_rot = 45,
      column_title = ifelse(eth == "South-Asian Surinamese", "South-Asian\nSurinamese", eth),
      column_title_gp = gpar(fontsize = 13, fontface = "bold"),
      show_heatmap_legend = show_legend,
      cell_fun = function(j, i, x, y, width, height, fill) {
        pval <- pm[i, j]
        if (!is.na(pval)) {
          sig <- ""
          if (pval < 0.001) sig <- "***"
          else if (pval < 0.01) sig <- "**"
          else if (pval < 0.05) sig <- "*"
          if (sig != "") grid.text(sig, x, y, gp = gpar(fontsize = 16), vjust = 0.75)
        }
      }
    )
  }

  ht_dutch <- make_strat_heatmap("Dutch")
  ht_sas   <- make_strat_heatmap("South-Asian Surinamese")

  pdf("results/3_species_change/2_species/heatmap_stratified_ethnicity.pdf",
      width = 10, height = 6)
  draw(
    ht_dutch + ht_sas,
    annotation_legend_list = list(lgd_packed),
    padding = unit(c(5, 25, 5, 5), "mm")
  )
  dev.off()

  #### Figure 3C — Alistipes putredinis × cardiometabolic change ####

  ali_base <- data.frame(
      sampleID = rownames(mb),
      Alistipes = log10(mb[, "Alistipes_putredinis"] + 0.01)
  ) %>%
      filter(str_detect(sampleID, "HELIBA")) %>%
      mutate(ID = str_c("S", str_remove(sampleID, "HELIBA_")))

  tot_ali <- inner_join(ali_base, df, by = "ID")

  p_bmi <- ggplot(tot_ali %>% filter(!is.na(BMI_delta)),
                  aes(x = Alistipes, y = BMI_delta)) +
      geom_point(alpha = 0.35, size = 1.5, color = "royalblue") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.15, linewidth = 0.9, color = "black") +
      stat_cor(size = 3, method = "spearman") +
      labs(x = "A. putredinis (log10)", y = "Delta BMI", title = "BMI") +
      theme_Publication()

  p_trig <- ggplot(tot_ali %>% filter(!is.na(Trig_delta)),
                  aes(x = Alistipes, y = Trig_delta)) +
      geom_point(alpha = 0.35, size = 1.5, color = "royalblue") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.15, linewidth = 0.9, color = "black") +
      stat_cor(size = 3, method = "spearman") +
      labs(x = "A. putredinis (log10)", y = "Delta Triglycerides", title = "Triglycerides") +
      theme_Publication()

  p_hba1c <- ggplot(tot_ali %>% filter(!is.na(HbA1c_delta)),
                    aes(x = Alistipes, y = HbA1c_delta)) +
      geom_point(alpha = 0.35, size = 1.5, color = "royalblue") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.15, linewidth = 0.9, color = "black") +
      stat_cor(size = 3, method = "spearman") +
      labs(x = "A. putredinis (log10)", y = "Delta HbA1c", title = "HbA1c") +
      theme_Publication()

  (pl_fig3_C <- ggarrange(p_bmi, p_trig, p_hba1c, ncol = 3))
  ggsave(pl_fig3_C, filename = "results/3_species_change/2_species/alistipes_clinical_scatter.pdf",
        width = 12, height = 4.5)
