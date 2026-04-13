# Clinical associations with 16S species/genera
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
mbs <- readRDS("data/16s/selectedspecies.RDS")
mb <- mbs
head(mb)

for (a in 1:(ncol(mb)-4)){
  microbe <- names(mb)[a]
  mb1 <- mb |> select(sampleID, all_of(a)) |>
              mutate(timepoint = case_when(str_detect(sampleID, "HELIBA") ~ "baseline",
                                          str_detect(sampleID, "HELIFU") ~ "follow-up"),
                      ID = str_c("S", str_remove(str_remove(sampleID, "HELIBA_"), "HELIFU_"))) |> 
              pivot_wider(id_cols = ID, names_from = c("timepoint"), values_from = microbe) |> 
              mutate(change = `follow-up` - baseline)
  summary(mb1$change)
  tot <- full_join(mb1, df, by = "ID")
  head(tot)
  print(names(tot))
  ## Change in these microbes and clinical changes
  pl <- ggplot(data = tot, aes(x = change, y = SBP_delta, color = EthnicityTot)) +
          scale_color_jco() +
          geom_point(alpha = 0.5) +
          stat_cor(method = "spearman", color = "darkgrey") +
          geom_smooth(formula = "y~x", method = "lm", color = "darkgrey") +
          theme_Publication() +
          labs(title = microbe)
  print(pl)

  pl <- ggplot(data = tot, aes(x = change, y = BMI_delta, color = EthnicityTot)) +
          scale_color_jco() +
          geom_point(alpha = 0.5) +
          stat_cor(method = "spearman", color = "darkgrey") +
          geom_smooth(formula = "y~x", method = "lm", color = "darkgrey") +
          theme_Publication() +
          labs(title = microbe)
  print(pl)

  pl <- ggplot(data = tot, aes(x = change, y = Trig_delta, color = EthnicityTot)) +
          scale_color_jco() +
          geom_point(alpha = 0.5) +
          geom_smooth(formula = "y~x", method = "lm", color = "darkgrey") +
          stat_cor(method = "spearman", color = "darkgrey") +
          theme_Publication() +
          labs(title = microbe)
  print(pl)

  pl <- ggplot(data = tot, aes(x = change, y = HbA1c_delta, color = EthnicityTot)) +
          scale_color_jco() +
          geom_point(alpha = 0.5) +
          stat_cor(method = "spearman", color = "darkgrey") +
          geom_smooth(formula = "y~x", method = "lm", color = "darkgrey") +
          theme_Publication() +
          labs(title = microbe)
  print(pl)
  
  ## Baseline
  pl <- ggplot(data = tot, aes(x = baseline, y = SBP_delta, color = EthnicityTot)) +
          scale_color_jco() +
          geom_point(alpha = 0.5) +
          stat_cor(method = "spearman", color = "darkgrey") +
          geom_smooth(formula = "y~x", method = "lm", color = "darkgrey") +
          theme_Publication() +
          labs(title = microbe)
  print(pl)

  pl <- ggplot(data = tot, aes(x = baseline, y = BMI_delta, color = EthnicityTot)) +
          scale_color_jco() +
          geom_point(alpha = 0.5) +
          stat_cor(method = "spearman", color = "darkgrey") +
          geom_smooth(formula = "y~x", method = "lm", color = "darkgrey") +
          theme_Publication() +
          labs(title = microbe)
  print(pl)

  pl <- ggplot(data = tot, aes(x = baseline, y = Trig_delta, color = EthnicityTot)) +
          scale_color_jco() +
          geom_point(alpha = 0.5) +
          geom_smooth(formula = "y~x", method = "lm", color = "darkgrey") +
          stat_cor(method = "spearman", color = "darkgrey") +
          theme_Publication() +
          labs(title = microbe)
  print(pl)

  pl <- ggplot(data = tot, aes(x = baseline, y = HbA1c_delta, color = EthnicityTot)) +
          scale_color_jco() +
          geom_point(alpha = 0.5) +
          stat_cor(method = "spearman", color = "darkgrey") +
          geom_smooth(formula = "y~x", method = "lm", color = "darkgrey") +
          theme_Publication() +
          labs(title = microbe)
  print(pl)

}

# heatmap of all correlations
# Initialize matrices
mb <- mbs
cor_mat <- matrix(nrow = ncol(mb)-4, ncol = 8)
pval_mat <- matrix(nrow = ncol(mb)-4, ncol = 8)
rownames(cor_mat) <- rownames(pval_mat) <- names(mb)[1:(ncol(mb)-4)]
colnames(cor_mat) <- colnames(pval_mat) <- c("change_SBP", "change_BMI", "change_Trig", "change_HbA1c",
                                             "baseline_SBP", "baseline_BMI", "baseline_Trig", "baseline_HbA1c")
#mb <- mb[rownames(mb) != "HELIFU_103370",]
for (a in 1:(ncol(mb)-4)){
  microbe <- names(mb)[a]
  mb1 <- mb |>
    select(sampleID, all_of(a)) |>
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

  # Functie om correlatie en p-waarde te berekenen
  cor_pval <- function(x, y) {
    idx <- complete.cases(x, y)
    if (sum(idx) == 0) return(c(estimate = NA, p.value = NA))
    
    res <- suppressWarnings(cor.test(x[idx], y[idx], method = "spearman"))
    return(c(estimate = res$estimate, p.value = res$p.value))
  }

  vars <- list(
    change_SBP = cor_pval(tot$change, tot$SBP_delta),
    change_BMI = cor_pval(tot$change, tot$BMI_delta),
    change_Trig = cor_pval(tot$change, tot$Trig_delta),
    change_HbA1c = cor_pval(tot$change, tot$HbA1c_delta),
    baseline_SBP = cor_pval(tot$baseline, tot$SBP_delta),
    baseline_BMI = cor_pval(tot$baseline, tot$BMI_delta),
    baseline_Trig = cor_pval(tot$baseline, tot$Trig_delta),
    baseline_HbA1c = cor_pval(tot$baseline, tot$HbA1c_delta)
  )
  
  cor_mat[microbe, ] <- sapply(vars, function(v) v[["estimate.rho"]])
  pval_mat[microbe, ] <- sapply(vars, function(v) v[["p.value"]])
}

col_fun <- circlize::colorRamp2(
  c(-0.5, 0, 0.5), 
  c(pal_nejm()(6)[6], "white", pal_nejm()(3)[3])
)

colnames(cor_mat) <- c(
  "dSBP", "dBMI", "dTriglycerides", "dHbA1c",
  "dSBP", "dBMI", "dTriglycerides", "dHbA1c"
)

group_labels <- rep(c("d microbe", "baseline microbe"), each = 4)

# Heatmap
heatmap_bugs <- Heatmap(
  as.matrix(cor_mat),
  name = "Spearman\nCorrelation",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 2),
  na_col = "grey95",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_split = factor(group_labels, levels = c("Δ microbe", "baseline microbe")),
  row_names_side = "left",
  show_row_dend = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  row_names_rot = 0,
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "Spearman\nCorrelation",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-0.50", "-0.25", "0", "0.25", "0.50")
  ),
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
heatmap_bugs

# Legenda voor significantie
lgd_sig <- Legend(
  pch = c("*", "**", "***"), type = "points",
  labels = c("p < 0.05", "p < 0.01", "p < 0.001"),
  legend_gp = gpar(fontsize = 10)
)

heatmap_grob <- grid.grabExpr(draw(
      heatmap_bugs, 
      annotation_legend_list = list(lgd_sig),
      padding = unit(c(15, 20, 10, 10), "mm"), # bottom, left, top, right padding
  )
)
ggsave("results/3_species_change/2_species/heatmap_spec_clinvars_16s.pdf", plot = as_ggplot(heatmap_grob),
        device = cairo_pdf, width = 9, height = 4.5)
