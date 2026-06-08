## Strain sharing — figure panels G, H, I, J, K
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggridges)

## theme
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
resultsfolder <- "results/3_species_change/5_strain_stability"
dir.create(resultsfolder, showWarnings = FALSE, recursive = TRUE)

#### Data ####
df <- rio::import("data/shotgun/strainsharing_merged.csv")
thres <- rio::import("data/shotgun/thresholds_merged.csv") %>%
    mutate(
        across(c("n_markers", "n_samples", "aln_length", "avg_gap_prop",
                 "threshold_value", "max_youden", "false_positive_rate", "false_negative_rate"),
               transfnum))
thres$n_markers <- NULL # bug in pipeline: n_samples = n_markers, n_samples not extracted from info
thres$n_markers <- thres$n_samples
thres$n_samples <- NULL
abundance <- rio::import("data/shotgun/combined_table_fixedlab.tsv")

#### Calculation strain sharing metric ####
colnames(df) <- str_remove(colnames(df), "sharing_")
sharing_sgb <- apply(df[,3:ncol(df)], 2, function(x) sum(x, na.rm = TRUE))
sharing_sum <- apply(df[,3:ncol(df)], 1, function(x) sum(x, na.rm = TRUE))
df$sharing_sum <- sharing_sum
strain_total <- apply(df[,3:ncol(df)], 1, function(x) sum(!is.na(x), na.rm = TRUE))
df$strain_total <- strain_total
df$sharing_perc <- (sharing_sum / strain_total) * 100

# Coverage by SGBs in strain sharing set
sgbs <- str_remove(names(df)[which(str_detect(names(df),"t__"))], "sharing_")
abundance <- abundance[,which(colnames(abundance) != "HELIBA_103370")] # only NAs
abundance <- abundance[,which(colnames(abundance) != "HELIFU_103370")] # not paired
abundance2 <- abundance %>% filter(str_detect(clade_name, "t__"))
colnames(abundance2) <- c(colnames(abundance2)[1], str_replace(colnames(abundance2)[2:ncol(abundance2)], "_T1", ""))
mean(colSums(abundance2[2:ncol(abundance2)])); sd(colSums(abundance2[2:ncol(abundance2)]))

# Tax
clade <- abundance2$clade_name
cladesplit <- str_split(clade, "\\|", n = 8, simplify = TRUE)
cladesplit <- as.data.frame(cladesplit)
colnames(cladesplit) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB")
cladesplit <- cladesplit %>% mutate(across(everything(.), ~str_remove_all(.x, "[a-z]__")))
cladesplit$rowname <- clade
rownames(abundance2) <- cladesplit$SGB[match(cladesplit$rowname, abundance2$clade_name)]
saveRDS(cladesplit, "data/shotgun/shotgun_taxtable_sgb.RDS")

# Select SGBs in strain sharing set
abundance3 <- abundance2 %>%
    mutate(clade_name = str_extract(clade_name, pattern = "(t__)[A-z]*[0-9]*[a-z_]*")) %>%
    filter(clade_name %in% sgbs)
mean(colSums(abundance3[2:ncol(abundance3)])); sd(colSums(abundance3[2:ncol(abundance3)]))
abundance3$average <- rowSums(abundance3[2:ncol(abundance3)])/(ncol(abundance3)-1)

# Paired vs different sample sets
dfsame <- df %>% filter(str_remove(sampleid_1, "HELIBA_") == str_remove(sampleid_2, "HELIFU_"))
dfdiff <- df %>% filter(str_remove(sampleid_1, "HELIBA_") != str_remove(sampleid_2, "HELIFU_"))

#### Merge with clinical data ####
clin <- readRDS("data/clinicaldata/clinicaldata_long.RDS")
dfsh <- dfsame %>%
    filter(str_detect(sampleid_1, "HELIBA_")) %>%
    dplyr::select(sampleID = sampleid_1, sharing_perc)
dftot <- inner_join(clin, dfsh, by = "sampleID")

bray_alphadiv <- readRDS("data/shotgun/alphabetadiversity_shotgun.RDS") %>%
    dplyr::select(sampleID, shannon, shannon_delta, richness, richness_delta, distance)
dftot <- dftot %>% left_join(., bray_alphadiv)

#### Figure panels ####

## pl_fig3_J: strain sharing by ethnicity
(pl_fig3_J <- ggplot(data = dftot, aes(x = EthnicityTot, y = sharing_perc)) +
    geom_violin(colour = NA, aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_manual(values = rev(pal_simpsons()(2)), guide = "none") +
    labs(y = "% of stable strains", title = "Strain stability\nbetween timepoints", x = "") +
    theme_Publication())
ggsave("results/3_species_change/5_strain_stability/ethnicities.pdf", width = 4.5, height = 5)

## p-value formatter for stat_cor: shows p < 0.0001 instead of full scientific notation
stat_cor_fmt <- function(...) {
    stat_cor(aes(label = after_stat(paste0(
        "R = ", round(r, 2),
        ", p", ifelse(p < 0.0001, " < 0.0001", paste0(" = ", round(p, 4)))
    ))), output.type = "text", ...)
}

## pl_fig3_H: follow-up time vs strain stability
(pl_fig3_H <- ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = FUtime, y = sharing_perc)) +
    geom_jitter(color = "#197EC0FF", alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_y_continuous(expand = expansion(add = c(2, 15))) +
    labs(y = "% of stable strains", x = "Follow-up time (years)", title = "Follow-up time and\nstrain stability") +
    stat_cor_fmt(label.y = 100) +
    theme_Publication())
ggsave("results/3_species_change/5_strain_stability/futime.pdf", width = 4.5, height = 5)

## pl_fig3_I: Bray-Curtis vs strain stability
(pl_fig3_I <- ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = distance, y = sharing_perc)) +
    geom_jitter(color = "#197EC0FF", alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_y_continuous(expand = expansion(add = c(2, 15))) +
    labs(y = "% of stable strains", x = "Bray-Curtis dissimilarity",
         title = "Strain stability and\nmicrobiome change") +
    stat_cor_fmt(label.y = 100) +
    theme_Publication())
ggsave("results/3_species_change/5_strain_stability/bray_strainsharing.pdf", width = 4.5, height = 5)

## pl_fig3_G: baseline Shannon vs strain stability
(pl_fig3_G <- ggplot(data = dftot %>% filter(!is.na(shannon)), aes(x = shannon, y = sharing_perc)) +
    geom_jitter(color = "#197EC0FF", alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_y_continuous(expand = expansion(add = c(2, 15))) +
    labs(y = "% of stable strains", x = "Shannon index (baseline)", title = "Baseline diversity and\nstrain stability") +
    stat_cor_fmt(label.y = 100) +
    theme_Publication())
ggsave("results/3_species_change/5_strain_stability/shannon_strainsharing.pdf", width = 4.5, height = 5)

#### Taxonomy and per-SGB strain sharing ####
tax <- cladesplit %>% filter(SGB %in% rownames(abundance3))
tax <- tax[match(tax$SGB, rownames(abundance3)),]
tax$abundance <- rowSums(abundance3[2:ncol(abundance3)])/(ncol(abundance3)-1)

othernames <- names(dfsame)[which(!str_detect(names(dfsame), "t__"))]
dfsame <- dfsame[,c(othernames, str_c("t__", tax$SGB))]
tax$sharing_sum <- apply(dfsame[,6:ncol(dfsame)], 2, function(x) sum(x, na.rm = TRUE))
dfsame2 <- dfsame[,c(othernames, str_c("t__", tax$SGB))]
tax$sharing_perc <- apply(dfsame2[,6:ncol(dfsame2)], 2, function(x) (sum(x, na.rm = TRUE) /
                                (sum(!is.na(x), na.rm = TRUE))) * 100)
tax$n <- apply(dfsame2[,6:ncol(dfsame2)], 2, function(x) (sum(!is.na(x))))
tax <- tax %>% right_join(., thres, by = "SGB") %>% filter(n > 50)

# Calculate per-ethnicity params
idsdutch <- dftot$sampleID[which(dftot$EthnicityTot == "Dutch")]
dfsamedutch <- dfsame %>% filter(sampleid_1 %in% idsdutch)
dfsamedutch <- dfsamedutch[,c(othernames, str_c("t__", tax$SGB))]
taxdutch <- tax
taxdutch$sharing_sum <- apply(dfsamedutch[,6:ncol(dfsamedutch)], 2, function(x) sum(x, na.rm = TRUE))
dfsamedutch2 <- dfsamedutch[,c(othernames, str_c("t__", taxdutch$SGB))]
taxdutch$sharing_perc <- apply(dfsamedutch2[,6:ncol(dfsamedutch2)], 2,
                             function(x) (sum(x, na.rm = TRUE) / sum(!is.na(x))) * 100)
taxdutch$n <- apply(dfsamedutch2[,6:ncol(dfsamedutch2)], 2, function(x) sum(!is.na(x)))
taxdutch$EthnicityTot <- "Dutch"

idssas <- dftot$sampleID[which(dftot$EthnicityTot == "South-Asian Surinamese")]
dfsamesas <- dfsame %>% filter(sampleid_1 %in% idssas)
dfsamesas <- dfsamesas[,c(othernames, str_c("t__", tax$SGB))]
taxsas <- tax
taxsas$sharing_sum <- apply(dfsamesas[,6:ncol(dfsamesas)], 2, function(x) sum(x, na.rm = TRUE))
dfsamesas2 <- dfsamesas[,c(othernames, str_c("t__", tax$SGB))]
taxsas$sharing_perc <- apply(dfsamesas2[,6:ncol(dfsamesas2)], 2,
                          function(x) (sum(x, na.rm = TRUE) / sum(!is.na(x))) * 100)
taxsas$n <- apply(dfsamesas2[,6:ncol(dfsamesas2)], 2, function(x) sum(!is.na(x)))
taxsas$EthnicityTot <- "South-Asian Surinamese"
taxtot <- rbind(taxdutch, taxsas)

# Select top 10 most significant SGBs from per-SGB logistic regression (strainsharing_persgb.R)
persgb_results <- read.csv2("results/3_species_change/5_strain_stability/covariates/persgb_ethnicity.csv",
                             stringsAsFactors = FALSE)
top10_sgbs <- persgb_results %>%
    filter(!is.na(p.value)) %>%
    arrange(qval, p.value) %>%
    slice_head(n = 10) %>%
    pull(SGB)

difftrue <- taxtot %>%
    filter(SGB %in% top10_sgbs) %>%
    pivot_wider(id_cols = c(SGB, Species), names_from = EthnicityTot, values_from = sharing_perc) %>%
    mutate(highest = case_when(Dutch >= `South-Asian Surinamese` ~ Dutch,
                               Dutch < `South-Asian Surinamese` ~ `South-Asian Surinamese`),
           Species = fct_reorder(Species, .x = highest)) %>%
    arrange(Species) %>%
    mutate(SGB = as.factor(SGB), SGB = fct_inorder(SGB))

## pl_fig3_K component: dumbbell plot of differential strains
leg <- c("Dutch" = pal_simpsons()(2)[2], "South-Asian Surinamese" = pal_simpsons()(1))
(pl <- ggplot(data = difftrue, aes(x = Species)) +
    geom_segment(aes(y = Dutch, yend = `South-Asian Surinamese`), color = "darkgrey") +
    geom_point(aes(y = `South-Asian Surinamese`, color = "South-Asian Surinamese"), size = 3) +
    geom_point(aes(y = Dutch, color = "Dutch"), size = 3) +
    scale_color_manual(values = leg) +
    labs(y = "% subjects with stable strain",
         x = "",
         title = "Strain stability", color = "") +
    theme_Publication() +
    coord_flip())
ggsave("results/3_species_change/5_strain_stability/diffsgbs_ethnicity_connect.pdf", width = 7, height = 7)

## pl_fig3_K component: abundance density ridgeline
ab <- abundance3 %>% filter(clade_name %in% str_c("t__", difftrue$SGB)) %>%
    dplyr::select(contains("HELI"))
ab <- as.data.frame(t(as.matrix(ab)))
ab <- ab[,levels(difftrue$SGB)]
colnames(ab) <- difftrue$Species[match(difftrue$SGB, colnames(ab))]
ab <- ab %>% rownames_to_column(var = "sampleID") %>%
    right_join(dftot, by = "sampleID") %>%
    pivot_longer(., cols = 2:(nlevels(difftrue$SGB)+1), names_to = "Species", values_to = "abundance") %>%
    mutate(Species = fct_inorder(as.factor(Species)))

(dens <- ggplot(data = ab, aes(x = abundance + 0.01,
                              y = Species,
                              fill = EthnicityTot)) +
    ggridges::geom_density_ridges(alpha = 0.5, rel_min_height = 0.01) +
    scale_x_log10(limits = c(0.1, 26)) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "", fill = "", x = "log10(abundance)") +
    theme_Publication() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank()))

## pl_fig3_K component: n subjects per SGB and ethnicity
taxtrue <- taxtot %>% filter(SGB %in% difftrue$SGB) %>%
    right_join(., difftrue, by = c("SGB", "Species")) %>%
    mutate(SGB = fct_reorder(SGB, highest), Species = fct_reorder(Species, highest))
saveRDS(taxtrue, "data/shotgun/sharing_tax.RDS")

(pl4 <- ggplot(data = taxtrue,
              aes(x = Species, y = n, fill = EthnicityTot)) +
    geom_bar(stat = "identity") +
    scale_fill_simpsons(guide = "none") +
    coord_flip() +
    theme_void() +
    labs(y = "subjects", x = "") +
    theme_Publication() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank()))

## pl_fig3_K: combined dumbbell + density + n
pl_fig3_K <- ggarrange(pl, dens, pl4,
                       ncol = 3, widths = c(3, 1.2, 1.2),
                       common.legend = TRUE, legend = "bottom")
ggsave("results/3_species_change/5_strain_stability/diffstrains_complete.pdf",
       plot = pl_fig3_K, width = 15, height = 10)
