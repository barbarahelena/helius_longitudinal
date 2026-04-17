## Strain sharing — exploratory / QC plots
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

#### Libraries ####
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggridges)

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
resultsfolder <- "results/3_species_change/3_strain_stability"
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

#### Strain sharing metric ####
colnames(df) <- str_remove(colnames(df), "sharing_")
sharing_sgb <- apply(df[,3:ncol(df)], 2, function(x) sum(x, na.rm = TRUE))
sharing_sum <- apply(df[,3:ncol(df)], 1, function(x) sum(x, na.rm = TRUE))
df$sharing_sum <- sharing_sum
strain_total <- apply(df[,3:ncol(df)], 1, function(x) sum(!is.na(x), na.rm = TRUE))
df$strain_total <- strain_total
df$sharing_perc <- (sharing_sum / strain_total) * 100

# SGBs in strain sharing set
sgbs <- str_remove(names(df)[which(str_detect(names(df),"t__"))], "sharing_")
abundance <- abundance[,which(colnames(abundance) != "HELIBA_103370")] # only NAs
abundance <- abundance[,which(colnames(abundance) != "HELIFU_103370")] # not paired
abundance2 <- abundance %>% filter(str_detect(clade_name, "t__"))
colnames(abundance2) <- c(colnames(abundance2)[1], str_replace(colnames(abundance2)[2:ncol(abundance2)], "_T1", ""))

# Taxonomy
clade <- abundance2$clade_name
cladesplit <- str_split(clade, "\\|", n = 8, simplify = TRUE)
cladesplit <- as.data.frame(cladesplit)
colnames(cladesplit) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB")
cladesplit <- cladesplit %>% mutate(across(everything(.), ~str_remove_all(.x, "[a-z]__")))
cladesplit$rowname <- clade
rownames(abundance2) <- cladesplit$SGB[match(cladesplit$rowname, abundance2$clade_name)]

# Select SGBs in strain sharing set
abundance3 <- abundance2 %>%
    mutate(clade_name = str_extract(clade_name, pattern = "(t__)[A-z]*[0-9]*[a-z_]*")) %>%
    filter(clade_name %in% sgbs)
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

# Per-ethnicity params
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

# Keep only SGBs with >20% prevalence in both ethnicities
n_dutch <- length(idsdutch)
n_sas   <- length(idssas)
sgbs_keep <- taxtot %>%
    pivot_wider(id_cols = SGB, names_from = EthnicityTot, values_from = n) %>%
    filter(Dutch / n_dutch > 0.2 & `South-Asian Surinamese` / n_sas > 0.2) %>%
    pull(SGB)
taxtot <- taxtot %>% filter(SGB %in% sgbs_keep)

# SGBs with >10% difference between ethnicities
difftax <- taxtot %>% pivot_wider(., id_cols = c(SGB, Species), names_from = c(EthnicityTot),
                                  values_from = sharing_perc) %>%
    mutate(diff = Dutch - `South-Asian Surinamese`,
           diff_bin = case_when(diff < -10 | diff > 10 ~ TRUE, .default = FALSE))
difftrue <- difftax %>% filter(diff_bin == TRUE) %>%
    mutate(highest = case_when(Dutch >= `South-Asian Surinamese` ~ Dutch,
                               Dutch < `South-Asian Surinamese` ~ `South-Asian Surinamese`),
           Species = fct_reorder(Species, .x = highest)) %>%
    arrange(Species) %>%
    mutate(SGB = as.factor(SGB), SGB = fct_inorder(SGB))

# Abundance data for differential SGBs
ab <- abundance3 %>% filter(clade_name %in% str_c("t__", difftrue$SGB)) %>%
    dplyr::select(contains("HELI"))
ab <- as.data.frame(t(as.matrix(ab)))
ab <- ab[,levels(difftrue$SGB)]
colnames(ab) <- difftrue$Species[match(difftrue$SGB, colnames(ab))]
ab <- ab %>% rownames_to_column(var = "sampleID") %>%
    right_join(dftot, by = "sampleID") %>%
    pivot_longer(., cols = 2:(nlevels(difftrue$SGB)+1), names_to = "Species", values_to = "abundance") %>%
    mutate(Species = fct_inorder(as.factor(Species)))

#### Distribution: paired vs different samples ####
gghistogram(dfsame$sharing_perc, fill = "royalblue", bins = 30) +
    geom_vline(aes(xintercept = median(dfsame$sharing_perc)), color = "firebrick", linewidth = 1) +
    labs(title = "Strain sharing - paired samples", x = "percentage of SGBs") +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/pairedsamples_hist.pdf", width = 4.5, height = 5)

gghistogram(dfdiff$sharing_perc, fill = "firebrick", bins = 30) +
    labs(title = "Strain sharing - different samples", x = "percentage of SGBs") +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/difsamples_hist.pdf", width = 4.5, height = 5)

mean(dfsame$sharing_perc); median(dfsame$sharing_perc)
mean(dfdiff$sharing_perc); median(dfdiff$sharing_perc)

#### Sex ####
ggplot(data = dftot, aes(x = fct_reorder(Sex, .x = sharing_perc, .fun = median), y = sharing_perc)) +
    geom_violin(colour = NA, aes(fill = Sex)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing between timepoints", x = "Sex") +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/sex.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(Sex)), aes(x = Sex, y = sharing_perc)) +
    geom_violin(colour = NA, aes(fill = Sex)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - sex differences", x = "Sex") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/sex_ethnicities.pdf", width = 6, height = 5)

#### Disease outcomes ####
dftot %>% group_by(DM_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc), .groups = "drop_last")
dftot %>% group_by(HT_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc), .groups = "drop_last")
dftot %>% group_by(MetSyn_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc), .groups = "drop_last")

ggplot(data = dftot %>% filter(!is.na(DM_new)), aes(x = DM_new, y = sharing_perc)) +
    geom_violin(colour = NA, aes(fill = DM_new)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - new DM", x = "") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/dmnew.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(MetSyn_new)), aes(x = MetSyn_new, y = sharing_perc)) +
    geom_violin(colour = NA, aes(fill = MetSyn_new)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - new MetSyn", x = "") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/metsynnew.pdf", width = 4.5, height = 5)

dftot %>% group_by(EthnicityTot, DM_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc), .groups = "drop_last")
ggplot(data = dftot %>% filter(!is.na(DM_new)), aes(x = DM_new, y = sharing_perc)) +
    geom_violin(colour = NA, aes(fill = DM_new)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - new DM", x = "New diabetes diagnosis") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/dmnew_ethnicities.pdf", width = 6, height = 5)

dftot %>% group_by(EthnicityTot, MetSyn_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc), .groups = "drop_last")
ggplot(data = dftot %>% filter(!is.na(MetSyn_new)), aes(x = MetSyn_new, y = sharing_perc)) +
    geom_violin(colour = NA, aes(fill = MetSyn_new)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - new MetSyn", x = "New MetSyn diagnosis") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/metsynnew_ethnicities.pdf", width = 6, height = 5)

#### Alpha diversity associations ####
ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = shannon_delta, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, height = 0) +
    geom_smooth(color = "black", method = "loess") +
    labs(y = "Percentage of stable strains", x = "Difference in Shannon index", title = "Shannon change and strain sharing") +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/shannonchange_strainsharing.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = richness, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, height = 0) +
    geom_smooth(color = "black", method = "lm") +
    stat_cor() +
    labs(y = "Percentage of stable strains", x = "Richness (baseline)", title = "Richness (baseline) and strain sharing") +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/richness_strainsharing.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = richness_delta, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, height = 0) +
    geom_smooth(color = "black", method = "loess") +
    labs(y = "Percentage of stable strains", x = "Difference in richness", title = "Richness change and strain sharing") +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/richnesschange_strainsharing.pdf", width = 4.5, height = 5)

#### FU time and age by ethnicity ####
ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = FUtime, y = sharing_perc)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    stat_cor() +
    labs(y = "Percentage of stable strains", x = "FU time (years)", title = "FU time and strain sharing") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/futime_ethnicity.pdf", width = 6, height = 5)

ggplot(data = dftot %>% filter(!is.na(Age)), aes(x = Age, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    stat_cor() +
    labs(y = "Percentage of stable strains", x = "Age (years)", title = "Baseline age and strain sharing") +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/age.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(Age)), aes(x = Age, y = sharing_perc)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    stat_cor() +
    labs(y = "Percentage of stable strains", x = "Age (years)", title = "Baseline age and strain sharing") +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave("results/3_species_change/3_strain_stability/age_ethnicity.pdf", width = 6, height = 5)

#### SGB-level QC ####
gghistogram(tax$sharing_perc, fill = "firebrick", bins = 30) +
    labs(title = "Strain sharing per SGB", x = "percentage of samples") +
    theme_Publication()

tax %>% arrange(sharing_perc) %>% dplyr::select(Family, Species, SGB, sharing_perc, n) %>% slice(1:20)
tax %>% arrange(-abundance) %>% dplyr::select(Species, SGB, abundance) %>% slice(1:10)

ggplot(data = tax, aes(x = n, y = sharing_perc)) +
    geom_jitter(color = "darkgreen", alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    stat_cor(label.x = 250) +
    labs(title = "Strain sharing percentage vs prevalence", x = "Number of subjects detected",
         y = "Sharing percentage (%)") +
    theme_Publication()

#### All SGBs by ethnicity ####
ggplot(data = taxtot %>% arrange(sharing_perc),
       aes(x = fct_reorder(Species, .x = sharing_perc), y = sharing_perc, group = EthnicityTot)) +
    geom_segment(aes(y = 0, yend = sharing_perc, color = EthnicityTot)) +
    geom_point(aes(color = EthnicityTot), stat = "identity") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of subjects with stable strains", x = "", title = "Strain sharing") +
    theme_Publication() +
    coord_flip()
ggsave("results/3_species_change/3_strain_stability/allsgbs_ethnicity.pdf", width = 10, height = 20)

#### Per-SGB density (differential strains) ####
plist <- c()
for(a in difftrue$Species) {
    plist[[a]] <- ggplot(data = ab %>% filter(Species == a),
                         aes(x = abundance + 0.01, fill = EthnicityTot)) +
        geom_density(alpha = 0.5) +
        scale_x_log10() +
        scale_fill_simpsons(guide = "none") +
        labs(y = "", fill = "", x = "") +
        theme_void()
}
pl_dens_qc <- ggarrange(plotlist = plist, ncol = 1, common.legend = TRUE)
ggsave("results/3_species_change/3_strain_stability/diffsgbs_eth_dens.pdf",
       plot = pl_dens_qc, width = 4, height = 20)

#### Threshold values per differential SGB ####
thres2 <- thres %>% filter(SGB %in% difftrue$SGB) %>% right_join(difftrue, ., by = "SGB")
ggplot(data = thres2, aes(x = Species, y = threshold_value)) +
    geom_bar(stat = "identity", fill = pal_simpsons()(3)[3]) +
    coord_flip() +
    labs(y = "threshold (ngd)", x = "") +
    theme_Publication()

#### nGD distributions per differential SGB ####
ngd <- rio::import("data/shotgun/ngd_merged.csv")
colnames(ngd) <- str_remove(colnames(ngd), "dist_t__")
ngd2 <- ngd[,c("sampleid_1", "sampleid_2", paste(difftrue$SGB))]
colnames(ngd2)[3:ncol(ngd2)] <- as.character(difftrue$Species[match(difftrue$SGB, colnames(ngd2)[3:ncol(ngd2)])])
ngd2 <- ngd2 %>% mutate(
    relation = case_when(
        str_remove(str_remove(sampleid_1, "HELIFU_"), "HELIBA_") ==
            str_remove(str_remove(sampleid_2, "HELIFU_"), "HELIBA_") ~ "same",
        .default = "different"
    )
)

plist2 <- c()
for(b in difftrue$Species) {
    ngd3 <- ngd2 %>% pivot_longer(., cols = all_of(b), names_to = "Species", values_to = "nGD") %>%
        mutate(nGD = transfnum(nGD)) %>%
        filter(Species == b)

    ngd4 <- ngd3 %>% filter(relation == "same") %>%
        mutate(sampleID = sampleid_1) %>%
        left_join(., dftot, by = "sampleID")

    plist2[[b]] <- ggplot(data = ngd4) +
        ggridges::geom_density_ridges(aes(x = nGD, y = Species, fill = EthnicityTot),
                                      alpha = 0.5, rel_min_height = 0.01, bandwidth = 0.005) +
        theme_Publication() +
        scale_fill_simpsons(guide = "none") +
        labs(x = "nGD", y = "", fill = "", title = b) +
        theme(axis.text.y = element_blank(),
              plot.title = element_text(size = rel(0.7)))
}

ngdpl <- ggarrange(plotlist = rev(plist2), ncol = 7, nrow = 6, common.legend = TRUE)
ggsave("results/3_species_change/3_strain_stability/ngd.pdf", plot = ngdpl, width = 18, height = 12)
