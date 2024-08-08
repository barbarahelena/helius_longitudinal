## Strain sharing
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)

## theme
theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
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

transfnum <- function(var) {
    var1 <- as.numeric(gsub(",", ".", gsub("\\.", "", var)))
    return(var1)
}

#### Output folder ####
resultsfolder <- "results/strainsharing"
dir.create(resultsfolder, showWarnings = FALSE)

#### Data ####
df <- rio::import("data/shotgun/strainsharing_merged.csv")
head(df)[1:5,1:5]
thres <- rio::import("data/shotgun/thresholds_merged.csv") %>% 
    mutate(
        across(c("n_markers", "n_samples", "aln_length", "avg_gap_prop",
                 "threshold_value", "max_youden", "false_positive_rate", "false_negative_rate"), 
               transfnum))
abundance <- rio::import("data/shotgun/combined_table_fixedlab.tsv")

#### Calculation strain sharing metric ####
colnames(df) <- str_remove(colnames(df), "sharing_")
sharing_sgb <- apply(df[,3:ncol(df)], 2, function(x) sum(x, na.rm = TRUE))
length(sharing_sgb[which(sharing_sgb == 0)])
# df <- df %>% dplyr::select(-which(sharing_sgb == 0))
sharing_sum <- apply(df[,3:ncol(df)], 1, function(x) sum(x, na.rm = TRUE))
sharing_sum[1:5]
df$sharing_sum <- sharing_sum
df$sharing_perc <- (sharing_sum / ncol(df[,which(str_detect(colnames(df), "t__"))])) * 100

# How much of microbiota composition is covered by these SGBs?
sgbs <- str_remove(names(df)[which(str_detect(names(df),"t__"))], "sharing_")
sgbs[1:5]
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
sgbs[which(!sgbs %in% abundance3$clade_name)]

#### Exploratory analyses ####
# How much strain sharing is there between same and different samples
dfsame <- df %>% filter(str_remove(sampleid_1, "HELIBA_") == str_remove(sampleid_2, "HELIFU_"))
gghistogram(dfsame$sharing_perc, fill = "royalblue") + 
    geom_vline(aes(xintercept = median(dfsame$sharing_perc)), color = "firebrick", size = 1) +
    labs(title = "Strain sharing - paired samples", x = "percentage of SGBs") +
    theme_Publication()
dfdiff <- df %>% filter(str_remove(sampleid_1, "HELIBA_") != str_remove(sampleid_2, "HELIFU_"))
gghistogram(dfdiff$sharing_perc, fill = "firebrick") + 
    labs(title = "Strain sharing - different samples", x = "percentage of SGBs") +
    theme_Publication()
mean(dfsame$sharing_perc); median(dfsame$sharing_perc)
mean(dfdiff$sharing_perc); median(dfdiff$sharing_perc)

#### Merge with clinical data ####
clin <- readRDS("data/clinicaldata_long.RDS")
dfsh <- dfsame %>% 
    filter(str_detect(sampleid_1, "HELIBA_")) %>% 
    dplyr::select(sampleID = sampleid_1, sharing_perc)
dim(dfsh)
dftot <- inner_join(clin, dfsh, by = "sampleID")
dim(dftot)

bray_alphadiv <- readRDS("data/shotgun/alphabetadiversity_shotgun.RDS") %>% 
    dplyr::select(sampleID, shannon, shannon_delta, richness, richness_delta, distance)
dftot <- dftot %>% left_join(., bray_alphadiv)

dftot %>% group_by(EthnicityTot) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc))
ggplot(data = dftot, aes(x = fct_reorder(EthnicityTot, .x = sharing_perc, .fun = median), y = sharing_perc)) +
    geom_violin(aes(fill = EthnicityTot)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing between timepoints", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() 
ggsave("results/strainsharing/ethnicities.pdf", width = 4.5, height = 5)

ggplot(data = dftot, aes(x = fct_reorder(Sex, .x = sharing_perc, .fun = median), y = sharing_perc)) +
    geom_violin(aes(fill = Sex)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing between timepoints", x = "Sex") +
    stat_compare_means(comparisons = list(c("Male", "Female")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() 
ggsave("results/strainsharing/sex.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(Sex)), aes(x = Sex, y = sharing_perc)) +
    geom_violin(aes(fill = Sex)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - sex differences", x = "Sex") +
    stat_compare_means(comparisons = list(c("Female", "Male")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    facet_wrap(~EthnicityTot) +
    theme_Publication() 
ggsave("results/strainsharing/sex_ethnicities.pdf", width = 6, height = 5)

dftot %>% group_by(DM_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc))
dftot %>% group_by(HT_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc))
dftot %>% group_by(MetSyn_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc))
ggplot(data = dftot %>% filter(!is.na(DM_new)), aes(x = DM_new, y = sharing_perc)) +
    geom_violin(aes(fill = DM_new)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - new DM", x = "") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() 
ggsave("results/strainsharing/dmnew.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(MetSyn_new)), aes(x = MetSyn_new, y = sharing_perc)) +
    geom_violin(aes(fill = MetSyn_new)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - new MetSyn", x = "") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    theme_Publication() 
ggsave("results/strainsharing/metsynnew.pdf", width = 4.5, height = 5)

dftot %>% group_by(EthnicityTot, DM_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc))
ggplot(data = dftot %>% filter(!is.na(DM_new)), aes(x = DM_new, y = sharing_perc)) +
    geom_violin(aes(fill = DM_new)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - new DM", x = "New diabetes diagnosis") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    facet_wrap(~EthnicityTot) +
    theme_Publication() 
ggsave("results/strainsharing/dmnew_ethnicities.pdf", width = 6, height = 5)

dftot %>% group_by(EthnicityTot, MetSyn_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc))
ggplot(data = dftot %>% filter(!is.na(MetSyn_new)), aes(x = MetSyn_new, y = sharing_perc)) +
    geom_violin(aes(fill = MetSyn_new)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing - new MetSyn", x = "New MetSyn diagnosis") +
    stat_compare_means(comparisons = list(c("Yes", "No")), tip.length = 0, hide.ns = TRUE,
                       label = "p.signif", method = "t.test") +
    facet_wrap(~EthnicityTot) +
    theme_Publication() 
ggsave("results/strainsharing/metsynnew_ethnicities.pdf", width = 6, height = 5)

ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = FUtime, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", x = "FU time (years)", title = "FU time and strain sharing") +
    stat_cor() +
    theme_Publication()
ggsave("results/strainsharing/futime.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = distance, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", x = "Bray-Curtis dissimilarity (paired)", title = "Bray dissimilarity and strain sharing") +
    stat_cor() +
    theme_Publication()
ggsave("results/strainsharing/bray_strainsharing.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = shannon, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", x = "Shannon index (baseline)", title = "Shannon at baseline and strain sharing") +
    stat_cor() +
    theme_Publication()
ggsave("results/strainsharing/shannon_strainsharing.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = shannon_delta, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, height = 0) +
    geom_smooth(color = "black", method = "loess") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", x = "Difference in Shannon index", title = "Shannon change and strain sharing") +
    # stat_cor() +
    theme_Publication()
ggsave("results/strainsharing/shannonchange_strainsharing.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = richness, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, height = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", x = "Richness (baseline)", title = "Richness (baseline) and strain sharing") +
    stat_cor() +
    theme_Publication()
ggsave("results/strainsharing/richness_strainsharing.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = richness_delta, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, height = 0) +
    geom_smooth(color = "black", method = "loess") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", x = "Difference in richness", title = "Richness change and strain sharing") +
    # stat_cor() +
    theme_Publication()
ggsave("results/strainsharing/richnesschange_strainsharing.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(FUtime)), aes(x = FUtime, y = sharing_perc)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", x = "FU time (years)", title = "FU time and strain sharing") +
    stat_cor() +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave("results/strainsharing/futime_ethnicity.pdf", width = 6, height = 5)

ggplot(data = dftot %>% filter(!is.na(Age)), aes(x = Age, y = sharing_perc)) +
    geom_jitter(color = "royalblue", alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", x = "Age (years)", title = "Baseline age and strain sharing") +
    stat_cor() +
    theme_Publication()
ggsave("results/strainsharing/age.pdf", width = 4.5, height = 5)

ggplot(data = dftot %>% filter(!is.na(Age)), aes(x = Age, y = sharing_perc)) +
    geom_jitter(aes(color = EthnicityTot), alpha = 0.5, width = 0) +
    geom_smooth(color = "black", method = "lm") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", x = "Age (years)", title = "Baseline age and strain sharing") +
    stat_cor() +
    facet_wrap(~EthnicityTot) +
    theme_Publication()
ggsave("results/strainsharing/age_ethnicity.pdf", width = 6, height = 5)


#### Taxonomy of strains ####
tax <- cladesplit %>% filter(SGB %in% rownames(abundance3))
tax <- tax[match(tax$SGB, rownames(abundance3)),]
tax$abundance <- rowSums(abundance3[2:ncol(abundance3)])/(ncol(abundance3)-1)
tax %>% arrange(-abundance) %>% dplyr::select(Species, SGB, abundance) %>%  slice(1:10)

othernames <- names(dfsame)[which(!str_detect(names(dfsame), "t__"))]
dfsame <- dfsame[,c(othernames, str_c("t__", tax$SGB))]
tax$sharing_sum <- apply(dfsame[,5:ncol(dfsame)], 2, function(x) sum(x, na.rm = TRUE) )
tax <- tax %>% filter(sharing_sum != 0)
nrow(tax) # 752 SGBs left that are non-zero
dfsame2 <- dfsame[,c(othernames, str_c("t__", tax$SGB))]
tax$sharing_perc <- apply(dfsame2[,5:ncol(dfsame2)], 2, function(x) (sum(x, na.rm = TRUE) / 
                                (sum(!x, na.rm = TRUE)+sum(x, na.rm = TRUE))) * 100)
tax$n <- apply(dfsame2[,5:ncol(dfsame2)], 2, function(x) sum(!x, na.rm = TRUE)+sum(x, na.rm = TRUE))
tax <- tax %>% right_join(., thres, by = "SGB") %>% 
    filter(n_samples > 100) %>% 
    filter(n > 50)

gghistogram(tax$sharing_perc, fill = "firebrick") + 
    labs(title = "Strain sharing per SGB", x = "percentage of samples") +
    theme_Publication()

tax %>% arrange(sharing_perc) %>% dplyr::select(Family, Species, SGB, sharing_perc, n) %>% slice(1:20)

ggplot(data = tax, aes(x = n, y = sharing_perc)) +
    geom_jitter(color = "darkgreen", alpha = 0.5) +
    geom_smooth(color = "black", method = "lm") +
    stat_cor(label.x = 250) +
    labs(title = "Strain sharing percentage vs prevalence", x = "Number of subjects detected",
         y = "Sharing percentage (%)") +
    theme_Publication()

# Calculate params for Dutch
idsdutch <- dftot$sampleID[which(dftot$EthnicityTot == "Dutch")]
dfsamedutch <- dfsame %>% filter(sampleid_1 %in% idsdutch) # select Dutch
dfsamedutch <- dfsamedutch[,c(othernames, str_c("t__", tax$SGB))] # put cols in same sequence as tax
taxdutch <- tax
taxdutch$sharing_sum <- apply(dfsamedutch[,5:ncol(dfsamedutch)], 2, function(x) sum(x, na.rm = TRUE) )
nrow(taxdutch) # 152 SGBs left that are non-zero
dfsamedutch2 <- dfsamedutch[,c(othernames, str_c("t__", taxdutch$SGB))]
taxdutch$sharing_perc <- apply(dfsamedutch2[,5:ncol(dfsamedutch2)], 2, 
                          function(x) (sum(x, na.rm = TRUE) / 
                                        (sum(!x, na.rm = TRUE)+sum(x, na.rm = TRUE))) * 100)
taxdutch$n <- apply(dfsamedutch2[,5:ncol(dfsamedutch2)], 2, 
               function(x) sum(!x, na.rm = TRUE)+sum(x, na.rm = TRUE))
taxdutch$EthnicityTot <- "Dutch"

# Calculate params for SAS
idssas <- dftot$sampleID[which(dftot$EthnicityTot == "South-Asian Surinamese")]
dfsamesas <- dfsame %>% filter(sampleid_1 %in% idssas) # select SAS
dfsamesas <- dfsamesas[,c(othernames, str_c("t__", tax$SGB))] # put cols in same sequence as tax
taxsas <- tax
taxsas$sharing_sum <- apply(dfsamesas[,5:ncol(dfsamesas)], 2, function(x) sum(x, na.rm = TRUE) )
nrow(taxsas) # 752 SGBs left that are non-zero
dfsamesas2 <- dfsamesas[,c(othernames, str_c("t__", tax$SGB))]
taxsas$sharing_perc <- apply(dfsamesas2[,5:ncol(dfsamesas2)], 2, 
                          function(x) (sum(x, na.rm = TRUE) / 
                                           (sum(!x, na.rm = TRUE)+sum(x, na.rm = TRUE))) * 100)
taxsas$n <- apply(dfsamesas2[,5:ncol(dfsamesas2)], 2, 
               function(x) sum(!x, na.rm = TRUE)+sum(x, na.rm = TRUE))
taxsas$EthnicityTot <- "South-Asian Surinamese"
taxtot <- rbind(taxdutch, taxsas)

ggplot(data = taxtot %>% arrange(sharing_perc), 
       aes(x = fct_reorder(SGB, .x = sharing_perc), y = sharing_perc, group = EthnicityTot)) +
    geom_segment(aes(y = 0, yend = sharing_perc, color = EthnicityTot)) +
    geom_point(aes(color = EthnicityTot), stat = "identity") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of subjects with stable strains",
         x = "",
         title = "Strain sharing") +
    theme_Publication() + 
    coord_flip()

ggsave("results/strainsharing/allsgbs_ethnicity.pdf", width = 10, height = 20)

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

leg <- c("Dutch" = pal_simpsons()(1), "South-Asian Surinamese" = pal_simpsons()(2)[2])
(pl <- ggplot(data = difftrue, aes(x = Species)) +
    geom_segment(aes(y = Dutch, yend = `South-Asian Surinamese`), color = "darkgrey") +
    geom_point(aes(y = `South-Asian Surinamese`, color = "South-Asian Surinamese"), size = 3) +
    geom_point(aes(y = Dutch, color = "Dutch"), size = 3) +
    scale_color_manual(values = leg) +
    labs(y = "Percentage of subjects with stable strains",
         x = "",
         title = "Strain sharing", color = "") +
    theme_Publication() + 
    coord_flip())
ggsave("results/strainsharing/diffsgbs_ethnicity_connect.pdf", width = 7, height = 7)

ab <- abundance3 %>% filter(clade_name %in% str_c("t__", difftrue$SGB)) %>% 
    dplyr::select(contains("HELI"))
ab <- as.data.frame(t(as.matrix(ab)))
ab <- ab[,levels(difftrue$SGB)]
colnames(ab) <- difftrue$Species[match(difftrue$SGB, colnames(ab))]
ab <- ab %>% rownames_to_column(var = "sampleID") %>% 
    right_join(dftot, by = "sampleID") %>% 
    pivot_longer(., cols = 2:35, names_to = "Species", values_to = "abundance") %>% 
    mutate(Species = fct_inorder(as.factor(Species))
           )

plist <- c()
for(a in difftrue$Species){
    print(a)
    dens <- ggplot(data = ab %>% filter(Species == a), 
               aes(x = abundance + 0.01, fill = EthnicityTot)) +
            geom_density(alpha = 0.5) +
            scale_x_log10(limits = c(0.1, 26)) +
            scale_fill_simpsons(guide = "none") +
            labs(y = "", fill = "", x = "") +
            theme_void()
    plist[[a]] <- dens
}
pl2 <- ggarrange(plotlist = plist, ncol = 1, common.legend = TRUE)
pl %>% aplot::insert_right(pl2, width = 0.2)
ggsave("results/strainsharing/diffsgbs_eth_dens.pdf", width = 10, height = 7)

(dens <- ggplot(data = ab, aes(x = abundance + 0.01, 
                              y = Species, 
                              fill = EthnicityTot)) +
    geom_density_ridges(alpha = 0.5, rel_min_height = 0.01) +
    scale_x_log10(limits = c(0.1, 26)) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "", fill = "", x = "log10(abundance)") + 
    theme_Publication() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank()))

# which strains have higher cut-offs (in general higher mutation rates?)
thres2 <- thres %>% filter(SGB %in% difftrue$SGB) %>% right_join(difftrue, ., by = "SGB")
pl3 <- ggplot(data = thres2, aes(x = Species, y = threshold_value)) +
    geom_bar(stat = "identity", fill = pal_simpsons()(3)[3]) +
    coord_flip() +
    theme_void() +
    labs(y = "threshold (ngd)", x = "") +
    theme_Publication() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank())

# on how many samples are the data based
pl4 <- ggplot(data = thres2, aes(x = Species, y = n_samples)) +
    geom_bar(stat = "identity", fill = pal_simpsons()(4)[4]) +
    coord_flip() +
    theme_void() +
    labs(y = "number of samples", x = "") +
    theme_Publication() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank())

# putting everything together
options("aplot_guides" = "keep")
pl %>% aplot::insert_right(dens, width = 0.3) %>% 
    aplot::insert_right(pl3, width = 0.3) %>% 
    aplot::insert_right(pl4, width = 0.3)
ggsave("results/strainsharing/diffstrains_complete.pdf", width = 15, height = 10)


