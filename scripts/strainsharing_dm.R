## Strain sharing
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
thres$n_markers <- NULL # bug in pipeline: n_samples = n_markers, n_samples not extracted from info..
thres$n_markers <- thres$n_samples
thres$n_samples <- NULL
abundance <- rio::import("data/shotgun/combined_table_fixedlab.tsv")

#### Calculation strain sharing metric ####
colnames(df) <- str_remove(colnames(df), "sharing_")
sharing_sgb <- apply(df[,3:ncol(df)], 2, function(x) sum(x, na.rm = TRUE))
length(sharing_sgb[which(sharing_sgb == 0)])
# df <- df %>% dplyr::select(-which(sharing_sgb == 0))
sharing_sum <- apply(df[,3:ncol(df)], 1, function(x) sum(x, na.rm = TRUE))
sharing_sum[1:5]
df$sharing_sum <- sharing_sum
strain_total <- apply(df[,3:ncol(df)], 1, function(x) sum(!is.na(x), na.rm = TRUE))
strain_total[1:5]
df$strain_total <- strain_total
# df$sharing_perc <- (sharing_sum / ncol(df[,which(str_detect(colnames(df), "t__"))])) * 100
df$sharing_perc <- (sharing_sum / strain_total) * 100

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

dftot %>% group_by(DM_new) %>% summarise(mean_sh = mean(sharing_perc, na.rm = TRUE), n_sh = length(sharing_perc))
dftot <- dftot %>% mutate(DM_new = case_when(DM_new == "Yes" ~ "New DM",
                                             DM_new == "No" ~ "No DM",
                                             is.na(DM_new) ~ "Baseline DM"
    
))
table(dftot$DM_new)
comp <- list(c("Baseline DM", "No DM"), c("Baseline DM", "New DM"), c("No DM", "New DM"))
ggplot(data = dftot, aes(x = fct_reorder(DM_new, .x = sharing_perc, .fun = median), y = sharing_perc)) +
    geom_violin(aes(fill = DM_new)) +
    geom_boxplot(fill = "white", width = 0.2) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "Percentage of stable strains", title = "Strain sharing between timepoints", x = "") +
    stat_compare_means(comparisons = comp, tip.length = 0, label.y = 95, step.increase = 0.07,
                       label = "p.format", method = 't.test') +
    theme_Publication() 
ggsave("results/strainsharing/dm_diagnoses.pdf", width = 4.5, height = 5)

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
    stat_compare_means(comparisons = comp, tip.length = 0,
                       label = "p.format", method = "t.test") +
    facet_wrap(~EthnicityTot) +
    theme_Publication() 
ggsave("results/strainsharing/dmnew_ethnicities_2.pdf", width = 6, height = 5)


#### Taxonomy of strains ####
tax <- cladesplit %>% filter(SGB %in% rownames(abundance3))
tax <- tax[match(tax$SGB, rownames(abundance3)),]
tax$abundance <- rowSums(abundance3[2:ncol(abundance3)])/(ncol(abundance3)-1)
tax %>% arrange(-abundance) %>% dplyr::select(Species, SGB, abundance) %>%  slice(1:10)
othernames <- names(dfsame)[which(!str_detect(names(dfsame), "t__"))]
dfsame <- dfsame[,c(othernames, str_c("t__", tax$SGB))]
tax$sharing_sum <- apply(dfsame[,6:ncol(dfsame)], 2, function(x) sum(x, na.rm = TRUE) )
# tax <- tax %>% filter(sharing_sum != 0)
nrow(tax) # 752 SGBs left that are non-zero
dfsame2 <- dfsame[,c(othernames, str_c("t__", tax$SGB))]
tax$sharing_perc <- apply(dfsame2[,6:ncol(dfsame2)], 2, function(x) (sum(x, na.rm = TRUE) / 
                                (sum(!is.na(x), na.rm = TRUE))) * 100)
tax$n <- apply(dfsame2[,6:ncol(dfsame2)], 2, function(x) (sum(!is.na(x))))
tax <- tax %>% right_join(., thres, by = "SGB") %>% 
    filter(n > 50)

# Calculate params for Dutch
idsDM <- dftot$sampleID[which(dftot$DM == "Yes")]
dfsameDM <- dfsame %>% filter(sampleid_1 %in% idsDM) # select DM
dfsameDM <- dfsameDM[,c(othernames, str_c("t__", tax$SGB))] # put cols in same sequence as tax
taxDM <- tax
taxDM$sharing_sum <- apply(dfsameDM[,6:ncol(dfsameDM)], 2, function(x) sum(x, na.rm = TRUE) )
nrow(taxDM) # 153 SGBs left that are non-zero
dfsameDM2 <- dfsameDM[,c(othernames, str_c("t__", taxDM$SGB))]
taxDM$sharing_perc <- apply(dfsameDM2[,6:ncol(dfsameDM2)], 2, 
                             function(x) (sum(x, na.rm = TRUE) / 
                                              sum(!is.na(x))) * 100)
taxDM$n <- apply(dfsameDM2[,6:ncol(dfsameDM2)], 2, 
                  function(x) sum(!is.na(x)))
taxDM$DM <- "DM"

# Calculate params for SAS
idsctrl <- dftot$sampleID[which(dftot$DM == "No")]
dfsamectrl <- dfsame %>% filter(sampleid_1 %in% idsctrl) # select ctrl
dfsamectrl <- dfsamectrl[,c(othernames, str_c("t__", tax$SGB))] # put cols in same sequence as tax
taxctrl <- tax
taxctrl$sharing_sum <- apply(dfsamectrl[,6:ncol(dfsamectrl)], 2, function(x) sum(x, na.rm = TRUE) )
nrow(taxctrl) # 153 SGBs left that are non-zero
dfsamectrl2 <- dfsamectrl[,c(othernames, str_c("t__", tax$SGB))]
taxctrl$sharing_perc <- apply(dfsamectrl2[,6:ncol(dfsamectrl2)], 2, 
                          function(x) (sum(x, na.rm = TRUE) / 
                                           sum(!is.na(x))) * 100)
taxctrl$n <- apply(dfsamectrl2[,6:ncol(dfsamectrl2)], 2, 
               function(x) sum(!is.na(x)))
taxctrl$DM <- "No DM"
taxtot <- rbind(taxDM, taxctrl)

ggplot(data = taxtot %>% arrange(sharing_perc), 
       aes(x = fct_reorder(Species, .x = sharing_perc), y = sharing_perc, group = DM)) +
    geom_segment(aes(y = 0, yend = sharing_perc, color = DM)) +
    geom_point(aes(color = DM), stat = "identity") +
    scale_color_simpsons(guide = "none") +
    labs(y = "Percentage of subjects with stable strains",
         x = "",
         title = "Strain sharing") +
    theme_Publication() + 
    coord_flip()

ggsave("results/strainsharing/allsgbs_dm.pdf", width = 10, height = 20)

difftax <- taxtot %>% pivot_wider(., id_cols = c(SGB, Species), names_from = c(DM), 
                                  values_from = sharing_perc) %>% 
    mutate(diff = DM - `No DM`,
           diff_bin = case_when(diff < -20 | diff > 20 ~ TRUE, .default = FALSE)) 
difftrue <- difftax %>% filter(diff_bin == TRUE) %>% 
    mutate(highest = case_when(DM >= `No DM` ~ DM,
                               DM < `No DM` ~ `No DM`),
           Species = fct_reorder(make.unique(Species), .x = highest)) %>% arrange(highest) %>% 
    mutate(SGB = as.factor(SGB), SGB = fct_inorder(SGB))

leg <- c("DM" = pal_simpsons()(1), "No DM" = pal_simpsons()(2)[2])
(pl <- ggplot(data = difftrue, aes(x = Species)) +
    geom_segment(aes(y = DM, yend = `No DM`), color = "darkgrey") +
    geom_point(aes(y = `No DM`, color = "No DM"), size = 3) +
    geom_point(aes(y = DM, color = "DM"), size = 3) +
    scale_color_manual(values = leg) +
    labs(y = "Percentage of subjects with stable strains",
         x = "",
         title = "Strain sharing", color = "") +
    theme_Publication() + 
    coord_flip())
ggsave("results/strainsharing/diffsgbs_dm_connect.pdf", width = 7, height = 7)

ab <- abundance3 %>% filter(clade_name %in% str_c("t__", difftrue$SGB)) %>% 
    dplyr::select(contains("HELI"))
ab <- as.data.frame(t(as.matrix(ab)))
ab <- ab[,levels(difftrue$SGB)]
colnames(ab) <- make.unique(as.character(difftrue$Species[match(difftrue$SGB, colnames(ab))]))
ab <- ab %>% rownames_to_column(var = "sampleID") %>% 
    right_join(dftot, by = "sampleID") %>% 
    pivot_longer(., cols = 2:(nlevels(difftrue$SGB)+1), names_to = "Species", values_to = "abundance") %>% 
    mutate(Species = fct_inorder(as.factor(Species)))

plist <- c()
for(a in difftrue$Species){
    print(a)
    dens <- ggplot(data = ab %>% filter(Species == a), 
               aes(x = abundance + 0.01, fill = DM)) +
            geom_density(alpha = 0.5) +
            scale_x_log10() + # limits = c(0.1, 26)
            scale_fill_simpsons(guide = "none") +
            labs(y = "", fill = "", x = "") +
            theme_void()
    plist[[a]] <- dens
}
pl2 <- ggarrange(plotlist = plist, ncol = 1, common.legend = TRUE)
pl %>% aplot::insert_right(pl2, width = 0.2)
ggsave("results/strainsharing/diffsgbs_dm_dens.pdf", width = 10, height = 7)

(dens <- ggplot(data = ab, aes(x = abundance + 0.01, 
                              y = Species, 
                              fill = DM)) +
    geom_density_ridges(alpha = 0.5, rel_min_height = 0.01) +
    scale_x_log10(limits = c(0.1, 26)) +
    scale_fill_simpsons(guide = "none") +
    labs(y = "", fill = "", x = "log10(abundance)") + 
    theme_Publication() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank()))

# which strains have higher cut-offs (in general higher mutation rates?)
thres2 <- thres %>% filter(SGB %in% difftrue$SGB) %>% right_join(difftrue, ., by = "SGB")
(pl3 <- ggplot(data = thres2, aes(x = Species, y = threshold_value)) +
    geom_bar(stat = "identity", fill = pal_simpsons()(3)[3]) +
    coord_flip() +
    theme_void() +
    labs(y = "threshold (ngd)", x = "") +
    theme_Publication() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank()
          ))

# on how many samples are the data based
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
    labs(y = "number of subjects", x = "") +
    theme_Publication() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank()
          ))

# putting everything together
options("aplot_guides" = "keep")
pl %>% aplot::insert_right(dens, width = 0.3) %>% 
    aplot::insert_right(pl3, width = 0.3) %>% 
    aplot::insert_right(pl4, width = 0.3)
ggsave("results/strainsharing/diffstrains_complete.pdf", width = 15, height = 10)

#### Distance plots of diff strains ####
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
    
    pln <- ggplot(data = ngd4) +
            geom_density_ridges(aes(x = nGD, y = Species, fill = EthnicityTot),
                                alpha = 0.5, rel_min_height = 0.01, bandwidth = 0.005) +
            theme_Publication() +
            scale_fill_simpsons(guide = "none") +
            # scale_x_continuous(limits = c(0, 0.06)) +
            labs(x = "nGD", y = "", fill = "", title = b) +
            theme(axis.text.y = element_blank(),
                  plot.title = element_text(size = rel(0.7)))
    # if(b != "Ruminococcus_bromii") {
    #     plnb <- pln + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    #                         axis.text.y = element_blank(), axis.text.x = element_blank())
    # } else{ plnb <- pln }
    plist2[[b]] <- pln
}

ngdpl <- ggarrange(plotlist = rev(plist2), ncol = 7, nrow = 6, common.legend = TRUE)
ggsave(ngdpl, "results/strainsharing/ngd.pdf", width = 18, height = 12)
