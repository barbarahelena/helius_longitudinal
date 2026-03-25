# LMMs on species from lmm script
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(afex)
library(phyloseq)

# Theme
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
meta <- readRDS("data/clinicaldata_wide.RDS")
mbs <- readRDS("data/16s/phyloseq_paired16s.RDS")
mbs
tax <- as.data.frame(mbs@tax_table)
head(tax)
alist <- tax |> filter(Genus == "Alistipes" & Species == "putredinis") # 1: ASV_66
odori <- tax |> filter(Genus == "Odoribacter" & Species == "splanchnicus") # 1: ASV_208
parabac <- tax |> filter(Genus == "Parabacteroides" & Species == "distasonis") # 6
gastra <- tax |> filter(Order == "Gastranaerophilales") # 62
actin <- tax |> filter(Genus == "Actinomyces") # 23
senegal <- tax |> filter(Genus == "Senegalimassilia") # 4
erysip <- tax |> filter(Genus == "Erysipelatoclostridium" & Species == "ramosum") # 1: ASV_364

taxa <- rbind(alist, odori, parabac, gastra, actin, senegal, erysip)
rownames(taxa)
mbs <- prune_taxa(taxa_names(mbs) %in% rownames(taxa), mbs)
mbs
sample_names(mbs)
sample_names(sample_data(mbs)) <- sample_names(mbs)
sample_names(otu_table(mbs))
mbs_genus <- tax_glom(mbs, taxrank = "Genus")
mbs_genus
mbs_genus@tax_table
tax <- as.data.frame(tax_table(mbs_genus))
taxname_map <- c(
    "Alistipes"      = "Alistipes putredinis",
    "Odoribacter"    = "Odoribacter splanchnicus",
    "Parabacteroides"= "Parabacteroides distasonis",
    "Actinomyces"    = "Actinomyces spp.",
    "Senegalimassilia" = "Senegalimassilia spp.",
    "CAG-196"        = "CAG-196",
    "Zag_111"        = "Zag_111"
)
tax$taxname <- taxname_map[tax$Genus]
tax$taxname[is.na(tax$taxname)] <- tax$Genus[is.na(tax$taxname)]  # fallback to genus
tax_table(mbs_genus) <- tax_table(as.matrix(tax))

mbdf <- as.data.frame(as(otu_table(mbs_genus), "matrix"))
rownames(mbdf)
head(mbdf)[1:5,1:5]
tax

meta <- meta |> 
  dplyr::select(ID, EthnicityTot) |> 
  filter(EthnicityTot %in% c("Dutch", "South-Asian Surinamese")) 
mbn <- as.data.frame(t(as.matrix(mbdf)))
colnames(mbn) <- tax$taxname[match(rownames(tax), colnames(mbn))]
head(mbn)
mbn <- mbn |> mutate_all(function(x) log10(x + 1 / 150))
mbn <- mbn |> mutate(
    sampleID = rownames(mbn),
    ID = str_c("S", str_remove(str_remove(sampleID, "HELIBA_"), "HELIFU_")),
    timepoint = case_when(
      str_detect(sampleID, "HELIBA") ~ "baseline",
      str_detect(sampleID, "HELIFU") ~ "follow-up"
    )
)
df_tot <- left_join(mbn, meta, by = "ID") |> filter(!is.na(EthnicityTot))
head(df_tot)

statres <- c()
for(a in 1:6){
  df_tot$mb <- df_tot[[a]]
  mbname <- colnames(df_tot)[a]
  model1 <- lmer(mb ~ EthnicityTot*timepoint + (1|ID), data = df_tot)
  res <- summary(model1)
  print(res)
  confint_model1 <- confint(model1)
  estimate <- as.numeric(format(round(res$coefficients[4,1], 3), nsmall = 3))
  conflow <- as.numeric(format(round(confint_model1[6,1], 3), nsmall = 3))
  confhigh <- as.numeric(format(round(confint_model1[6,2], 3), nsmall = 3))
  pval <- format(round(res$coefficients[4,5], 3), nsmall = 5)
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
statres

plist <- list()
for(i in 1:nrow(statres)){
    nm <- statres$mbname[i]
    print(nm)
    pval <- statres$pval[i]
    df_tot$mb <- df_tot[,statres$mbname[i]]
    df_means <- df_tot |> group_by(EthnicityTot, timepoint) |> 
        summarise(mean = mean(mb), sd = sd(mb), n = length(mb), .groups = "drop_last")
    print(df_means)
    res_lmm <- statres |> filter(mbname == nm) |> dplyr::select(-mbname)
    if(max(df_tot$mb) < 0) mbmax <- max(df_tot$mb*0.8) else mbmax <- max(df_tot$mb*1.2)
    if(max(df_tot$mb) < 0) mbstat <- max(df_tot$mb*0.9) else mbstat <- max(df_tot$mb*0.7)
    mbmin <- min(df_tot$mb)
    print(mbmax); print(mbmin)
    pl2 <- ggplot() +
        geom_line(data = df_tot, aes(x = timepoint, y = mb,
                  color = EthnicityTot, group = ID), alpha = 0.01, linewidth = 0.5) +
        geom_point(data = df_tot, aes(x = timepoint, y = mb,
                  color = EthnicityTot, group = EthnicityTot), alpha = 0.01, size = 0.8) +
        geom_line(data = df_means, aes(x = timepoint, y = mean, 
                  color = EthnicityTot, group = EthnicityTot), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = timepoint, y = mean, 
                  color = EthnicityTot, group = EthnicityTot), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = mean - (sd/sqrt(n)),
                          ymax = mean + (sd/sqrt(n)),
                          x = timepoint,
                          color = EthnicityTot), width=0.1) +
        stat_pvalue_manual(res_lmm, y.position = mbstat, label = "{sigq}", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5) +
        scale_color_jco() + 
        coord_cartesian(ylim = c(mbmin,mbmax)) +
        theme_Publication() +
        labs(x = "Timepoint", y = "log10(Relative abundance (%))", title = nm, color = "")
        plist[[i]] <- pl2
}

(plots <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:11],
          nrow = 2, ncol = 3))

dir.create("results/3_species_change/2_species/lmer", recursive = TRUE, showWarnings = FALSE)
ggsave(plots, filename = "results/3_species_change/2_species/lmer/16s_lmer_plots.pdf", width = 10, height = 7)
write.csv2(statres, "results/3_species_change/2_species/lmer/16s_lmm_results.csv")
df_tot$mb <- NULL
saveRDS(df_tot, "data/16s/selectedspecies.RDS")
