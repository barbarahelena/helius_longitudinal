## Compositional plots and diversity metrics

## Libraries
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(vegan)
library(circlize)
library(ggsci)

## Functions
# cols <- colorRampPalette(c(pal_flatui()(8)))
cols <- colorRampPalette(c(pal_simpsons()(16)[c(1,2,4,6,8:16)]))

theme_composition <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text.x =  element_text(angle = 45, hjust = 1),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                legend.key.size= unit(0.4, "cm"),
                legend.spacing  = unit(0, "cm"),
                legend.text = element_text(size = rel(0.7)),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
}

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
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

#### Output folder ####
resultsfolder <- "results/composition"
dir.create(resultsfolder, showWarnings = FALSE)

## Open data
phyloseqrare <- readRDS("data/phyloseq_sampledata.RDS")
tax <- readRDS("data/16s/phyloseq/rarefied/taxtable_rarefied.RDS")

tab <- as(phyloseqrare@otu_table, 'matrix')
counts <- sample_sums(phyloseqrare@otu_table)
tab <- as.data.frame(t(tab/sample_sums(phyloseqrare))*100)
rowSums(tab) # samples should all sum up to 100%
clindata <- as(phyloseqrare@sam_data, "data.frame")

#### Species-level ####
N <- 20
# convert to long format
d <- tab %>% 
    rownames_to_column(var = 'Sample') %>% 
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance') %>% 
    mutate(sampleID = Sample)
d$Tax <- tax$Tax[match(d$ASV, tax$ASV)]
d$Species <- tax$Species[match(d$ASV, tax$ASV)]

top_taxa <- d %>%
    filter(!is.na(Species)) %>% 
    group_by(Species, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>% 
    group_by(Species) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    arrange(-Abundance) %>% 
    dplyr::select(Species) %>% 
    head(N)
toptax <- d$Tax[match(top_taxa$Species, d$Species)]

dx <- d %>% 
    mutate(
        Species2 = case_when(
            is.na(Species) ~ "Unknown",
            Tax %in% toptax ~ paste(Tax),
            !(Tax %in% toptax) & !(is.na(Species)) ~ "Other species"
        ),
        Species2 = as.factor(Species2)
    ) %>% 
    group_by(Species2, Sample) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(timepoint = clindata$timepoint.y[match(Sample, clindata$sampleID)]) %>% 
    group_by(Species2, timepoint) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    ungroup(.) %>% 
    mutate(
        Species2 = fct_reorder(Species2, Abundance),
        Species2 = fct_relevel(Species2, "Other species", after = 0L),
        Species2 = fct_relevel(Species2, "Unknown", after = 0L)
    )

lev <- levels(dx$Species2)
# colsp <- c(cols[1:(length(lev)-2)], cols[21:22])

set.seed(1234)
(comp_species <- dx %>% 
        ggplot(aes(x = timepoint, y = Abundance, fill = Species2)) +
        geom_bar(stat = "identity", color = "black") +
        scale_fill_manual(values = rev(c(sample(cols(24),20), "grey70", "grey90")), labels = lev) +
        guides(fill = guide_legend(ncol = 1)) +
        labs(y="Composition (%)", x = "", title = "Species", fill = "") +
        scale_y_continuous(expand = c(0, 0)) +
        theme_composition())
ggsave(comp_species, filename = "results/composition/composition_species.pdf", width = 7, height = 5)

#### Genus-level ####
N <- 20
d <- tab %>% 
    rownames_to_column(var = 'Sample') %>% 
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance') %>% 
    mutate(sampleID = Sample)
d$Genus <- tax$Genus[match(d$ASV, tax$ASV)]

top_gen <- d %>%
    group_by(Genus, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Genus) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    arrange(-Abundance) %>%
    dplyr::select(Genus) %>%
    filter(Genus != 'ambiguous') %>%
    head(N) %>%
    unlist()
top_gen

dx_genus <- d %>% mutate(
    Genus2 = case_when(
        Genus %in% top_gen ~ paste(Genus),
        is.na(Genus) ~ paste("Unknown"),
        !(Genus %in% top_gen) ~ paste("Other genera")
    ),
    Genus2 = as.factor(Genus2)
) %>% 
    group_by(Genus2, Sample) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(timepoint = clindata$timepoint.y[match(Sample, clindata$sampleID)]) %>% 
    group_by(Genus2, timepoint) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    ungroup(.) %>% 
    mutate(group = "all samples",
           Genus2 = fct_reorder(Genus2, Abundance),
           Genus2 = fct_relevel(Genus2, "Other genera", after = 0L),
           Genus2 = fct_relevel(Genus2, "Unknown", after = 0L)
    )

lev <- levels(dx_genus$Genus2)

set.seed(18)
(comp_genus <- dx_genus %>% 
    ggplot(aes(x = timepoint, y = Abundance, fill = Genus2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = rev(c(sample(cols(20)), "grey70", "grey90")), labels = lev) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(y="Composition (%)", x = "", title = "Genus", fill = "") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition())
ggsave(comp_genus, filename = "results/composition/composition_genus.pdf", width = 6, height = 5)

#### Family level ####
N <- 20
# convert to long format
d <- tab %>% 
    rownames_to_column(var = 'Sample') %>% 
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance') %>% 
    mutate(sampleID = Sample)
# add species taxonomy (including ambiguous)
d$Family <- tax$Family[match(d$ASV, tax$ASV)]

top_fam <- d %>%
    group_by(Family, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Family) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    arrange(-Abundance) %>%
    dplyr::select(Family) %>%
    filter(Family != 'ambiguous') %>%
    head(N) %>% 
    unlist(.)

dx_family <- d %>% mutate(
    Family2 = case_when(
        Family %in% top_fam ~ paste(Family),
        is.na(Family) ~ paste("Unknown"),
        !(Family %in% top_taxa) ~ paste("Other families")
    ),
    Family2 = as.factor(Family2)
) %>% 
    group_by(Family2, Sample) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(timepoint = clindata$timepoint.y[match(Sample, clindata$sampleID)]) %>% 
    group_by(Family2, timepoint) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    ungroup(.) %>% 
    mutate(Family2 = fct_reorder(Family2, Abundance),
           Family2 = fct_relevel(Family2, "Other families", after = 0L),
           Family2 = fct_relevel(Family2, "Unknown", after = 0L),
           
    )

lev <- levels(dx_family$Family2)

set.seed(1234)
(comp_family <- dx_family %>% 
    ggplot(aes(x = timepoint, y = Abundance, fill = Family2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = rev(c(sample(sample(cols(20),20)), "grey70", "grey90")), labels = lev) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(y="Composition (%)", x = "", title = "Family", fill = "") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition())
ggsave(comp_family, filename = "results/composition/composition_family.pdf", width = 6, height = 5)


#### Phylum level ####
N <- 6
# convert to long format
d <- tab %>% 
    rownames_to_column(var = 'Sample') %>% 
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance') %>% 
    mutate(sampleID = Sample)
d$Phylum <- tax$Phylum[match(d$ASV, tax$ASV)]

top_taxa <- d %>%
    group_by(Phylum, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Phylum) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    arrange(-Abundance) %>%
    dplyr::select(Phylum) %>%
    filter(Phylum != 'ambiguous') %>%
    head(N) %>%
    unlist()
top_taxa

dx <- d %>% mutate(
    Phylum2 = case_when(
        Phylum %in% top_taxa ~ paste(Phylum),
        is.na(Phylum) ~ paste("Unknown"),
        !(Phylum %in% top_taxa) ~ paste("Other phyla")
    ),
    Phylum2 = as.factor(Phylum2)
) %>% 
    group_by(Phylum2, Sample) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(timepoint = clindata$timepoint.y[match(Sample, clindata$sampleID)]) %>% 
    group_by(Phylum2, timepoint) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    ungroup(.) %>% 
    mutate(group = "all samples",
           Phylum2 = fct_reorder(Phylum2, Abundance),
           Phylum2 = fct_relevel(Phylum2, "Other phyla", after = 0L),
           Phylum2 = fct_relevel(Phylum2, "Unknown", after = 0L)
    )

lev <- levels(dx$Phylum2)
# colphyl <- c(cols[1:N], cols[21:22])

set.seed(1234)
(comp_phylum <- dx %>% 
    ggplot(aes(x = timepoint, y = Abundance, fill = Phylum2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = rev(c(sample(cols(10),6), "grey70", "grey90")), labels = lev) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(y="Composition (%)", x = "", title = "Phylum", fill = "") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition())
ggsave(comp_phylum, filename = "results/composition/composition_phylum.pdf", width = 6, height = 5)

pl_total <- ggarrange(comp_species, comp_genus, comp_family, comp_phylum,
                      ncol = 4, labels = c("A", "B", "C", "D"),
                      widths = c(1.4, 1.2, 1.3, 1.0))
pl_comp <- ggarrange(comp_phylum, comp_family, comp_genus,
                     nrow = 1, widths = c(1.0, 1.3, 1.2))
ggsave(pl_total, filename = "results/composition/composition_total.pdf", width = 30, height = 8)

