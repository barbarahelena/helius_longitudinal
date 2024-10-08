#### Alpha diversity indices

## Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggsci)
library(phyloseq)

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

#### alpha diversity data ####
df <- readRDS("data/16s/clin_alphadiversity.RDS")
# counts <- sample_sums(phydata@otu_table)
# counts # samples should all sum up to 14932

## Output folder
resultsfolder <- "results/alphadiversity"
dir.create(resultsfolder, showWarnings = FALSE)

## Diversity metrics
# Shannon plots
(plshan <- ggplot(data = df, aes(x = timepoint, y = shannon, fill = timepoint)) +
    geom_violin() +
    scale_fill_manual(values = rev(pal_simpsons()(7)[c(1,7)]), guide = "none") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    stat_compare_means(label.y = 5.5) +
    labs(title = "Shannon index", y = "Shannon index", x="") + 
    theme_Publication())
ggsave(plshan, filename = "results/alphadiversity/shannon.svg", width = 4, height = 5)
ggsave(plshan, filename = "results/alphadiversity/shannon.pdf", width = 4, height = 5)

# Richness plots
(plrich <- ggplot(data = df, aes(x = timepoint, y = SR, fill = timepoint)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_manual(values = rev(pal_simpsons()(7)[c(1,7)]), guide = "none") + 
    labs(title = "Species richness", y = "Number of species", x = "") +
    stat_compare_means(method = "wilcox.test"))
ggsave(plrich, filename = "results/alphadiversity/richness.pdf", width = 4, height = 5)
ggsave(plrich, filename = "results/alphadiversity/richness.svg", width = 4, height = 5)

## Faith's PD
(plfaith <- ggplot(data = df, aes(x = timepoint, y = PD, fill = timepoint)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_manual(values = rev(pal_simpsons()(7)[c(1,7)]), guide = "none") + 
    labs(title = "Faith's PD", y = "Faith's phylogenetic diversity", x = "") +
    stat_compare_means(method = "wilcox.test"))
ggsave(plfaith, filename = "results/alphadiversity/faiths.pdf", device = "pdf", width = 4, height = 5)
ggsave(plfaith, filename = "results/alphadiversity/faiths.svg", device = "svg", width = 4, height = 5)

## Ggarrange
pl_total <- ggarrange(plshan, plrich, plfaith, labels = c("A", "B", "C"), nrow =1)
ggsave(pl_total, filename = "results/alphadiversity/alphadivplots.pdf", width = 11, height = 5.5)
ggsave(pl_total, filename = "results/alphadiversity/alphadivplots.svg", width = 11, height = 5.5)


linearmixed_arg <- function(data, var){
    data$var <- data[[var]]
    model1_v4 <- lmer(var ~ Ethnicity*timepoint + (1|ID),
                      data = data)
    res_v4 <- summary(model1_v4)
    pval <- format(round(res_v4$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres <- cbind(gene = var, group1 = 0, group2 = 1, pval)
    statres <- tibble::as_tibble(statres)
    statres$p_signif <- case_when(
        statres$pval < 0.05 ~paste0("*"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.001 ~paste0("***"),
        .default = ""
    )
    return(statres)
}

plot_lmm <- function(data, var, lmm){
    data$var <- log10(data[[var]]+1)
    df_means <- data %>% 
        group_by(Ethnicity, timepoint) %>% 
        summarise(mean_var=mean(var, na.rm = TRUE), sd_var = sd(var, na.rm = TRUE),
                  n_var = length(var))
    pl <- ggplot() +
        geom_line(data = df_means, aes(x = timepoint, y = mean_var, 
                                       color = Ethnicity, group = Ethnicity), alpha = 1, linewidth = 0.75) +
        geom_line(data = data, aes(x = timepoint, y = var, color = Ethnicity, group = ID), alpha = 0.1) +
        geom_errorbar(data = df_means,
                      aes(ymin = mean_var - (sd_var/sqrt(n_var)),
                          ymax = mean_var + (sd_var/sqrt(n_var)),
                          x = timepoint,
                          color = Ethnicity), width=0.05, linewidth = 0.75) +
        ggpubr::stat_pvalue_manual(lmm, y.position = max(data$var), label = "p_signif", 
                                   tip.length = 0, bracket.shorten = 0.1, size = 5) +
        scale_color_manual(values = rev(pal_simpsons()(2))) + 
        # scale_y_continuous(limits = c(0,220), breaks = seq(from = 0, to = 220, by = 20)) +
        theme_Publication() +
        labs(x = "Timepoint", y = "log10(cpm+1)", title = var,
             color = "")
    print(pl)
}


names(df)
lmm_var <- c("PD", "SR", "shannon")
df2 <- df %>% 
    filter(Ethnicity != "Other" & Ethnicity != "Unknown") %>% 
    filter(sampleID %in% mbids$sampleID) %>% droplevels(.) %>% 
    mutate(timepoint = case_when(
        timepoint == "baseline" ~ 0,
        timepoint == "follow-up" ~ 1
    ),
    timepoint = as.numeric(timepoint),
    ID = str_c("S", str_remove(str_remove(ID, "HELIFU_"), "HELIBA_")),
    ID = as.factor(ID)
    ) %>% droplevels(.)
dim(df2)

names(df2)
for(a in lmm_var){
    print(a)
    pl <- plot_lmm(data = df2, var = a, lmm = statsig)
    plist[[a]] <- pl
}



#### Shotgun data ####
# ## Load data
# df_new <- rio:: import("data/clinicaldata_shotgun.RDS")
# 
# ## Diversity metrics
# # Shannon plots
# shannon <- rio::import("data/metaphlan/diversity/combined_table_shannon.tsv") 
# shannon <- shannon %>% select(ID = V1, shannon = diversity_shannon) %>% 
#     mutate(ID = str_remove(ID, "_T1"))
# df_shan <- left_join(shannon, df_new, by = "ID")
# plshan <- ggplot(data = df_shan, aes(x = Sex, y = shannon, fill = Sex)) +
#     geom_violin() +
#     scale_fill_manual(values = rev(pal_nejm()(2)), guide = "none") +
#     geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
#     stat_compare_means(label.y = 5.5) +
#     labs(title = "Shannon index", y = "Shannon index", x="") + 
#     theme_Publication()
# ggsave(plshan, filename = "results/alphadiversity/shannon_shotgun.svg", width = 4, height = 5)
# ggsave(plshan, filename = "results/alphadiversity/shannon_shotgun.pdf", width = 4, height = 5)
# 
# 
# ## Species richness
# richness <- rio::import("data/metaphlan/diversity/combined_table_richness.tsv")
# richness <- richness %>% select(ID = V1, richness = observed) %>% 
#     mutate(ID = str_remove(ID, "_T1"))
# dfspec <- left_join(richness, df_new, by = "ID")
# 
# # Male-female
# plrich <- ggplot(data = dfspec, aes(x = Sex, y = richness, fill = Sex)) +
#     geom_violin()+
#     geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
#     theme_Publication() + 
#     scale_fill_manual(values = rev(pal_nejm()(2)), guide = "none") + 
#     labs(title = "Species richness", y = "Number of species", x = "") +
#     stat_compare_means(method = "wilcox.test")
# ggsave(plrich, filename = "results/alphadiversity/richness_shotgun_sex.pdf", width = 4, height = 5)
# ggsave(plrich, filename = "results/alphadiversity/richness_shotgun_sex.svg", width = 4, height = 5)
#
# 
# ## Ggarrange male-female
# pl_total <- ggarrange(plshan, plrich, labels = c("A", "B"), nrow =1)
# ggsave(pl_total, filename = "results/alphadiversity/alphadivplots_shotgun.pdf", width = 8, height = 5.5)
# ggsave(pl_total, filename = "results/alphadiversity/alphadivplots_shotgun.svg", width = 8, height = 5.5)
# 

