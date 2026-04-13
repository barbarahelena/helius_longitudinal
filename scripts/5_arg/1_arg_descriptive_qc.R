## ARG descriptives and QC plots
## Exploratory plots, not part of the assembled figure

library(tidyverse)
library(ggsci)
library(ggpubr)

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family) +
     theme(plot.title = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA, fill = NA),
           plot.background = element_rect(colour = NA, fill = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold", size = rel(0.8)),
           axis.title.y = element_text(angle=90, vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(size = rel(0.7)),
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
           strip.text = element_text(face="bold")))
}

dir.create("results/5_arg", showWarnings = FALSE, recursive = TRUE)

df_raw <- rio::import("data/shotgun/arg/all_samples.merged_arg_counts.tsv") %>%
  mutate(Class = str_to_title(Class),
         Subclass = str_to_title(Subclass), 
          timepoint = case_when(str_detect(Sample, "HELIBA") ~ "Baseline",
                                str_detect(Sample, "HELIFU") ~ "Follow-up"),
          timepoint = as.factor(timepoint)) |> 
  filter(! Sample %in% c("HELIBA_"))

# Class-level prevalence
class_prevalence <- df_raw |>
  filter(Prevalence > 0) %>%
  group_by(Class, timepoint) %>%
  summarise(n_genes = n_distinct(Gene_Symbol),
            n_samples = n_distinct(Sample),
            total_reads = sum(Mapped_Reads),
            mean_coverage = mean(Coverage, na.rm = TRUE), .groups = "drop_last") %>%
  arrange(-n_samples) %>%
  mutate(Class = fct_reorder(Class, n_samples, .fun = mean))
write.csv2(class_prevalence, "results/5_arg/qc/class_prevalence.csv", row.names = FALSE)

# Subclass-level prevalence
subclass_prevalence <- df_raw %>%
  filter(Prevalence == 1) %>%
  group_by(Class, Subclass, timepoint) %>%
  summarise(n_genes = n_distinct(Gene_Symbol),
            n_samples = n_distinct(Sample),
            total_reads = sum(Mapped_Reads),
            mean_coverage = mean(Coverage, na.rm = TRUE), .groups = "drop") %>%
  arrange(-n_samples)
write.csv2(subclass_prevalence, "results/5_arg/qc/subclass_prevalence.csv", row.names = FALSE)

# Gene-level prevalence
gene_prevalence <- df_raw %>%
  group_by(Gene_Symbol, Class, Subclass, timepoint) %>%
  summarise(n_samples = n_distinct(Sample),
            total_reads = sum(Mapped_Reads),
            mean_coverage = mean(Coverage, na.rm = TRUE),
            prevalence_pct = (n_distinct(Sample) / n_distinct(df_raw$Sample)) * 100,
            .groups = "drop") %>%
  arrange(-n_samples)
write.csv2(gene_prevalence, "results/5_arg/qc/gene_prevalence.csv", row.names = FALSE)

# Class distribution
(p1 <- ggplot(class_prevalence, aes(x = Class, y = n_samples, fill = timepoint)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_simpsons() +
  theme_Publication() +
  labs(x = "ARG Class", y = "Number of samples", title = "ARG Class Distribution", fill = "Timepoint"))
ggsave("results/5_arg/qc/class_distribution.pdf", p1, width = 8, height = 6)

# QC Plot 2: Top 20 subclasses
top_subclasses <- subclass_prevalence %>%
  slice_max(n_samples, n = 20) %>%
  mutate(Subclass = fct_reorder(Subclass, n_samples))

n_classes <- n_distinct(top_subclasses$Class)
p2 <- ggplot(top_subclasses, aes(x = Subclass, y = n_samples, fill = Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  {if(n_classes <= 7) scale_fill_bmj() else scale_fill_viridis_d(option = "turbo")} +
  theme_Publication() +
  facet_wrap(~timepoint) +
  labs(x = "ARG Subclass", y = "Number of samples", title = "Top 20 ARG Subclasses")
ggsave("results/5_arg/qc/top20_subclass_distribution.pdf", p2, width = 10, height = 7)

# QC Plot 3: Sample-level ARG diversity
sample_diversity <- df_raw %>%
  filter(Prevalence > 0) %>%
  group_by(Sample) %>%
  summarise(n_genes = n_distinct(Gene_Symbol),
            n_classes = n_distinct(Class),
            total_reads = sum(Mapped_Reads), .groups = "drop_last") %>%
  mutate(timepoint = case_when(str_detect(Sample, "HELIBA") ~ "baseline",
                               str_detect(Sample, "HELIFU") ~ "follow-up",
                               TRUE ~ NA_character_))

p3 <- ggplot(sample_diversity, aes(x = timepoint, y = n_genes, fill = timepoint)) +
  geom_violin(colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white") +
  stat_compare_means() +
  scale_fill_simpsons() +
  theme_Publication() +
  labs(x = "Timepoint", y = "Number of ARGs per sample", title = "ARG Diversity per Sample") +
  theme(legend.position = "none")
ggsave("results/5_arg/qc/sample_arg_diversity.pdf", p3, width = 5, height = 6)

# QC Plot 4: Class diversity per sample
p4 <- ggplot(sample_diversity, aes(x = timepoint, y = n_classes, fill = timepoint)) +
  geom_violin(colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white") +
  stat_compare_means() +
  scale_fill_simpsons() +
  theme_Publication() +
  labs(x = "Timepoint", y = "Number of ARG classes per sample", title = "ARG Class Diversity per Sample") +
  theme(legend.position = "none")
ggsave("results/5_arg/qc/sample_class_diversity.pdf", p4, width = 5, height = 6)

# QC Plot 5: Gene prevalence distribution
gene_prev_dist <- gene_prevalence %>%
  mutate(prev_category = cut(prevalence_pct,
                             breaks = c(0, 5, 10, 25, 50, 75, 100),
                             labels = c("<5%", "5-10%", "10-25%", "25-50%", "50-75%", ">75%")))

p5 <- ggplot(gene_prev_dist, aes(x = prev_category, fill = prev_category)) +
  geom_bar() +
  scale_fill_viridis_d() +
  theme_Publication() +
  labs(x = "Prevalence across samples", y = "Number of genes", title = "Gene Prevalence Distribution") +
  theme(legend.position = "none")
ggsave("results/5_arg/qc/gene_prevalence_distribution.pdf", p5, width = 6, height = 5)

# QC Plot 6: Coverage distribution by class
coverage_by_class <- df_raw %>%
  filter(Coverage > 0) %>%
  group_by(Class) %>%
  filter(n() > 50) %>%
  ungroup() %>%
  mutate(Class = fct_reorder(Class, Coverage, .fun = median))

p6 <- ggplot(coverage_by_class, aes(x = Class, y = log10(Coverage + 1), fill = Class)) +
  geom_violin(colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA) +
  coord_flip() +
  scale_fill_viridis_d() +
  theme_Publication() +
  facet_wrap(~timepoint) +
  labs(x = "", y = "log10(Coverage + 1)", title = "Coverage Distribution by ARG Class") +
  theme(legend.position = "none")
ggsave("results/5_arg/qc/coverage_by_class.pdf", p6, width = 12, height = 7)

# QC Plot 7: Sequencing depth analysis
# Check if total ARG burden is biased by sequencing depth differences
sample_depth <- df_raw %>%
  group_by(Sample) %>%
  summarise(sequencing_depth = mean(Total_Reads, na.rm = TRUE),
            total_arg_rpkm = sum(RPKM, na.rm = TRUE),
            n_arg_genes = sum(Prevalence > 0),
            .groups = "drop") %>%
  mutate(timepoint = case_when(str_detect(Sample, "HELIBA") ~ "baseline",
                               str_detect(Sample, "HELIFU") ~ "follow-up",
                               TRUE ~ NA_character_))

# Depth comparison between timepoints
depth_summary <- sample_depth %>%
  group_by(timepoint) %>%
  summarise(mean_depth = mean(sequencing_depth),
            median_depth = median(sequencing_depth),
            sd_depth = sd(sequencing_depth),
            .groups = "drop")
write.csv2(depth_summary, "results/5_arg/qc/depth_by_timepoint.csv", row.names = FALSE)

(p7 <- ggplot(sample_depth, aes(x = timepoint, y = sequencing_depth, fill = timepoint)) +
  geom_violin(colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA) +
  stat_compare_means() +
  scale_fill_simpsons() +
  scale_y_log10() +
  theme_Publication() +
  labs(x = "", y = "Sequencing Depth - Total Reads (log scale)",
       title = "Sequencing Depth by Timepoint") +
  theme(legend.position = "none"))
ggsave("results/5_arg/qc/depth_by_timepoint.pdf", p7, width = 6, height = 6)

# Correlation between depth and ARG burden
cor_result <- cor.test(sample_depth$sequencing_depth, sample_depth$total_arg_rpkm,
                       method = "spearman")

(p8 <- ggplot(sample_depth, aes(x = sequencing_depth, y = total_arg_rpkm, color = timepoint)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black", linetype = "dashed") +
  scale_color_simpsons() +
  scale_x_log10() +
  scale_y_log10() +
  theme_Publication() +
  labs(x = "Sequencing Depth - Total Reads (log scale)",
       y = "Total ARG Burden - RPKM (log scale)",
       title = "Correlation: Sequencing Depth vs ARG Burden",
       subtitle = paste0("Spearman's rho = ", round(cor_result$estimate, 3),
                        ", p = ", formatC(cor_result$p.value, format = "e", digits = 2)),
       color = "Timepoint"))
ggsave("results/5_arg/qc/depth_vs_burden_correlation.pdf", p8, width = 8, height = 6)

# Faceted by timepoint
(p9 <- ggplot(sample_depth, aes(x = sequencing_depth, y = total_arg_rpkm)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~timepoint) +
  theme_Publication() +
  labs(x = "Sequencing Depth - Total Reads (log scale)",
       y = "Total ARG Burden - RPKM (log scale)",
       title = "Depth vs Burden by Timepoint"))
ggsave("results/5_arg/qc/depth_vs_burden_by_timepoint.pdf", p9, width = 10, height = 5)

# Save correlation statistics
cor_baseline <- sample_depth %>%
  filter(timepoint == "baseline") %>%
  with(cor.test(sequencing_depth, total_arg_rpkm, method = "spearman"))

cor_followup <- sample_depth %>%
  filter(timepoint == "follow-up") %>%
  with(cor.test(sequencing_depth, total_arg_rpkm, method = "spearman"))

depth_burden_cor <- data.frame(
  group = c("Overall", "Baseline", "Follow-up"),
  spearman_rho = c(cor_result$estimate, cor_baseline$estimate, cor_followup$estimate),
  p_value = c(cor_result$p.value, cor_baseline$p.value, cor_followup$p.value)
)
write.csv2(depth_burden_cor, "results/5_arg/qc/depth_burden_correlation.csv", row.names = FALSE)

# QC Plot 10: Cleveland dot plot — ARG class prevalence at baseline vs follow-up
n_total_tp <- df_raw %>%
  filter(! Sample %in% c("HELIBA_")) %>%
  mutate(tp = case_when(str_detect(Sample, "HELIBA") ~ "Baseline",
                        str_detect(Sample, "HELIFU") ~ "Follow-up")) %>%
  group_by(tp) %>%
  summarise(n_total = n_distinct(Sample), .groups = "drop")

class_prev_dot <- class_prevalence %>%
  left_join(n_total_tp, by = c("timepoint" = "tp")) %>%
  mutate(prevalence_pct = (n_samples / n_total) * 100) %>%
  group_by(Class) %>%
  mutate(base_prev = prevalence_pct[timepoint == "Baseline"][1]) %>%
  ungroup() %>%
  mutate(Class = fct_reorder(Class, base_prev),
         Class = gsub("_", " ", Class),
         timepoint = factor(timepoint, levels = c("Baseline", "Follow-up")))

n_arg_classes <- n_distinct(class_prev_dot$Class)

(p10 <- ggplot(class_prev_dot, aes(y = Class, x = prevalence_pct, colour = Class)) +
  geom_line(aes(group = Class), colour = "gray75", linewidth = 0.6) +
  geom_point(aes(shape = timepoint), size = 3) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "gray55") +
  scale_shape_manual(values = c("Baseline" = 16, "Follow-up" = 1), name = "") +
  scale_x_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%"),
                     breaks = seq(0, 100, 25)) +
  {if (n_arg_classes <= 8) scale_colour_brewer(palette = "Dark2", guide = "none")
   else scale_colour_viridis_d(option = "turbo", guide = "none")} +
  theme_Publication() +
  labs(x = "Prevalence (% of samples)", y = "",
       title = "ARG Class Prevalence: Baseline vs Follow-up") +
  theme(axis.text.y = element_text(size = 9)))
ggsave("results/5_arg/qc/class_prevalence_dotplot.pdf", p10, width = 8, height = 7)
