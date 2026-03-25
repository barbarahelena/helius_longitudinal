## Test plot to check for extra lines
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

df_tot <- readRDS("data/arg_prepared_for_lmm.RDS")
statres <- read.csv2("results/5_arg/longitudinal/lmm_ethnicity_timepoint_results.csv")

# Test with first significant gene
statres_sig <- statres %>%
  filter(pval < 0.05) %>%
  arrange(pval)

nm <- statres_sig$mbname[1]
cat("Testing gene:", nm, "\n")

df_tot$mb <- log10(df_tot[[nm]] + 1)

df_means <- df_tot %>%
  group_by(EthnicityTot, timepoint) %>%
  summarise(mean = mean(mb, na.rm = TRUE), sd = sd(mb, na.rm = TRUE),
            n = n(), .groups = "drop")

cat("\ndf_means (should have 4 rows):\n")
print(df_means)

cat("\nNumber of unique ethnic groups in df_tot:", n_distinct(df_tot$EthnicityTot), "\n")
cat("Number of unique ethnic groups in df_means:", n_distinct(df_means$EthnicityTot), "\n\n")

# Create the plot
p <- ggplot() +
  geom_line(data = df_tot, aes(x = timepoint, y = mb, color = EthnicityTot, group = ID),
            alpha = 0.01, linewidth = 0.5) +
  geom_point(data = df_tot, aes(x = timepoint, y = mb, color = EthnicityTot),
             alpha = 0.01, size = 0.8) +
  geom_line(data = df_means, aes(x = timepoint, y = mean, color = EthnicityTot, group = EthnicityTot),
            alpha = 1, linewidth = 0.8) +
  geom_point(data = df_means, aes(x = timepoint, y = mean, color = EthnicityTot, group = EthnicityTot),
             alpha = 1, size = 1.3) +
  geom_errorbar(data = df_means,
                aes(ymin = mean - (sd/sqrt(n)), ymax = mean + (sd/sqrt(n)),
                    x = timepoint, color = EthnicityTot), width = 0.1) +
  scale_color_jco() +
  theme_Publication() +
  labs(x = "Timepoint", y = "log10(RPKM + 1)",
       title = paste0("Test Plot: ", nm),
       color = "")

ggsave("results/5_arg/longitudinal/TEST_PLOT.pdf", p, width = 8, height = 6)
cat("\nPlot saved to: results/5_arg/longitudinal/TEST_PLOT.pdf\n")
cat("This plot should show exactly 2 bold colored lines (one per ethnicity group)\n")
cat("plus many faint gray individual lines in the background.\n")
