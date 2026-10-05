## Supplementary Figure 11 — A. putredinis strain retention & microdiversity
## Assembles panels from two already-existing analysis scripts by sourcing
## each and capturing its ggplot objects (same pattern as
## scripts/4_functional_change/assemble_figure.R), rather than duplicating
## their plotting code. Run order matches the pipeline dependency: script 1
## writes instrain_strain_retention.csv, which script 3 reads.
##   A  Genetic distance (popANI) by ethnicity           <- 1_instrain_strain_retention.R (pl_popani)
##   B  Divergence distribution vs. between-person floor <- 1_instrain_strain_retention.R (pl_snps_hist)
##   C  Microdiversity by ethnicity, cross-sectional      <- 3_instrain_microdiversity.R (pl_crosssectional)
##   D  Microdiversity trajectory, retained strains       <- 3_instrain_microdiversity.R (pl_trajectory)
##   E  Microdiversity change by ethnicity                <- 3_instrain_microdiversity.R (pl_delta)
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

library(ggpubr)
library(tidyverse)

out_dir <- "results/3_species_change/5_instrain"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### 1. Source strain retention — capture panels A, B ####
source("scripts/3_species_change/5_instrain/1_instrain_strain_retention.R")
pl_A <- pl_popani
pl_B <- pl_snps_hist

#### 2. Source microdiversity — capture panels C, D, E ####
source("scripts/3_species_change/5_instrain/3_instrain_microdiversity.R")
pl_C <- pl_crosssectional
pl_D <- pl_trajectory
# Original title carries a "B. " prefix meant for pl_delta's own standalone
# A/B figure (instrain_microdiversity.pdf) — drop it here, this is panel E.
pl_E <- pl_delta + ggtitle("Change by ethnicity")

#### 3. Assemble — 3 rows (2/1/2 panels) for a portrait page ####
row1 <- ggarrange(pl_A, pl_B, ncol = 2, labels = c("A", "B"))
row2 <- ggarrange(pl_C, ncol = 1, labels = "C")
row3 <- ggarrange(pl_D, pl_E, ncol = 2, widths = c(1.3, 1), labels = c("D", "E"))

suppl11 <- ggarrange(row1, row2, row3, nrow = 3, heights = c(1, 1, 1))

ggsave(file.path(out_dir, "suppl_figure11_microdiversity.pdf"), suppl11, width = 10, height = 14)
cat("\nSupplementary Figure 11 saved to:", file.path(out_dir, "suppl_figure11_microdiversity.pdf"), "\n")
