## Figure 5 — Antimicrobial resistance gene (ARG) dynamics
##
## Panels:
##   A — Total ARG burden over time: paired violin, global decline (script 3)
##   B — ARG burden by ethnicity × timepoint: no difference at baseline or follow-up (script 3)
##   C — ARG class prevalence dumbbell: aminoglycoside/tetracycline asymmetry, stability (script 2)
##   D — Gene prevalence dumbbell: specific genes with ethnic differences at baseline (script 2)
##   E — Gene abundance dumbbell: trajectories, vancomycin convergence (script 2)
##   F — ARG Shannon diversity by ethnicity: NL higher at baseline, convergence over time (script 3)

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggsci)

source("scripts/5_arg/1_arg_descriptive_qc.R")
source("scripts/5_arg/2_arg_cross_comparisons.R")   # produces pl_E, pl_E_abund, pl_class_prev
source("scripts/5_arg/3_arg_longitudinal_lmm.R")    # produces pl_A, pl_B, pl_B_shan

## ── Assemble Figure 5 ────────────────────────────────────────────────────────

# Row 1: ARG burden summary
row1 <- ggarrange(
  pl_A, pl_B,
  ncol          = 2,
  widths        = c(0.8, 1.2),
  labels        = c("A", "B"),
  font.label    = list(size = 14, face = "bold"),
  common.legend = FALSE
)

# Row 2: class-level patterns and diversity
row2 <- ggarrange(
  pl_class_prev, pl_B_shan,
  ncol          = 2,
  widths        = c(1, 1),
  labels        = c("C", "F"),
  font.label    = list(size = 14, face = "bold"),
  common.legend = FALSE
)

# Row 3: gene-level patterns (needs more height for 20-gene dumbbells)
row3 <- ggarrange(
  pl_E, pl_E_abund,
  ncol          = 2,
  widths        = c(1, 1),
  labels        = c("D", "E"),
  font.label    = list(size = 14, face = "bold"),
  common.legend = TRUE,
  legend        = "bottom"
)

fig5 <- ggarrange(
  row1, row2, row3,
  nrow    = 3,
  heights = c(1, 1, 1.3)   # bottom row taller for 20-gene dumbbells
)

dir.create("results/5_arg", showWarnings = FALSE, recursive = TRUE)
ggsave(fig5, filename = "results/5_arg/figure5.pdf",
       width = 9, height = 12, device = cairo_pdf)

