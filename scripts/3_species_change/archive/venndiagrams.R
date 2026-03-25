


library(VennDiagram)
library(tidyverse)

# Data
feat_ba <- read.csv('mlmodels/eth_base_16s/output_XGB_class_eth_baseline_16s_2024_08_19__20-32-34/feature_importance.txt', sep = "\t") # baseline
feat_fu <- read.csv('mlmodels/eth_fu_16s/output_XGB_class_eth_followup_16s_2024_08_19__21-53-15/feature_importance.txt', sep = "\t") # follow-up

head(feat_ba)
head(feat_fu)
featba20 <- feat_ba |> slice(1:25) |> pull(FeatName)
featfu20 <- feat_fu |> slice(1:25) |> pull(FeatName)

# Union of all elements
all_elems <- union(featba20, featfu20)

# Create binary matrix
mat <- data.frame(
  baseline = as.integer(all_elems %in% featba20),
  followup = as.integer(all_elems %in% featfu20),
  row.names = all_elems
)
mat
dim(mat)

upset(mat, sets = c("baseline", "followup"))

venn.plot <- venn.diagram(x = list(featba20, featfu20),
                            category.names = c("baseline", "followup"),
                            filename = "results/ml_figures/venndiagram_16s.png",         # keep in R instead of saving to file
                            imagetype = "png",
                            fill = c("royalblue", "firebrick"),
                            alpha = 0.8
                          )


# Data
feat_ba <- read.csv('mlmodels/eth_base/output_XGB_class_eth_baseline_2024_08_19__21-02-46/feature_importance.txt', sep = "\t") # baseline
feat_fu <- read.csv('mlmodels/eth_fu/output_XGB_class_eth_followup_2024_08_20__09-08-32/feature_importance.txt', sep = "\t") # follow-up

head(feat_ba)
head(feat_fu)
featba20 <- feat_ba |> slice(1:25) |> pull(FeatName)
featfu20 <- feat_fu |> slice(1:25) |> pull(FeatName)

# Union of all elements
all_elems <- union(featba20, featfu20)

# Create binary matrix
mat <- data.frame(
  baseline = as.integer(all_elems %in% featba20),
  followup = as.integer(all_elems %in% featfu20),
  row.names = all_elems
)
mat
dim(mat)

upset(mat, sets = c("baseline", "followup"))

venn.plot <- venn.diagram(x = list(featba20, featfu20),
                            category.names = c("baseline", "followup"),
                            filename = "results/ml_figures/venndiagram_sg.png",         # keep in R instead of saving to file
                            imagetype = "png",
                            fill = c("royalblue", "firebrick"),
                            alpha = 0.8
                          )

# Data
feat_ba <- read.csv('mlmodels/eth_base_pathways/output_XGB_class_eth_baseline_2024_08_30__21-45-31/feature_importance.txt', sep = "\t") # baseline
feat_fu <- read.csv('mlmodels/eth_fu_pathways/output_XGB_class_eth_followup_2024_08_30__23-09-43/feature_importance.txt', sep = "\t") # follow-up

head(feat_ba)
head(feat_fu)
featba20 <- feat_ba |> slice(1:25) |> pull(FeatName)
featfu20 <- feat_fu |> slice(1:25) |> pull(FeatName)

# Union of all elements
all_elems <- union(featba20, featfu20)

# Create binary matrix
mat <- data.frame(
  baseline = as.integer(all_elems %in% featba20),
  followup = as.integer(all_elems %in% featfu20),
  row.names = all_elems
)
mat
dim(mat)

upset(mat, sets = c("baseline", "followup"))

venn.plot <- venn.diagram(x = list(featba20, featfu20),
                            category.names = c("baseline", "followup"),
                            filename = "results/ml_figures/venndiagram_sg_pw.png",         # keep in R instead of saving to file
                            imagetype = "png",
                            fill = c("royalblue", "firebrick"),
                            alpha = 0.8
                          )
