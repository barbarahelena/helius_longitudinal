## Color Palette Options for Many Categories (>10)

library(ggsci)
library(viridis)
library(RColorBrewer)

# Summary of palette capacities:
# ggsci palettes (LIMITED):
#   - scale_fill_jco()      : ~10 colors
#   - scale_fill_bmj()      : ~7 colors
#   - scale_fill_aaas()     : ~10 colors
#   - scale_fill_nejm()     : ~8 colors
#   - scale_fill_lancet()   : ~9 colors
#   - scale_fill_npg()      : ~10 colors
#
# Viridis palettes (UNLIMITED discrete values):
#   - scale_fill_viridis_d(option = "viridis")  # default, purple-blue-green-yellow
#   - scale_fill_viridis_d(option = "magma")    # dark purple-magenta-yellow
#   - scale_fill_viridis_d(option = "inferno")  # dark purple-red-yellow
#   - scale_fill_viridis_d(option = "plasma")   # purple-red-yellow
#   - scale_fill_viridis_d(option = "cividis")  # blue-yellow (colorblind friendly)
#   - scale_fill_viridis_d(option = "turbo")    # rainbow-like (good discrimination)
#
# RColorBrewer (VARIABLE limits):
#   - scale_fill_brewer(palette = "Set1")      : max 9 colors
#   - scale_fill_brewer(palette = "Set2")      : max 8 colors
#   - scale_fill_brewer(palette = "Set3")      : max 12 colors
#   - scale_fill_brewer(palette = "Paired")    : max 12 colors
#
# Manual/Custom palettes (UNLIMITED):
#   - scale_fill_manual(values = rainbow(n))
#   - scale_fill_manual(values = hcl.colors(n, palette = "Dynamic"))

## Recommended palettes for ARG analysis (23 classes)

# Option 1: Viridis turbo (best color discrimination for many categories)
scale_fill_viridis_d(option = "turbo")

# Option 2: Manual rainbow
n_classes <- 23
scale_fill_manual(values = rainbow(n_classes))

# Option 3: HCL colors (perceptually uniform)
scale_fill_manual(values = hcl.colors(n_classes, palette = "Dynamic"))

# Option 4: Custom palette combining multiple sets
custom_palette <- c(
  ggsci::pal_jco()(10),           # First 10 colors
  ggsci::pal_lancet()(9),         # Next 9 colors
  ggsci::pal_nejm()(8)            # Additional colors
)
scale_fill_manual(values = custom_palette[1:n_classes])

## Example usage in a plot:

# Test with your ARG data
if(FALSE){
  # Load data
  df_raw <- rio::import("data/all_samples.merged_arg_counts.tsv")

  class_prevalence <- df_raw %>%
    group_by(Class) %>%
    summarise(n_observations = n(), .groups = "drop_last") %>%
    arrange(-n_observations) %>%
    mutate(Class = fct_reorder(Class, n_observations))

  # Compare different palettes

  # Viridis turbo (recommended)
  p1 <- ggplot(class_prevalence, aes(x = Class, y = n_observations, fill = Class)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_viridis_d(option = "turbo") +
    theme_minimal() +
    labs(title = "Viridis Turbo (good discrimination)") +
    theme(legend.position = "none")

  # Viridis plasma
  p2 <- ggplot(class_prevalence, aes(x = Class, y = n_observations, fill = Class)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_viridis_d(option = "plasma") +
    theme_minimal() +
    labs(title = "Viridis Plasma") +
    theme(legend.position = "none")

  # HCL Dynamic
  p3 <- ggplot(class_prevalence, aes(x = Class, y = n_observations, fill = Class)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = hcl.colors(23, palette = "Dynamic")) +
    theme_minimal() +
    labs(title = "HCL Dynamic") +
    theme(legend.position = "none")

  # Rainbow
  p4 <- ggplot(class_prevalence, aes(x = Class, y = n_observations, fill = Class)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = rainbow(23)) +
    theme_minimal() +
    labs(title = "Rainbow") +
    theme(legend.position = "none")

  library(gridExtra)
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}

## Tip: For too many categories, consider grouping
# Instead of showing all 23 classes, group rare ones into "Other"
# Example:
group_rare_classes <- function(df, min_count = 100){
  df %>%
    mutate(
      Class_grouped = case_when(
        n_observations >= min_count ~ Class,
        TRUE ~ "Other (rare)"
      )
    )
}
