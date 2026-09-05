library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(sparseLDA)

# Load global prep script for palettes and themes
source(here("Scripts/00_visualization_prep.R"))
load(here("Outputs", "photophys_preprocessed_data.RData"))


# Ensure Cluster_Label is properly mapped on your metadata
meta_lda <- meta_filteredK %>%
  mutate(
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels)
  )

# --- 1. LDA for Species ---
# Fits sparse discriminant analysis
sda_fit <- sda(
  x = as.matrix(permanova_matrix_filtered), 
  y = as.factor(meta_lda$Species)
)
# Extract LD1 coefficients (loadings)
loadings_spp <- as.data.frame(lda_spp$scaling) %>%
  tibble::rownames_to_column(var = "Metric") %>%
  pivot_longer(cols = starts_with("LD"), names_to = "Axis", values_to = "Loading")

# --- 2. LDA for Clusters (K_3) ---
lda_cluster <- lda(permanova_matrix_filtered, grouping = meta_lda$Cluster_Label)

# Extract LD1 coefficients
loadings_cluster <- as.data.frame(lda_cluster$scaling) %>%
  tibble::rownames_to_column(var = "Metric") %>%
  pivot_longer(cols = starts_with("LD"), names_to = "Axis", values_to = "Loading")