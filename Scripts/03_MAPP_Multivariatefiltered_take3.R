rm(list = ls())
library(tidyverse)
library(vegan)
library(here)
library(ggplot2)
library(ggnewscale)
library(pheatmap)
library(viridis)
library(caret)

# =========================================================================
# Aesthetics & Color Setup
# =========================================================================
species_colors <- c(
  "AABR" = "#046E8F", 
  "AGLO" = "lightgreen", 
  "AHYA" = "#D44D5C",
  "ICRA" = "#462255"  
)

# Define explicit named vector for site shapes
custom_shapes <- c(
  "1"  = 23, # Diamond (border + fill)
  "3"  = 22, # Square (border + fill)
  "4"  = 15, # Solid Square
  "5"  = 16, # Solid Circle
  "7"  = 24, # Open / Bordered Triangle
  "8"  = 18, # Solid Diamond
  "9"  = 21, # Circle (border + fill)
  "10" = 1,  # Open Circle
  "11" = 17  # Solid Closed Triangle
)

category_colors <- c(
  "ABQ (Absorption)"                       = "#1F77B4", 
  "Quant (Quantum Yield)"                   = "#2CA02C", 
  "Sigma (Antenna Size)"                    = "#FF7F0E", 
  "qm (Max / Non-Photochemical Quenching)"  = "#D62728", 
  "qP (Photochemical Quenching)"            = "#9467BD", 
  "Tau1 (Transport Kinetics 1)"             = "#8C564B", 
  "Tau2 (Transport Kinetics 2)"             = "#E377C2", 
  "Connect (Connectivity)"                  = "#17BECF", 
  "Other Metric"                            = "#7F7F7F"  
)

# =========================================================================
# 1. Read & Format Raw Data
# =========================================================================
load(here("Outputs", "photophys_preprocessed_data.RData"))

# Compute Pearson correlation matrix across metrics
cor_matrix <- cor(permanova_matrix, use = "pairwise.complete.obs", method = "pearson")
cor_matrix[is.na(cor_matrix) | is.nan(cor_matrix)] <- 0

# Run findCorrelation at desired cutoff threshold
cor_cutoff <- 0.99
high_cor_idx <- findCorrelation(cor_matrix, cutoff = cor_cutoff)

# Safely identify traits to retain (handles edge case where high_cor_idx is empty)
if (length(high_cor_idx) > 0) {
  retained_trait_names <- colnames(permanova_matrix)[-high_cor_idx]
} else {
  retained_trait_names <- colnames(permanova_matrix)
}

cat("Total initial traits:", ncol(permanova_matrix), "\n")
cat("Traits retained after correlation (r =", cor_cutoff, ") filtering:", length(retained_trait_names), "\n")

# Subset filtered traits (Samples x Retained Metrics)
permanova_matrix_filtered <- permanova_matrix[, retained_trait_names, drop = FALSE]

# Ensure metadata strictly matches the exact sample order of the filtered matrix
meta_filtered <- ed_aligned %>%
  filter(SampleID_clean %in% rownames(permanova_matrix_filtered)) %>%
  arrange(match(SampleID_clean, rownames(permanova_matrix_filtered)))

# Verify row order is identical (CRITICAL for vegan functions)
stopifnot(identical(rownames(permanova_matrix_filtered), meta_filtered$SampleID_clean))

# =========================================================================
# 4. nMDS & Trait Vector Fitting
# =========================================================================
set.seed(123)
nmds_filtered <- metaMDS(
  permanova_matrix_filtered, 
  distance      = "euclidean", 
  k             = 2, 
  trymax        = 200, 
  maxit         = 500, 
  autotransform = FALSE,
  trace         = 0
)

cat("\nFiltered nMDS Stress Value:", nmds_filtered$stress, "\n")

# Fit vectors for retained traits
ef_filtered <- envfit(nmds_filtered, permanova_matrix_filtered, permutations = 999)
saveRDS(ef_filtered, here("Outputs", "ef_filtered.rds"))
#readRDS(here("Outputs/ef_filtered.rds"))


vector_scores <- as.data.frame(scores(ef_filtered, display = "vectors"))
vector_scores$Metric <- rownames(vector_scores)
vector_scores$r2     <- ef_filtered$vectors$r
vector_scores$p_val  <- ef_filtered$vectors$pvals

vector_scores <- vector_scores %>%
  left_join(metric_families, by = "Metric")

# Filter top metric drivers
top_vectors <- vector_scores %>%
  filter(p_val < 0.002) %>%
  arrange(desc(r2)) %>%
  head(10)

cat("\n--- TOP METRIC DRIVERS OF FILTERED ORDINATION ---\n")
print(top_vectors)

nmds_scores <- as.data.frame(scores(nmds_filtered, display = "sites")) %>%
  mutate(SampleID_clean = rownames(permanova_matrix_filtered)) %>%
  left_join(meta_filtered, by = "SampleID_clean")

max_site_coord   <- max(abs(c(nmds_scores$NMDS1, nmds_scores$NMDS2)))
max_vector_coord <- max(sqrt(top_vectors$NMDS1^2 + top_vectors$NMDS2^2))
arrow_mult       <- (max_site_coord * 0.8) / max_vector_coord

top_vectors_clean <- top_vectors %>%
  mutate(
    Metric_Category = Family, # Use upstream classification directly
    Trait_Family    = str_extract(Metric, "(?i)(mQuant|qqP|mqP|ppq|mpq|qP|rqm|qm|rABQ|ABQ|Sigma|Tau1|Tau2|NPQ|nCon)"),
    Phase           = str_extract(Metric, "^[A-Za-z0-9]+"),
    Clean_Label     = case_when(
      !is.na(Trait_Family) & !is.na(Phase) ~ paste0(Trait_Family, " (", Phase, ")"),
      !is.na(Trait_Family)                 ~ Trait_Family,
      TRUE                                 ~ Metric
    ),
    NMDS1_scaled = NMDS1 * arrow_mult,
    NMDS2_scaled = NMDS2 * arrow_mult
  )


# =========================================================================
# 5. Plot A spp elipses
# =========================================================================
spp <- ggplot() +
  geom_point(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = Species, shape = factor(Site)), 
    size = 2.8, alpha = 0.75
  ) + 
  stat_ellipse(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = Species, fill = Species), 
    geom = "polygon", alpha = 0.15, level = 0.95
  ) +
  scale_shape_manual(values = custom_shapes, name = "Site") +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_fill_manual(values = species_colors, name = "Species") +
  new_scale_color() +
  geom_segment(
    data = top_vectors_clean, 
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
    arrow = arrow(length = unit(0.20, "cm")), linewidth = 0.85
  ) +
  geom_text(
    data = top_vectors_clean,
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Clean_Label),
    size = 3
  ) +
  scale_color_manual(values = category_colors, name = "Trait Family") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(plot.margin = margin(15, 25, 15, 25, "pt"), legend.position = "right") +
  labs(title = "Species Ellipses", x = "nMDS Dimension 1", y = "nMDS Dimension 2")
# =========================================================================
# 4. PLOT B: HCLUST (K=3) ELLIPSES
# =========================================================================

# Ensure K_3 factor exists in meta_filtered
cluster_colors <- c("1" = "orange", "2" = "violet", "3" = "grey")

kplot <- ggplot() +
  geom_point(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = Species, shape = factor(Site)), 
    size = 2.8, alpha = 0.75
  ) + 
  scale_shape_manual(values = custom_shapes, name = "Site") +
  scale_color_manual(values = species_colors, name = "Species") +
  new_scale_color() +
  new_scale_fill() +
  stat_ellipse(
    data = nmds_scores %>% filter(!is.na(K_3)), 
    aes(x = NMDS1, y = NMDS2, color = factor(K_3), fill = factor(K_3)), 
    geom = "polygon", alpha = 0.15, level = 0.95
  ) +
  scale_color_manual(values = cluster_colors, name = "Hclust Group (K_3)") +
  scale_fill_manual(values = cluster_colors, name = "Hclust Group (K_3)") +
  new_scale_color() +
  geom_segment(
    data = top_vectors_clean, 
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
    arrow = arrow(length = unit(0.20, "cm")), linewidth = 0.85
  ) +
  geom_text(
    data = top_vectors_clean,
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Clean_Label),
    size = 3
  ) +
  scale_color_manual(values = category_colors, name = "Trait Family") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(plot.margin = margin(15, 25, 15, 25, "pt"), legend.position = "right") +
  labs(title = "Hclust Groupings (K=3)", x = "nMDS Dimension 1", y = NULL)

combined_plot <- (spp + kplot) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = paste0("Filtered (r = ", cor_cutoff, ", nMDS Stress = ", round(nmds_filtered$stress, 3), ")"),
    tag_levels = 'A'
  ) &
  theme(legend.position = "right")

ggsave(
  filename = here("Plots", paste0("nMDS_combined_r_10", cor_cutoff * 100, ".png")),
  plot     = combined_plot,
  width    = 16.0,
  height   = 9.0,
  dpi      = 300,
  units    = "in"
)

##########################
#hclust
##########################
set.seed(123)
fit_cbb <- pvclust(
  t(permanova_matrix_filtered), 
  method.hclust = "ward.D",
  method.dist   = "canberra",
  nboot         = 1000,
  parallel      = TRUE
)
# Save fit object
saveRDS(fit_cbb, here("Outputs", "fit_pvclust_99.rds"))

#or load
readrds(here("Outputs", "fit_pvclust_99.rds"))
# Grab hclust tree directly from pvclust fit
my_hclust <- fit_cbb$hclust
original_labels <- my_hclust$labels # SampleIDs (e.g., "11_AABR_1")

# Subset matrix to retain only collinearity-filtered metrics
# Transpose so metrics are rows and samples are columns
heatmap_matrix <- t(permanova_matrix_filtered)

# Align columns to match leaf order of the dendrogram
heatmap_matrix <- heatmap_matrix[, match(original_labels, colnames(heatmap_matrix))]

# Z-score scale across metrics (rows) and cap outliers [-3, 3]
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
heatmap_matrix_scaled[heatmap_matrix_scaled > 3]  <- 3
heatmap_matrix_scaled[heatmap_matrix_scaled < -3] <- -3

# Drop raw baseline metrics (Fm/Fo) if present
keep_metrics <- !rownames(heatmap_matrix_scaled) %in% c("Fm", "Fo", "fm", "fo", "F0", "FM")
heatmap_matrix_scaled <- heatmap_matrix_scaled[keep_metrics, , drop = FALSE]

# =========================================================================
# 3. PREPARE ROW ANNOTATIONS (METRIC FAMILIES)
# =========================================================================

# Filter metric_families table to match active matrix rows
annotation_row <- metric_families %>%
  filter(Metric %in% rownames(heatmap_matrix_scaled)) %>%
  select(Metric, Family) %>%
  column_to_rownames("Metric")

# Sort Matrix and Row Annotations alphabetically by Family
sorted_order <- order(annotation_row$Family)
heatmap_matrix_scaled <- heatmap_matrix_scaled[sorted_order, ]
annotation_row         <- annotation_row[sorted_order, , drop = FALSE]

# Define palette for metric families
metric_colors <- c(
  "Quant (Quantum Yield)"             = "#1F77B4", # Strong Blue
  "Sigma (Antenna Size)"              = "#FF7F0E", # Vivid Orange
  "Connect (Connectivity)"            = "#2CA02C", # Leaf Green
  "Connect (Normalized Curves)"       = "#2CA02C", 
  "Tau1 (Transport Kinetics 1)"       = "#D62728", # Crimson Red
  "Tau2 (Transport Kinetics 2)"       = "#9467BD", # Purple
  "NPQ (Non-Photochemical Quenching)" = "#8C564B", # Warm Brown
  "qP (Photochemical Quenching)"      = "#E377C2", # Pink/Rose
  "ABQ (Absorption)"                  = "#17BECF", # Cyan
  "qm (Max Quenching)"                = "#BCBD22", # Olive Yellow
  "OJIP Kinetics (Derivatives)"       = "#98DF8A", # Light Green
  "OJIP Ratios & Area"                = "#FFBB78", # Light Orange
  "Other Metric"                      = "#7F7F7F"  # Slate Grey
)

# =========================================================================
# 4. PREPARE COLUMN ANNOTATIONS (SAMPLES)
# =========================================================================

# Match ed_aligned metadata directly to tree labels
annotation_col <- ed_aligned %>%
  filter(SampleID_clean %in% colnames(heatmap_matrix_scaled)) %>%
  arrange(match(SampleID_clean, colnames(heatmap_matrix_scaled))) %>%
  column_to_rownames("SampleID_clean") %>%
  select(Species, any_of(c("K_3")), ED50) %>%
  mutate(across(c(Species, any_of(c("K_3"))), as.factor))

# Combine color maps for pheatmap
ann_colors <- list(
  Family  = metric_colors,
  Species = species_colors
)

# =========================================================================
# 5. RENDER & EXPORT HEATMAP
# =========================================================================

my_palette <- colorRampPalette(c("darkblue", "white", "darkred"))(100)
my_breaks  <- seq(-3, 3, length.out = 101)

jpeg(
  file   = here("Plots", "Sorted_Grouped_Physiological_Heatmap_Filtered1.jpg"), 
  width  = 12, 
  height = 10, 
  units  = "in", 
  res    = 300
)

pheatmap(
  heatmap_matrix_scaled,
  annotation_colors = ann_colors,
  annotation_col    = annotation_col,
  annotation_row    = annotation_row,
  show_colnames     = FALSE,         
  show_rownames     = FALSE,         
  cluster_rows      = FALSE, # false= alphabetical Family sorting
  cluster_cols      = my_hclust, # Direct pvclust tree      
  color             = my_palette,            
  breaks            = my_breaks,            
  main              = "Collinearity-Filtered MAPP Fluorescence Metrics vs. Clusters",
  margins           = c(5, 20)       
)

dev.off()

