rm(list = ls())
library(tidyverse)
library(vegan)
library(here)
library(ggplot2)
library(ggnewscale)
library(pheatmap)
library(viridis)
library(caret)
library(pvclust)
library(patchwork)

source(here("Scripts/00_visualization_prep.R"))
# =========================================================================
# 1. Read in and format raw data
# =========================================================================
load(here("Outputs", "photophys_preprocessed_data.RData"))

# calc Pearson correlation matrix across metrics
cor_matrix <- cor(permanova_matrix, use = "pairwise.complete.obs", method = "pearson")
cor_matrix[is.na(cor_matrix) | is.nan(cor_matrix)] <- 0

# run findCorrelation at 0.99
cor_cutoff <- 0.99
high_cor_idx <- findCorrelation(cor_matrix, cutoff = cor_cutoff)

# id traits to retain (handles edge case where high_cor_idx is empty)
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


# =================================
# hclust
# =================================
# set.seed(123)
# fit_cbb <- pvclust(
#   t(permanova_matrix_filtered),
#   method.hclust = "ward.D",
#   method.dist   = "canberra",
#   nboot         = 1000,
#   parallel      = TRUE
# )
# #Save fit object
# saveRDS(fit_cbb, here("Outputs", "fit_pvclust_99.rds"))

#or load
fit_cbb<-readRDS(here("Outputs", "fit_pvclust_99.rds"))
# get hclust tree from pvclust fit
my_hclust <- fit_cbb$hclust
original_labels <- my_hclust$labels # SampleIDs (e.g., "11_AABR_1")

# extract cluster across K = 2 through K = 15
k_matrix <- cutree(my_hclust, k = 2:15)

# convert into a clean lookup dataframe with "K_2", "K_3" columns
k_lookup <- as.data.frame(k_matrix) %>%
  rename_with(~ paste0("K_", 2:15)) %>%
  rownames_to_column("SampleID_clean") %>%
  mutate(across(starts_with("K_"), ~ factor(paste0("Cluster_", .x))))

# join all K_2 through K_15 cluster columns into meta_filtered
meta_filtered <- meta_filtered %>%
  left_join(k_lookup, by = "SampleID_clean")

# verify row alignment
stopifnot(identical(rownames(permanova_matrix_filtered), meta_filtered$SampleID_clean))

cat("\n--- K_3 CLUSTER DISTRIBUTION ACROSS SPECIES ---\n")
print(table(meta_filtered$K_3, meta_filtered$Species))

#save in output
saveRDS(meta_filtered, here("Outputs", "master_cluster_metadata_K2_K15.rds"))
write.csv(meta_filtered, here("Outputs", "master_cluster_metadata_K2_K15.csv"), 
  row.names = FALSE)

# transpose so metrics are rows and samples are columns
heatmap_matrix <- t(permanova_matrix_filtered)

# align columns to match leaf order of the dendrogram
heatmap_matrix <- heatmap_matrix[, match(original_labels, colnames(heatmap_matrix))]

# Z-score scale across metrics (rows) and cap outliers [-3, 3]
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
heatmap_matrix_scaled[heatmap_matrix_scaled > 3]  <- 3
heatmap_matrix_scaled[heatmap_matrix_scaled < -3] <- -3

# drop raw baseline metrics (Fm/Fo) if present (triple check)
keep_metrics <- !rownames(heatmap_matrix_scaled) %in% c("Fm", "Fo", "fm", "fo", "F0", "FM")
heatmap_matrix_scaled <- heatmap_matrix_scaled[keep_metrics, , drop = FALSE]

# =========================================================================
# 3. PREPARE ROW ANNOTATIONS (METRIC FAMILIES)
# =========================================================================

# Filter metric_families table to match active matrix rows
annotation_row <- metric_families %>%
  filter(Metric %in% rownames(heatmap_matrix_scaled)) %>%
  dplyr::select(Metric, Family) %>%
  column_to_rownames("Metric")

# Sort Matrix and Row Annotations alphabetically by Family
sorted_order <- order(annotation_row$Family)
heatmap_matrix_scaled <- heatmap_matrix_scaled[sorted_order, ]
annotation_row         <- annotation_row[sorted_order, , drop = FALSE]


# =========================================================================
# 4. PREPARE COLUMN ANNOTATIONS
# =========================================================================

# Match ed_aligned metadata directly to tree labels
annotation_col <- meta_filtered %>%
  filter(SampleID_clean %in% colnames(heatmap_matrix_scaled)) %>%
  arrange(match(SampleID_clean, colnames(heatmap_matrix_scaled))) %>%
  column_to_rownames("SampleID_clean") %>%
  dplyr::select(Species, any_of(c("K_3")), ED50) %>%
  mutate(across(c(Species, any_of(c("K_3"))), as.factor))

# =========================================================================
# 5. HEATMAP
# =========================================================================

my_palette <- colorRampPalette(c("darkblue", "white", "darkred"))(100)
my_breaks  <- seq(-3, 3, length.out = 101)

jpeg(
  file   = here("Plots", "Sorted_Physiological_Heatmap_Filtered1.jpg"), 
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
  cluster_rows      = TRUE, # false= alphabetical Family sorting
  cluster_cols      = my_hclust, # Direct pvclust tree      
  color             = my_palette,            
  breaks            = my_breaks,            
  main              = "Collinearity-Filtered MAPP Fluorescence Metrics vs. Clusters",
  margins           = c(5, 20)       
)

dev.off()

##############################################################################
#### MULTIVARIATE ANALYSIS####################################################
##############################################################################
# =========================================================================
# PERMANOVA & Dispersion Tests
# =========================================================================
# # spp PERMANOVA
# permanova_filtered <- adonis2(
#   permanova_matrix_filtered ~ Species, 
#   data         = meta_filtered, 
#   method       = "euclidean", 
#   permutations = 999
# )
# cat("\n--- FILTERED PERMANOVA RESULTS ---\n")
# print(permanova_filtered)

# Species and site global PERMANOVA
permanova_filtered <- adonis2(
  permanova_matrix_filtered ~ factor(Site) + Species, 
  data         = meta_filtered, 
  method       = "euclidean", 
  by = "terms", #site and spp
  permutations = 999
)
cat("\n--- FILTERED PERMANOVA RESULTS w spp and site ---\n")
print(permanova_filtered)

# Homogeneity test of Multivariate Dispersion (BETADISPER)
dist_filtered <- vegdist(permanova_matrix_filtered, method = "euclidean")
bd_filtered   <- betadisper(dist_filtered, meta_filtered$Species)
bd_filtered_site  <- betadisper(dist_filtered, meta_filtered$Site)
anova_bd      <- anova(bd_filtered)
anova_bd_site <- anova(bd_filtered_site)

cat("\n--- FILTERED species DISPERSION TEST p-value:", anova_bd[1, "Pr(>F)"], "---\n")

cat("\n--- FILTERED site DISPERSION TEST p-value:", anova_bd_site[1, "Pr(>F)"], "---\n")
#both pass

#######################
# Pairwise PERMANOVA ##
#######################
#set up function
pairwise.adonis2 <- function(x, data, strata = NULL, nperm = 999, ... ) {
  ststri <- ifelse(is.null(strata), 'Null', strata)
  fostri <- as.character(x)
  
  x1 <- x
  lhs <- eval(x1[[2]], environment(x1), globalenv())
  environment(x1) <- environment()
  rhs <- x1[[3]]
  
  x1[[2]] <- NULL
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE)
  co <- combn(unique(as.character(rhs.frame[, 1])), 2)
  
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres, paste(co[1, elem], co[2, elem], sep = '_vs_'))
  }
  
  res <- vector(mode = "list", length = length(nameres))
  names(res) <- nameres
  res['parent_call'] <- list(paste(fostri[2], fostri[1], fostri[3], ', strata =', ststri, ', permutations', nperm))
  
  for(elem in 1:ncol(co)){
    if(inherits(eval(lhs), 'dist')){
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[, 1] %in% c(co[1, elem], co[2, elem]),
                                           rhs.frame[, 1] %in% c(co[1, elem], co[2, elem])])
    } else {
      xred <- eval(lhs)[rhs.frame[, 1] %in% c(co[1, elem], co[2, elem]), ]
    }
    
    mdat1 <- data[rhs.frame[, 1] %in% c(co[1, elem], co[2, elem]), ]
    
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred', as.character(rhs), sep = '~'))
    } else {
      xnew <- as.formula(paste('xred', paste(rhs[-1], collapse = as.character(rhs[1])), sep = '~'))
    }
    
    if(is.null(strata)){
      ad <- adonis2(xnew, data = mdat1, ... )
    } else {
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[, ststri])
      ad <- adonis2(xnew, data = mdat1, permutations = perm, ... )
    }
    
    res[nameres[elem + 1]] <- list(ad[1:5])
  }
  class(res) <- c("pwadstrata", "list")
  return(res)
}
# 1. Run Pairwise PERMANOVA for species only (dont want bc site explains 17% of variance per global permanova)
# pw_filtered <- pairwise.adonis2(permanova_matrix_filtered ~ Species, data = meta_filtered, method = "euclidean")
# cat("\n--- FILTERED PAIRWISE PERMANOVA ---\n")
# print(pw_filtered)

# 2. Run Pairwise PERMANOVA for spp, stratified by site
pw_species_strata <- pairwise.adonis2(
  permanova_matrix_filtered ~ Species, 
  data    = meta_filtered, 
  strata  = "Site", 
  method  = "euclidean",
  nperm   = 999
)
# 
cat("\n--- PAIRWISE PERMANOVA (STRATIFIED BY SITE) ---\n")
print(pw_species_strata)

#3. Run pairwise comparison of site
# pw_filtered_formula <- pairwise.adonis2(
#   permanova_matrix_filtered ~  factor(Site), 
#   data   = meta_filtered,
#   strata  = "Species",
#   method = "euclidean",
#   nperm  = 999
# )
# 
# cat("\n--- PAIRWISE PERMANOVA (PARTIALING OUT SITE) ---\n")
# print(pw_filtered_formula)
pw_filtered_strata<-pw_species_strata

pw_summary_df <- map_dfr(
  names(pw_filtered_strata)[names(pw_filtered_strata) != "parent_call"], 
  function(pair_name) {
    res <- pw_filtered_strata[[pair_name]]
    data.frame(
      Comparison = pair_name,
      Df         = res$Df[1],
      SumOfSqs   = round(res$SumOfSqs[1], 2),
      R2         = round(res$R2[1], 4),
      F_stat     = round(res$F[1], 2),
      p_val      = res$`Pr(>F)`[1]
    )
  }
) %>%
  filter(!is.na(p_val)) %>%
  mutate(
    p_adj = round(p.adjust(p_val, method = "BH"), 4)
  )

cat("\n--- PAIRWISE PERMANOVA SUMMARY TABLE ---\n")
print(pw_summary_df)


# =========================================================================
# 3. nMDS & Trait Vector Fitting
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
#ef_filtered <- envfit(nmds_filtered, permanova_matrix_filtered, permutations = 999)
#saveRDS(ef_filtered, here("Outputs", "ef_filtered.rds"))
ef_filtered <- readRDS(here("Outputs", "ef_filtered.rds"))

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
    Metric_Category = Family, 
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

############################################
#plot w no ellipses
p <- ggplot() +
  geom_point(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = Species, shape = factor(Site)), 
    size = 2.8, alpha = 0.75
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
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Metric),
    size = 3
  ) +
  scale_color_manual(values = metric_colors, name = "Trait Family") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(plot.margin = unit(c(15, 25, 15, 25), "pt"))+
  labs( x = "nMDS Dimension 1", y = "nMDS Dimension 2")
p

# =========================================================================
# plot spp ellipses
# =========================================================================
#calc centroids
spp_centroids <- nmds_scores %>%
  filter(!is.na(Species)) %>%
  group_by(Species) %>%
  summarize(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  )
spp <- ggplot() +
  geom_point(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = Species, shape = factor(Site)), 
    size = 2.8, alpha = 0.75
  ) + 
  stat_ellipse(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = Species, fill = Species), 
    geom = "polygon", alpha = 0.1, level = 0.95
  ) +
  #centroid :
  geom_point(
    data = spp_centroids,
    aes(x = NMDS1, y = NMDS2, color = Species, fill = Species),
    shape = 4,        
     #color = "black",  
    size = 5,        
    stroke = 4,    
    show.legend = FALSE
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
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Metric),
    size = 3
  ) +
  scale_color_manual(values = metric_colors, name = "Trait Family") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(plot.margin = unit(c(15, 25, 15, 25), "pt"))+
  labs(title = "Species Ellipses", x = "nMDS Dimension 1", y = "nMDS Dimension 2")

spp
# =========================================================================
# plot by K-3 ellipses
# =========================================================================
#prep labels
nmds_scores <- nmds_scores %>%
  mutate(
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels),
    Site          = factor(as.character(Site), levels = c("1", "3", "4", "5", "7", "9", "10", "11"))
  )
#calc centroids
k3_centroids <- nmds_scores %>%
  filter(!is.na(Cluster_Label)) %>%
  group_by(Cluster_Label) %>%
  summarize(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  )

kplot <- ggplot() +
  # Points: Color by Species, Shape by Site
  geom_point(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = Species, shape = Site), 
    size = 2.8, 
    alpha = 0.75
  ) + 
  scale_shape_manual(values = custom_shapes, name = "Site") +
  scale_color_manual(values = species_colors, limits = species_order, labels = species_labels, name = "Species") +
  
  # Scale layer for Ellipses (K_3 functional clusters)
  new_scale_color() +
  new_scale_fill() +
  
  # 95% Confidence Ellipses
  stat_ellipse(
    data = nmds_scores %>% filter(!is.na(Cluster_Label)),
    aes(x = NMDS1, y = NMDS2, color = Cluster_Label, fill = Cluster_Label),
    geom = "polygon", 
    alpha = 0.1, 
    level = 0.95
  ) +
  
  # Centroid Markers
  geom_point(
    data = k3_centroids,
    aes(x = NMDS1, y = NMDS2, color = Cluster_Label),
    shape = 4,        
    size = 5,          
    stroke = 4,      
    show.legend = FALSE
  ) +
  scale_color_manual(values = cluster_colors, limits = target_levels, name = "Photophysiological Cluster") +
  scale_fill_manual(values = cluster_colors, limits = target_levels, name = "Photophysiological Cluster") +
  
  # Scale layer for envfit vectors
  new_scale_color() +
  geom_segment(
    data = top_vectors_clean, 
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
    arrow = arrow(length = unit(0.20, "cm")), 
    linewidth = 0.85
  ) +
  geom_text(
    data = top_vectors_clean,
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Metric),
    size = 3
  ) +
  scale_color_manual(values = metric_colors, name = "Trait Family", drop = FALSE) +
  
  # Formatting
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    plot.margin  = unit(c(15, 25, 15, 25), "pt"),
    legend.text  = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "Hclust Groupings (K = 3)", 
    x     = "nMDS Dimension 1", 
    y     = "nMDS Dimension 2"
  )

print(kplot)

combined_plot <- (spp + kplot) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = paste0("Filtered (r = ", cor_cutoff, ", nMDS Stress = ", round(nmds_filtered$stress, 3), ")"),
    tag_levels = 'A'
  ) &
  theme(legend.position = "right")

print(combined_plot)

ggsave(
  filename = here("Plots", paste0("Fig4_nMDS_filtered_top10",".png")),
  plot     = combined_plot,
  width    = 16.0,
  height   = 9.0,
  dpi      = 300,
  units    = "in"
)

#all 3
# combined_plot3 <- p + spp + kplot + 
#   plot_layout(
#     ncol   = 3, 
#     guides = "collect" # Keeps specific legends per panel; use "collect" to combine shared legends
#   ) + 
#   plot_annotation(
#     title    = "Multivariate Photophysiological Separation across Coral Species and Clusters",
#     tag_levels = "A" # Automatically labels panels A, B, C
#   )
# 
# print(combined_plot3)
# 
# ggsave(
#   filename = here("Plots", paste0("nMDS_3combined",".png")),
#   plot     = combined_plot3,
#   width    = 16.0,
#   height   = 9.0,
#   dpi      = 300,
#   units    = "in"
# )

# # 1. Clear stuck graphics devices
while(!is.null(dev.list())) dev.off()

# 2. Re-open a clean graphics device
dev.new()

# 3. Print your plot object directly
print(spp)