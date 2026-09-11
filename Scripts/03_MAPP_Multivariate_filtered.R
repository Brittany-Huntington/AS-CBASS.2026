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
library(factoextra)

source(here("Scripts/00_visualization_prep.R"))
load(here("Outputs", "photophys_preprocessed_data.RData"))


##############################################################################
#### MULTIVARIATE ANALYSIS####################################################
##############################################################################
# =========================================================================
# PERMANOVA & Dispersion Tests
# =========================================================================
# Homogeneity test of Multivariate Dispersion (BETADISPER)
dist_filtered <- vegdist(permanova_matrix_filtered, method = "euclidean")
bd_filtered   <- betadisper(dist_filtered, meta_filteredK$Species)
bd_filtered_site  <- betadisper(dist_filtered, meta_filteredK$Site)
anova_bd      <- anova(bd_filtered)
anova_bd_site <- anova(bd_filtered_site)

cat("\n--- FILTERED species DISPERSION TEST p-value:", anova_bd[1, "Pr(>F)"], "---\n")
#pass
cat("\n--- FILTERED site DISPERSION TEST p-value:", anova_bd_site[1, "Pr(>F)"], "---\n")
#does NOT pass

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
# permanova_spp_filtered <- adonis2(
#   permanova_matrix_filtered ~ Species, 
#   data         = meta_filtered, 
#   method       = "euclidean", 
#   by = "terms", #site and spp
#   permutations = 999
# )
# cat("\n--- FILTERED PERMANOVA RESULTS w spp  ---\n")
# print(permanova_spp_filtered)

#BH notes:  vegan can't do random effects AND the order of the factors matter when running by = "terms". So this model design asks:
#how much multivariate variation in the photophys data is explained first by Site?, and then
#how much additional variation is explained by Species, after Site has already be accounted for?
#in other words: "Does a species effect persist after accounting for site-to-site variation?"


#If instead you what to ask: Does photophys differences exists among species, after accounting for differences among sites?
#and vice versa for look at site; both factors get their own 'adjusted' test
#then us by 'margin' to look at partial-effects

#round 2: using by 'margin' to get partial effects of each term rather than sequential sum of squares
permanova_spp_site <- adonis2(
  permanova_matrix_filtered ~ Species + factor(Site),
  data = meta_filtered,
  method = "euclidean",
  by = "margin", #since spp are unbalanced
  #strata       = meta_filtered$Site, #permutations constrained within sites
  permutations = 999
)

permanova_spp_site
#Species explains 15% of the multivariate photophysiological variation after accounting for Site.
#Site not signigicant

#######################
# Pairwise PERMANOVA ##
#######################

#fr BH round 2: using by 'margin' to get partial effects of each term rather than sequential sum of squares
permanova_species <- adonis2(
  permanova_matrix_filtered ~ Species + factor(Site),
  data = meta_filtered,
  method = "euclidean",
  by = "margin", #need this bc. spp unbalances across sites.measures what Species explains after accounting for Site bias.
  strata       = meta_filtered$Site,  # Controls for site dispersion heterogeneity
  #excluding strata results in site not being significant!
  permutations = 999
)

# Since site is nt sig, only need posthoc between spp
# Run pairwise comparison of site
# 1. Get all unique species pairs
sp_pairs <- combn(unique(as.character(meta_filtered$Species)), 2, simplify = FALSE)

# 2. Iterate over pairs matching global adonis2 structure
pw_species_margin <- map_dfr(sp_pairs, function(pair) {
  
  # Subset data for the current pair
  sub_idx  <- meta_filtered$Species %in% pair
  sub_meta <- meta_filtered[sub_idx, ]
  sub_mat  <- permanova_matrix_filtered[sub_idx, ]
  
  # Ensure clean factor drops
  sub_meta$Species <- factor(sub_meta$Species, levels = pair)
  sub_meta$Site    <- factor(sub_meta$Site)
  
  # Run identical marginal model with site strata
  res <- adonis2(
    sub_mat ~ Species + factor(Site),
    data         = sub_meta,
    method       = "euclidean",
    by           = "margin",
    strata       = sub_meta$Site,
    permutations = 999
  )
  
  # Extract Species term (row 1)
  data.frame(
    Comparison = paste(pair[1], "vs", pair[2]),
    Df         = res$Df[1],
    SumOfSqs   = round(res$SumOfSqs[1], 2),
    R2         = round(res$R2[1], 4),
    F_stat     = round(res$F[1], 2),
    p_val      = res$`Pr(>F)`[1]
  )
}) %>%
  filter(!is.na(p_val)) %>%
  mutate(
    p_adj  = round(p.adjust(p_val, method = "BH"), 4),
    Signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

cat("\n--- PAIRWISE SPECIES PERMANOVA (MARGINAL & STRATIFIED) ---\n")
pw_species_margin<-pw_species_margin%>% 
  arrange(desc(R2)) 


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
ef_filtered <- envfit(nmds_filtered, permanova_matrix_filtered, permutations = 999)
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
  left_join(meta_filteredK, by = "SampleID_clean")

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
# p <- ggplot() +
#   geom_point(
#     data = nmds_scores, 
#     aes(x = NMDS1, y = NMDS2, color = Species, shape = factor(Site)), 
#     size = 2.8, alpha = 0.75
#   ) + 
#   scale_shape_manual(values = custom_shapes, name = "Site") +
#   scale_color_manual(values = species_colors, name = "Species") +
#   scale_fill_manual(values = species_colors, name = "Species") +
#   new_scale_color() +
#   geom_segment(
#     data = top_vectors_clean, 
#     aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
#     arrow = arrow(length = unit(0.20, "cm")), linewidth = 0.85
#   ) +
#   geom_text(
#     data = top_vectors_clean,
#     aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Metric),
#     size = 3
#   ) +
#   scale_color_manual(values = metric_colors, name = "Trait Family") +
#   coord_cartesian(clip = "off") +
#   theme_bw() +
#   theme(plot.margin = unit(c(15, 25, 15, 25), "pt"))+
#   labs( x = "nMDS Dimension 1", y = "nMDS Dimension 2")
# p

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
    aes(x = NMDS1, y = NMDS2, color = Species), #shape = factor(Site)), 
    size = 2.8, alpha = 0.75
  ) + 
  stat_ellipse(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = Species, fill = Species), 
    geom = "polygon", alpha = 0.1, level = 0.90
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
  #scale_shape_manual(values = custom_shapes, name = "Site") +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_fill_manual(values = species_colors, name = "Species") +
  new_scale_color() +
  # geom_segment(
  #   data = top_vectors_clean, 
  #   aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
  #   arrow = arrow(length = unit(0.20, "cm")), linewidth = 0.85
  # ) +
  # geom_text(
  #   data = top_vectors_clean,
  #   aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Metric),
  #   size = 3
  # ) +
  scale_color_manual(values = metric_colors, name = "Trait Family") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(plot.margin = unit(c(15, 25, 15, 25), "pt"))+
  labs( x = "nMDS Dimension 1", y = "nMDS Dimension 2")

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
    title = paste0("Filtered ( nMDS Stress = ", round(nmds_filtered$stress, 3), ")"),
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
print(p)
