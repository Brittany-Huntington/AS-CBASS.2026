library(tidyverse)
library(vegan)
library(here)
library(ggplot2)
library(ggnewscale)

species_colors <- c(
  "aabr" = "#046E8F", 
  "aglo" = "lightgreen", 
  "ahya" = "#D44D5C",
  "icra" = "#462255"  
)
# Define explicit named vector for site shapes
custom_shapes = c(
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

# =========================================================================
# 1. Read & Format Unfiltered Data
# =========================================================================
pp           = read.csv(here("data", "permanova_matrix_clean.csv")) 
annot_master = read.csv(here("Data", "cluster_and_ed50_data.csv")) # phenotype_heatmap script
ed           = annot_master

pp_formatted = pp %>%
  mutate(
    sample_file = X,
    site        = str_extract(sample_file, "(?<=^t)\\d+"),
    species     = str_extract(sample_file, "(?<=_)[A-Za-z]+"),
    genotype    = str_extract(sample_file, "\\d+$")
  )

meta_cols  = c("sample_file", "site", "species", "genotype")
trait_cols = setdiff(names(pp_formatted), c(meta_cols, "X"))
pp_formatted = pp_formatted[, c(meta_cols, trait_cols)]

# Matching Key Function to join ED50/Species metadata
create_matching_key = function(df) {
  df %>%
    mutate(
      tp_num           = str_match(sample_file, "^t(\\d+)")[,2],
      geno_name        = toupper(str_match(sample_file, "^t\\d+_([a-zA-Z]+)")[,2]),
      rep_num          = str_match(sample_file, "^t\\d+_[a-zA-Z]+(\\d+)")[,2],
      GroupingProperty = paste(tp_num, geno_name, rep_num, sep = "_")
    ) %>%
    dplyr::select(-tp_num, -geno_name, -rep_num)
}

# Apply key creation to full unfiltered dataset
pp_unfiltered = create_matching_key(pp_formatted)

# Join Species and metadata
pp_unfiltered_species = pp_unfiltered %>%
  left_join(
    ed %>% dplyr::select(SampleID_clean), 
    by = c("sample_file" = "SampleID_clean")
  ) %>%
  filter(!is.na(species))

# =========================================================================
# 2. Extract & Scale Unfiltered Traits (Clean Zero-Variance Traits)
# =========================================================================
# 1. Extract raw numeric trait matrix
trait_matrix_raw = pp_unfiltered_species %>% dplyr::select(all_of(trait_cols))

# 2. Drop zero-variance or all-NA columns that break scaling/dist
valid_traits_mask = sapply(trait_matrix_raw, function(x) {
  v = var(x, na.rm = TRUE)
  !is.na(v) && v > 0
})

trait_matrix_clean = trait_matrix_raw[, valid_traits_mask]

cat("Unfiltered Total Traits:", length(trait_cols), "\n")
cat("Retained Non-Zero-Variance Traits:", ncol(trait_matrix_clean), "\n")

# 3. Z-score scale the matrix
matrix_unfiltered_scaled = scale(trait_matrix_clean)
matrix_unfiltered_scaled[is.na(matrix_unfiltered_scaled)] = 0

# Extract metadata corresponding to rows
meta_unfiltered = pp_unfiltered_species %>% 
  dplyr::select(sample_file, genotype, sample_file, species, site)

# =========================================================================
# 3. PERMANOVA & Dispersion Tests (Unfiltered)
# =========================================================================
# Overall PERMANOVA
permanova_unfiltered = adonis2(
  matrix_unfiltered_scaled ~ species, 
  data         = meta_unfiltered, 
  method       = "euclidean", 
  permutations = 999
)

cat("\n--- UNFILTERED PERMANOVA RESULTS ---\n")
print(permanova_unfiltered)

# Homogeneity of Multivariate Dispersion (BETADISPER)
dist_unfiltered = vegdist(matrix_unfiltered_scaled, method = "euclidean")
bd_unfiltered   = betadisper(dist_unfiltered, meta_unfiltered$species)
anova_bd        = anova(bd_unfiltered)

cat("\n--- UNFILTERED DISPERSION TEST p-value:", anova_bd[1, "Pr(>F)"], "---\n")

# Pairwise PERMANOVA Function
pairwise.adonis2 = function(x, data, strata = NULL, nperm = 999, ... ) {
  ststri = ifelse(is.null(strata), 'Null', strata)
  fostri = as.character(x)
  
  x1 = x
  lhs = eval(x1[[2]], environment(x1), globalenv())
  environment(x1) = environment()
  rhs = x1[[3]]
  
  x1[[2]] = NULL
  rhs.frame = model.frame(x1, data, drop.unused.levels = TRUE)
  co = combn(unique(as.character(rhs.frame[,1])), 2)
  
  nameres = c('parent_call')
  for (elem in 1:ncol(co)){
    nameres = c(nameres, paste(co[1,elem], co[2,elem], sep = '_vs_'))
  }
  
  res = vector(mode = "list", length = length(nameres))
  names(res) = nameres
  res['parent_call'] = list(paste(fostri[2], fostri[1], fostri[3], ', strata =', ststri, ', permutations', nperm))
  
  for(elem in 1:ncol(co)){
    if(inherits(eval(lhs), 'dist')){
      xred = as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem], co[2,elem]),
                                          rhs.frame[,1] %in% c(co[1,elem], co[2,elem])])
    } else {
      xred = eval(lhs)[rhs.frame[,1] %in% c(co[1,elem], co[2,elem]), ]
    }
    
    mdat1 = data[rhs.frame[,1] %in% c(co[1,elem], co[2,elem]), ]
    
    if(length(rhs) == 1){
      xnew = as.formula(paste('xred', as.character(rhs), sep = '~'))
    } else {
      xnew = as.formula(paste('xred', paste(rhs[-1], collapse = as.character(rhs[1])), sep = '~'))
    }
    
    if(is.null(strata)){
      ad = adonis2(xnew, data = mdat1, ... )
    } else {
      perm = how(nperm = nperm)
      setBlocks(perm) = with(mdat1, mdat1[, ststri])
      ad = adonis2(xnew, data = mdat1, permutations = perm, ... )
    }
    
    res[nameres[elem+1]] = list(ad[1:5])
  }
  class(res) = c("pwadstrata", "list")
  return(res)
}

pw_unfiltered = pairwise.adonis2(matrix_unfiltered_scaled ~ species, data = meta_unfiltered, method = "euclidean")
cat("\n--- UNFILTERED PAIRWISE PERMANOVA ---\n")
print(pw_unfiltered)

# =========================================================================
# 4. nMDS & Trait Vector Fitting (Unfiltered)
# =========================================================================
set.seed(123)
nmds_unfiltered = metaMDS(
  matrix_unfiltered_scaled, 
  distance      = "euclidean", 
  k             = 2, 
  trymax        = 200, 
  maxit         = 500, 
  autotransform = FALSE,
  trace         = 0
)

cat("\nUnfiltered nMDS Stress Value:", nmds_unfiltered$stress, "\n")

# Fit vectors for all traits
ef_unfiltered = envfit(nmds_unfiltered, matrix_unfiltered_scaled, permutations = 999)

saveRDS(ef_unfiltered, here("Outputs", "ef_unfiltered_nmds.rds"))

vector_scores = as.data.frame(scores(ef_unfiltered, display = "vectors"))
vector_scores$Metric = rownames(vector_scores)
vector_scores$r2     = ef_unfiltered$vectors$r
vector_scores$p_val  = ef_unfiltered$vectors$pvals

# Filter top 10 metric drivers (these are collinear!)
top_vectors = vector_scores %>%
  filter(p_val < 0.01) %>%
  arrange(desc(r2)) %>%
  head(15)

cat("\n--- TOP METRIC DRIVERS OF UNFILTERED ORDINATION ---\n")
print(top_vectors)
# Keep all statistically significant candidate traits
sig_vectors = vector_scores %>%
  filter(p_val < 0.05) %>%
  arrange(desc(r2))

# =========================================================================
# 2. Compute Correlation Matrix & Variable Clustering (|r| >= 0.90threshold)
# =========================================================================
# Extract scaled raw trait data for candidate vectors
sig_trait_matrix = matrix_unfiltered_scaled[, sig_vectors$Metric]

# Pearson correlation matrix
cor_mat = cor(sig_trait_matrix, use = "pairwise.complete.obs")

# Convert correlation to distance matrix (1 - |r|)
trait_dist = as.dist(1 - abs(cor_mat))

# Hierarchical clustering of variables
trait_tree = hclust(trait_dist, method = "complete")

# Cut tree at distance = 0.15 (equivalent to |r| = 0.90 collinearity threshold)
collinear_cluster_id = cutree(trait_tree, h = 0.10)

# =========================================================================
# 3. Build Reference Dataframe of Collinear Clusters
# =========================================================================
trait_cluster_ref = data.frame(
  Metric       = names(collinear_cluster_id),
  Cluster_ID   = factor(collinear_cluster_id)
) %>%
  left_join(sig_vectors, by = "Metric") %>%
  group_by(Cluster_ID) %>%
  # Identify the top representative trait per cluster (highest envfit R2)
  mutate(
    Is_Cluster_Exemplar = ifelse(r2 == max(r2), TRUE, FALSE)
  ) %>%
  ungroup() %>%
  arrange(Cluster_ID, desc(r2))

# --- VIEW REFERENCE DATAFRAME ---
# This dataframe lists every trait, its collinear cluster, and whether it's the exemplar
cat("--- TRAIT COLLINEARITY CLUSTER REFERENCE TABLE ---\n")
print(head(trait_cluster_ref, 20))

# =========================================================================
# 4. Extract Top 10 NON-COLLINEAR Drivers
# =========================================================================
top_10_non_collinear = trait_cluster_ref %>%
  filter(Is_Cluster_Exemplar == TRUE) %>%
  distinct(Cluster_ID, .keep_all = TRUE) %>% # Break any rare R2 ties
  arrange(desc(r2)) %>%
  head(10)

#pull associated drivers for the top 10
top_10_cluster_ids = trait_cluster_ref %>%
  filter(Is_Cluster_Exemplar == TRUE) %>%
  distinct(Cluster_ID, .keep_all = TRUE) %>%
  arrange(desc(r2)) %>%
  head(10) %>%
  pull(Cluster_ID)


top_10_with_all_collinears = trait_cluster_ref %>%
  # Filter for only clusters belonging to the top 10
  filter(Cluster_ID %in% top_10_cluster_ids) %>%
  # Arrange so the cluster exemplar is always listed first, followed by collinear traits by R2
  arrange(match(Cluster_ID, top_10_cluster_ids), desc(Is_Cluster_Exemplar), desc(r2)) %>%
  dplyr::select(
    Cluster_ID, 
    Is_Cluster_Exemplar, 
    Metric, 
    r2, 
    p_val
  )

write.csv(
  top_10_with_all_collinears, 
  here("Outputs", "top_10_non_collinear_clusters_with_all_traits.csv"), 
  row.names = FALSE
)

top_10_summary_wide = top_10_with_all_collinears %>%
  group_by(Cluster_ID) %>%
  summarise(
    Exemplar_Trait      = Metric[Is_Cluster_Exemplar == TRUE][1],
    Exemplar_R2         = r2[Is_Cluster_Exemplar == TRUE][1],
    Total_Traits_In_Group = n(),
    Collinear_Traits    = paste(Metric[Is_Cluster_Exemplar == FALSE], collapse = ", ")
  ) %>%
  arrange(desc(Exemplar_R2))

print(top_10_summary_wide)

#visualize correlations per group
library(pheatmap)
library(viridis)

#LOOP THROUGH and make correlation matrices for each trait
# 1. Create output directory if it doesn't exist
output_dir <- here("Plots", "cor")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Get unique Cluster IDs from the top 10 non-collinear dataset
unique_clusters <- unique(top_10_with_all_collinears$Cluster_ID)

# 3. Loop through each cluster and generate/save the correlation heatmap
for (cid in unique_clusters) {
  
  # Select traits belonging to the current cluster
  cluster_traits <- top_10_with_all_collinears %>%
    filter(Cluster_ID == cid) %>%
    pull(Metric)
  
  # Skip clusters with fewer than 2 traits (can't build a matrix/heatmap for 1 trait)
  if (length(cluster_traits) < 2) {
    cat(sprintf("Skipping Cluster %s: contains only 1 trait.\n", cid))
    next
  }
  
  # Extract correlation sub-matrix
  sub_cor_mat <- cor_mat[cluster_traits, cluster_traits, drop = FALSE]
  
  # Clean row and column names for clean display
  # Fallback to full metric name if regex pattern doesn't match
  short_names <- str_extract(cluster_traits, "^[^.]+\\.[^.]+")
  short_names[is.na(short_names)] <- cluster_traits[is.na(short_names)]
  
  rownames(sub_cor_mat) <- short_names
  colnames(sub_cor_mat) <- short_names
  
  # Define output file path
  file_path <- file.path(output_dir, sprintf("Cluster_%s_correlation_heatmap.png", cid))
  
  # Dynamic dimension sizing based on the number of traits
  n_traits <- length(cluster_traits)
  img_width <- max(6, n_traits * 0.8)
  img_height <- max(5, n_traits * 0.8)
  
  # Generate and save heatmap via pheatmap filename argument
  pheatmap(
    sub_cor_mat,
    color            = viridis(100),
    display_numbers  = TRUE,
    number_format    = "%.2f",
    fontsize_number  = max(6, 10 - n_traits * 0.2), # Adjust font size dynamically
    main             = sprintf("Cluster %s Trait Correlation Matrix (r)", cid),
    filename         = file_path,
    width            = img_width,
    height           = img_height
  )
  
  cat(sprintf("Saved heatmap for Cluster %s -> %s\n", cid, file_path))
}

################

top_vectors_clean = top_10_non_collinear %>%
  mutate(
    # Match raw string to standard metric palette categories
    Metric_Category = case_when(
      grepl("Quant", Metric, ignore.case = TRUE) ~ "Quant (Quantum Yield)",
      grepl("qP",    Metric, ignore.case = TRUE) ~ "qP (Photochemical Quenching)",
      grepl("qm",    Metric, ignore.case = TRUE) ~ "qm (Max Quenching)",
      grepl("ABQ",   Metric, ignore.case = TRUE) ~ "ABQ (Absorption)",
      grepl("Sigma", Metric, ignore.case = TRUE) ~ "Sigma (Antenna Size)",
      grepl("Connect", Metric, ignore.case = TRUE) ~ "Connect (Connectivity)",
      grepl("Tau1",  Metric, ignore.case = TRUE) ~ "Tau1 (Transport Kinetics 1)",
      grepl("Tau2",  Metric, ignore.case = TRUE) ~ "Tau2 (Transport Kinetics 2)",
      grepl("NPQ",   Metric, ignore.case = TRUE) ~ "NPQ (Non-Photochemical Quenching)",
      TRUE                                       ~ "Other Metric"
    ),
    
    # Extract short trait family
    Trait_Family = str_extract(Metric, "(mQuant|qqP|mqP|qP|rqm|qm|rABQ|Sigma|Tau1|Tau2|NPQ)"),
    
    # Extract condition/lighting phase prefix (e.g. L3CampR1, DCompL12, DCampL1, mean)
    Phase = str_extract(Metric, "^[A-Za-z0-9]+"),
    
    # Construct clean publication label (e.g. "mQuant (L3CampR1)")
    Clean_Label = case_when(
      !is.na(Trait_Family) & !is.na(Phase) ~ paste0(Trait_Family, " (", Phase, ")"),
      !is.na(Trait_Family)                 ~ Trait_Family,
      TRUE                                 ~ Metric
    ),
    
    # Scale vector arrows for nMDS canvas
    NMDS1_scaled = NMDS1 * sqrt(r2) * arrow_mult,
    NMDS2_scaled = NMDS2 * sqrt(r2) * arrow_mult
  )

# Verify clean label extraction
print(top_vectors_clean %>% dplyr::select(Metric, Metric_Category, Clean_Label, Trait_Family))

# Site scores for plotting
# nmds_scores = as.data.frame(scores(nmds_unfiltered, display = "sites")) %>%
#   bind_cols(meta_unfiltered %>% dplyr::select(species, genotype, site))

nmds_scores = as.data.frame(scores(nmds_unfiltered, display = "sites")) %>%
  bind_cols(meta_unfiltered %>% dplyr::select(site, sample_file, species, genotype))

# join K_3 from `ed` using GroupingProperty (or sample_file)
nmds_scores = nmds_scores %>%
  left_join(
    ed %>% dplyr::select(SampleID_clean, K_3), 
    by = c("sample_file" = "SampleID_clean")
  ) %>%
  mutate(K_3 = factor(K_3))

# Define clean color palette for K=3 Clusters
cluster_colors = c(
  "1" = "orange",  
  "2" = "purple",  
  "3" = "grey"   
)

# Auto-calculate arrow scaling multiplier to fit canvas perfectly
max_site_coord   = max(abs(c(nmds_scores$NMDS1, nmds_scores$NMDS2)))
max_vector_coord = max(sqrt(top_vectors_clean$NMDS1^2 + top_vectors_clean$NMDS2^2))
arrow_mult       = (max_site_coord * 0.8) / max_vector_coord

top_vectors_clean = top_vectors_clean %>%
  mutate(
    NMDS1_scaled = NMDS1 * sqrt(r2) * arrow_mult,
    NMDS2_scaled = NMDS2 * sqrt(r2) * arrow_mult
  )

# =========================================================================
# 5. Plot Unfiltered nMDS
# =========================================================================
spp<-ggplot() +
geom_point(
  data = nmds_scores, 
  aes(x = NMDS1, y = NMDS2, color = species, shape = site), 
  size = 2.8, 
  alpha = 0.75
) +
  stat_ellipse(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = species, fill = species), 
    geom = "polygon", 
    alpha = 0.12, 
    level = 0.50
  ) +
  scale_shape_manual(values = custom_shapes, name = "Site") +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_fill_manual(values = species_colors, name = "Species") +
  new_scale_color() +
geom_segment(
  data = top_10_vectors_clean, 
  aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Trait_Family),
  arrow = arrow(length = unit(0.20, "cm")), 
  linewidth = 0.85
) +
  geom_text(
    data = top_10_vectors_clean, 
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Metric, color = Trait_Family),
    fontface = "bold", 
    size = 3.5
  ) +
  scale_color_brewer(palette = "Set1", name = "Trait Family", na.value = "grey70") +
theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box      = "vertical"
  ) +
  labs(
    title = paste0("Unfiltered (nMDS, Stress = ", round(nmds_unfiltered$stress, 3), ")"),
    x = "nMDS Dimension 1",
    y = "nMDS Dimension 2"
  )
ggsave(
  filename = here("Plots", "nMDS_spp_unfiltered.png"),
  plot     = spp,
  width    = 9.5,
  height   = 8.0,
  dpi      = 300,
  units    = "in"
)

#####plot by cluster, site

# Ensure site is a factor for discrete shape mapping
nmds_scores$site = as.factor(nmds_scores$site)


p_nmds<- ggplot() +
  # --- Layer 1: Species points & Site shapes ---
  geom_point(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = species, shape = site), 
    size = 3, 
    alpha = 0.85
  ) +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_shape_manual(values = custom_shapes, name = "Site") +
  
  new_scale_color() +
  new_scale_fill() +
  
  # --- Layer 2: K_3 Cluster Ellipses ---
  stat_ellipse(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = K_3, fill = K_3), 
    geom = "polygon", 
    alpha = 0.15, 
    level = 0.95
  ) +
  scale_color_manual(values = cluster_colors, name = "Hclust Group (K_3)") +
  scale_fill_manual(values = cluster_colors, name = "Hclust Group (K_3)") +
  
  new_scale_color() +
  
  # --- Layer 3: Vectors & Labels COLORED BY TRAIT_FAMILY ---
  geom_segment(
    data = top_10_vectors_clean, 
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Trait_Family),
    arrow = arrow(length = unit(0.20, "cm")), 
    linewidth = 0.85
  ) +
  geom_text(
    data = top_10_vectors_clean, 
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Metric, color = Trait_Family),
    fontface = "bold", 
    size = 3.5
  ) +
  # Option A: Distinct qualitative palette automatically assigned per family
  scale_color_brewer(palette = "Set1", name = "Trait Family",  na.value = "grey70") +
  
  # Option B: (Alternative) Viridis palette for high-contrast accessibility
  # scale_color_viridis_d(option = "Dark1", name = "Trait Family") +
  
  # --- Formatting ---
  theme_classic(base_size = 14) +
  theme(legend.position = "right") +
  labs(
    x = "nMDS Dimension 1",
    y = "nMDS Dimension 2"
  )
ggsave(
  filename = here("Plots", "nMDS_K3_unfiltered.png"),
  plot     = p_nmds,
  width    = 9.5,
  height   = 8.0,
  dpi      = 300,
  units    = "in"
)


#not taking out collinear!
collinear_metric_plot <- ggplot() +
  # --- 1. Site Points & Ellipses (Use unscaled NMDS1 & NMDS2) ---
  geom_point(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = species, shape = site), 
    size = 2.8, 
    alpha = 0.75
  ) +
  stat_ellipse(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = species, fill = species), 
    geom = "polygon", 
    alpha = 0.12, 
    level = 0.50
  ) +
  scale_shape_manual(values = custom_shapes, name = "Site") +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_fill_manual(values = species_colors, name = "Species") +
  
  # --- Reset color scale for vectors ---
  new_scale_color() +
  
  # --- 2. Filtered Non-Collinear Vectors (Use top_10_vectors_clean & Scaled coords) ---
  geom_segment(
    data = top_vectors, 
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Trait_Family),
    arrow = arrow(length = unit(0.20, "cm")), 
    linewidth = 0.85
  ) +
  geom_text(
    data = top_vectors, 
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Clean_Label, color = Trait_Family),
    fontface = "bold", 
    size = 3.5
  ) +
  scale_color_brewer(palette = "Set1", name = "Trait Family", na.value = "grey70") +
  
  # --- 3. Formatting ---
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box      = "vertical"
  ) +
  labs(
    title = paste0("Unfiltered (nMDS, Stress = ", round(nmds_unfiltered$stress, 3), ")"),
    x = "nMDS Dimension 1",
    y = "nMDS Dimension 2"
  )

# Save
ggsave(
  filename = here("Plots", "nMDS_spp_unfiltered_COR.png"),
  plot     = collinear_metric_plot,
  width    = 9.5,
  height   = 8.0,
  dpi      = 300,
  units    = "in"
)
