library(tidyverse)
library(vegan)
library(here)

species_colors <- c(
  "aabr" = "#046E8F", 
  "aglo" = "lightgreen", 
  "ahya" = "#D44D5C",
  "icra" = "#462255"  
)
custom_shapes <- c(23, 22, 15, 16, 17, 18, 21, 24, 1)

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

# Filter top 10 metric drivers
top_vectors = vector_scores %>%
  filter(p_val < 0.01) %>%
  arrange(desc(r2)) %>%
  head(15)

cat("\n--- TOP METRIC DRIVERS OF UNFILTERED ORDINATION ---\n")
print(top_vectors)

top_vectors_clean = top_vectors %>%
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
    
    # Extract short trait family (mQuant, qqP, mqP, qm, rABQ)
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
print(top_vectors_clean %>% dplyr::select(Metric, Metric_Category, Clean_Label))

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
  "1" = "yellow",  
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
ggplot() +
  geom_point(data = nmds_scores, aes(x = NMDS1, y = NMDS2, color = species), 
             size = 2.8, alpha = 0.75) +
  stat_ellipse(data = nmds_scores, aes(x = NMDS1, y = NMDS2, color = species, fill = species), 
               geom = "polygon", alpha = 0.12, level = 0.50) +
  geom_segment(data = top_vectors_clean, 
               aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled),
               arrow = arrow(length = unit(0.22, "cm")), 
               color = "grey10", linewidth = 0.8) +
  geom_text(data = top_vectors_clean, 
            aes(x = NMDS1_scaled * 1.12, y = NMDS2_scaled * 1.12, label = Clean_Label),
            color = "black", fontface = "bold", size = 3.5) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right") +
  labs(
    title = paste0("Unfiltered Photophysiological Ordination (nMDS, Stress = ", round(nmds_unfiltered$stress, 3), ")"),
    x = "nMDS Dimension 1",
    y = "nMDS Dimension 2"
  )

#plot by cluster, site
ggplot() +
  # --- Points: Shape by Species, Color by K_3 Cluster from `ed` ---
  geom_point(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = species, shape = site), 
    size = 3, 
    alpha = 0.85
  ) +
  
  # --- Enclose K_3 Clusters (Ellipses) ---
  stat_ellipse(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = K_3, fill = K_3), 
    geom = "polygon", 
    alpha = 0.15, 
    level = 0.95
  ) +
  
  # Set Cluster Fill & Color Scales
  scale_color_manual(values = cluster_colors, name = "Hclust Group (K_3)") +
  scale_fill_manual(values = cluster_colors, name = "Hclust Group (K_3)") +
  
  # --- Switch Color Scale for Trait Vector Arrows ---
  new_scale_color() +
  
  # --- Trait Vectors & Labels Matched to Metric_Colors ---
  geom_segment(
    data = top_vectors_clean, 
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled), #, color = Metric_Category),
    arrow = arrow(length = unit(0.20, "cm")), 
    linewidth = 0.8
  ) +
  geom_text(
    data = top_vectors_clean, 
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Clean_Label),# color = species_colors),
    fontface = "bold", 
    size = 3.5
  ) +
  scale_color_manual(values = metric_colors, name = "Photophys Metric") +
  
  # --- Theme Formatting ---
  theme_classic(base_size = 14) +
  theme(legend.position = "right") +
  labs(
    #title = paste0("nMDS with K_3 Clusters & Top Trait Vectors (Stress = ", round(nmds_unfiltered$stress, 3), ")"),
    x = "nMDS Dimension 1",
    y = "nMDS Dimension 2"
  )+
  scale_shape_manual(values = custom_shapes) 

