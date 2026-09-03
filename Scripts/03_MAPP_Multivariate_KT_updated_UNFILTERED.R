rm(list = ls())
library(tidyverse)
library(vegan)
library(here)
library(ggplot2)
library(ggnewscale)
library(pheatmap)
library(viridis)

# species_colors <- c(
#   "aabr" = "#046E8F", 
#   "aglo" = "lightgreen", 
#   "ahya" = "#D44D5C",
#   "icra" = "#462255"  
# )
species_colors <- c(
  "AABR" = "#046E8F", 
  "AGLO" = "lightgreen", 
  "AHYA" = "#D44D5C",
  "ICRA" = "#462255"  
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
category_colors <- c(
  "ABQ (Absorption)"                       = "#1F77B4", 
  "Quant (Quantum Yield)"                  = "#2CA02C", 
  "Sigma (Antenna Size)"                   = "#FF7F0E", 
  "qm (Max / Non-Photochemical Quenching)"  = "#D62728", 
  "qP (Photochemical Quenching)"           = "#9467BD", 
  "Tau1 (Transport Kinetics 1)"            = "#8C564B", 
  "Tau2 (Transport Kinetics 2)"            = "#E377C2", 
  "Connect (Connectivity)"                 = "#17BECF", 
  "Other Metric"                           = "#7F7F7F"  
)

pp = read.csv(here("Outputs", "permanova_matrix_clean1.csv")) 
ed = read.csv(here("Outputs", "ed50_colorslope_cluster_metadata.csv")) #from data compiler
load(here("Outputs", "photophys_preprocessed_data.RData"))


#run hclust
cbb<-pvclust(t(Alls.cor_clean), method.dist="canberra", parallel=TRUE, method.hclust="ward.D", nboot=1000)
save(cbb, file = here("Data", "cbb_AS26_expt_1000bs_updatedw9icra3.RData"))
#or if you dont need to rerun
#load(file = here("Data", "cbb_AS26_expt_1000bs_updated.RData")) #loads cbb
#load(file = here("Data", "cbb_AS26_expt_1000bs_updatedw9icra3.RData"))


load(file = here("Data", "cbb_AS26_expt_1000bs_updatedw9icra3.RData"))
#Cut the hclust tree for K = 2 through 15
cluster_levels <- lapply(2:15, function(k) {
  groups <- cutree(cbb$hclust, k = k)
  data.frame(
    SampleID = names(groups),
    Cluster = as.integer(groups),
    k_val = k
  )
})

#Pivot the cluster assignments wide  to match cols
cluster_matrix_wide <- bind_rows(cluster_levels) %>%
  pivot_wider(
    names_from = k_val,
    values_from = Cluster,
    names_prefix = "K_"
  )

#Clean the SampleID names and parse site, species, and genotype
cluster_summary <- cluster_matrix_wide %>%
  mutate(
    SampleID_clean = SampleID %>%
      toupper() %>%
      trimws() %>%
      str_replace("[-_\\.]?\\d{2,3}R?[-_\\.]?PROCESSED\\.CSV$", "") %>%
      str_replace("[-_\\.]?\\d{2,3}R$", "") %>%
      str_replace("^T(?=\\d)", "") %>%
      str_replace("([A-Z]{3,4})(\\d+)$", "\\1_\\2"),
    # site     = str_extract(SampleID_clean, "^\\d+"),                      
    # species  = str_extract(SampleID_clean, "(?<=_)[A-Z]{3,4}(?=_)"),      
    # genotype = str_extract(SampleID_clean, "\\d+$")                        
  ) %>%
  distinct(SampleID_clean, .keep_all = TRUE)
# 
# #join ed50 and color slopes to cluster summary and save
# cluster_summary_merged <- cluster_summary %>%
#   left_join(
#     color,
#       #dplyr::select(-c(Image, Folder, Species, SITE)), 
#     #by = "SampleID_clean")
#     by = "SampleID" = "SampleID_clean")
#   
# 
# cluster_summary_merged <- cluster_summary_merged %>%
#   left_join(
#     cbass_metadata, 
#     by = c("SampleID_clean")
#   )
# 
# write.csv(cluster_summary_merged, file = here("Outputs", "fvfm_color_cluster_metadata.csv"), row.names = FALSE)

cluster_summary_ed50_slope<-cluster_summary %>%
  left_join(
    ed50,
    # dplyr::select(-c(site, species, Site, Species)), 
    by = "SampleID_clean")

cluster_summary_ed50_slope<-cluster_summary_ed50_slope%>%
  left_join(
    color%>%
      dplyr::select(SampleID, int, slope, r2, p), 
    by = c("SampleID_clean" = "SampleID")
  )

write.csv(cluster_summary_ed50_slope, file = here("Outputs", "ed50_colorslope_cluster_metadata.csv"), row.names = FALSE)
###############
ed_cleaned = ed %>%
  mutate(
    SampleID_clean = str_trim(as.character(SampleID_clean))
  )

pp_formatted = pp %>%
  mutate(
    SampleID_clean = as.character(X),
    site        = str_extract(SampleID_clean, "^\\d+"),
    species     = str_extract(SampleID_clean, "(?<=_)[A-Za-z]+"),
    genotype    = str_extract(SampleID_clean, "\\d+$")
  )

meta_cols  = c("SampleID_clean", "site", "species", "genotype")
trait_cols = setdiff(names(pp_formatted), c(meta_cols, "X"))
pp_formatted = pp_formatted[, c(meta_cols, trait_cols)]


# Matching Key Function to join ED50/Species metadata
create_matching_key = function(df) {
  df %>%
    mutate(
      clean_file = str_remove(SampleID_clean, "\\.[a-zA-Z0-9]+$"),
      tp_num     = str_match(clean_file, "^t(\\d+)")[,2],
      geno_name  = toupper(str_match(clean_file, "^t\\d+_([a-zA-Z]+)")[,2]),
      rep_num    = str_match(clean_file, "^t\\d+_[a-zA-Z]+(\\d+)")[,2],
      GroupingProperty = paste(tp_num, geno_name, rep_num, sep = "_")
    ) %>%
    dplyr::select(-clean_file, -tp_num, -geno_name, -rep_num)
}

# Apply key creation to full unfiltered dataset
pp_unfiltered = create_matching_key(pp_formatted)

# Join Species and metadata
pp_unfiltered_species <- pp_unfiltered %>%
  mutate(SampleID_clean = trimws(SampleID_clean)) %>%
  left_join(
    ed_cleaned %>%
      mutate(SampleID_clean = trimws(SampleID_clean)) %>%
      dplyr::select(SampleID_clean, Species),
    by =  "SampleID_clean"
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
  dplyr::select(all_of(meta_cols), SampleID_clean) %>%
  left_join(
    ed_cleaned %>% dplyr::select(SampleID_clean, K_3, ED50, slope), 
    by = "SampleID_clean"
  ) %>%
  mutate(K_3 = factor(K_3))

# =========================================================================
# 3. PERMANOVA & Dispersion Tests 
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
#ef_unfiltered = envfit(nmds_unfiltered, matrix_unfiltered_scaled, permutations = 999)

#saveRDS(ef_unfiltered, here("Outputs", "ef_unfiltered_nmds.rds"))
#ReadRDS(here("Outputs/ef_unfiltered_nmds.rds"))

vector_scores = as.data.frame(scores(ef_unfiltered, display = "vectors"))
vector_scores$Metric = rownames(vector_scores)
vector_scores$r2     = ef_unfiltered$vectors$r
vector_scores$p_val  = ef_unfiltered$vectors$pvals

vector_scores = vector_scores %>%
  left_join(metric_families, by = "Metric")

# Filter top 10 metric drivers (these are collinear!)
top_vectors = vector_scores %>%
  filter(p_val < 0.01) %>%
  arrange(desc(r2)) %>%
  head(25)

cat("\n--- TOP METRIC DRIVERS OF UNFILTERED ORDINATION ---\n")
print(top_vectors)
# Keep all statistically significant candidate traits
sig_vectors = vector_scores %>%
  filter(p_val < 0.0019) %>%
  arrange(desc(r2))

nmds_scores = as.data.frame(scores(nmds_unfiltered, display = "sites")) %>%
  bind_cols(meta_unfiltered %>% dplyr::select(site, SampleID_clean, species, genotype))

# join K_3 from `ed` using GroupingProperty (or sample_file)
nmds_scores = nmds_scores %>%
  left_join(
    ed %>% dplyr::select(SampleID_clean, K_3, Site), 
    by =  "SampleID_clean"
  ) %>%
  mutate(K_3 = factor(K_3))

# Define clean color palette for K=3 Clusters
cluster_colors = c(
  "1" = "orange",  
  "2" = "violet",  
  "3" = "grey"   
)
max_site_coord   = max(abs(c(nmds_scores$NMDS1, nmds_scores$NMDS2)))
max_vector_coord = max(sqrt(top_vectors$NMDS1^2 + top_vectors$NMDS2^2))
arrow_mult       = (max_site_coord * 0.8) / max_vector_coord


top_vectors_clean = top_vectors %>%
  mutate(
    # Master Family Category Mapping
    Metric_Category = case_when(
      # Quenching - Max / Non-Photochemical (qm, rqm)
      grepl("r?qm|npq", Metric, ignore.case = TRUE) ~ "qm (Max / Non-Photochemical Quenching)",
      
      # Quenching - Photochemical (qp, mpq, mqP, ppq, qqP)
      grepl("m?qp|p?qp|qqp|mpq|ppq", Metric, ignore.case = TRUE) ~ "qP (Photochemical Quenching)",
      
      # Absorption (ABQ, rABQ)
      grepl("r?abq", Metric, ignore.case = TRUE) ~ "ABQ (Absorption)",
      
      # Quantum Yield
      grepl("quant|fvfm|yield|fv_fm", Metric, ignore.case = TRUE) ~ "Quant (Quantum Yield)",
      
      # Kinetics & Antenna
      grepl("sigma|sig", Metric, ignore.case = TRUE) ~ "Sigma (Antenna Size)",
      grepl("tau1|tau_1|\\bt1\\b", Metric, ignore.case = TRUE) ~ "Tau1 (Transport Kinetics 1)",
      grepl("tau2|tau_2|\\bt2\\b", Metric, ignore.case = TRUE) ~ "Tau2 (Transport Kinetics 2)",
      grepl("connect|conn|ncon", Metric, ignore.case = TRUE) ~ "Connect (Connectivity)",
      
      TRUE ~ "Other Metric"
    ),
    
    # Extract clean 3-letter/4-letter short tag for in-plot labels
    Trait_Family = str_extract(Metric, "(?i)(mQuant|qqP|mqP|ppq|mpq|qP|rqm|qm|rABQ|ABQ|Sigma|Tau1|Tau2|NPQ|nCon)"),
    
    # Extract condition/lighting phase prefix (e.g. L3CampR1)
    Phase = str_extract(Metric, "^[A-Za-z0-9]+"),
    
    # Construct clean publication label
    Clean_Label = case_when(
      !is.na(Trait_Family) & !is.na(Phase) ~ paste0(Trait_Family, " (", Phase, ")"),
      !is.na(Trait_Family)                 ~ Trait_Family,
      TRUE                                 ~ Metric
    ),
    
    # Apply arrow_mult to scale vector endpoints for NMDS display
    NMDS1_scaled = NMDS1 * arrow_mult,
    NMDS2_scaled = NMDS2 * arrow_mult
  )

# Verify clean label extraction
print(top_vectors_clean %>% dplyr::select(Metric, Metric_Category, Clean_Label, Trait_Family))

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
  data = top_vectors_clean, 
  aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
  arrow = arrow(length = unit(0.20, "cm")), 
  linewidth = 0.85
) +
  geom_text(
    data = top_vectors_clean,
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Trait_Family),
    #fontface = "bold",
    size = 3
  ) +
  scale_color_manual(values = category_colors, name = "Trait Family") +
 # scale_color_viridis_d(option = "turbo", name = "Trait Family") +
  
  # --- 3. BOUNDS & SPACING ---
  coord_cartesian(clip = "off") +
  expand_limits(
    x = c(min(nmds_scores$NMDS1) * 1.25, max(nmds_scores$NMDS1) * 1.25),
    y = c(min(nmds_scores$NMDS2) * 1.25, max(nmds_scores$NMDS2) * 1.25)
  ) +
  theme_bw() +
  theme(
    plot.margin = margin(15, 25, 15, 25, "pt"),
    legend.position = "right"
  )+
#   scale_color_brewer(palette = "Set1", name = "Trait Family", na.value = "grey70") +
# theme_classic(base_size = 14) +
#   theme(
#     legend.position = "right",
#     legend.box      = "vertical"
#   ) +
coord_cartesian(clip = "off") +
  expand_limits(
    x = c(min(nmds_scores$NMDS1) * 1.25, max(nmds_scores$NMDS1) * 1.25),
    y = c(min(nmds_scores$NMDS2) * 1.25, max(nmds_scores$NMDS2) * 1.25)
  ) +
  theme_bw() +
  theme(
    plot.margin = margin(15, 25, 15, 25, "pt"),
    legend.position = "right"
  ) +
  labs(
    title = paste0("Unfiltered (nMDS, Stress = ", round(nmds_unfiltered$stress, 3), ")"),
    x = "nMDS Dimension 1",
    y = "nMDS Dimension 2"
  )
ggsave(
  filename = here("Plots", "nMDS_spp_unfiltered_ez_top20.png"),
  plot     = spp,
  width    = 9.5,
  height   = 8.0,
  dpi      = 300,
  units    = "in"
)

# # Fit ED50 directly onto the full NMDS space
# ef_ed50 <- envfit(nmds_unfiltered, meta_unfiltered$ED50, permutations = 999, na.rm = TRUE)
# 
# # View direction and R2 of ED50
# print(ef_ed50)

#####plot by cluster, site

# Ensure site is a factor for discrete shape mapping
nmds_scores$site = as.factor(nmds_scores$site)

p_nmds<-ggplot() +
  geom_point(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = species, shape = site), 
    size = 2.8, 
    alpha = 0.75
  ) +
 
  scale_shape_manual(values = custom_shapes, name = "Site") +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_fill_manual(values = species_colors, name = "Species") +
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
  geom_segment(
    data = top_vectors_clean, 
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
    arrow = arrow(length = unit(0.20, "cm")), 
    linewidth = 0.85
  ) +
  geom_text(
    data = top_vectors_clean,
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Trait_Family),
    #fontface = "bold",
    size = 3
  ) +
  scale_color_manual(values = category_colors, name = "Trait Family") +
  # scale_color_viridis_d(option = "turbo", name = "Trait Family") +
  
  # --- 3. BOUNDS & SPACING ---
  coord_cartesian(clip = "off") +
  expand_limits(
    x = c(min(nmds_scores$NMDS1) * 1.25, max(nmds_scores$NMDS1) * 1.25),
    y = c(min(nmds_scores$NMDS2) * 1.25, max(nmds_scores$NMDS2) * 1.25)
  ) +
  theme_bw() +
  theme(
    plot.margin = margin(15, 25, 15, 25, "pt"),
    legend.position = "right"
  )+
  #   scale_color_brewer(palette = "Set1", name = "Trait Family", na.value = "grey70") +
  # theme_classic(base_size = 14) +
  #   theme(
  #     legend.position = "right",
  #     legend.box      = "vertical"
  #   ) +
  coord_cartesian(clip = "off") +
  expand_limits(
    x = c(min(nmds_scores$NMDS1) * 1.25, max(nmds_scores$NMDS1) * 1.25),
    y = c(min(nmds_scores$NMDS2) * 1.25, max(nmds_scores$NMDS2) * 1.25)
  ) +
  theme_bw() +
  theme(
    plot.margin = margin(15, 25, 15, 25, "pt"),
    legend.position = "right"
  ) +
  labs(
    title = paste0("Unfiltered (nMDS, Stress = ", round(nmds_unfiltered$stress, 3), ")"),
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

###################################
#now taking out collinear metrics!#
###################################
library(vegan)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(here)



# =========================================================================
# 2. Compute Correlation Matrix & Variable Clustering (|r| >= 0.99threshold)
# =========================================================================
# Extract scaled raw trait data for candidate vectors
sig_trait_matrix = matrix_unfiltered_scaled[, sig_vectors$Metric]

# Pearson correlation matrix
cor_mat = cor(sig_trait_matrix, use = "pairwise.complete.obs")

# Convert correlation to distance matrix (1 - |r|)
trait_dist = as.dist(1 - abs(cor_mat))

# Hierarchical clustering of variables
trait_tree = hclust(trait_dist, method = "complete")

# Cut tree at distance = 0.01 (equivalent to |r| = 0.99 collinearity threshold)
collinear_cluster_id = cutree(trait_tree, h = 0.01)

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
top_10_non_collinear <- trait_cluster_ref %>%
  filter(Is_Cluster_Exemplar == TRUE) %>%
  distinct(Cluster_ID, .keep_all = TRUE) %>% # Break any rare R2 ties
  arrange(desc(r2)) %>%
  head(20)

# #pull associated drivers for the top 10
# top_10_cluster_ids = trait_cluster_ref %>%
#   filter(Is_Cluster_Exemplar == TRUE) %>%
#   distinct(Cluster_ID, .keep_all = TRUE) %>%
#   arrange(desc(r2)) %>%
#   head(10)


top_10_with_all_collinears = trait_cluster_ref %>%
  # Filter for only clusters belonging to the top 10
  filter(Cluster_ID %in% top_10_non_collinear) %>%
  # Arrange so the cluster exemplar is always listed first, followed by collinear traits by R2
  arrange(match(Cluster_ID, top_10_non_collinear), desc(Is_Cluster_Exemplar), desc(r2)) %>%
  dplyr::select(
    Cluster_ID, 
    Is_Cluster_Exemplar, 
    Metric, 
    r2, 
    p_val
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

# write.csv(
#   top_10_summary_wide, 
#   here("Outputs", "top_10_non_collinear_clusters_with_all_traits.csv"), 
#   row.names = FALSE
# )

# =========================================================================
# 1. Site Scores & Ellipses from FULL UNFILTERED NMDS
# =========================================================================
nmds_scores_unfiltered <- as.data.frame(scores(nmds_unfiltered, display = "sites")) %>%
  bind_cols(meta_unfiltered %>% dplyr::select(site, SampleID_clean, species, genotype)) %>%
  left_join(
    ed %>% dplyr::select(SampleID_clean, K_3), 
    by =  "SampleID_clean"
  ) %>%
  mutate(
    K_3  = factor(K_3),
    site = factor(site),
    HH   = paste0(toupper(species), "_", K_3)
  )

# =========================================================================
# 2. Extract Vectors from FULL envfit, filtered to Top 10 Exemplars
# =========================================================================
# Vector of exemplar metric names from your collinearity clustering
exemplar_metrics <- top_10_non_collinear %>% pull(Metric)

# Extract ALL vector scores from the full unfiltered envfit model
vector_scores_unfiltered <- as.data.frame(scores(ef_unfiltered, display = "vectors")) %>%
  mutate(
    Metric = rownames(.),
    r2     = ef_unfiltered$vectors$r,
    p_val  = ef_unfiltered$vectors$pvals
  ) %>%
  # FILTER STEP: Keep ONLY the 10 exemplar drivers!
  filter(Metric %in% exemplar_metrics)

# Calculate arrow multiplier relative to full NMDS sample coordinates
max_site_coord   <- max(abs(c(nmds_scores_unfiltered$NMDS1, nmds_scores_unfiltered$NMDS2)))
max_vector_coord <- max(sqrt(vector_scores_unfiltered$NMDS1^2 + vector_scores_unfiltered$NMDS2^2))
arrow_mult       <- (max_site_coord * 0.8) / max_vector_coord

# Apply Master Trait Family regex mapping & coordinate scaling
top_vectors_clean_10 <- vector_scores_unfiltered %>%
  mutate(
    # Master Family Category Mapping
    Metric_Category = case_when(
      # Quenching - Max / Non-Photochemical (qm, rqm)
      grepl("r?qm|npq", Metric, ignore.case = TRUE) ~ "qm (Max / Non-Photochemical Quenching)",
      
      # Quenching - Photochemical (qp, mpq, mqP, ppq, qqP)
      grepl("m?qp|p?qp|qqp|mpq|ppq", Metric, ignore.case = TRUE) ~ "qP (Photochemical Quenching)",
      
      # Absorption (ABQ, rABQ)
      grepl("r?abq", Metric, ignore.case = TRUE) ~ "ABQ (Absorption)",
      
      # Quantum Yield
      grepl("quant|fvfm|yield|fv_fm", Metric, ignore.case = TRUE) ~ "Quant (Quantum Yield)",
      
      # Kinetics & Antenna
      grepl("sigma|sig", Metric, ignore.case = TRUE) ~ "Sigma (Antenna Size)",
      grepl("tau1|tau_1|\\bt1\\b", Metric, ignore.case = TRUE) ~ "Tau1 (Transport Kinetics 1)",
      grepl("tau2|tau_2|\\bt2\\b", Metric, ignore.case = TRUE) ~ "Tau2 (Transport Kinetics 2)",
      grepl("connect|conn|ncon", Metric, ignore.case = TRUE) ~ "Connect (Connectivity)",
      
      TRUE ~ "Other Metric"
    ),
    
    # Extract clean 3-letter/4-letter short tag for in-plot labels
    Trait_Family = str_extract(Metric, "(?i)(mQuant|qqP|mqP|ppq|mpq|qP|rqm|qm|rABQ|ABQ|Sigma|Tau1|Tau2|NPQ|nCon)"),
    
    # Extract condition/lighting phase prefix (e.g. L3CampR1)
    Phase = str_extract(Metric, "^[A-Za-z0-9]+"),
    
    # Construct clean publication label
    Clean_Label = case_when(
      !is.na(Trait_Family) & !is.na(Phase) ~ paste0(Trait_Family, " (", Phase, ")"),
      !is.na(Trait_Family)                 ~ Trait_Family,
      TRUE                                 ~ Metric
    ),
    
    # Apply arrow_mult to scale vector endpoints for NMDS display
    NMDS1_scaled = NMDS1 * arrow_mult,
    NMDS2_scaled = NMDS2 * arrow_mult
  )

# =========================================================================
# 3. Plot Full Space with Non-Collinear Vector Overlay
# =========================================================================
p_nmds_noncollinear<-
  ggplot() +
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
    data = top_vectors_clean_10, 
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
    arrow = arrow(length = unit(0.20, "cm")), 
    linewidth = 0.85
  ) +
  geom_text(
    data = top_vectors_clean_10,
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Trait_Family),
    #fontface = "bold",
    size = 3
  ) +
  scale_color_manual(values = category_colors, name = "Trait Family") +
  # scale_color_viridis_d(option = "turbo", name = "Trait Family") +
  
  # --- 3. BOUNDS & SPACING ---
  coord_cartesian(clip = "off") +
  expand_limits(
    x = c(min(nmds_scores$NMDS1) * 1.25, max(nmds_scores$NMDS1) * 1.25),
    y = c(min(nmds_scores$NMDS2) * 1.25, max(nmds_scores$NMDS2) * 1.25)
  ) +
  theme_bw() +
  theme(
    plot.margin = margin(15, 25, 15, 25, "pt"),
    legend.position = "right"
  )+
  #   scale_color_brewer(palette = "Set1", name = "Trait Family", na.value = "grey70") +
  # theme_classic(base_size = 14) +
  #   theme(
  #     legend.position = "right",
  #     legend.box      = "vertical"
  #   ) +
  coord_cartesian(clip = "off") +
  expand_limits(
    x = c(min(nmds_scores$NMDS1) * 1.25, max(nmds_scores$NMDS1) * 1.25),
    y = c(min(nmds_scores$NMDS2) * 1.25, max(nmds_scores$NMDS2) * 1.25)
  ) +
  theme_bw() +
  theme(
    plot.margin = margin(15, 25, 15, 25, "pt"),
    legend.position = "right"
  ) +
  labs(
    title = paste0("Unfiltered Ordination with Top 10 Non-Collinear Vectors (Stress = ", round(nmds_unfiltered$stress, 3), ")"),
    x = "nMDS Dimension 1",
    y = "nMDS Dimension 2"
  )

ggsave(
  filename = here("Plots", "nMDS_K3_unfiltered_NONCOLLINEAR20.png"),
  plot     = p_nmds_noncollinear,
  width    = 9.5,
  height   = 8.0,
  dpi      = 300,
  units    = "in"
)

