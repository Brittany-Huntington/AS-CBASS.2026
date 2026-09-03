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
pp <- read.csv(here("data", "permanova_matrix_clean.csv")) 
ed <- read.csv(here("Outputs", "ed50_colorslope_cluster_metadata.csv")) # from data compiler
load(here("Data", "metric_families.Rdata"))

ed_cleaned <- ed %>%
  mutate(SampleID_clean = str_trim(as.character(SampleID_clean)))

pp_formatted <- pp %>%
  mutate(
    SampleID_clean = as.character(X),
    site           = str_extract(SampleID_clean, "^\\d+"),
    species        = str_extract(SampleID_clean, "(?<=_)[A-Za-z]+"),
    genotype       = str_extract(SampleID_clean, "\\d+$")
  )

meta_cols  <- c("SampleID_clean", "site", "species", "genotype")
trait_cols <- setdiff(names(pp_formatted), c(meta_cols, "X"))
pp_formatted <- pp_formatted[, c(meta_cols, trait_cols)]

# Matching Key Function to join ED50/Species metadata
create_matching_key <- function(df) {
  df %>%
    mutate(
      clean_file       = str_remove(SampleID_clean, "\\.[a-zA-Z0-9]+$"),
      tp_num           = str_match(clean_file, "^t(\\d+)")[, 2],
      geno_name        = toupper(str_match(clean_file, "^t\\d+_([a-zA-Z]+)")[, 2]),
      rep_num          = str_match(clean_file, "^t\\d+_[a-zA-Z]+(\\d+)")[, 2],
      GroupingProperty = paste(tp_num, geno_name, rep_num, sep = "_")
    ) %>%
    dplyr::select(-clean_file, -tp_num, -geno_name, -rep_num)
}

# Apply key creation to full dataset
pp_with_key <- create_matching_key(pp_formatted)

# Join Species and metadata
pp_species <- pp_with_key %>%
  mutate(SampleID_clean = trimws(SampleID_clean)) %>%
  left_join(
    ed_cleaned %>%
      mutate(SampleID_clean = trimws(SampleID_clean)) %>%
      dplyr::select(SampleID_clean, Species),
    by = "SampleID_clean"
  ) %>%
  filter(!is.na(species))

# =========================================================================
# 2. Filtering Process (from Script 2)
# =========================================================================
# Extract raw trait matrix
trait_data <- pp_species %>%
  dplyr::select(all_of(trait_cols)) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

# Drop zero-variance columns and columns with all NAs
non_zero_var_mask <- sapply(trait_data, function(x) {
  v <- var(x, na.rm = TRUE)
  !is.na(v) && v > 0
})

trait_data_clean <- trait_data[, non_zero_var_mask]

# Compute Pearson correlation matrix
cor_matrix <- cor(trait_data_clean, use = "pairwise.complete.obs", method = "pearson")
cor_matrix[is.na(cor_matrix)]  <- 0
cor_matrix[is.nan(cor_matrix)] <- 0

# Run findCorrelation at desired cutoff threshold (e.g., 0.90)
cor_cutoff <- 0.99
high_cor_idx <- findCorrelation(cor_matrix, cutoff = cor_cutoff)

# Identify traits to retain
retained_trait_names <- colnames(trait_data_clean)[-high_cor_idx]

cat("Total initial traits:", length(trait_cols), "\n")
cat("Traits retained after non-zero variance & correlation (r =", cor_cutoff, ") filtering:", length(retained_trait_names), "\n")

# Subset filtered traits
trait_matrix_filtered <- trait_data_clean[, retained_trait_names]

# Z-score scale the matrix
matrix_filtered_scaled <- scale(trait_matrix_filtered)
matrix_filtered_scaled[is.na(matrix_filtered_scaled)] <- 0

# Extract metadata corresponding to rows
meta_filtered <- pp_species %>%
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
permanova_filtered <- adonis2(
  matrix_filtered_scaled ~ species, 
  data         = meta_filtered, 
  method       = "euclidean", 
  permutations = 999
)

cat("\n--- FILTERED PERMANOVA RESULTS ---\n")
print(permanova_filtered)

# Homogeneity of Multivariate Dispersion (BETADISPER)
dist_filtered <- vegdist(matrix_filtered_scaled, method = "euclidean")
bd_filtered   <- betadisper(dist_filtered, meta_filtered$species)
anova_bd      <- anova(bd_filtered)

cat("\n--- FILTERED DISPERSION TEST p-value:", anova_bd[1, "Pr(>F)"], "---\n")

# Pairwise PERMANOVA Function
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

pw_filtered <- pairwise.adonis2(matrix_filtered_scaled ~ species, data = meta_filtered, method = "euclidean")
cat("\n--- FILTERED PAIRWISE PERMANOVA ---\n")
print(pw_filtered)

# =========================================================================
# 4. nMDS & Trait Vector Fitting
# =========================================================================
set.seed(123)
nmds_filtered <- metaMDS(
  matrix_filtered_scaled, 
  distance      = "euclidean", 
  k             = 2, 
  trymax        = 200, 
  maxit         = 500, 
  autotransform = FALSE,
  trace         = 0
)

cat("\nFiltered nMDS Stress Value:", nmds_filtered$stress, "\n")

# Fit vectors for retained traits
#ef_filtered <- envfit(nmds_filtered, matrix_filtered_scaled, permutations = 999)
#saveRDS(ef_filtered, here("Outputs", "ef_filtered_nmds.rds"))
ReadRDS(here("Outputs/ef_filtered_nmds.rds"))


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
  head(30)

cat("\n--- TOP METRIC DRIVERS OF FILTERED ORDINATION ---\n")
print(top_vectors)

nmds_scores <- as.data.frame(scores(nmds_filtered, display = "sites")) %>%
  bind_cols(meta_filtered %>% dplyr::select(site, SampleID_clean, species, genotype))

nmds_scores <- nmds_scores %>%
  left_join(
    ed %>% dplyr::select(SampleID_clean, K_3, Site), 
    by = "SampleID_clean"
  ) %>%
  mutate(K_3 = factor(K_3))

max_site_coord   <- max(abs(c(nmds_scores$NMDS1, nmds_scores$NMDS2)))
max_vector_coord <- max(sqrt(top_vectors$NMDS1^2 + top_vectors$NMDS2^2))
arrow_mult       <- (max_site_coord * 0.8) / max_vector_coord

top_vectors_clean <- top_vectors %>%
  mutate(
    # Master Family Category Mapping
    Metric_Category = case_when(
      grepl("r?qm|npq", Metric, ignore.case = TRUE) ~ "qm (Max / Non-Photochemical Quenching)",
      grepl("m?qp|p?qp|qqp|mpq|ppq", Metric, ignore.case = TRUE) ~ "qP (Photochemical Quenching)",
      grepl("r?abq", Metric, ignore.case = TRUE) ~ "ABQ (Absorption)",
      grepl("quant|fvfm|yield|fv_fm", Metric, ignore.case = TRUE) ~ "Quant (Quantum Yield)",
      grepl("sigma|sig", Metric, ignore.case = TRUE) ~ "Sigma (Antenna Size)",
      grepl("tau1|tau_1|\\bt1\\b", Metric, ignore.case = TRUE) ~ "Tau1 (Transport Kinetics 1)",
      grepl("tau2|tau_2|\\bt2\\b", Metric, ignore.case = TRUE) ~ "Tau2 (Transport Kinetics 2)",
      grepl("connect|conn|ncon", Metric, ignore.case = TRUE) ~ "Connect (Connectivity)",
      TRUE ~ "Other Metric"
    ),
    
    # Extract clean short tag for in-plot labels
    Trait_Family = str_extract(Metric, "(?i)(mQuant|qqP|mqP|ppq|mpq|qP|rqm|qm|rABQ|ABQ|Sigma|Tau1|Tau2|NPQ|nCon)"),
    Phase        = str_extract(Metric, "^[A-Za-z0-9]+"),
    Clean_Label  = case_when(
      !is.na(Trait_Family) & !is.na(Phase) ~ paste0(Trait_Family, " (", Phase, ")"),
      !is.na(Trait_Family)                 ~ Trait_Family,
      TRUE                                 ~ Metric
    ),
    
    NMDS1_scaled = NMDS1 * arrow_mult,
    NMDS2_scaled = NMDS2 * arrow_mult
  )


# Define clean color palette for K=3 Clusters
cluster_colors = c(
  "1" = "orange",  
  "2" = "violet",  
  "3" = "grey"   
)

# =========================================================================
# 5. Plot Filtered nMDS
# =========================================================================
spp <- ggplot() +
  geom_point(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = species, shape = site), 
    size = 2.8, 
    alpha = 0.75
   ) + 
  # scale_shape_manual(values = custom_shapes, name = "Site") +
  # scale_color_manual(values = species_colors, name = "Species") +
  # scale_fill_manual(values = species_colors, name = "Species") +
  #new_scale_color() +
  #new_scale_fill() +
  stat_ellipse(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = species, fill = species), 
    geom = "polygon", 
    alpha = 0.15, 
    level = 0.95
  ) +
    scale_shape_manual(values = custom_shapes, name = "Site") +
    scale_color_manual(values = species_colors, name = "Species") +
    scale_fill_manual(values = species_colors, name = "Species") +
  # scale_color_manual(values = cluster_colors, name = "Hclust Group (K_3)") +
  # scale_fill_manual(values = cluster_colors, name = "Hclust Group (K_3)") +
  new_scale_color() +
  geom_segment(
    data = top_vectors_clean, 
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled, color = Metric_Category),
    arrow = arrow(length = unit(0.20, "cm")), 
    linewidth = 0.85
  ) +
  geom_text(
    data = top_vectors_clean,
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Clean_Label),
    size = 3
  ) +
  scale_color_manual(values = category_colors, name = "Trait Family") +
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
    #title = paste0("Filtered (r = ", cor_cutoff, ", nMDS, Stress = ", round(nmds_filtered$stress, 3), ")"),
    title = "Ellipses are species",
    x = "nMDS Dimension 1",
    y = "nMDS Dimension 2"
  )

# ggsave(
#   filename = here("Plots", paste0("nMDS_spp_filtered_r_10", cor_cutoff * 100, ".png")),
#   plot     = spp,
#   width    = 9.5,
#   height   = 8.0,
#   dpi      = 300,
#   units    = "in"
# )

# =========================================================================
# 6. Plot Filtered nMDS by K
# =========================================================================
kplot <- ggplot() +
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
  stat_ellipse(
    data = nmds_scores, 
    aes(x = NMDS1, y = NMDS2, color = K_3, fill = K_3), 
    geom = "polygon", 
    alpha = 0.15, 
    level = 0.95
  ) +
  # scale_shape_manual(values = custom_shapes, name = "Site") +
  # scale_color_manual(values = species_colors, name = "Species") +
  # scale_fill_manual(values = species_colors, name = "Species") +
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
    aes(x = NMDS1_scaled * 1.10, y = NMDS2_scaled * 1.10, label = Clean_Label),
    size = 3
  ) +
  scale_color_manual(values = category_colors, name = "Trait Family") +
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
   # title = paste0("Filtered (r = ", cor_cutoff, ", nMDS, Stress = ", round(nmds_filtered$stress, 3), ")"),
    title = "Ellipses are groupings from hclust (K=3)",
    x = "nMDS Dimension 1",
    y = "nMDS Dimension 2"
  )

# ggsave(
#   filename = here("Plots", paste0("nMDS_K_filtered_r_10", cor_cutoff * 100, ".png")),
#   plot     = kplot,
#   width    = 9.5,
#   height   = 8.0,
#   dpi      = 300,
#   units    = "in"
# )

library(patchwork)

# 1. Clean up individual plot titles/labels for a side-by-side layout
# Remove the y-axis label from the right plot to avoid redundancy
kplot_clean <- kplot + 
  labs(y = NULL)

# 2. Combine side-by-side and collect shared legends
combined_plot <- (spp + kplot_clean) +
  plot_layout(
    ncol = 2, 
    guides = "collect"  # Combines identical scale legends (Site, Species, Trait Family)
  ) +
  plot_annotation(
    title = paste0("Filtered (r = ", cor_cutoff, ", nMDS, Stress = ", round(nmds_filtered$stress, 3), ")"),
    tag_levels = 'A'  # Adds 'A' and 'B' tags to the subplots
  ) &
  theme(legend.position = "right")

# 3. Save the combined plot with an expanded width
ggsave(
  filename = here("Plots", paste0("nMDS_combined_r_10", cor_cutoff * 100, ".png")),
  plot     = combined_plot,
  width    = 16.0,  # Widen to comfortably fit both panels and legends
  height   = 9.0,
  dpi      = 300,
  units    = "in"
)

