library(tidyverse)
library(here)
library(stringr)
library(pvclust)
library(dendextend)
library(parallel)
library(broom)
library(vegan)
library(caret)

rm(list = ls())
# species_colors <- c(
#   "AABR" = "#EE8080FF", 
#   "AGLO" = "#83BA75FF", 
#   "AHYA" = "#5F8DAAFF", 
#   "ICRA" = "#6E5B8AFF"  
# )"AGLO" = "#9DD9D2", # Light Sage Teal (#9DD9D2FF)
# "AHYA" = "#046E8F", # Deep Ocean Blue (#046E8FFF)
# 
# species_colors <- c(
#   "AABR" = "#D44D5C", # Muted Crimson/Rose Red (#D44D5CFF)
#   "AGLO" = "#83BA75FF", # Light Sage Teal (#9DD9D2FF)
#   "AHYA" = "#462255", # Deep Ocean Blue (#046E8FFF)
#   "ICRA" = "#FF8811"  # Bright Amber/Orange (#FF8811FF) -- (or "#462255" for Deep Purple)
# )
species_colors <- c(
  "AABR" = "#046E8F", # Muted Crimson/Rose Red (#D44D5CFF)
  "AGLO" = "#83BA75FF", # Light Sage Teal (#9DD9D2FF)
  "AHYA" = "#D44D5C", # Deep Ocean Blue (#046E8FFF)
  "ICRA" = "#462255"  # Bright Amber/Orange (#FF8811FF) -- (or "#462255" for Deep Purple)
)

species_colors <- c(
  "AABR" = "#046E8F", 
  "AGLO" = "lightgreen", 
  "AHYA" = "#D44D5C",
  "ICRA" = "#462255"  
)

# =========================================================================
# 1. Read & Format Data
# =========================================================================
pp <- read.csv(here("data", "permanova_matrix_clean.csv")) 
ed <- read.csv(here("data", "EDsdf-all.csv"))

pp_formatted <- pp %>%
  mutate(
    sample_file = X,
    timepoint   = str_extract(sample_file, "^[^_]+"),
    genotype    = str_extract(sample_file, "(?<=_)[A-Za-z0-9-]+")
  )

meta_cols <- c("sample_file", "timepoint", "genotype")
trait_cols <- setdiff(names(pp_formatted), c(meta_cols, "X"))
pp_formatted <- pp_formatted[, c(meta_cols, trait_cols)]

# =========================================================================
# 2. Robust Correlation Filtering (caret::findCorrelation)
# =========================================================================

# 1. Extract pure numeric trait data
trait_data <- pp_formatted[, trait_cols]

# 2. Drop zero-variance columns (sd == 0 or all NAs)
non_zero_var_mask <- sapply(trait_data, function(x) {
  v <- var(x, na.rm = TRUE)
  !is.na(v) && v > 0
})
trait_data_clean <- trait_data[, non_zero_var_mask]

# Compute correlation matrix
cor_matrix <- cor(trait_data_clean, use = "pairwise.complete.obs", method = "pearson")

# Replace any residual NAs/NaNs in cor_matrix with 0
cor_matrix[is.na(cor_matrix)] <- 0
cor_matrix[is.nan(cor_matrix)] <- 0

# Run findCorrelation safely
high_cor_idx <- findCorrelation(cor_matrix, cutoff = 0.99)
med_cor_idx  <- findCorrelation(cor_matrix, cutoff = 0.90)

# Get column names to remove
cols_to_remove_99 <- c(
  colnames(trait_data)[!non_zero_var_mask], 
  colnames(trait_data_clean)[high_cor_idx]  
)

cols_to_remove_90 <- c(
  colnames(trait_data)[!non_zero_var_mask], 
  colnames(trait_data_clean)[med_cor_idx]   
)

# Subset formatted dataframes cleanly
pp_99 <- pp_formatted[, !(colnames(pp_formatted) %in% cols_to_remove_99)]
pp_90 <- pp_formatted[, !(colnames(pp_formatted) %in% cols_to_remove_90)]

cat("Number of traits retained at r = 0.99 threshold:", ncol(pp_99) - length(meta_cols), "\n")
cat("Number of traits retained at r = 0.90 threshold:", ncol(pp_90) - length(meta_cols), "\n")

#Save pre-filtered datasets
saveRDS(pp_99, here("Outputs", "pp_99_filtered.rds"))
saveRDS(pp_90, here("Outputs", "pp_90_filtered.rds"))

# Matching Key Function
create_matching_key <- function(df) {
  df %>%
    mutate(
      tp_num = str_match(sample_file, "^t(\\d+)")[,2],
      geno_name = toupper(str_match(sample_file, "^t\\d+_([a-zA-Z]+)")[,2]),
      rep_num = str_match(sample_file, "^t\\d+_[a-zA-Z]+(\\d+)")[,2],
      GroupingProperty = paste(tp_num, geno_name, rep_num, sep = "_")
    ) %>%
    dplyr::select(-tp_num, -geno_name, -rep_num)
}

pp_99 <- create_matching_key(pp_99)
pp_90 <- create_matching_key(pp_90)

# =========================================================================
# 3. Hierarchical Clustering & Multiverse (K_3 to K_15)
# =========================================================================
prep_for_pvclust <- function(df) {
  traits_only <- df[, !(names(df) %in% c("sample_file", "timepoint", "genotype", "GroupingProperty"))]
  pv_matrix <- t(traits_only)
  colnames(pv_matrix) <- df$sample_file
  return(as.data.frame(pv_matrix))
}

pv_data_99 <- prep_for_pvclust(pp_99)
pv_data_90 <- prep_for_pvclust(pp_90)

cat("Running bootstraps for pp_99 (Canberra + Ward.D)...\n")
fit_99 <- pvclust(
  pv_data_99, 
  method.dist   = "canberra", 
  method.hclust = "ward.D", 
  nboot         = 1000, 
  parallel      = TRUE, 
  iseed         = 123
)

cat("Running bootstraps for pp_90 (Canberra + Ward.D)...\n")
fit_90 <- pvclust(
  pv_data_90, 
  method.dist   = "canberra", 
  method.hclust = "ward.D", 
  nboot         = 1000, 
  parallel      = TRUE, 
  iseed         = 123
)

# Save heavy bootstrap outputs (pvclust objects)
saveRDS(fit_99, here("Outputs", "fit_pvclust_99.rds"))
saveRDS(fit_90, here("Outputs", "fit_pvclust_90.rds"))

# Build distance trees directly on scaled trait matrices to ensure balanced k cuts
dist_99 <- vegdist(t(pv_data_99), method = "canberra", na.rm = TRUE)
tree_99 <- hclust(dist_99, method = "ward.D")

dist_90 <- vegdist(t(pv_data_90), method = "canberra", na.rm = TRUE)
tree_90 <- hclust(dist_90, method = "ward.D")

saveRDS(tree_99, here("Outputs", "ed_99_tree.rds"))
saveRDS(tree_90, here("Outputs", "ed_90_tree.rds"))



# Multiverse ANOVA across K_3 to K_15
ed_99_multiverse <- pp_99 %>% dplyr::select(GroupingProperty)
ed_90_multiverse <- pp_90 %>% dplyr::select(GroupingProperty)

# Name columns K_3, K_4, ..., K_15
for (k in 3:15) {
  col_name <- paste0("K_", k)
  ed_99_multiverse[[col_name]] <- factor(cutree(tree_99, k = k))
  ed_90_multiverse[[col_name]] <- factor(cutree(tree_90, k = k))
}

# Join ED50 metadata
ed_99_multiverse <- ed_99_multiverse %>% left_join(ed, by = "GroupingProperty")
ed_90_multiverse <- ed_90_multiverse %>% left_join(ed, by = "GroupingProperty")

# Verify K_3 groups
cat("Unique groups in ed_99_multiverse K_3:", length(unique(ed_99_multiverse$K_3)), "\n")
cat("Unique groups in ed_90_multiverse K_3:", length(unique(ed_90_multiverse$K_3)), "\n")

# Save multiverse files
saveRDS(ed_99_multiverse, here("Outputs", "ed_99_multiverse.rds"))
saveRDS(ed_90_multiverse, here("Outputs", "ed_90_multiverse.rds"))


screen_cluster_resolutions <- function(multiverse_df) {
  results <- list()
  for (k in 3:15) {
    col_name <- paste0("K_", k)
    fit <- aov(reformulate(col_name, response = "ED50"), data = multiverse_df)
    tidy_fit <- tidy(fit)
    
    ss_group <- tidy_fit$sumsq[1]
    ss_total <- sum(tidy_fit$sumsq)
    
    results[[col_name]] <- data.frame(
      Resolution = paste0("k = ", k),
      p_value = tidy_fit$p.value[1],
      R_squared = ss_group / ss_total
    )
  }
  bind_rows(results)
}

cat("--- Evaluation of Cluster Sizes for pp_99 ---\n")
print(screen_cluster_resolutions(ed_99_multiverse))

cat("--- Evaluation of Cluster Sizes for pp_90 ---\n")
print(screen_cluster_resolutions(ed_90_multiverse))

# =========================================================================
# 4. Species-Level Multivariate Analysis (PERMANOVA & SIMPER)
# =========================================================================
pp_99_species <- pp_99 %>%
  left_join(ed %>% dplyr::select(GroupingProperty, Species), by = "GroupingProperty") %>%
  filter(!is.na(Species))

pp_90_species <- pp_90 %>%
  left_join(ed %>% dplyr::select(GroupingProperty, Species), by = "GroupingProperty") %>%
  filter(!is.na(Species))

matrix_99_scaled <- pp_99_species %>% dplyr::select(where(is.numeric)) %>% scale()
matrix_99_scaled[is.na(matrix_99_scaled)] <- 0

matrix_90_scaled <- pp_90_species %>% dplyr::select(where(is.numeric)) %>% scale()
matrix_90_scaled[is.na(matrix_90_scaled)] <- 0

meta_99 <- pp_99_species %>% dplyr::select(sample_file, timepoint, genotype, GroupingProperty, Species)
meta_90 <- pp_90_species %>% dplyr::select(sample_file, timepoint, genotype, GroupingProperty, Species)

# Main PERMANOVAs
permanova_99 <- adonis2(matrix_99_scaled ~ Species, data = meta_99, method = "euclidean", permutations = 999)
permanova_90 <- adonis2(matrix_90_scaled ~ Species, data = meta_90, method = "euclidean", permutations = 999)

cat("--- PERMANOVA (0.99 Threshold) ---\n")
print(permanova_99)

cat("--- PERMANOVA (0.90 Threshold) ---\n")
print(permanova_90)

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
  
  co <- combn(unique(as.character(rhs.frame[,1])), 2)
  
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres, paste(co[1,elem], co[2,elem], sep = '_vs_'))
  }
  
  res <- vector(mode = "list", length = length(nameres))
  names(res) <- nameres
  res['parent_call'] <- list(paste(fostri[2], fostri[1], fostri[3], ', strata =', ststri, ', permutations', nperm))
  
  for(elem in 1:ncol(co)){
    if(inherits(eval(lhs), 'dist')){
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem], co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem], co[2,elem])])
    } else {
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem], co[2,elem]), ]
    }
    
    mdat1 <- data[rhs.frame[,1] %in% c(co[1,elem], co[2,elem]), ]
    
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
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  class(res) <- c("pwadstrata", "list")
  return(res)
}

pw99_SP <- pairwise.adonis2(matrix_99_scaled ~ Species, data = meta_99, method = "euclidean", p.adjust = "BH")
pw90_SP <- pairwise.adonis2(matrix_90_scaled ~ Species, data = meta_90, method = "euclidean", p.adjust = "BH")

# =========================================================================
# Homogeneity of Dispersion Tests & Ordination Plots
# =========================================================================
taxa_distmat_99 <- vegdist(matrix_99_scaled, method = "euclidean")
bd_99           <- betadisper(taxa_distmat_99, meta_99$Species)
anova_bd_99     <- anova(bd_99)

cat("--- Dispersion Test (0.99) ANOVA p-value:", anova_bd_99[1, "Pr(>F)"], "---\n")

taxa_distmat_90 <- vegdist(matrix_90_scaled, method = "euclidean")
bd_90           <- betadisper(taxa_distmat_90, meta_90$Species)
anova_bd_90     <- anova(bd_90)

cat("--- Dispersion Test (0.90) ANOVA p-value:", anova_bd_90[1, "Pr(>F)"], "---\n")

# PCA Visualization
pca_fit <- prcomp(matrix_99_scaled, center = FALSE, scale. = FALSE)
pca_scores <- as.data.frame(pca_fit$x) %>% bind_cols(meta_99)
var_exp <- round(100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2)), 1)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = Species, fill = Species)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95) +
  theme_classic(base_size = 14) +
  labs(
    #title = "Multivariate Physiological Divergence Across Coral Species",
    x = paste0("PC1 (", var_exp[1], "% Variance)"),
    y = paste0("PC2 (", var_exp[2], "% Variance)")
  )+
  # Apply muted manual palette
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) 

# =========================================================================
# nMDS & Trait Vector Fitting
# =========================================================================
set.seed(123)
nmds_fit <- metaMDS(
  matrix_90_scaled, 
  distance = "euclidean", 
  k = 2, 
  trymax = 200,          # Gives vegan more attempts to repeat the best solution
  maxit = 500,           # Fixes the 'no. of iterations >= maxit' warning
  autotransform = FALSE,
  trace = 1              # Set to 1 if you want to watch convergence in terminal
)


cat("nMDS Stress Value:", nmds_fit$stress, "\n")

ef <- envfit(nmds_fit, matrix_90_scaled, permutations = 999)

saveRDS(ef, here("Outputs", "ef_nmds.rds"))

vector_scores <- as.data.frame(scores(ef, display = "vectors"))
vector_scores$Metric <- rownames(vector_scores)
vector_scores$r2     <- ef$vectors$r
vector_scores$p_val  <- ef$vectors$pvals

top_vectors <- vector_scores %>%
  filter(p_val < 0.002) %>%
  arrange(desc(r2)) %>%
  head(11)

cat("--- TOP METRIC DRIVERS OF SPECIES ORDINATION ---\n")
print(top_vectors)

nmds_scores <- as.data.frame(scores(nmds_fit, display = "sites")) %>%
  bind_cols(meta_90 %>% dplyr::select(Species, genotype))

top_vectors_clean <- top_vectors %>%
  mutate(
    Trait_Family = str_extract(Metric, "(mQuant|qQuant|rABQ|rqm)"),
    Phase        = str_extract(Metric, "(L1Camp|L3Camp|DCamp|DCompL[0-9]+|DCompL3|L12Comp)"),
    Clean_Label  = case_when(
      !is.na(Phase) & !is.na(Trait_Family) ~ paste0(Trait_Family, " (", Phase, ")"),
      str_detect(Metric, "rABQ")            ~ "rABQ (Mean)",
      TRUE                                 ~ Trait_Family
    )
  )

arrow_mult <- 50.2 


top_vectors_clean <- top_vectors %>%
  mutate(
    # 1. Extract functional trait family
    Trait_Family = str_extract(Metric, "(qQuant|Quant|mTau2|rTau1|ABQ)"),
    
    # 2. Extract lighting phase/condition
    Phase = str_extract(Metric, "(L32Comp|L22Comp|L3Comp|L3Camp|L12Comp|L1)"),
    
    # 3. Build a clean, publication-ready label cleanly without NAs
    Clean_Label = case_when(
      !is.na(Trait_Family) & !is.na(Phase) ~ paste0(Trait_Family, " (", Phase, ")"),
      !is.na(Trait_Family)                  ~ Trait_Family,
      TRUE                                 ~ str_extract(Metric, "^[^.]+") # Fallback to first word
    )
  )

# Verify no NAs remain in Clean_Label
cat("Missing labels:", sum(is.na(top_vectors_clean$Clean_Label)), "\n")
print(top_vectors_clean %>% dplyr::select(Metric, Trait_Family, Phase, Clean_Label))

top_vectors_clean <- top_vectors_clean %>%
  mutate(
    NMDS1_scaled = NMDS1 * sqrt(r2) * arrow_mult,
    NMDS2_scaled = NMDS2 * sqrt(r2) * arrow_mult
  )

ggplot() +
  geom_point(data = nmds_scores, aes(x = NMDS1, y = NMDS2, color = Species), 
             size = 2.8, alpha = 0.75) +
  stat_ellipse(data = nmds_scores, aes(x = NMDS1, y = NMDS2, color = Species, fill = Species), 
               geom = "polygon", alpha = 0.12, level = 0.50) +
  geom_segment(data = top_vectors_clean, 
               aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled),
               arrow = arrow(length = unit(0.22, "cm")), 
               color = "grey10", linewidth = 0.8) +
  geom_text(data = top_vectors_clean, 
            aes(x = NMDS1_scaled * 1.12, y = NMDS2_scaled * 1.12, label = Clean_Label),
            color = "black", fontface = "bold", size = 3.8) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right") +
  labs(
    #title = "Interspecific Photophysiological Ordination (nMDS)",
    #subtitle = paste0("Stress = ", round(nmds_fit$stress, 3), 
   #                   " | Arrows indicate top metric drivers (R² > 0.87, p < 0.002)"),
    x = "nMDS Dimension 1",
    y = "nMDS Dimension 2"
  )+
  # Apply muted manual palette
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) 
