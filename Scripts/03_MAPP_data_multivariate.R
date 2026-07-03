library(tidyverse)
library(here)
library(stringr)
library(pvclust)
library(dendextend)
library(parallel)
library(broom)

rm(list = ls())

####Read Data-------------
pp <- read.csv(here("data", "permanova_matrix_clean.csv")) 
ed <- read.csv(here("data", "EDsdf-all.csv"))

str(pp)


pp_formatted <- pp
pp_formatted$sample_file <- pp_formatted$X
pp_formatted$timepoint <- str_extract(pp_formatted$sample_file, "^[^_]+")
pp_formatted$genotype <- str_extract(pp_formatted$sample_file, "(?<=_)[A-Za-z0-9-]+")

meta_cols <- c("sample_file", "timepoint", "genotype")
trait_cols <- setdiff(names(pp_formatted), c(meta_cols, "X"))

pp_formatted <- pp_formatted[, c(meta_cols, trait_cols)]


#run correlation filter on the new formatted MAPP df and drop redundant MAPP traits-------------------
metadata_cols <- c("X", "sample_file", "timepoint", "genotype")
trait_data <- pp_formatted[, !(names(pp_formatted) %in% c("sample_file", "timepoint", "genotype"))]
cor_matrix <- cor(trait_data, use = "pairwise.complete.obs", method = "pearson")
abs_cor_matrix <- abs(cor_matrix) #Clear out the lower triangle and diagonal to avoid double-dropping
abs_cor_matrix[lower.tri(abs_cor_matrix, diag = TRUE)] <- 0

# Identify columns with zero variance (which cause NaNs in correlations) and columns that have a correlation > 0.99 with another trait (following Hoadley et al. 2024)
zero_var_cols <- colnames(trait_data)[apply(trait_data, 2, var, na.rm = TRUE) == 0]
high_cor_cols <- colnames(abs_cor_matrix)[apply(abs_cor_matrix, 2, function(x) any(x > 0.99, na.rm = TRUE))] #following Hoadley
med_cor_cols <- colnames(abs_cor_matrix)[apply(abs_cor_matrix, 2, function(x) any(x > 0.90, na.rm = TRUE))] #exploratory to be slightly less conservative

cols_to_remove <- unique(c(zero_var_cols, high_cor_cols)) #n=150
cols_to_remove.90 <- unique(c(zero_var_cols, med_cor_cols)) #n=619

pp_99 <- pp_formatted[, !(colnames(pp_formatted) %in% cols_to_remove)]
pp_90 <- pp_formatted[, !(colnames(pp_formatted) %in% cols_to_remove.90)] #half the number of metrics



#create a column for linking this pp dataframe back to the ED50 df
create_matching_key <- function(df) {
  df %>%
    mutate(
      tp_num = str_match(sample_file, "^t(\\d+)")[,2],
      geno_name = toupper(str_match(sample_file, "^t\\d+_([a-zA-Z]+)")[,2]),
      rep_num = str_match(sample_file, "^t\\d+_[a-zA-Z]+(\\d+)")[,2],
      GroupingProperty = paste(tp_num, geno_name, rep_num, sep = "_")
    ) %>%
    select(-tp_num, -geno_name, -rep_num)
}

pp_99 <- create_matching_key(pp_99)
pp_90 <- create_matching_key(pp_90)



#building the phenotypic dendrograms and cutting them into 4 distinct clusters-------
#transpose df so coral genotypes are columns
prep_for_pvclust <- function(df) {
  traits_only <- df[, !(names(df) %in% c("sample_file", "timepoint", "genotype"))]
  pv_matrix <- t(traits_only)
  colnames(pv_matrix) <- df$sample_file
  return(as.data.frame(pv_matrix))
}

pv_data_99 <- prep_for_pvclust(pp_99)
pv_data_90 <- prep_for_pvclust(pp_90)


#Run the 10,000 Bootstrap Clusterings
num_cores <- detectCores() - 1  # Detect available cores for parallel processing

# Run pvclust for the 0.99 threshold dataset
cat("Running 10,000 bootstraps for pp_99 (this may take a few minutes)... \n")
fit_99 <- pvclust(pv_data_99, 
                  method.hclust = "average", 
                  method.dist = "correlation", 
                  nboot = 500, 
                  parallel = TRUE,
                  iseed = 123) # Set seed for reproducible bootstrap swaps

# Run pvclust for the 0.90 threshold dataset
cat("Running 10,000 bootstraps for pp_90... \n")
fit_90 <- pvclust(pv_data_90, 
                  method.hclust = "average", 
                  method.dist = "correlation", 
                  nboot = 500, 
                  parallel = TRUE,
                  iseed = 123)


ed_99_multiverse <- pp_99 %>% select(GroupingProperty)
ed_90_multiverse <- pp_90 %>% select(GroupingProperty)

#Run the loop to create clusters for k = 4 through k = 8
for (k in 4:8) {
  # Create a dynamic column name (e.g., "cluster_k4", "cluster_k5")
  col_name <- paste0("cluster_k", k)
  
  # Cut the trees and assign as factors
  ed_99_multiverse[[col_name]] <- factor(cutree(fit_99$hclust, k = k))
  ed_90_multiverse[[col_name]] <- factor(cutree(fit_90$hclust, k = k))
}


# join with ed df
ed_99_multiverse <- ed_99_multiverse %>% left_join(ed, by = "GroupingProperty")
ed_90_multiverse <- ed_90_multiverse %>% left_join(ed, by = "GroupingProperty")


#ANOVA for ED50 to clusters-----------------------------
#Q: Does a coral's multi-trait physiological profile predict its ultimate thermal threshold (ED_50) better than its taxonomic identity?----

# Function to screen ANOVAs across all cluster sizes
screen_cluster_resolutions <- function(multiverse_df) {
  results <- list()
  
  for (k in 4:8) {
    col_name <- paste0("cluster_k", k)
    
    # Run the ANOVA dynamically using reformulate()
    fit <- aov(reformulate(col_name, response = "ED50"), data = multiverse_df)
    tidy_fit <- tidy(fit)
    
    # Calculate R-squared (Sum of Squares Group / Total Sum of Squares)
    ss_group <- tidy_fit$sumsq[1]
    ss_total <- sum(tidy_fit$sumsq)
    r_squared <- ss_group / ss_total
    
    results[[col_name]] <- data.frame(
      Resolution = paste0("k = ", k),
      p_value = tidy_fit$p.value[1],
      R_squared = r_squared
    )
  }
  
  bind_rows(results)
}

# Check the results for both thresholds
print("--- Evaluation of Cluster Sizes for pp_99 ---")
print(screen_cluster_resolutions(ed_99_multiverse))

print("--- Evaluation of Cluster Sizes for pp_90 ---")
print(screen_cluster_resolutions(ed_90_multiverse))
print(abs_tukey_table)



#Multivariate analysis---------------------------------
#Q: Do different coral species have distinct, signature physiological profiles under acute heat stress?----

library(vegan)
# 1. Bring 'Species' into the physiological dataframes
pp_99_species <- pp_99%>%
  left_join(ed %>% select(GroupingProperty, Species), by = "GroupingProperty") %>%
  filter(!is.na(Species)) # Remove any samples that don't have matching species data

pp_90_species <- pp_90 %>%
  left_join(ed %>% select(GroupingProperty, Species), by = "GroupingProperty") %>%
  filter(!is.na(Species))


# 2. Extract clean, standardized physiological matrices (Y)
matrix_99_scaled <- pp_99_species %>%
  select(where(is.numeric)) %>%  # This automatically drops ALL non-numeric columns
  scale()
matrix_99_scaled[is.na(matrix_99_scaled)] <- 0

matrix_90_scaled <- pp_90_species %>%
  select(where(is.numeric)) %>%
  scale()
matrix_90_scaled[is.na(matrix_90_scaled)] <- 0

# 3. Extract matching metadata frames (X)
meta_99 <- pp_99_species %>% select(sample_file, timepoint, genotype, GroupingProperty, Species)
meta_90 <- pp_90_species %>% select(sample_file, timepoint, genotype, GroupingProperty, Species)


permanova_99 <- adonis2(matrix_99_scaled ~ Species, data = meta_99, method = "euclidean", permutations = 999)
print(permanova_99)

permanova_90 <- adonis2(matrix_90_scaled ~ Species,data = meta_90, method = "euclidean", permutations = 999)
print(permanova_90)

#pairwise post-hoc tests to see which years differed from one another
#load function
pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- eval(x1[[2]], environment(x1), globalenv())
  environment(x1) <- environment()
  # extract factors on right hand side of formula
  rhs <- x1[[3]]
  # create model.frame matrix
  x1[[2]] <- NULL
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE)
  
  # create unique pairwise combination of factors
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors
  for(elem in 1:ncol(co)){
    
    #reduce model elements
    if(inherits(eval(lhs),'dist')){
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))
    }else{
      xnew <- as.formula(paste('xred' ,
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names
  class(res) <- c("pwadstrata", "list")
  return(res)
}


pw99_SP <- pairwise.adonis2(matrix_99_scaled ~ Species, data = meta_99, method = "euclidean", p.adjust = "BH")
pw99_SP 

pw90_SP <- pairwise.adonis2(matrix_90_scaled ~ Species, data = meta_99, method = "euclidean", p.adjust = "BH")
pw90_SP #no species differences

#explore which specific taxa are differing among OBS_YEAR
sim <- simper(matrix_99_scaled, group = meta_99$Species, permutations = 999)
sim
summary(sim) 



sim.90 <- simper(matrix_90_scaled, group = meta_90$Species, permutations = 999)
sim.90
summary(sim.90) 


#check PERMANOVA assumption of equal dispersion among groups/levels for each factor
taxa_distmat <- vegdist((matrix_99_scaled), method = "euclidian")
bd <-  betadisper(taxa_distmat, meta_99$Species)
boxplot(bd)
anova(bd) #passes assumption

taxa_distmat <- vegdist((matrix_90_scaled), method = "euclidian")
bd <-  betadisper(taxa_distmat, meta_90$Species)
boxplot(bd)
anova(bd) #passes assumption
#TukeyHSD(bd, ordered = FALSE, conf.level = 0.95) #optional post-hoc to see where dispersion may differ among factor levels



#plotting---PCA
pca_fit <- prcomp(matrix_99_scaled, center = FALSE, scale. = FALSE) # Data is already scaled

pca_scores <- as.data.frame(pca_fit$x) %>%
  bind_cols(pp_99_species %>% select(Species, genotype))

var_explained <- round(100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2)), 1)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = Species), geom = "polygon", alpha = 0.1, level = 0.95) +
  labs(
    title = "Multivariate Physiological Profiles by Coral Species",
    x = paste0("PC1 (", var_explained[1], "% variance explained)"),
    y = paste0("PC2 (", var_explained[2], "% variance explained)")
  ) +
  theme_minimal() +
  theme(legend.position = "right")


#plotting----nMDS
# 1. Run the nMDS using Euclidean distance on your scaled, patched matrix
# autotransform = FALSE prevents metaMDS from trying to log-transform your negative scaled numbers
nmds_fit <- metaMDS(matrix_99_scaled, 
                    distance = "euclidean", 
                    k = 2, 
                    trymax = 100, 
                    autotransform = FALSE)

# 2. CRITICAL: Check the stress value in your console
cat("nMDS Stress Value:", nmds_fit$stress, "\n")

# 3. Extract the coordinates for plotting
nmds_scores <- as.data.frame(scores(nmds_fit, display = "sites")) %>%
  bind_cols(pp_99_species %>% select(Species, genotype))

# 4. Plot the nMDS
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Species)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = Species), geom = "polygon", alpha = 0.1, level = 0.95) +
  labs(
    title = "nMDS Ordination of Physiological Profiles (Euclidean)",
    subtitle = paste0("Stress = ", round(nmds_fit$stress, 3))
  ) +
  theme_minimal()


library(vegan)
library(dplyr)
library(ggplot2)

