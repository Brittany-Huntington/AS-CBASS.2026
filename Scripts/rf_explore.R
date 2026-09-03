
##############################################################################
#### RANDOM FOREST FEATURE SELECTION (SPECIES IDENTITY) ###################
##############################################################################
library(randomForest)
# Build dataset for Random Forest Classification
rf_data <- as.data.frame(permanova_matrix_filtered)
rf_data$Species <- factor(meta_filtered$Species)

set.seed(123)
rf_species <- randomForest(
  Species ~ ., 
  data        = rf_data, 
  ntree       = 1000, 
  importance  = TRUE,
  proximity   = TRUE
)

cat("\n--- RANDOM FOREST SPECIES CLASSIFICATION CONFUSION MATRIX ---\n")
print(rf_species$confusion)

# Extract top 15 most important features driving species classification
imp_spp <- as.data.frame(importance(rf_species)) %>%
  rownames_to_column("Metric") %>%
  left_join(metric_families, by = "Metric") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  head(10)

cat("\n--- TOP 10 METRICS DISTINGUISHING SPECIES (RANDOM FOREST) ---\n")
print(imp_spp %>% select(Metric, Family))


# =========================================================================
# A. RANDOM FOREST REGRESSION: PREDICTING ED50 (THERMAL TOLERANCE)
# =========================================================================

# Add ED50 and filter missing values
rf_ed50_df <- rf_data
rf_ed50_df$ED50 <- meta_filtered$ED50
rf_ed50_df <- rf_ed50_df %>% filter(!is.na(ED50))

set.seed(123)
rf_ed50 <- randomForest(
  ED50 ~ ., 
  data       = rf_ed50_df, 
  ntree      = 1000, 
  importance = TRUE
)

cat("\n--- RANDOM FOREST ED50 REGRESSION SUMMARY ---\n")
print(rf_ed50)

# Extract top 15 features predicting ED50
imp_ed50 <- as.data.frame(importance(rf_ed50)) %>%
  rownames_to_column("Metric") %>%
  left_join(metric_families, by = "Metric") %>%
  arrange(desc(`%IncMSE`)) %>%
  head(10)

cat("\n--- TOP 10 METRICS PREDICTING ED50 (%IncMSE) ---\n")
print(imp_ed50 %>% select(Metric, Family))


# =========================================================================
# B. RANDOM FOREST CLASSIFICATION: PREDICTING K_3 CLUSTERS
# =========================================================================

# Add K_3 clusters as target factor
rf_k3_df <- rf_data
rf_k3_df$K_3 <- factor(meta_filtered$K_3)

set.seed(123)
rf_k3 <- randomForest(
  K_3 ~ ., 
  data       = rf_k3_df, 
  ntree      = 1000, 
  importance = TRUE
)

cat("\n--- RANDOM FOREST K_3 CONFUSION MATRIX ---\n")
print(rf_k3$confusion)

# Extract top 15 features distinguishing functional clusters
imp_k3 <- as.data.frame(importance(rf_k3)) %>%
  rownames_to_column("Metric") %>%
  left_join(metric_families, by = "Metric") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  head(15)

cat("\n--- TOP 15 METRICS DISTINGUISHING K_3 CLUSTERS ---\n")
print(imp_k3 %>% select(Metric, Family, MeanDecreaseGini, MeanDecreaseAccuracy))


# =========================================================================
# C. VARIABLE IMPORTANCE PLOT (EXAMPLE FOR ED50)
# =========================================================================

p_imp_ed50 <- ggplot(imp_ed50, aes(x = reorder(Metric, `%IncMSE`), y = `%IncMSE`, fill = Family)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Top Photophysiological Drivers of Thermal Tolerance (ED50)",
    x     = "Photophysiology Metric",
    y     = "% Increase in MSE (Feature Importance)"
  )

print(p_imp_ed50)