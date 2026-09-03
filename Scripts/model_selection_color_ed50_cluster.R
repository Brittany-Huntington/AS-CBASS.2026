library(lme4)
library(lmerTest)
library(tidyverse)
library(performance)
library(ggrepel)
library(here)
meta_filtered<-readRDS(here("Outputs", "master_cluster_metadata_K2_K15.rds")) #from phenotype_heatmap script
df_modeling <- meta_filtered %>%
  mutate(across(c(starts_with("K_"), Site), as.factor)) %>%
  filter(!is.na(ED50), !is.na(Site))

# Set outcome variable ("ED50" or "slope")
Ychoice <- "ED50"

# Clean dataset
df_modeling <- df_modeling %>%
  filter(!is.na(.data[[Ychoice]]), !is.na(Species), !is.na(Site)) %>%
  mutate(
    YY      = .data[[Ychoice]],
    Species = droplevels(as.factor(Species)),
    Site    = droplevels(as.factor(Site))
  )

# ==============================================================================
# 1. Fit Baseline Species-Only Model
# ==============================================================================
all_models <- list()
all_models[["spp_m"]] <- lmer(YY ~ Species + (1 | Site), data = df_modeling, REML = FALSE)

# ==============================================================================
# 2. Dynamically Loop Across K = 2 to 15
# ==============================================================================
k_range <- 2:15

for (k in k_range) {
  k_col <- paste0("K_", k)
  
  # Ensure K column exists in dataset
  if (k_col %in% names(df_modeling)) {
    
    # Ensure K is formatted as a factor
    df_modeling[[k_col]] <- droplevels(as.factor(df_modeling[[k_col]]))
    
    # Dynamic formula construction
    f_k_model      <- as.formula(paste0("YY ~ ", k_col, " + (1 | Site)"))
    f_kms          <- as.formula(paste0("YY ~ Species + ", k_col, " + (1 | Site)"))
    f_kmsi         <- as.formula(paste0("YY ~ Species * ", k_col, " + (1 | Site)"))
    #f_k_spp_random <- as.formula(paste0("YY ~ ", k_col, " + (1 | Site) + (1 | Species)"))
    
    # Fit and store models (safely wrapped to catch potential singular/rank error fits)
    all_models[[paste0("k", k, "_model")]] <- tryCatch(
      lmer(f_k_model, data = df_modeling, REML = FALSE),
      error = function(e) NULL
    )
    
    all_models[[paste0("k", k, "ms")]] <- tryCatch(
      lmer(f_kms, data = df_modeling, REML = FALSE),
      error = function(e) NULL
    )
    
    all_models[[paste0("k", k, "msi")]] <- tryCatch(
      lmer(f_kmsi, data = df_modeling, REML = FALSE),
      error = function(e) NULL
    )
    
    all_models[[paste0("k", k, "_spp_random")]] <- tryCatch(
      lmer(f_k_spp_random, data = df_modeling, REML = FALSE),
      error = function(e) NULL
    )
  }
}

# Filter out any models that failed to converge or returned NULL
all_models <- compact(all_models)

# ==============================================================================
# 3. Build Master Model Performance Table
# ==============================================================================
master_model_table <- imap_dfr(all_models, function(model_obj, model_name) {
  
  # Information criteria with ML estimator
  perf <- model_performance(model_obj, estimator = "ML")
  
  # Extract variance partition metrics safely
  r2_vals <- tryCatch(r2(model_obj), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
  
  # Extract cluster number K if applicable
  k_val <- ifelse(model_name == "spp_m", "None", str_extract(model_name, "\\d+"))
  
  # Model family categorization
  model_type <- case_when(
    model_name == "spp_m" ~ "Species Only",
    str_detect(model_name, "msi$") ~ "Species * K (Interaction)",
    str_detect(model_name, "ms$") ~ "Species + K (Additive)",
    #str_detect(model_name, "_spp_random$") ~ "K + (1|Species)",
    TRUE ~ "K Only"
  )
  
  data.frame(
    Model_Name      = model_name,
    K               = k_val,
    Model_Type      = model_type,
    AIC             = round(perf$AIC, 1),
    BIC             = round(perf$BIC, 1),
    Marginal_R2     = round(r2_vals$R2_marginal, 3),
    Conditional_R2  = round(r2_vals$R2_conditional, 3),
    stringsAsFactors = FALSE
  )
})

# View sorted table
master_model_table <- master_model_table %>% arrange(BIC)
top_models <- head(master_model_table, 6)

# ==============================================================================
# 4. Compute Pareto Efficiency Frontier Data
# ==============================================================================
# Using Marginal R2 as target (capturing fixed effect explanatory power)
# Change to Conditional_R2 if preferred
plot_data <- master_model_table %>%
  filter(!is.na(BIC), !is.na(Marginal_R2)) %>%
  rename(R2_target = Marginal_R2)

# Calculate Pareto frontier (minimum BIC for given or increasing R2)
frontier_data <- plot_data %>%
  arrange(BIC, desc(R2_target)) %>%
  mutate(cum_max_R2 = cummax(R2_target)) %>%
  filter(R2_target == cum_max_R2) %>%
  distinct(R2_target, .keep_all = TRUE)

# Tag frontier models
plot_data <- plot_data %>%
  mutate(is_frontier = Model_Name %in% frontier_data$Model_Name)

# ==============================================================================
# 5. Generate Pareto Efficiency Frontier Plot
# ==============================================================================
frontier_plot <- ggplot(plot_data, aes(x = BIC, y = R2_target)) +
  # Step line connecting frontier models
  geom_step(
    data = frontier_data,
    aes(x = BIC, y = R2_target),
    direction = "hv",
    color = "#1f77b4",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  # Suboptimal models
  geom_point(
    data = filter(plot_data, !is_frontier),
    aes(shape = Model_Type),
    color = "gray60",
    size = 2.5,
    alpha = 0.7
  ) +
  # Optimal frontier models
  geom_point(
    data = filter(plot_data, is_frontier),
    aes(shape = Model_Type),
    color = "#1f77b4",
    size = 4.2
  ) +
  # Non-overlapping labels for frontier models + top candidates
  geom_text_repel(
    data = filter(plot_data, is_frontier | rank(BIC) <= 10),
    aes(label = Model_Name, color = is_frontier),
    size = 3.5,
    fontface = ifelse(filter(plot_data, is_frontier | rank(BIC) <= 10)$is_frontier, "bold", "plain"),
    max.overlaps = 30,
    box.padding = 0.4
  ) +
  scale_color_manual(values = c("TRUE" = "#1f77b4", "FALSE" = "gray30")) +
  theme_classic(base_size = 14) +
  labs(
    title = expression("Model Efficiency Frontier: Model Complexity vs. Marginal " * R^2),
    subtitle = "Dashed line indicates Pareto-optimal models (highest R² for given BIC penalty)",
    x = "BIC (Lower = Better Fit/Parsimony)",
    y = expression("Marginal " * R^2 * " (Fixed Effect Variance Explained)"),
    shape = "Model Structure"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(frontier_plot)


######################
# Best model is k3ms
#######################
spp_m    <- lmer(ED50 ~ Species + (1 | Site), data = df_modeling, REML = FALSE)
k3ms     <- lmer(ED50 ~ Species + K_3 + (1 | Site), data = df_modeling, REML = FALSE)
k3_model <- lmer(ED50 ~ K_3 + (1 | Site), data = df_modeling, REML = FALSE)
k4ms     <- lmer(ED50 ~ Species + K_4 + (1 | Site), data = df_modeling, REML = FALSE)
k2_model <- lmer(ED50 ~ K_2 + (1 | Site), data = df_modeling, REML = FALSE)

# 2. Extract Model Performance (explicitly using estimator = "ML" to eliminate warning)
p_spp <- model_performance(spp_m, estimator = "ML")
p_k3ms <- model_performance(k3ms, estimator = "ML")
p_k3  <- model_performance(k3_model, estimator = "ML")
p_k4ms <- model_performance(k4ms, estimator = "ML")
p_k2  <- model_performance(k2_model, estimator = "ML")

# 3. Helper Function to Format Fixed Effects Coefficients
get_coef_str <- function(m) {
  df <- tidy(m, effects = "fixed")
  paste(
    sprintf("%s = %.2f ± %.2f (p = %.3f)", df$term, df$estimate, df$std.error, df$p.value),
    collapse = " | "
  )
}

coef_spp  <- get_coef_str(spp_m)
coef_k3ms <- get_coef_str(k3ms)
coef_k3   <- get_coef_str(k3_model)
coef_k4ms <- get_coef_str(k4ms)
coef_k2   <- get_coef_str(k2_model)

# 4. Helper Function for Variance Component Percentages
get_var_pct <- function(m) {
  vc <- as.data.frame(VarCorr(m))
  site_var  <- ifelse("Site" %in% vc$grp, vc$vcov[vc$grp == "Site"], 0)
  res_var   <- vc$vcov[vc$grp == "Residual"]
  total_var <- site_var + res_var
  c(site = (site_var / total_var) * 100, res = (res_var / total_var) * 100)
}

v_spp  <- get_var_pct(spp_m)
v_k3ms <- get_var_pct(k3ms)
v_k3   <- get_var_pct(k3_model)
v_k4ms <- get_var_pct(k4ms)
v_k2   <- get_var_pct(k2_model)

# 5. Build Complete Table 2 with Coefficients
table_2_full <- tibble(
  `Model Name` = c("spp_m (Baseline)", "k3ms", "k3_model", "k4ms", "k2_model"),
  `Fixed Predictors` = c("Species Only", "Species + K3 Cluster", "K3 Cluster Only", "Species + K4 Cluster", "K2 Cluster Only"),
  `AIC` = round(c(p_spp$AIC, p_k3ms$AIC, p_k3$AIC, p_k4ms$AIC, p_k2$AIC), 1),
  `BIC` = round(c(p_spp$BIC, p_k3ms$BIC, p_k3$BIC, p_k4ms$BIC, p_k2$BIC), 1),
  `Marginal R²` = round(c(p_spp$R2_marginal, p_k3ms$R2_marginal, p_k3$R2_marginal, p_k4ms$R2_marginal, p_k2$R2_marginal), 3),
  `Conditional R²` = round(c(p_spp$R2_conditional, p_k3ms$R2_conditional, p_k3$R2_conditional, p_k4ms$R2_conditional, p_k2$R2_conditional), 3),
  `Site Variance (%)` = round(c(v_spp["site"], v_k3ms["site"], v_k3["site"], v_k4ms["site"], v_k2["site"]), 1),
  `Residual Variance (%)` = round(c(v_spp["res"], v_k3ms["res"], v_k3["res"], v_k4ms["res"], v_k2["res"]), 1),
  #`Fixed Effect Coefficients (Estimate ± SE, p-val)` = c(coef_spp, coef_k3ms, coef_k3, coef_k4ms, coef_k2)
)

print(table_2_full)

# Save to CSV for manuscript inclusion
write.csv(table_2_full, here::here("Outputs", "Table2_Top_Models_With_Coefficients.csv"), row.names = FALSE)

# Combine unnested fixed-effect tables across models
top_models_coef_long <- bind_rows(
  tidy(spp_m, effects = "fixed") %>% mutate(Model = "spp_m"),
  tidy(k3ms, effects = "fixed") %>% mutate(Model = "k3ms"),
  tidy(k3_model, effects = "fixed") %>% mutate(Model = "k3_model"),
  tidy(k4ms, effects = "fixed") %>% mutate(Model = "k4ms"),
  tidy(k2_model, effects = "fixed") %>% mutate(Model = "k2_model")
) %>%
  transmute(
    Model = Model,
    Term = term,
    Estimate = round(estimate, 3),
    `Std. Error` = round(std.error, 3),
    `t-statistic` = round(statistic, 2),
    `p-value` = round(p.value, 4)
  )

print(top_models_coef_long)

write.csv(top_models_coef_long, here::here("Outputs", "Supplementary_Top_Models_Coefficients.csv"), row.names = FALSE)
