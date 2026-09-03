library(lme4)
library(lmerTest)
library(tidyverse)
library(performance)
library(ggrepel)
library(here)
library(broom.mixed)

custom_shapes <- c(23, 22, 15, 16, 17, 18, 21, 24)
species_colors <- c(
  "AABR" = "#046E8F",
  "AGLO" = "#83BA75FF",
  "AHYA" = "#D44D5C",
  "ICRA" = "#462255"
)
species_order <- c("ICRA", "AGLO", "AABR", "AHYA")
species_labels <- c(
  "AGLO" = "A. globiceps",
  "AABR" = "A. abrotanoides",
  "AHYA" = "A. hyacinthus",
  "ICRA" = "I. crateriformis"
)

roman_map <- c(
  "Cluster_2" = "Cluster I",
  "Cluster_1" = "Cluster II",
  "Cluster_3" = "Cluster III"
)
target_levels <- c("Cluster I", "Cluster II", "Cluster III")

# Load metadata
meta_filtered <- readRDS(here("Outputs", "master_cluster_metadata_K2_K15.rds"))

df_modeling <- meta_filtered %>%
  filter(!is.na(ED50), !is.na(Site), !is.na(Species)) %>%
  mutate(
    YY      = ED50,
    Species = droplevels(as.factor(Species)),
    Site    = droplevels(as.factor(Site)),
    across(starts_with("K_"), ~ droplevels(as.factor(.)))
  )

# ==============================================================================
# STEP 1: Baseline Species Model (Does ED50 differ by Species?)
# ==============================================================================
spp_m_reml <- lmer(YY ~ Species + (1 | Site), data = df_modeling, REML = TRUE)

cat("\n--- STEP 1: SPECIES BASELINE MODEL (REML = TRUE) ---\n")
print(anova(spp_m_reml))

# ==============================================================================
# STEP 2: Table 1 - Compare Additive Cluster Models (k2ms to k15ms) to find Best K
# ==============================================================================
k_ms_models <- list()
k_range <- 2:15

for (k in k_range) {
  k_col <- paste0("K_", k)
  if (k_col %in% names(df_modeling)) {
    f_kms <- as.formula(paste0("YY ~ Species + ", k_col, " + (1 | Site)"))
    
    k_ms_models[[paste0("k", k, "ms")]] <- tryCatch(
      lmer(f_kms, data = df_modeling, REML = FALSE),
      error = function(e) NULL
    )
  }
}

k_ms_models <- compact(k_ms_models)

# Build Table 1 (Evaluating optimal K)
table_1_kms <- imap_dfr(k_ms_models, function(model_obj, model_name) {
  perf    <- model_performance(model_obj, estimator = "ML")
  r2_vals <- tryCatch(r2(model_obj), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
  k_val   <- str_extract(model_name, "\\d+")
  
  data.frame(
    Model_Name      = model_name,
    K               = k_val,
    Model_Type      = "Species + K (Additive)",
    AIC             = round(perf$AIC, 1),
    BIC             = round(perf$BIC, 1),
    Marginal_R2     = round(r2_vals$R2_marginal, 3),
    Conditional_R2  = round(r2_vals$R2_conditional, 3),
    stringsAsFactors = FALSE
  )
}) %>% arrange(BIC)

cat("\n--- STEP 2: TABLE 1 - ADDITIVE MODEL SELECTION (BEST K) ---\n")
print(head(table_1_kms))
write.csv(table_1_kms, here::here("Outputs", "Table1_K_Additive_Selection.csv"), row.names = FALSE)

# ==============================================================================
# STEP 3: Table 2 - Two-Model Head-to-Head Comparison (k3ms vs spp_m)
# ==============================================================================
# 1. Fit models via ML for model comparison
spp_m_ml <- lmer(YY ~ Species + (1 | Site), data = df_modeling, REML = FALSE)
k3ms_ml  <- lmer(YY ~ Species + K_3 + (1 | Site), data = df_modeling, REML = FALSE)

# 2. Extract Type II ANOVA F-statistics and p-values
an_spp  <- anova(spp_m_ml, type = 2)
an_k3ms <- anova(k3ms_ml, type = 2)

spp_term_spp  <- sprintf("F(%.0f, %.1f) = %.2f, p = %.4f", 
                         an_spp["Species", "NumDF"], an_spp["Species", "DenDF"], 
                         an_spp["Species", "F value"], an_spp["Species", "Pr(>F)"])

spp_term_k3ms <- sprintf("F(%.0f, %.1f) = %.2f, p = %.4f", 
                         an_k3ms["Species", "NumDF"], an_k3ms["Species", "DenDF"], 
                         an_k3ms["Species", "F value"], an_k3ms["Species", "Pr(>F)"])

k3_term_k3ms  <- sprintf("F(%.0f, %.1f) = %.2f, p = %.4f", 
                         an_k3ms["K_3", "NumDF"], an_k3ms["K_3", "DenDF"], 
                         an_k3ms["K_3", "F value"], an_k3ms["K_3", "Pr(>F)"])

# 3. Performance & Variance Metrics
p_spp  <- model_performance(spp_m_ml, estimator = "ML")
p_k3ms <- model_performance(k3ms_ml, estimator = "ML")

get_var_pct <- function(m) {
  vc <- as.data.frame(VarCorr(m))
  site_var  <- ifelse("Site" %in% vc$grp, vc$vcov[vc$grp == "Site"], 0)
  res_var   <- vc$vcov[vc$grp == "Residual"]
  total_var <- site_var + res_var
  c(site = (site_var / total_var) * 100, res = (res_var / total_var) * 100)
}

v_spp  <- get_var_pct(spp_m_ml)
v_k3ms <- get_var_pct(k3ms_ml)

lrt_result <- anova(spp_m_ml, k3ms_ml)

# 4. Build Complete Table 2
table_2_comparison <- tibble(
  `Model Performance Metric` = c(
    "Fixed Predictors",
    "AIC",
    "BIC",
    "Marginal R² (Fixed Effects)",
    "Conditional R² (Fixed + Site)",
    "Site Variance (%)",
    "Residual Error Variance (%)",
    "Species Effect (F, p-value)",
    "K3 Cluster Effect (F, p-value)",
    "LRT vs. Species Only"
  ),
  `Species Only (spp_m)` = c(
    "Species",
    sprintf("%.1f", p_spp$AIC),
    sprintf("%.1f", p_spp$BIC),
    sprintf("%.3f (%.1f%%)", p_spp$R2_marginal, p_spp$R2_marginal * 100),
    sprintf("%.3f (%.1f%%)", p_spp$R2_conditional, p_spp$R2_conditional * 100),
    sprintf("%.1f%%", v_spp["site"]),
    sprintf("%.1f%%", v_spp["res"]),
    spp_term_spp,
    "—",
    "Reference Model"
  ),
  `Species + K3 Phenotype (k3ms)` = c(
    "Species + K3 Cluster",
    sprintf("%.1f", p_k3ms$AIC),
    sprintf("%.1f", p_k3ms$BIC),
    sprintf("%.3f (%.1f%%)", p_k3ms$R2_marginal, p_k3ms$R2_marginal * 100),
    sprintf("%.3f (%.1f%%)", p_k3ms$R2_conditional, p_k3ms$R2_conditional * 100),
    sprintf("%.1f%%", v_k3ms["site"]),
    sprintf("%.1f%%", v_k3ms["res"]),
    spp_term_k3ms,
    k3_term_k3ms,
    sprintf("χ² = %.2f, p = %.4f", lrt_result$Chisq[2], lrt_result$`Pr(>Chisq)`[2])
  )
)

print(table_2_comparison)
write.csv(table_2_comparison, here::here("Outputs", "Table2_spp_m_vs_k3ms_Full.csv"), row.names = FALSE)
# ==============================================================================
# FINAL: Refit Winning Model (k3ms) via REML = TRUE for Final Coefficients
# ==============================================================================
k3ms_reml <- lmer(YY ~ Species + K_3 + (1 | Site), data = df_modeling, REML = TRUE)

# Overall ANOVA (Type II Satterthwaite)
cat("--- Type II ANOVA for k3ms (REML = TRUE) ---\n")
print(anova(k3ms_reml, type = 2))


supp_coefs <- tidy(k3ms_reml, effects = "fixed", conf.int = TRUE) %>%
  transmute(
    Term = term,
    Estimate = round(estimate, 3),
    `Std. Error` = round(std.error, 3),
    `2.5% CI` = round(conf.low, 3),
    `97.5% CI` = round(conf.high, 3),
    `p-value` = round(p.value, 4)
  )

write.csv(supp_coefs, here::here("Outputs", "Supplementary_k3ms_REML_Coefficients.csv"), row.names = FALSE)

#plot
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(tidyverse)
library(ggplot2)

# ==============================================================================
# Diagnostics
# ==============================================================================
# --- Residual Diagnostics ---
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(resid(k3ms_reml) ~ fitted(k3ms_reml), main = "Residuals vs Fitted (k3ms)", xlab = "Fitted", ylab = "Residuals")
qqnorm(resid(k3ms_reml), main = "Normal Q-Q Plot (k3ms)")
qqline(resid(k3ms_reml), col = "red")
par(mfrow = c(1, 1))

# ==============================================================================
# 2. Extract EMMs & Compact Letter Display (CLD)
# ==============================================================================
k3_emmeans_obj <- emmeans(k3ms_reml, pairwise ~ K_3, adjust = "tukey")
k3_emm_df       <- as.data.frame(k3_emmeans_obj$emmeans)
k3_contrasts_df <- as.data.frame(k3_emmeans_obj$contrasts)

cat("\n--- Estimated Marginal Means (ED50) for K_3 Clusters ---\n")
print(k3_emm_df %>% dplyr::select(K_3, emmean, SE, df, lower.CL, upper.CL))

cat("\n--- Pairwise Tukey Contrasts Between K_3 Clusters ---\n")
print(k3_contrasts_df %>% dplyr::select(contrast, estimate, SE, p.value))

# Compute Tukey CLD Letters
cld_clusters   <- cld(k3_emmeans_obj$emmeans, adjust = "tukey", Letters = letters) %>%
  as.data.frame() %>%
  mutate(.group = trimws(.group))

# # ==============================================================================
# # 3. Align Cluster Labels & Factor Ordering
# # ==============================================================================
df_modeling <- df_modeling %>%
  mutate(
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels)
  )

cld_clusters <- cld_clusters %>%
  mutate(
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels)
  )

# Verify counts in console
cat("Raw Data Counts per Cluster Label:\n")
print(table(df_modeling$Cluster_Label, useNA = "always"))

cat("\nEMM Summary Counts per Cluster Label:\n")
print(table(cld_clusters$Cluster_Label, useNA = "always"))
# ==============================================================================
# 4. Publication-Ready Plot
# ==============================================================================
y_min_val <- min(df_modeling$ED50, na.rm = TRUE) - 0.5
y_max_val <- max(cld_clusters$upper.CL, na.rm = TRUE) + 0.8

plot_k3ms <- ggplot() +
  # Layer A: Raw data points
  geom_jitter(
    data = df_modeling, 
    aes(x = Cluster_Label, y = ED50, color = Species, shape = Site),
    width = 0.12, 
    alpha = 0.5, 
    size = 3
  ) +
  # Layer B: Model Estimated Marginal Means & 95% CIs
  geom_errorbar(
    data = cld_clusters,
    aes(x = Cluster_Label, ymin = lower.CL, ymax = upper.CL),
    width = 0.15, 
    linewidth = 0.8, 
    color = "black"
  ) +
  geom_point(
    data = cld_clusters,
    aes(x = Cluster_Label, y = emmean),
    size = 5,
    fill = "black",
    color = "black"
  ) +
  # Layer C: Significance Letters (CLD)
  geom_text(
    data = cld_clusters,
    aes(x = Cluster_Label, y = upper.CL + 0.5, label = .group),
    fontface = "bold", 
    hjust = 0.5,
    size = 6,
    color = "black",
    inherit.aes = FALSE
  ) +
  # Custom Aesthetics & Formatting
  scale_color_manual(values = species_colors, limits = species_order, labels = species_labels) +
  scale_shape_manual(values = custom_shapes) +
  coord_cartesian(ylim = c(y_min_val, y_max_val)) +
  labs(
    x = "Photophysiological Cluster (K = 3)",
    y = expression("Estimated Marginal Mean ED"[50] ~ "(°C)"),
    color = "Species",
    shape = "Site"
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title      = element_text(face = "bold", size = 13, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    axis.text.x     = element_text(face = "bold", size = 11, color = "black"),
    axis.line       = element_line(linewidth = 0.6, color = "black"),
    axis.ticks      = element_line(linewidth = 0.6, color = "black"),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 11, face = "italic"),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
  )


print(plot_k3ms)
# Save High-Res PDF/PNG Exports
ggsave(here::here("Plots", "ED50_K3_emmeans_clusters.pdf"), plot = plot_k3ms, width = 8, height = 6, device = cairo_pdf)
ggsave(here::here("Plots", "ED50_K3_emmeans_clusters.png"), plot = plot_k3ms, width = 8, height = 6, dpi = 600)
