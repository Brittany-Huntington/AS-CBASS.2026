#######################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)
library(lme4)
library(performance)
library(lmerTest)
library(emmeans)
library(multcomp)
library(patchwork)

# 1. LOAD AND CLEAN DATA
master_dataset_fvfm <- read.csv(here("Data/fvfm_color_coral_dataset_merged.csv"))

cluster <- read.csv(here("Data/cluster_and_ed50_data.csv")) %>%
  mutate(SampleID_clean = as.character(SampleID_clean)) %>%
  mutate(SampleID_clean = toupper(trimws(as.character(SampleID_clean)))) %>%
  mutate(SampleID_clean = sub("^T", "", trimws(as.character(SampleID_clean)))) %>%
  mutate(SampleID_clean = sub("([A-Za-z]+)(\\d+)$", "\\1_\\2", SampleID_clean)) %>%
  mutate(K_3 = factor(K_3)) %>%
  dplyr::select(SampleID_clean, K_3) %>%
  distinct(SampleID_clean, .keep_all = TRUE)

plot_data_clean <- master_dataset_fvfm %>%
  filter(!is.na(Pam_value), !is.na(Mean_inverted)) %>%
  mutate(
    Pam_value      = as.numeric(Pam_value),
    Color_mean     = as.numeric(Color_mean),
    Mean_inverted  = as.numeric(Mean_inverted),
    Species        = as.factor(Species),
    Site           = as.factor(Site),
    SampleID_clean = as.character(SampleID_clean)
  )

plot_data_clean <- left_join(plot_data_clean, cluster, by = "SampleID_clean")
plot_data_clean$Rel_Temperature <- plot_data_clean$Temperature - 29

# 2. SPECIES SUBSETS AND EXPLORATORY PLOTS
ICRA <- plot_data_clean %>% filter(Species == "ICRA")
ICRA <- ICRA %>% filter(!(Site %in% c("11", 11) & Mean_inverted > 6))

AGLO <- plot_data_clean %>% filter(Species == "AGLO")
AABR <- plot_data_clean %>% filter(Species == "AABR")
AHYA <- plot_data_clean %>% filter(Species == "AHYA")

ggplot(AHYA, aes(x = Color_mean, y = Pam_value, color = factor(Site), fill = factor(Site))) +
  geom_point(alpha = 0.25, size = 1.8, stroke = 0) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, size = 1.2) +
  facet_wrap(~ Site, labeller = label_both) +
  labs(
    x = "Calibrated Color Score",
    y = expression(F[v]/F[m])
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

plot_data_clean %>% 
  group_by(Species, Temperature) %>%
  summarize(
    N   = length(Color_mean),
    MnC = mean(Color_mean),
    MnP = mean(Pam_value),
    sdC = sd(Color_mean)/sqrt(N),
    sdP = sd(Pam_value)/sqrt(N)
  ) %>%
  ggplot(aes(x = MnC, y = MnP, color = Temperature)) +
  geom_point() +
  geom_errorbar(aes(ymin = MnP - sdP, ymax = MnP + sdP)) +
  geom_errorbarh(aes(xmin = MnC - sdC, xmax = MnC + sdC)) +
  facet_grid(Species ~ .) +
  stat_smooth(method = "lm", color = "gold")

plot_data_clean %>% 
  ggplot(aes(x = Color_mean, y = Pam_value, color = Temperature)) +
  geom_point() +
  facet_grid(Temperature ~ Species)

plot_data_clean %>% 
  ggplot(aes(x = Temperature, y = Color_mean/10, color = factor(Genotype), group = Genotype)) +
  geom_point() +
  stat_smooth(method = "lm", span = .25, se = FALSE) +
  facet_grid(Species ~ Site)

plot_data_clean %>% 
  ggplot(aes(x = Rel_Temperature, y = Pam_value, color = factor(Genotype), group = Genotype)) +
  geom_point() +
  stat_smooth(method = "lm", span = .75, se = FALSE, alpha = .1) +
  facet_grid(Species ~ Site)

# 3. INDIVIDUAL LOOP MODELLING (EXTRACTING FRAGMENT SLOPES)
uG <- unique(plot_data_clean$SampleID_clean)

lmresults <- NULL
for (ig in 1:length(uG)) {
  thisgenotype <- subset(plot_data_clean, SampleID_clean == uG[ig])
  thismod      <- lm(Color_mean ~ Rel_Temperature, data = thisgenotype)
  stm          <- summary(thismod)
  thismodresults <- data.frame(
    SampleID = uG[ig],
    Species  = thisgenotype$Species[1],
    Site     = thisgenotype$Site[1],
    int      = thismod$coefficients[1],
    slope    = thismod$coefficients[2],
    r2       = stm$r.squared,
    p        = stm$coefficients[2, 4]
  )
  lmresults <- rbind(lmresults, thismodresults)
}

lmresults %>%
  ggplot(aes(x = Species, y = slope)) +
  geom_boxplot()

ggmod <- lmer(slope ~ Species + (1 | Site), data = lmresults)
summary(ggmod)
emmeans(ggmod, specs = "Species")

# Repeat loop including K_3 factor for lmresults1
lmresults1 <- NULL
for (ig in 1:length(uG)) {
  thisgenotype <- subset(plot_data_clean, SampleID_clean == uG[ig])
  if (nrow(thisgenotype) >= 3) {
    thismod <- lm(Color_mean ~ Rel_Temperature, data = thisgenotype)
    stm     <- summary(thismod)
    thismodresults <- data.frame(
      SampleID = uG[ig],
      Species  = thisgenotype$Species[1],
      Site     = thisgenotype$Site[1],
      K_3      = as.factor(thisgenotype$K_3[1]),
      int      = thismod$coefficients[1],
      slope    = thismod$coefficients[2],
      r2       = stm$r.squared,
      p        = stm$coefficients[2, 4]
    )
    lmresults1 <- rbind(lmresults1, thismodresults)
  }
}

ggmod1 <- lmer(slope ~ K_3 * Species + (1 | Site), data = lmresults1)
summary(ggmod1)
anova(ggmod1)
performance(ggmod1)

em_k3_1 <- emmeans(ggmod1, specs = ~ K_3 | Species)
em_k3_1
pairs(em_k3_1)

em_df1 <- as.data.frame(em_k3_1) %>% filter(!is.na(SE))

ggplot(em_df1, aes(x = K_3, y = emmean, color = K_3)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, linewidth = 0.8) +
  facet_wrap(~ Species, scales = "free_x") +
  scale_color_brewer(palette = "Set1", name = "Cluster (K_3)") +
  labs(
    x = "Photophysiological Cluster (K_3)",
    y = expression("Model Estimated Slope (" * Delta * "Color / " * Delta * "Temp)")
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

lmresults1_clean <- lmresults1 %>% filter(!is.na(K_3))

ggplot(lmresults1_clean, aes(x = K_3, y = slope)) +
  geom_violin(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2, aes(shape = Site)) +
  facet_grid(~ Species, scales = "free_x") +
  scale_fill_brewer(palette = "Set1", name = "Cluster (K_3)") +
  scale_color_brewer(palette = "Set1", name = "Cluster (K_3)") +
  labs(
    x = "Photophysiological Cluster (K_3)",
    y = expression("Extracted Color Slope (" * Delta * "Color / " * Delta * "Temp)")
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

# 4. PRIMARY MIXED-EFFECTS MODELS & STATISTICAL TESTING
linmod <- lm(Color_mean ~ Rel_Temperature * SampleID_clean, data = plot_data_clean)

lmmmmod1  <- lmer(Color_mean ~ Rel_Temperature + (1 + Rel_Temperature | Site/SampleID_clean), data = plot_data_clean)
lmmmodi2  <- lmer(Color_mean ~ Rel_Temperature + K_3 + (1 | Site) + (1 + Rel_Temperature | SampleID_clean), data = plot_data_clean)
lmmmodi2b <- lmer(
  Color_mean ~ Rel_Temperature + K_3 + (1 | Site) + (1 + Rel_Temperature | SampleID_clean), 
  data = plot_data_clean,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
lmmmodi3a <- lmer(Color_mean ~ Rel_Temperature * K_3 + (1 | Site) + (1 + Rel_Temperature | SampleID_clean), data = plot_data_clean)
lmmmodi3b <- lmer(
  Color_mean ~ Rel_Temperature * K_3 + (1 | Site) + (1 + Rel_Temperature | SampleID_clean), 
  data = plot_data_clean,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

lmmmodi3_clean <- lmer(Color_mean ~ Rel_Temperature * K_3 + (1 | Site) + (1 | SampleID_clean), data = plot_data_clean)
lmmmod3_clean  <- lmer(Color_mean ~ Rel_Temperature + K_3 + (1 | Site) + (1 | SampleID_clean), data = plot_data_clean)

lmmmodi_species <- lmer(
  Color_mean ~ Rel_Temperature * K_3 * Species + (1 | Site) + (1 | SampleID_clean),
  data = plot_data_clean
)

lmmmodi_final <- lmer(
  Color_mean ~ Rel_Temperature * K_3 + Rel_Temperature * Species + (1 | Site) + (1 | SampleID_clean),
  data = plot_data_clean
)

anova(lmmmodi_final, type = "III")
anova(lmmmodi_species, type = "III")
summary(lmmmodi_final)
performance(lmmmodi3_clean)
performance(lmmmodi_species)
anova(lmmmodi3_clean)

slope_emms <- emtrends(lmmmodi_final, specs = "K_3", var = "Rel_Temperature")
slope_emms

slope_emms_byspecies <- emtrends(lmmmodi_final, specs = ~ K_3 | Species, var = "Rel_Temperature")
slope_emms_byspecies
pairs(slope_emms_byspecies)

# 5. FIGURES & VISUALIZATIONS
em_species_df <- as.data.frame(slope_emms_byspecies)

ggplot() +
  geom_jitter(data = lmresults1_clean, aes(x = K_3, y = slope, color = K_3), width = 0.15, alpha = 0.4) +
  geom_point(data = em_species_df, aes(x = K_3, y = Rel_Temperature.trend), size = 3.5) +
  geom_errorbar(data = em_species_df, aes(x = K_3, ymin = lower.CL, ymax = upper.CL), width = 0.15) +
  facet_wrap(~ Species, nrow = 1)

slope_emms_global <- emtrends(lmmmodi_final, specs = "K_3", var = "Rel_Temperature")
em_global_df <- as.data.frame(slope_emms_global)

lmresults1_clean <- lmresults1 %>% 
  filter(!is.na(K_3)) %>% 
  mutate(
    K_3     = factor(K_3),
    Species = factor(Species),
    Site    = factor(Site)
  )

p_k3_global <- ggplot() +
  geom_jitter(
    data = lmresults1_clean, 
    aes(x = K_3, y = slope, color = Species, shape = Site),
    width = 0.22, 
    alpha = 0.65, 
    size = 2.8
  ) +
  geom_errorbar(
    data = em_global_df, 
    aes(x = K_3, ymin = lower.CL, ymax = upper.CL),
    width = 0.12, 
    linewidth = 1.1, 
    color = "black"
  ) +
  geom_point(
    data = em_global_df, 
    aes(x = K_3, y = Rel_Temperature.trend),
    size = 4.5, 
    shape = 18, 
    color = "black"
  ) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8, 25), name = "Reef Site") +
  labs(
    x = "Photophysiological Cluster (K_3)",
    y = expression("Bleaching Slope (" * Delta * "Color / " * Delta * "°C)"),
    title = "Thermal Bleaching Sensitivity by Photophysiological Cluster",
    subtitle = "Points = raw fragment slopes (colored by species, shaped by site); Black diamonds = lmer model mean ± 95% CI"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.box       = "horizontal",
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3.5, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 3.5, alpha = 1))
  )

print(p_k3_global)

em_trends_df <- as.data.frame(slope_emms)

ggplot(em_trends_df, aes(x = K_3, y = Rel_Temperature.trend, color = K_3)) +
  geom_point(size = 4) +
  geom_jitter() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, linewidth = 1) +
  labs(
    x = "Photophysiological Cluster (K_3)",
    y = expression("Bleaching Rate (" * Delta * "Color Score / " * Delta * "°C)")
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

lmresults_clean <- lmresults1 %>% 
  mutate(K_3 = factor(K_3))

ggplot() +
  geom_jitter(
    data = lmresults_clean, 
    aes(x = K_3, y = slope, color = K_3),
    width = 0.15, 
    alpha = 0.35, 
    size = 2.2
  ) +
  geom_point(
    data = em_trends_df, 
    aes(x = K_3, y = Rel_Temperature.trend),
    size = 4.5, 
    color = "black"
  ) +
  geom_errorbar(
    data = em_trends_df, 
    aes(x = K_3, ymin = lower.CL, ymax = upper.CL),
    width = 0.12, 
    linewidth = 1.1, 
    color = "black"
  ) +
  scale_color_manual(
    values = c("1" = "#D95F02", "2" = "#1B9E77", "3" = "#E6AB02"),
    name = "Photophys Cluster"
  ) +
  labs(
    x = "Photophysiological Cluster (K_3)",
    y = expression("Bleaching Rate (" * Delta * "Color Score / " * Delta * "°C)"),
    title = "Thermal Bleaching Rates by Photophysiological Cluster",
    subtitle = "Faint points = raw fragment slopes; Black markers = lmer estimated mean ± 95% CI"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# Normalized Retention Loops
species_list <- unique(plot_data_clean$Species)
species_plots <- list()

for (spp in species_list) {
  plot_data_spp <- plot_data_clean %>%
    filter(
      Species == spp,
      !is.na(Temperature), 
      !is.na(Pam_value), 
      !is.na(Color_mean)
    ) %>%
    mutate(
      Temperature = as.numeric(Temperature),
      Pam_value   = as.numeric(Pam_value),
      Color_mean  = as.numeric(Color_mean)
    ) %>%
    group_by(Site) %>% 
    mutate(
      pam_min   = min(Pam_value, na.rm = TRUE),
      pam_range = max(Pam_value, na.rm = TRUE) - pam_min,
      col_min   = min(Color_mean, na.rm = TRUE),
      col_range = max(Color_mean, na.rm = TRUE) - col_min,
      `Fv/Fm` = if_else(pam_range > 0, (Pam_value - pam_min) / pam_range, 1),
      `Color` = if_else(col_range > 0, (Color_mean - col_min) / col_range, 1)
    ) %>%
    dplyr::select(-pam_min, -pam_range, -col_min, -col_range) %>%
    ungroup() %>%
    pivot_longer(
      cols = c(`Fv/Fm`, `Color`),
      names_to = "Indicator",
      values_to = "Retention"
    )
  
  if (nrow(plot_data_spp) == 0) next
  
  p <- ggplot(plot_data_spp, aes(x = Temperature, y = Retention, color = Indicator, linetype = Indicator)) +
    geom_point(alpha = 0.3, size = 1.8, position = position_jitter(width = 0.2)) +
    geom_smooth(method = "loess", span = 0.85, se = TRUE, alpha = 0.15, size = 1.2) +
    facet_wrap(~ Site, labeller = label_both) +
    scale_color_manual(values = c("Fv/Fm" = "#D95F02", "Color" = "#1B9E77")) +
    scale_linetype_manual(values = c("Fv/Fm" = "solid", "Color" = "dashed")) +
    scale_x_continuous(breaks = seq(26, 40, by = 2)) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(
      x = "Temperature (°C)",
      y = "Relative Retention",
      #title = paste0("Thermal Sensitivity Divergence: ", spp),
      #subtitle = "Fv/Fm vs. Tissue Color Loss across Temperature Steps"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank()
    )
  
  print(p)
  species_plots[[spp]] <- p
}

# 6. COMPOSITE PUBLICATION FIGURE WITH SIGNIFICANCE LETTERS
slope_cld <- cld(slope_emms_global, Letters = letters, adjust = "tukey") %>% 
  as.data.frame() %>% 
  mutate(.group = trimws(.group))

p_color_decay <- ggplot(
  plot_data_clean, 
  aes(x = Rel_Temperature, y = Color_mean, color = factor(K_3), fill = factor(K_3))
) +
  geom_point(alpha = 0.2, size = 1.2, position = position_jitter(width = 0.1)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, alpha = 0.2) +
  scale_color_manual(
    values = c("1" = "#D95F02", "2" = "#1B9E77", "3" = "#E6AB02"),
    name = "Photophys Cluster (K_3)"
  ) +
  scale_fill_manual(
    values = c("1" = "#D95F02", "2" = "#1B9E77", "3" = "#E6AB02"),
    name = "Photophys Cluster (K_3)"
  ) +
  labs(
    x = "Relative Temperature (°C)",
    y = "Color Score (Mean)",
   # title = "Thermal Color Loss Trajectories"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

p_slopes <- ggplot() +
  geom_jitter(
    data = lmresults1_clean, 
    aes(x = K_3, y = slope, color = Species, shape = Site),
    width = 0.2, 
    alpha = 0.6, 
    size = 2.5
  ) +
  geom_errorbar(
    data = slope_cld, 
    aes(x = K_3, ymin = lower.CL, ymax = upper.CL),
    width = 0.12, 
    linewidth = 1.1, 
    color = "black"
  ) +
  geom_point(
    data = slope_cld, 
    aes(x = K_3, y = Rel_Temperature.trend),
    size = 4.0, 
    shape = 18, 
    color = "black"
  ) +
  geom_text(
    data = slope_cld,
    aes(x = K_3, y = upper.CL + 0.015, label = .group),
    fontface = "bold",
    size = 5,
    color = "black"
  ) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8, 25), name = "Site") +
  labs(
    x = "Photophysiological Cluster (K_3)",
    y = expression("Bleaching Slope (" * Delta * "Color / " * Delta * "°C)"),
   # title = "Extracted Thermal Sensitivity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "bottom",
    legend.box       = "vertical",
    panel.grid.minor = element_blank()
  )

figure_color_complete <- p_color_decay + p_slopes + 
  plot_annotation(tag_levels = 'A')

print(figure_color_complete)

# Raw Data Dose-Response Visualization
ggplot(plot_data_clean, aes(x = Temperature, y = Color_mean, color = Genotype)) +
  geom_jitter(width = 0.3, size = 2.5, alpha = 0.75) +
  geom_smooth(method = "loess", se = FALSE, size = 0.8) +
  facet_grid(Species ~ Site, labeller = label_both) +
  scale_color_viridis_d(option = "turbo", name = "Genotype") +
  scale_x_continuous(breaks = seq(26, 40, by = 2)) +
  labs(
    x = "Temperature (°C)",
    y = "Calibrated Mean Color Score",
    title = "Raw CBASS Data: Temperature vs. Coral Color Score",
    subtitle = "Raw measurements jittered by Genotype; trendlines fitted via LOESS"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

# Fv/Fm vs. Color Score Visualization
plot_data <- master_dataset_fvfm %>%
  filter(!is.na(Pam_value), !is.na(Color_mean)) %>%
  mutate(
    Pam_value  = as.numeric(Pam_value),
    Color_mean = as.numeric(Color_mean),
    Genotype   = as.factor(Genotype),
    Species    = as.factor(Species),
    Site       = as.factor(Site)
  )

ggplot(plot_data, aes(x = Color_mean, y = Pam_value, color = Genotype)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 0.7) +
  facet_grid(Species ~ Site, labeller = label_both) +
  scale_color_viridis_d(option = "turbo", name = "Genotype") +
  labs(
    x = "Calibrated Color Score",
    y = expression(F[v]/F[m] ~ "(Photochemical Efficiency)"),
    title = expression(F[v]/F[m] ~ "vs. Tissue Color Score")
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

# Multi-metric Long Format Visualization
data_long <- master_dataset_fvfm %>%
  filter(!is.na(Temperature)) %>%
  mutate(
    `Fv/Fm`       = as.numeric(Pam_value),
    `Color Score` = as.numeric(Mean_inverted),
    Temperature   = as.numeric(Temperature),
    Species       = as.factor(Species),
    Site          = as.factor(Site)
  ) %>%
  tidyr::pivot_longer(
    cols = c(`Fv/Fm`, `Color Score`),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

ggplot(data_long, aes(x = Temperature, y = Value, color = Species, fill = Species)) +
  geom_point(alpha = 0.25, size = 1.8, position = position_jitter(width = 0.2)) +
  geom_smooth(method = "loess", span = 0.85, se = TRUE, alpha = 0.15, size = 1.1) +
  facet_grid(Metric ~ Site, scales = "free_y", switch = "y") +
  scale_color_brewer(palette = "Set1", name = "Species") +
  scale_fill_brewer(palette = "Set1", name = "Species") +
  scale_x_continuous(breaks = seq(26, 40, by = 2)) +
  labs(
    x = "Temperature (°C)",
    y = NULL,
    title = "Thermal Stress Thresholds: Photochemical Efficiency vs. Color Bleaching",
    subtitle = "Fitted with smooth non-linear curves (LOESS, 95% CI)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    strip.placement  = "outside",
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )