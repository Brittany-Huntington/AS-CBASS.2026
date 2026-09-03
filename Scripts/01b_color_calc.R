#Code for processing ED50 values

#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(rstudioapi)
library(RColorBrewer)
library(CBASSED50) #Voolstra R package
library(here)
library(readr)
library(purrr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(broom)


mandatory_columns()
rm(list = ls())

custom_shapes <- c(23, 22, 15, 16, 17, 18, 21, 24)
species_colors <- c(
  "AABR" = "#046E8F", # Muted Crimson/Rose Red (#D44D5CFF)
  "AGLO" = "#83BA75FF", # Light Sage Teal (#9DD9D2FF)
  "AHYA" = "#D44D5C", # Deep Ocean Blue (#046E8FFF)
  "ICRA" = "#462255"  # Bright Amber/Orange (#FF8811FF) -- (or "#462255" for Deep Purple)
)
#custom_shapes <- c(23, 22, 15, 16, 17, 18, 21, 24)
#species order
species_order <- c("ICRA", "AGLO", "AABR", "AHYA")

# Custom shapes & high-contrast colors
custom_shapes <- c(16, 15, 17, 18, 8, 11, 9, 10)

species_labels <- c(
  "ICRA" = "I. crateriformis",
  "AGLO" = "A. globiceps",
  "AABR" = "A. abrotanoides",
  "AHYA" = "A. hyacinthus"
)
#####Load and Prep Data-----------------------------

#Load metadata CSV file
metadata_img <- read_csv(here("Data/color/METADATA_color_scores.csv") )%>%
  mutate(
    Image = as.character(Image),
    Site  = as.character(Site)
  )

metadata2 <- read_csv(here("Data/TUT_combined_fvfm.csv") )%>%
  mutate(
    Cookie_no = as.numeric(Cookie_no), # match cookie numeric type
    Site      = as.character(Site)
  ) %>%
  dplyr::select(Cookie_no, Site, Bag_number, Tank_color, SampleID_clean)

fvfm <- read_csv(here("Data/TUT_CBASS_raw_PAM_fvfm.csv"))

fvfm <- fvfm %>%
  mutate(Site_number = as.character(Site_number))%>%
        # SampleID_clean = as.character(SampleID_clean)) %>%
  dplyr::select(-Date)


file_paths <- list.files(
  path = (here("Data/color") ), 
  pattern = "^color_scores.*\\.csv$", 
  full.names = TRUE
)

#get list of all files 
master_dataset <- file_paths %>%
  map_dfr(read_csv) %>%
  mutate(
    Image  = as.character(Image),
    cookie = as.numeric(cookie)
  ) %>%
  filter(Identifier != "scale" | is.na(Identifier))
  

master_dataset<-master_dataset%>%
  left_join(metadata_img, by = "Image")
  

# Join to get bag number
master_dataset_bag <- master_dataset %>%
  # Join metadata2 using Cookie_no and Site
  left_join(metadata2, by = c("cookie" = "Cookie_no", "Site" = "Site", "Tank_color" = "Tank_color"))

#join to get fvfm
fvfm_and_color <- master_dataset_bag %>%
  left_join(fvfm, by = c("SampleID_clean", "Bag_number" = "Bag_number", "Site" = "Site_number", "Temperature" = "Temperature", "Tank_color" = "Tank_color"))%>%
  dplyr::select(-Click_No, -Identifier)%>%
  dplyr::rename("Color_mean" = "Mean", "Color_SD" = "SD" )%>%
  mutate(
    Temperature = as.numeric(Temperature),
    Color_mean        = as.numeric(Color_mean))

fvfm_and_color <- fvfm_and_color %>%
  filter(!is.na(Pam_value), !is.na(Color_mean)) %>%
  mutate(
    Pam_value     = as.numeric(Pam_value),
    Color_mean = as.numeric(Color_mean),
    Species       = as.factor(Species),
    Site          = as.factor(Site),
    SampleID_clean = as.character(SampleID_clean),
    Rel_Temperature=Temperature-29
  )

write_csv(fvfm_and_color, "Outputs/fvfm_color_df.csv")

uG=unique(fvfm_and_color$SampleID_clean)
lmresults=NULL
#ig=1
for(ig in 1:length(uG)){
  thisgenotype=subset(fvfm_and_color,SampleID_clean==uG[ig])
  thismod=lm(Color_mean~Rel_Temperature,data=thisgenotype)
  stm=summary(thismod)
  thismodresults=data.frame(SampleID=uG[ig],Species=thisgenotype$Species[1],
                            Site=thisgenotype$Site[1],
                            int=thismod$coefficients[1],slope=thismod$coefficients[2],r2=stm$r.squared,p=stm$coefficients[2,4])
  lmresults=rbind(lmresults,thismodresults)
  print(ig)
}
write_csv(lmresults, "Outputs/color_slope.csv")

lmresults <- lmresults %>%
  mutate(
    Species = factor(Species, levels = species_order),
    Site = as.factor(Site)
  )

# Fit mixed model and calculate emmeans
m_color <- lmer(slope ~ Species + (1|Site), data = lmresults)
species_emms_color <- emmeans(m_color, ~ Species)
species_letters_color <- cld(species_emms_color, Letters = letters, adjust = "tukey")

# Clean emmeans summary data frame for plotting
color_emm_df <- as.data.frame(species_letters_color) %>%
  dplyr::rename(
    Mean = emmean,
    CI_lower = lower.CL,
    CI_upper = upper.CL,
    .group = .group
  ) %>%
  dplyr::mutate(
    .group = trimws(.group),
    Species = factor(Species, levels = species_order)
  )

# Plotting Color Slope
g_color <- ggplot(lmresults, aes(x = Species, y = slope)) +
  geom_jitter(
    aes(color = Species, shape = as.factor(Site)),
    width = 0.08, alpha = 0.45, size = 2.5, show.legend = FALSE
  ) +
  geom_point(
    data = color_emm_df,
    aes(x = Species, y = Mean, color = Species),
    size = 5, inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_errorbar(
    data = color_emm_df,
    aes(x = Species, ymin = CI_lower, ymax = CI_upper),
    width = 0.1, linewidth = 0.8, color = "black", inherit.aes = FALSE
  ) +
  geom_text(
    data = color_emm_df,
    aes(x = Species, y = CI_upper + 0.05, label = .group),
    hjust = 0.5, size = 6, fontface = "bold", color = "black", inherit.aes = FALSE
  ) +
  labs(
    x = "Species",
    y = expression("Color Slope (" * Delta * "Color / °C)")
  ) +
  scale_color_manual(values = species_colors, limits = species_order) +
  scale_fill_manual(values = species_colors, limits = species_order) + 
  scale_x_discrete(
    labels = species_labels, 
    limits = species_order, 
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_shape_manual(values = custom_shapes) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold", size = 13, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black", angle = 30, hjust = 1, face = "italic"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
  )

g_color

# Emmeans Plot for Color Slope
g_color_emm <- ggplot(lmresults, aes(x = Species, y = slope)) +
  # Raw individual data points jittered by Site shape
  geom_jitter(
    aes(color = Species, shape = as.factor(Site)),
    width = 0.08, 
    alpha = 0.45, 
    size = 2.5, 
    show.legend = FALSE
  ) +
  # Model Estimated Marginal Means (emmean)
  geom_point(
    data = color_emm_df,
    aes(x = Species, y = Mean, color = Species),
    size = 5,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  # 95% Confidence Intervals from lmer model emmeans
  geom_errorbar(
    data = color_emm_df,
    aes(x = Species, ymin = CI_lower, ymax = CI_upper),
    width = 0.1,
    linewidth = 0.8,
    color = "black",
    inherit.aes = FALSE
  ) +
  # Post-hoc significance letters directly from emmeans cld
  geom_text(
    data = color_emm_df,
    aes(x = Species, y = CI_upper + 0.2, label = .group),
    hjust = 0.5,
    size = 6,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  ) +
  # Axis Labels
  labs(
    x = "Species",
    y = expression("Emmeans Color Slope (" * Delta * "Color / °C)")
  ) +
  # Formatting & Aesthetics
  scale_color_manual(values = species_colors, limits = species_order) +
  scale_fill_manual(values = species_colors, limits = species_order) + 
  scale_x_discrete(
    labels = species_labels, 
    limits = species_order, 
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_shape_manual(values = custom_shapes) +
  # Theme Customizations
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold", size = 13, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black", angle = 30, hjust = 1, face = "italic"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
  )

g_color_emm

# 1. Generate predicted regression lines for each sample using lmresults
pred_lines <- fvfm_and_color %>%
  left_join(lmresults %>% dplyr::select(SampleID, int, slope), by = c("SampleID_clean" = "SampleID")) %>%
  mutate(
    pred_color = int + slope * Rel_Temperature,
    Species = factor(Species, levels = species_order)
  )

# 2. Extract species-level average slopes and intercepts from emmeans
m_slope <- lmer(slope ~ Species + (1|Site), data = lmresults)
m_int   <- lmer(int ~ Species + (1|Site), data = lmresults)

avg_slopes <- as.data.frame(emmeans(m_slope, ~ Species)) %>% dplyr::select(Species, avg_slope = emmean)
avg_ints   <- as.data.frame(emmeans(m_int, ~ Species))   %>% dplyr::select(Species, avg_int = emmean)

species_avg_lines <- left_join(avg_slopes, avg_ints, by = "Species") %>%
  mutate(Species = factor(Species, levels = species_order))

# 3. Plot individual slopes + average species slope lines
g_slopes <- ggplot() +
  # Individual sample lines (thin, semi-transparent)
  geom_line(
    data = pred_lines,
    aes(x = Rel_Temperature, y = pred_color, group = SampleID_clean, color = Species),
    alpha = 0.25, linewidth = 0.5
  ) +
  # Species-level average slope lines (bold)
  geom_abline(
    data = species_avg_lines,
    aes(intercept = avg_int, slope = avg_slope, color = Species),
    linewidth = 1.5
  ) +
  labs(
    x = "Relative Temperature (°C)",
    y = "Color Mean",
    color = "Species"
  ) +
  scale_color_manual(values = species_colors, limits = species_order, labels = species_labels) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    legend.text = element_text(face = "italic"),
    legend.position = "right"
  )

g_slopes

# 1. Generate predictions across temperature range with SE/CI from emmeans
temp_seq <- seq(min(fvfm_and_color$Rel_Temperature, na.rm = TRUE),
                max(fvfm_and_color$Rel_Temperature, na.rm = TRUE), 
                length.out = 100)

# Grid of Species x Temperatures
pred_grid <- expand.grid(
  Species = species_order,
  Rel_Temperature = temp_seq
)

# Extract predictions + SEs for slopes and intercepts using emmeans
emm_slope_df <- as.data.frame(emmeans(m_slope, ~ Species)) %>% dplyr::select(Species, slope = emmean, slope_se = SE)
emm_int_df   <- as.data.frame(emmeans(m_int, ~ Species))   %>% dplyr::select(Species, int = emmean, int_se = SE)

ci_ribbon_df <- left_join(pred_grid, emm_slope_df, by = "Species") %>%
  left_join(emm_int_df, by = "Species") %>%
  mutate(
    pred_color = int + slope * Rel_Temperature,
    # Approximate 95% CI bounds for linear response
    se_fit = sqrt(int_se^2 + (Rel_Temperature * slope_se)^2),
    ymin = pred_color - 1.96 * se_fit,
    ymax = pred_color + 1.96 * se_fit,
    Species = factor(Species, levels = species_order)
  )

# 2. Plot with ribbon CI
g_slopes_ribbon <- ggplot(pred_lines, aes(x = Rel_Temperature, y = Color_mean, color = Species)) +
  # Raw data points
  geom_point(alpha = 0.25, size = 1.2) +
  # Individual sample lines (thin/faint)
  geom_line(aes(y = pred_color, group = SampleID_clean), alpha = 0.2, linewidth = 0.4) +
  # 95% CI Ribbon around Species Average
  geom_ribbon(
    data = ci_ribbon_df,
    aes(x = Rel_Temperature, ymin = ymin, ymax = ymax, fill = Species),
    alpha = 0.25, color = NA, inherit.aes = FALSE
  ) +
  # Species Average Line (bold)
  geom_line(
    data = ci_ribbon_df,
    aes(x = Rel_Temperature, y = pred_color, color = Species),
    linewidth = 1.2
  ) +
  facet_wrap(~ Species, labeller = as_labeller(species_labels)) +
  labs(
    x = "Relative Temperature (°C)",
    y = "Color Mean"
  ) +
  scale_color_manual(values = species_colors, limits = species_order) +
  scale_fill_manual(values = species_colors, limits = species_order) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold.italic", size = 11),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black")
  )

g_slopes_ribbon
