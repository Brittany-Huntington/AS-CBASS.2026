rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(purrr)
library(sf)
library(viridis)
library(cowplot)
library(svglite)
library(here)

source(here("Scripts/00_visualization_prep.R"))
load(here("Outputs", "photophys_preprocessed_data.RData"))
hab       <- read.csv(here("Data", "fixed_sites_ensemble_suitability_2025.csv"))
density   <- read.csv(here("Data", "ESA_CORALS_SITE_LEVEL_SURVEY_TUT_2025.csv"))
candidate <- read.csv(here("Data", "Candidate Fixed Sites_final.csv"))

# --- 1. CALCULATE DENSITIES ---
density_clean <- density %>%
  filter(!is.na(LATITUDE) & !is.na(LONGITUDE)) %>%
  group_by(SITE, LATITUDE, LONGITUDE) %>%
  summarize(
    # Total counts divided by total area surveyed across all divers at that site
    AGLO_density = sum(AGLO_COUNT, na.rm = TRUE) / sum(AREA_SURVEYED_SQUARE_METERS, na.rm = TRUE),
    ICRA_density = sum(ICRA_COUNT, na.rm = TRUE) / sum(AREA_SURVEYED_SQUARE_METERS, na.rm = TRUE),
    Total_Area_m2 = sum(AREA_SURVEYED_SQUARE_METERS, na.rm = TRUE),
    Divers_Count  = n(),
    .groups = "drop"
  )

tutuila_shape <- st_read(here("Data/Tut_shapefiles/TUT.shp"))

# --- 2. PREPARE FIXED SITES (HABITAT SUITABILITY SQUARES) ---
classify_suitability <- function(x) {
  case_when(
    x < 0.33 ~ "Low (<0.33)",
    x < 0.66 ~ "Medium (0.33-0.66)",
    TRUE      ~ "High (>0.66)"
  )
}

fixed_sf <- hab %>%
  filter(!is.na(Lat) & !is.na(Long)) %>%
  mutate(
    AGLO_cat = factor(classify_suitability(aglo_suitability_mean), 
                      levels = c("Low (<0.33)", "Medium (0.33-0.66)", "High (>0.66)")),
    ICRA_cat = factor(classify_suitability(icra_suitability_mean), 
                      levels = c("Low (<0.33)", "Medium (0.33-0.66)", "High (>0.66)"))
  ) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

# Stack AGLO (top) and ICRA (bottom) squares for fixed sites
fixed_aglo <- fixed_sf %>% mutate(Species = "AGLO", Suitability = AGLO_cat)
fixed_icra <- fixed_sf %>% mutate(Species = "ICRA", Suitability = ICRA_cat)

st_geometry(fixed_aglo) <- (st_geometry(fixed_aglo) + c(0, 0.003)) %>% st_set_crs(4326)
st_geometry(fixed_icra) <- (st_geometry(fixed_icra) - c(0, 0.003)) %>% st_set_crs(4326)

fixed_stacked <- bind_rows(fixed_aglo, fixed_icra)

# --- 3. PREPARE RANDOM SURVEY SITES (CALCULATED DENSITIES) ---
# Using the updated 'density' dataframe (LATITUDE / LONGITUDE columns)
random_aglo <- density_clean %>%
  select(SITE, LATITUDE, LONGITUDE, Density = AGLO_density) %>%
  mutate(Species = "AGLO", Label = sprintf("%.2f", Density)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

random_icra <- density_clean %>%
  select(SITE, LATITUDE, LONGITUDE, Density = ICRA_density) %>%
  mutate(Species = "ICRA", Label = sprintf("%.2f", Density)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

# Apply spatial offsets so AGLO is nudged north and ICRA south
st_geometry(random_aglo) <- (st_geometry(random_aglo) + c(0, 0.002)) %>% st_set_crs(4326)
st_geometry(random_icra) <- (st_geometry(random_icra) - c(0, 0.002)) %>% st_set_crs(4326)

random_stacked <- bind_rows(random_aglo, random_icra)

# Ensure spatial projections align with island shapefile
random_stacked <- st_transform(random_stacked, st_crs(tutuila_shape))
fixed_stacked  <- st_transform(fixed_stacked, st_crs(tutuila_shape))
fixed_sf       <- st_transform(fixed_sf, st_crs(tutuila_shape))


# --- 1. SEPARATE ZERO AND NON-ZERO DENSITY DATA ---
random_nonzero <- random_stacked %>% filter(Density > 0)
random_zero    <- random_stacked %>% filter(Density == 0.00000000)

# --- 2. GENERATE COMBINED MAP ---
ggplot() +
  # Island Base Layer
  geom_sf(data = tutuila_shape, fill = "grey90", color = "grey50", linewidth = 0.3) +
  
  # Layer 1a: Zero Density Points (Rendered as 'X')
  geom_sf(
    data = random_zero, 
    aes(color = Species), 
    shape = 4, stroke = 1.2, size = 3, alpha = 0.8 # shape 4 is an 'X'
  ) +
  
  # Layer 1b: Non-Zero Density Points (Scaled Circles)
  geom_sf(
    data = random_nonzero, 
    aes(color = Species, size = Density), 
    shape = 16, alpha = 0.6 # shape 16 is a solid filled circle
  ) +
  
  # Layer 2: Raw Density Numeric Labels
  geom_sf_text(
    data = random_stacked %>% filter(Density > 0),
    aes(label = Label, color = Species),
    fontface = "bold",
    size = 2.8,
    nudge_x = 0.0035,
    show.legend = FALSE
  ) +
  
  # Layer 3: Fixed Site Suitability Squares (Stacked)
  geom_sf(
    data = fixed_stacked, 
    aes(fill = Suitability), 
    shape = 22, size = 4.5, color = "black", stroke = 0.8
  ) +
  
  # Layer 4: Fixed Site Names
  geom_sf_text(
    data = fixed_sf, 
    aes(label = Site), 
    nudge_y = 0.008, size = 3, fontface = "bold"
  ) +
  
  # Aesthetics & Formatting
  scale_color_manual(
    values = c("AGLO" = "#83BA75FF", "ICRA" = "#462255"),
    name = "Species"
  ) +
  scale_fill_brewer(palette = "YlOrRd", name = "Habitat Suitability") +
  scale_size_continuous(name = "Colony Density\n(colonies/m²)", range = c(2, 6)) +
  labs(
    title = "Tutuila Coral Survey: Calculated Densities & Site Suitability",
    subtitle = "Circles: Density > 0 | 'X': Zero Density | Top: AGLO | Bottom: ICRA",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "white"),
    legend.position = "right"
  )
############
# Process Fixed Sites (Blue Squares)
fixed_ref <- hab %>%
  filter(!is.na(Lat) & !is.na(Long)) %>%
  select(Site_Name = Site, Lat, Long) %>%
  mutate(Type = "Fixed Site (Chosen)") %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

# Process Random Sites (Orange Circles)
random_ref <- density_clean %>%
  filter(!is.na(LATITUDE) & !is.na(LONGITUDE)) %>%
  select(Site_Name = SITE, Lat = LATITUDE, Long = LONGITUDE) %>%
  mutate(Type = "Random Survey Site") %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

# Align projections with shapefile
fixed_ref  <- st_transform(fixed_ref, st_crs(tutuila_shape))
random_ref <- st_transform(random_ref, st_crs(tutuila_shape))

# Combine for easier labeling
all_sites <- bind_rows(
  fixed_ref %>% mutate(Coords = st_coordinates(.)),
  random_ref %>% mutate(Coords = st_coordinates(.))
) %>%
  mutate(
    X = Coords[, 1],
    Y = Coords[, 2]
  )

# --- 2. GENERATE CROSS-REFERENCE MAP ---

ggplot() +
  # Island Base Layer
  geom_sf(data = tutuila_shape, fill = "grey92", color = "grey60", linewidth = 0.3) +
  
  # Random Sites (Orange Circles)
  geom_sf(
    data = random_ref, 
    aes(color = Type), 
    shape = 16, size = 3, alpha = 0.8
  ) +
  
  # Fixed Sites (Blue Squares)
  geom_sf(
    data = fixed_ref, 
    aes(color = Type), 
    shape = 15, size = 4
  ) +
  
  # Non-Overlapping Labels for All Sites
  geom_text_repel(
    data = all_sites,
    aes(x = X, y = Y, label = Site_Name, color = Type),
    fontface = "bold",
    size = 2.8,
    max.overlaps = 50,      # Pushes text away to force visibility
    box.padding = 0.4,       # Distance around label box
    point.padding = 0.3,     # Distance from point
    segment.color = "grey40",# Line connecting label to point
    segment.size = 0.3
  ) +
  
  # Aesthetics & Formatting
  scale_color_manual(
    values = c(
      "Fixed Site (Chosen)" = "#0072B2",  # Blue
      "Random Survey Site"  = "#D55E00"   # Orange
    ),
    name = "Site Type"
  ) +
  labs(
    title = "Tutuila Site Cross-Reference Map",
    subtitle = "Compare Random Survey Site Names (Orange) vs. Fixed Site Names (Blue)",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "white"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )
