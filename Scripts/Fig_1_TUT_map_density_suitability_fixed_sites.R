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
library(ggrepel)

aglo_palette <- c("Low (<0.33)" = "#E8F5E9", "Medium (0.33-0.66)" = "#A5D6A7", "High (>0.66)" = "#83BA75FF")

# ICRA: Light to Dark Purples (High = #462255)
icra_palette <- c("Low (<0.33)" = "#E1D5E7", "Medium (0.33-0.66)" = "#9B7BB0", "High (>0.66)" = "#462255")

source(here("Scripts/00_visualization_prep.R"))
tutuila_shape <- st_read(here("Data/Tut_shapefiles/TUT.shp"))
load(here("Outputs", "photophys_preprocessed_data.RData"))
hab       <- read.csv(here("Data", "fixed_sites_ensemble_suitability_2025.csv"))
density   <- read.csv(here("Data", "ESA_CORALS_SITE_LEVEL_SURVEY_TUT_2025.csv"))
candidate <- read.csv(here("Data", "Candidate Fixed Sites_final.csv"))



# --- 1. CALCULATE SITE-LEVEL AGGREGATED DENSITIES ---
density_clean <- density %>%
  filter(!is.na(LATITUDE) & !is.na(LONGITUDE)) %>%
  group_by(SITE) %>%
  summarize(
    AGLO_density = sum(AGLO_COUNT, na.rm = TRUE) / sum(AREA_SURVEYED_SQUARE_METERS, na.rm = TRUE),
    ICRA_density = sum(ICRA_COUNT, na.rm = TRUE) / sum(AREA_SURVEYED_SQUARE_METERS, na.rm = TRUE),
    .groups = "drop"
  )

# --- 2. JOIN MATCHED DENSITIES DIRECTLY TO FIXED SITES (HAB) ---
classify_suitability <- function(x) {
  case_when(
    x < 0.33 ~ "Low (<0.33)",
    x < 0.66 ~ "Medium (0.33-0.66)",
    TRUE      ~ "High (>0.66)"
  )
}

hab_matched <- hab %>%
  filter(!is.na(Lat) & !is.na(Long)) %>%
  # Join calculated random densities via Random_site column
  left_join(density_clean, by = c("Random_site" = "SITE")) %>%
  mutate(
    AGLO_cat = factor(classify_suitability(aglo_suitability_mean), 
                      levels = c("Low (<0.33)", "Medium (0.33-0.66)", "High (>0.66)")),
    ICRA_cat = factor(classify_suitability(icra_suitability_mean), 
                      levels = c("Low (<0.33)", "Medium (0.33-0.66)", "High (>0.66)")),
    AGLO_label = ifelse(!is.na(AGLO_density) & AGLO_density > 0, sprintf("%.2f", AGLO_density), ""),
    ICRA_label = ifelse(!is.na(ICRA_density) & ICRA_density > 0, sprintf("%.2f", ICRA_density), "")
  ) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

# --- 3. STACK AND APPLY SPATIAL OFFSETS FOR PLOTTING ---
fixed_aglo <- hab_matched %>% 
  mutate(Species = "AGLO", Suitability = AGLO_cat, Density_Text = AGLO_label, Density = AGLO_density)

fixed_icra <- hab_matched %>% 
  mutate(Species = "ICRA", Suitability = ICRA_cat, Density_Text = ICRA_label, Density = ICRA_density)

# Apply spatial offset (AGLO shifted north +0.003 deg, ICRA south -0.003 deg)
st_geometry(fixed_aglo) <- (st_geometry(fixed_aglo) + c(0, 0.003)) %>% st_set_crs(4326)
st_geometry(fixed_icra) <- (st_geometry(fixed_icra) - c(0, 0.003)) %>% st_set_crs(4326)

fixed_stacked <- bind_rows(fixed_aglo, fixed_icra) %>%
  st_transform(st_crs(tutuila_shape))

fixed_sf <- st_transform(hab_matched, st_crs(tutuila_shape))


# Create interaction levels for dual palette mapping
fixed_stacked <- fixed_stacked %>%
  mutate(Spec_Suit = factor(paste(Species, Suitability, sep = "_")))

# Combined color palette mapping
combined_palette <- c(
  "AGLO_Low (<0.33)"       = "#E8F5E9",
  "AGLO_Medium (0.33-0.66)" = "#A5D6A7",
  "AGLO_High (>0.66)"      = "#83BA75FF",
  "ICRA_Low (<0.33)"       = "#E1D5E7",
  "ICRA_Medium (0.33-0.66)" = "#9B7BB0",
  "ICRA_High (>0.66)"      = "#462255"
)

ggplot() +
  geom_sf(data = tutuila_shape, fill = "grey92", color = "grey60", linewidth = 0.3) +
  geom_sf(
    data = fixed_stacked, 
    aes(fill = Spec_Suit), 
    shape = 22, size = 8, color = "black", stroke = 0.8
  ) +
  geom_sf_text(
    data = fixed_stacked,
    aes(label = Density_Text, color = Species),
    fontface = "bold", size = 2.4
  ) +
  geom_sf_text(
    data = fixed_sf, 
    aes(label = Site), 
    nudge_x = 0.001, 
    nudge_y = -0.009,
    size = 3, fontface = "bold"
  ) +
  scale_color_manual(values = c("AGLO" = "black", "ICRA" = "white"), guide = "none") +
  scale_fill_manual(
    values = combined_palette,
    name = "Species & Habitat Suitability Tier"
  ) +
  labs(
    title = "Tutuila Fixed Sites: Suitability & Density",
    subtitle = "Top Square: AGLO (Green) | Bottom Square: ICRA (Purple)"
  ) +
  theme_minimal()

# --- 4. GENERATE MAP WITH DENSITIES INSIDE SQUARES (ONE SCALE) ---
ggplot() +
  # Island Base Layer
  geom_sf(data = tutuila_shape, fill = "grey92", color = "grey60", linewidth = 0.3) +
  
  # Layer 1: Large Fixed Site Squares Filled by Habitat Suitability
  geom_sf(
    data = fixed_stacked, 
    aes(fill = Suitability), 
    shape = 22, size = 8, color = "black", stroke = 0.8
  ) +
  
  # Layer 2: Density Value Printed Directly INSIDE the Squares
  geom_sf_text(
    data = fixed_stacked,
    aes(label = Density_Text),
    fontface = "bold",
    size = 2.4,
    color = "black"
  ) +
  
  # Layer 3: Site Name Labels
  geom_sf_text(
    data = fixed_sf, 
    aes(label = Site), 
    nudge_y = 0.008, size = 3, fontface = "bold"
  ) +
  
  # Aesthetics & Formatting
  scale_fill_brewer(palette = "YlOrRd", name = "Habitat Suitability") +
  labs(
    title = "Tutuila Fixed Sites: Habitat Suitability & Matched Densities",
    subtitle = "Stacked Squares — Top: AGLO | Bottom: ICRA (Numbers inside = Density in colonies/m²)",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "white"),
    legend.position = "right"
  )