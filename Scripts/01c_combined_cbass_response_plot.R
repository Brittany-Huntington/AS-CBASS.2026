library(ggplot2)
library(patchwork)
library(here)

# Load shared theme/scales
source(here("Scripts/00_visualization_prep.R"))

# Load pre-made plots
ed50_emm  <- readRDS(here("Outputs", "plot_ed50.rds"))
color_emm <- readRDS(here("Outputs", "plot_color.rds"))

#removing spp from ed50 since share
ed50_emm <- ed50_emm +  
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
    #axis.ticks.x = element_blank()
  )
# ------------------------------------------------------------------------------
# Combine with Patchwork
# ------------------------------------------------------------------------------
# '+' places plots side-by-side; '/' places them stacked vertically
combined_fig <- (ed50_emm / color_emm) + 
  scale_color_manual(
    values = species_colors, 
    limits = species_order, 
    labels = species_labels
  ) +
  
  # Applies full names to x-axis tick labels
  scale_x_discrete(
    limits = species_order, 
    labels = species_labels
  )+
  plot_layout(
    guides = "collect", # Collects shared species/site legends into a single legend
    axis_titles = "collect") + 
  plot_annotation(
    tag_levels = "A" # Adds publication tags (A, B)
  ) & 
  theme(legend.position = "right")

print(combined_fig)


# ------------------------------------------------------------------------------
# Export Multi-Panel Figure
# ------------------------------------------------------------------------------
# For a 2-panel vertical plot (ED50 on top, Color on bottom):
ggsave(
  filename = here("Plots", "Fig2_ED50_and_Color_Combined.pdf"),
  plot     = combined_fig,
  width    = 6.5,       # Full page width (standard two-column journal width)
  height   = 8.5,       # Proportional height for a 2-panel vertical stack
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = here("Plots", "Fig2_ED50_and_Color_Combined.png"),
  plot     = combined_fig,
  width    = 6.5,
  height   = 8.5,
  units    = "in",
  dpi      = 600
)
