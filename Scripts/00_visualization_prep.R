library(ggplot2)

# --- Shapes & Species Palettes ---
custom_shapes  <- c(23, 22, 15, 16, 17, 18, 21, 24)
all_sites <- c("1", "3", "4", "5", "7", "9", "10", "11")

species_colors <- c(
  "AABR" = "#046E8F",
  "AGLO" = "#83BA75FF",
  "AHYA" = "#D44D5C",
  "ICRA" = "#462255"
)

species_order  <- c("ICRA", "AGLO", "AABR", "AHYA")

species_labels <- c(
  "AGLO" = "A. globiceps",
  "AABR" = "A. abrotanoides",
  "AHYA" = "A. hyacinthus",
  "ICRA" = "I. crateriformis"
)

# --- Cluster & Interaction Mappings ---
K_names       <- c(
  "2" = "Group I", 
  "1" = "Group II", 
  "3" = "Group III"
)

K_names_order <- c("Group I", "Group II", "Group III")

# Define Roman Numeral Map & Preferred Factor Order
roman_map <- c(
  "Cluster_1" = "Cluster II",
  "Cluster_2" = "Cluster I",
  "Cluster_3" = "Cluster III",
  "1"         = "Cluster II",
  "2"         = "Cluster I",
  "3"         = "Cluster III"
)

target_levels <- c("Cluster I", "Cluster II", "Cluster III")

# Updated Cluster Colors mapping directly to the new Roman Numeral labels
cluster_colors <- c(
  "Cluster I"   = "violet", # Group I (K=2)
  "Cluster II"  = "orange", # Group II (K=1)
  "Cluster III" = "cyan"    # Group III (K=3)
)



# --- Metric & Annotation Palettes (pheatmap / ComplexHeatmap) ---
metric_colors <- c(
  "Quant (Quantum Yield)"             = "green",
  "Sigma (Antenna Size)"              = "#FCF340",
  "Connect (Connectivity)"            = "#1F77B4", 
  "Connect (Normalized Curves)"       = "blue", 
  "Tau1 (Transport Kinetics 1)"       = "#17BECF", 
  "Tau2 (Transport Kinetics 2)"       = "#9467BD", 
  "NPQ (Non-Photochemical Quenching)" = "#8C564B",
  "qP (Photochemical Quenching)"      = "#E377C2", 
  "ABQ (Antenna Bed Quenching)"                  = "#FF7F0E",
  "qm (Max Quenching)"                = "#FF13F0", 
  "OJIP Kinetics (Derivatives)"       = "#98DF8A", 
  "OJIP Ratios & Area"                = "#FFBB78", 
  "Other Metric"                      = "#7F7F7F"  
)

ann_colors <- list(
  Family  = metric_colors,
  Species = species_colors,
  K_3     = cluster_colors
)

# --- Shared Publication Theme ---
theme_coral_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = 11),
      legend.text     = element_text(size = 10, face = "italic"),
      axis.title      = element_text(face = "bold", size = 12, color = "black"),
      axis.text.y     = element_text(size = 11, color = "black"),
      axis.text.x     = element_text(size = 11, color = "black", face = "bold"),
      axis.line       = element_line(linewidth = 0.6, color = "black"),
      axis.ticks      = element_line(linewidth = 0.6, color = "black"),
      plot.margin     = margin(t = 8, r = 8, b = 5, l = 8, unit = "pt")
    )
}