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

mandatory_columns()
rm(list = ls())

custom_shapes <- c(23, 22, 15, 16, 17, 18, 21, 24)
species_colors <- c(
  "AABR" = "#046E8F", # Muted Crimson/Rose Red (#D44D5CFF)
  "AGLO" = "#83BA75FF", # Light Sage Teal (#9DD9D2FF)
  "AHYA" = "#D44D5C", # Deep Ocean Blue (#046E8FFF)
  "ICRA" = "#462255"  # Bright Amber/Orange (#FF8811FF) -- (or "#462255" for Deep Purple)
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
  select(Cookie_no, Site, Bag_number, Tank_color)

fvfm <- read_csv(here("Data/TUT_CBASS_raw_PAM_fvfm.csv") )%>%
  mutate(
    Site_number  = as.character(Site_number)
  )%>%
  select(-Date)


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
master_dataset_fvfm <- master_dataset_bag %>%
  # Joinfvfm by 
  left_join(fvfm, by = c( "Bag_number" = "Bag_number", "Site" = "Site_number", "Temperature" = "Temperature", "Tank_color" = "Tank_color"))%>%
  select(-Click_No, -Identifier)%>%
  rename("Color_mean" = "Mean", "Color_SD" = "SD" )%>%
  mutate(
    Temperature = as.numeric(Temperature),
    Color_mean        = as.numeric(Color_mean),
    # If using raw 8-bit gray values (0 = dark, 255 = light):
    Mean_inverted = 10 - Color_mean)
  # )%>%
  # filter(!is.na(Temperature), !is.na(Mean))

# Save your master dataset
write_csv(master_dataset_fvfm, "master_coral_dataset_merged.csv")


# # Read data based on file format
# cbass_dataset <- read_data(input_data_path)
# cbass_dataset <-  cbass_dataset %>% select(-PAM1, -PAM2, -comments) #remove columns that contain blank rows
# View(cbass_dataset)

# To specify the prefix for output files
output_prefix <- tools::file_path_sans_ext(input_data_path)
output_plot <- here("Plots")

rlog::log_info(paste("Your current directory is", getwd()))
rlog::log_info(paste("Your input filename is", basename(input_data_path)))
rlog::log_info(paste("The output files will be written into", output_prefix))

#r process-and-validate-cbass-dataset


cbass_dataset <- preprocess_dataset(master_dataset_fvfm)
validate_cbass_dataset(cbass_dataset)
convert_columns(cbass_dataset)
dataset_has_mandatory_columns(cbass_dataset)



####Explore ED5s, ED50s, and ED95s---------------------
#create models
grouping_properties <- c("Site", "Species")
drm_formula <- "Color_mean ~ Temperature"
models <- fit_drms(cbass_dataset, grouping_properties, drm_formula, is_curveid = TRUE)

#get-eds
eds <- get_all_eds_by_grouping_property(models)
View(eds)
eds$GroupingProperty[eds$GroupingProperty == "1_AHYA"] <- "1_AHYA_1" #manual fix for when there is only one replicate of a species in the run

cbass_dataset <- define_grouping_property(cbass_dataset, grouping_properties) %>%
  mutate(GroupingProperty = paste(GroupingProperty, Genotype, sep = "_"))

eds_df <- 
  left_join(eds, cbass_dataset, by = "GroupingProperty") %>%
  select(names(eds), all_of(grouping_properties)) %>%
  distinct()

head(eds_df)
# write.csv(
#   eds_df,
#   paste0(output_prefix, "_EDsdf.csv"),
#   row.names = FALSE
# )


####Plotting-------------------------
#ED50 boxplot
eds_boxplot <- eds_df %>% ggplot(
  aes(x = Species, y = ED50, color = Species)) +
  geom_boxplot() + 
  stat_summary(
    fun = mean, 
    geom = "text", 
    aes(label = round(after_stat(y), 2)), show.legend = F,
    position = position_dodge(width = 0.75),
    vjust = -1
  ) +
  facet_grid(~ Site,labeller = as_labeller(function(x) paste("Site", sprintf("%03d", as.numeric(x))))) +
  ylab("ED50s - Temperatures [C°]")+
  scale_color_brewer(palette = "Set2")

eds_boxplot
#update site name
# ggsave(
#   here("Plots", "Site-001_ED50.pdf"),
#   eds_boxplot,
#   width = 16,
#   height = 9,
#   device = "pdf"
# )


exploratory_curve <- ggplot(data = cbass_dataset,
                            aes(x = Temperature, y = Pam_value,
                                group = GroupingProperty, # You can play around with the group value (e.g., Species, Site, Condition)
                                color = Genotype)) +
  geom_smooth(
    method = drc::drm,
    method.args = list(
      fct = drc::LL.3()),
    se = FALSE,
    size = 0.7
  ) +
  geom_point(size = 1.5) +
  facet_grid(Species ~ Site) +
  scale_color_brewer(palette = "Set2")


exploratory_curve
#update site name
# save_path <- file.path(output_plot, "Site-011_exploratory_curve.pdf")
# ggsave(save_path, exploratory_curve,  width = 16, height = 9, device = "pdf")



#Predict PAM values for assayed temperature range-------------
#Curves display the predicted PAM values, the 95% confidence intervals, and mean ED5s, ED50s, and ED95s for groupings (vertical line).

# First fit models with curveid = FALSE and with LL.4 = FALSE; If you get error messages, try LL.4 = TRUE
models <- fit_drms(cbass_dataset, grouping_properties, drm_formula, is_curveid = FALSE, LL.4 = FALSE)
# The default number of values for range of temperatures is 100
temp_ranges <- define_temperature_ranges(cbass_dataset$Temperature, n=100)
predictions <- get_predicted_pam_values(models, temp_ranges)

predictions_df <- 
  left_join(predictions,
            define_grouping_property(cbass_dataset, grouping_properties) %>% 
              select(c(all_of(grouping_properties), GroupingProperty)),
            by = "GroupingProperty",
            relationship = "many-to-many") %>%
  distinct()


summary_eds_df <- eds_df %>%
  group_by(Site, Species) %>%
  summarise(Mean_ED5 = mean(ED5),
            SD_ED5 = sd(ED5),
            SE_ED5 = sd(ED5) / sqrt(n()),
            Conf_Int_5 = qt(0.975, df = n() - 1) * SE_ED5,
            Mean_ED50 = mean(ED50),
            SD_ED50 = sd(ED50),
            SE_ED50 = sd(ED50) / sqrt(n()),
            Conf_Int_50 = qt(0.975, df = n() - 1) * SE_ED50,
            Mean_ED95 = mean(ED95),
            SD_ED95 = sd(ED95),
            SE_ED95 = sd(ED95) / sqrt(n()),
            # The value 0.975 corresponds to the upper tail probability
            # for a two-tailed t-distribution with a 95% 
            Conf_Int_95 = qt(0.975, df = n() - 1) * SE_ED95) %>%
  mutate(across(c(Mean_ED50, SD_ED50, SE_ED50,
                  Mean_ED5, SD_ED5, SE_ED5,
                  Mean_ED95, SD_ED95, SE_ED95,
                  Conf_Int_5,Conf_Int_50,Conf_Int_95), ~round(., 2)))

summary_eds_df
write.csv(summary_eds_df,paste(output_prefix, "summaryEDs_df.csv", sep = '_'), row.names = FALSE)

result_df <- predictions_df %>%
  left_join(summary_eds_df, by = c("Site", "Species"))


tempresp_curve <- ggplot(result_df,
                         aes(x = Temperature,
                             y = PredictedPAM,
                             group = GroupingProperty,
                             color = Species)) +
  geom_line() +
  geom_ribbon(aes(ymin = Upper,
                  ymax = Lower,
                  fill = Species),
              alpha = 0.2,
              linetype = "dashed") +
  geom_segment(aes(x = Mean_ED5,
                   y = 0,
                   xend = Mean_ED5,
                   yend = max(Upper)),
               linetype = 3) +
  geom_text(mapping=aes(x = Mean_ED5,
                        y = max(Upper) + 0.12,
                        label = round(Mean_ED5, 2)),
            size = 3, angle = 90, check_overlap = T) +
  geom_segment(aes(x = Mean_ED50,
                   y = 0,
                   xend = Mean_ED50,
                   yend = max(Upper)),
               linetype = 3) +
  geom_text(mapping=aes(x = Mean_ED50,
                        y = max(Upper) + 0.12,
                        label = round(Mean_ED50, 2)),
            size = 3, angle = 90, check_overlap = T) +
  geom_segment(aes(x = Mean_ED95,
                   y = 0,
                   xend = Mean_ED95,
                   yend = max(Upper)),
               linetype = 3) +
  geom_text(mapping=aes(x = Mean_ED95,
                        y = max(Upper) + 0.12,
                        label = round(Mean_ED95, 2)),
            size = 3, angle = 90, check_overlap = T) +
  facet_grid(Species ~ Site) +
  # To add the real PAM and compare with predicted values
  geom_point(data = cbass_dataset,
             aes(x = Temperature,
                 y = Pam_value)) +
  xlab("Temperature [C°]")+
  scale_y_continuous(expand = c(0, 0.2))+
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors)

tempresp_curve
#update site name
# save_path <- file.path(output_plot, "Site-011_tempresp_curve.pdf")
# ggsave(save_path, tempresp_curve,  width = 16, height = 9, device = "pdf")
