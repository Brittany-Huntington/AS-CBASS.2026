#Code for processing ED50 values

#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(rstudioapi)
library(RColorBrewer)
library(CBASSED50) #Voolstra R package

mandatory_columns()
rm(list = ls())

#####Load and Prep Data-----------------------------
# Get the input file path; click on which site file you want to work with
input_data_path <- selectFile(
  caption = "C:/github/CBASS/AS-CBASS.2026/Data")

# Read data based on file format
cbass_dataset <- read_data(input_data_path)
cbass_dataset <-  cbass_dataset %>% select(-PAM1, -PAM2, -comments) #remove columns that contain blank rows
View(cbass_dataset)

# To specify the prefix for output files
output_prefix <- tools::file_path_sans_ext(input_data_path)
output_plot <- "C:/github/CBASS/AS-CBASS.2026/Plots"

rlog::log_info(paste("Your current directory is", getwd()))
rlog::log_info(paste("Your input filename is", basename(input_data_path)))
rlog::log_info(paste("The output files will be written into", output_prefix))

#r process-and-validate-cbass-dataset
cbass_dataset <- preprocess_dataset(cbass_dataset)
validate_cbass_dataset(cbass_dataset)



####Explore ED5s, ED50s, and ED95s---------------------
#create models
grouping_properties <- c("Site", "Species")
drm_formula <- "Pam_value ~ Temperature"
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
write.csv(eds_df,
          paste(output_prefix, "EDsdf.csv", sep = '_'),
          row.names = FALSE)


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
save_path <- file.path(output_plot, "Site-001_ED50.pdf")
ggsave(save_path, eds_boxplot,  width = 16, height = 9, device = "pdf")


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
save_path <- file.path(output_plot, "Site-011_exploratory_curve.pdf")
ggsave(save_path, exploratory_curve,  width = 16, height = 9, device = "pdf")



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
  scale_y_continuous(expand = c(0, 0.2))

tempresp_curve
#update site name
save_path <- file.path(output_plot, "Site-011_tempresp_curve.pdf")
ggsave(save_path, tempresp_curve,  width = 16, height = 9, device = "pdf")
