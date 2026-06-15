#Code for exploring how ED values vary depending on which temp data is included.

#load packages
library(tidyverse)
library(readxl)
library(rstudioapi)
library(RColorBrewer)
library(CBASSED50) #Voolstra R package
library(lme4)
library(here)


rm(list = ls())
#####Load raw CBASS run data-----------------------------
S1 <- read.csv(here("Data", "RAW_ESA-001_20260214_PAM.csv"))
S3 <- read.csv(here("Data", "RAW_ESA-003_20260219_PAM.csv"))
S4 <- read.csv(here("Data", "RAW_ESA-004_20260220_PAM.csv"))
S5 <- read.csv(here("Data", "RAW_ESA-005_20260223_PAM.csv"))
S7 <- read.csv(here("Data", "RAW_ESA-007_20260213_PAM.csv"))
S9 <- read.csv(here("Data", "RAW_ESA-009_20260216_PAM.csv"))
S10 <- read.csv(here("Data", "RAW_ESA-010_20260217_PAM.csv"))
S11 <- read.csv(here("Data", "RAW_ESA-011_20260210_PAM.csv"))


dat <- bind_rows(S1,S3, S4, S5, S7, S9, S10, S11)
dat <- dat %>% mutate(across(where(is.character), as.factor))
dat$Site <- as.character(dat$Site)
dat$Genotype <- as.character(dat$Genotype )
colnames(dat)

#define ramp temp datasets------------------------
temp_sets <- list(
  Classic = c(29,36,38,39),
  Refined = c(29,35,36,37,37.5,38,38.5,39),
  Max_common = c(29,35,36,37,38,39),
  All = sort(unique(dat$Temperature))
)
tmin <- min(temp_sets$Classic, na.rm = TRUE)
tmax <- max(temp_sets$Classic, na.rm = TRUE)

#create function to process a ramp dataset
process_cbass_dataset <- function(dat, temps, dataset_name){
  
  # subset ramp temps and ensure full ramp exists for each site
  cbass_subset <- dat %>%
    filter(Temperature %in% temps)
  
  # Only enforce full ramp for specific ramp designs
  if(dataset_name != "All"){
    cbass_subset <- cbass_subset %>%
      group_by(Site) %>%
      filter(all(temps %in% unique(Temperature))) %>%
      ungroup()
  }
  
  cbass_subset <- cbass_subset %>%
    select(-PAM1, -PAM2, -comments)
  
  cbass_subset %>%
    group_by(Site) %>%
    summarise(n_temps = n_distinct(Temperature))
  
  # preprocess CBASS dataset
  cbass_subset <- preprocess_dataset(cbass_subset)
  validate_cbass_dataset(cbass_subset)
  
  # fit DRM models
  grouping_properties <- c("Site","Species")
  drm_formula <- "Pam_value ~ Temperature"
  
  models <- fit_drms(cbass_subset,
                     grouping_properties,
                     drm_formula,
                     is_curveid = TRUE)
  
  # extract ED values
  eds <- get_all_eds_by_grouping_property(models)
  
  # separate grouping property safely
  eds <- eds %>%
    separate(GroupingProperty,
             into = c("Site","Species","Genotype"),
             sep="_",
             fill="right") %>%
    mutate(Genotype = ifelse(is.na(Genotype),"1",Genotype))
  
  # keep relevant columns
  eds_df <- eds %>%
    select(ED5,ED50,ED95,Site,Species,Genotype) %>%
    mutate(Dataset = dataset_name)
  
  return(eds_df)
}

#run fuction across all ramp designs
eds_all <- bind_rows(lapply(names(temp_sets), function(name){
    process_cbass_dataset(dat,
                          temps = temp_sets[[name]],
                          dataset_name = name)
     })
)

#save combined dataset
write.csv(eds_all, here ("Data", "EDsdf_all_ramps.csv"), row.names = FALSE)


#run model comparing ED50 estimate among temp ramp designs------------------------------

#explore for outliers (getting those odd values under the "classic" temp ramp)
View(eds_all)
eds_all <- eds_all%>% mutate(across(where(is.character), as.factor))
#To ensure results were not influenced by unequal species representation, analyses were repeated for the two most commonly represented species (AGLO and ICRA).
eds_two <- eds_all %>% filter(Species %in% c("AGLO", "ICRA"), ED50 >= tmin, ED50 <=tmax)

#linear model
mod_full <- lmer(ED50 ~ Dataset + Species + (1|Site), data = eds_all)

#Remove  clear model outliers (e.g. standardized residuals > 3); n = 4 from the classic model
eds_all$resid <- resid(mod_full)
eds_all$z <- scale(eds_all$resid)
outliers <- eds_all %>% filter(abs(z) > 3)
eds_clean <- eds_all %>% filter(abs(z) <= 3)

#re-run linear model with outliers removed
mod_clean <- lmer(ED50 ~ Dataset + Species + (1|Site), data = eds_clean)
mod_two <- lmer(ED50 ~ Dataset + Species + (1|Site), data = eds_two)

anova(mod_full)
anova(mod_clean)
anova(mod_two)


g1 <- ggplot(eds_clean, aes(x = Dataset, y = ED50, fill = Species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35) +
  geom_jitter(aes(color = Species),  width = 0.15, alpha = 0.35, size = 1.5, show.legend = FALSE) +
  facet_wrap(~ Species, scales = "fixed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Temp Ramp",
    y = expression(paste("ED50 (", degree, "C)")),
    title = "ED50 Comparisons Across Temp Ramps"
  ) +
  theme(strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  theme(legend.position = "none")

g1
ggsave(filename = here ("Plots", "ED50 Comparisons Across Temp Ramps.jpeg"), plot = g1, 
       width = 10, height = 6,  dpi = 300 )  



#run model comparing ED50 estimate from same sites calculated under classic vs refined model------------------------------

shared_sites <- c(3,4,5,9,10)

dat_shared <- eds_two%>%
  filter(Site %in% shared_sites,
         Dataset %in% c("Classic","Refined"),
         Species %in% c("AGLO","ICRA"))


#paired test using site as blocking factor to test is refined ramp changes ED values
mod_paired <- lmer(ED50 ~ Dataset * Species + (1 | Site), data = dat_shared)
anova(mod_paired) #dataset still not significant



#plot
g2 <- ggplot(dat_shared, aes(x = Dataset, y = ED50, fill = Species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35) +
  geom_jitter(aes(color = Species),
              width = 0.15,
              alpha = 0.8,
              size = 2,
              show.legend = FALSE) +
  facet_wrap(~ Species, scales = "fixed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Temp Ramp",
    y = expression(paste("ED50 (", degree, "C)")),
    title = "ED50 Comparisons Across Temp Ramps",
    subtitle = "Paired Sites"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  theme(legend.position = "none")

g2
ggsave(filename = here ("Plots", "ED50 Comparisons Across Temp Ramps--paired sites.jpeg"), plot = g2, 
  width = 10, height = 6, dpi = 300)
