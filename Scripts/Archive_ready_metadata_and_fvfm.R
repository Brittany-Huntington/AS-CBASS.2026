library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(here)
library(lubridate)
library(readxl)
load(here("Outputs", "photophys_preprocessed_data.RData"))
dat<-read.csv(here ("Data/data_cbass.csv"))
dat2<-read.csv(here ("Outputs/fvfm_color_df.csv"))

slope_lookup <- dat2 %>%
  dplyr::select(SampleID_clean, Color_mean, Color_SD, Rel_Temperature)

# Merge color_slope into dat
dat <- dat %>%
  left_join(slope_lookup, by = c("SampleID" = "SampleID_clean"))

write.csv(dat, here("Outputs/cbass_archive.csv"))

head(dat)
head(meta_filtered)

dat<-dat%>%
  filter_out(dat$`CBASS / Sophie`== "Sophie")

dat2<-read.csv(here ("Outputs/fvfm_color_df_for_archive.csv"))

dat_clean <- dat %>%
  mutate(
    Bag_number = trimws(as.character(Bag_number)),
    Species    = trimws(tolower(as.character(Species)))
  )%>%
  dplyr::select(-Site)

dat2_clean <- dat2 %>%
  mutate(
    Bag_number = trimws(as.character(Bag_number)),
    Species    = trimws(tolower(as.character(Species)))
  )


dat_clean <- dat_clean %>%
  mutate(Date = ymd(Date))

dat2_clean <- dat2_clean %>%
  mutate(Date = ymd(Date))

dat3 <- inner_join(
  dat_clean, 
  dat2_clean, 
  by = c("Date", "Bag_number", "Species")
)

missing_rows <- anti_join(dat, dat2, by = c("Date", "Bag_number", "Species"))

# View which sites are disappearing
missing_rows %>% count(Site)

unmatched_bags <- anti_join(
  dat, 
  dat2, 
  by = c("Date", "Bag_number", "Species")
)

# View the missing bags grouped by Site
unmatched_bags %>% 
  select(Site, Date, Bag_number, Species)



write.csv(dat3, here("Data/TUT_experimental_data.csv"))
