library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(here)
library(lubridate)

dat<-read_excel(here ("Data/TUT_CBASS_colony_metadata.xlsx"))

dat<-dat%>%
  filter_out(dat$`CBASS / Sophie`== "Sophie")

dat2<-read.csv(here ("Data/TUT_CBASS_raw_PAM_fvfm.csv"))

dat_clean <- dat %>%
  mutate(
    Bag_number = trimws(as.character(Bag_number)),
    Species    = trimws(tolower(as.character(Species)))
  )

dat2_clean <- dat2 %>%
  mutate(
    Bag_number = trimws(as.character(Bag_number)),
    Species    = trimws(tolower(as.character(Species)))
  )


dat_clean <- dat %>%
  mutate(Date = ymd(Date))

dat2_clean <- dat2 %>%
  mutate(Date = ymd(Date))

dat3 <- inner_join(
  dat, 
  dat2, 
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
