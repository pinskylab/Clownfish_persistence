# Script to calculate persistence metrics
# For now, using fake data + parameters where necessary to get structure right - will note clearly, format might shift as real stuff comes in

#################### Set-up: ####################
#Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(here)

# Load some data files while deciding how various scripts should interact/communicate
# Load output from AnemLocations.R (or source that file) - mean gps coordinates across multiple observations of anemones
#load(file=here("Data",'AnemAllInfowLatLon2.RData')) #file with anems, after 2018 anems matched
#load(file=here("Data", "AnemLatLonObsbyAnem.RData")) #file with lat-lon info appended

# Load encounter history data frame from ClownfishMarkRecap.R
#load(file=here("Data", "encounters_all.RData"))

# Set a few things
sample_years = c(2015,2016,2017,2018)

# FAKE DATA FOR NOW THAT WILL GET REPLACED BY REAL DATA
eggs_per_clutch = 200
clutch_per_year = 5
eggs_intercept = 100
eggs_slope = 50

#################### Functions: ####################
# Functions and constants from my GitHub function/constant collection
# script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQAXExRItr2tnepmkRr7NRt4hylZrks5aciBtwA%3D%3D", ssl.verifypeer = FALSE)
# eval(parse(text = script))
script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Functions from Michelle's GitHub helpers script
#field_helpers (similar idea to helpers, in field repository) - this might be the newer version of the helpers collection?
script <- getURL("https://raw.githubusercontent.com/pinskylab/field/master/scripts/field_helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Eggs by female size function (from Adam Y's work), could change to have default intercept or slope in there once have those better estimated
eggsBySize <- function(intercept, slope, fish_size) {
  eggs <- intercept + slope*fish_size
  return(eggs)
}

# Calculate LEP - consider updating later to include some of the calculating of survivals and such in here instead of as inputs? Probably better as separate function, though
findLEP <- function(l_vec, f_vec) {
  out <- sum(l_vec*f_vec)
  return(out)
}

#################### Running things: ####################
########## Pull info from database
leyte <- read_db("Leyte")

allfish_fish <- leyte %>% 
  tbl("clownfish") %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, cap_id, anem_table_id, recap, tag_id, color, size) %>%
  collect() 

allfish_anems <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
  collect() %>%
  filter(anem_table_id %in% allfish_fish$anem_table_id)

allfish_dives <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, dive_type, date, site, gps) %>%
  collect() %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id)

# pull out just the year and put that in a separate column
allfish_dives$year <- as.integer(substring(allfish_dives$date,1,4))

#join together
allfish <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
allfish <- left_join(allfish, allfish_dives, by="dive_table_id")

allfish$size <- as.numeric(allfish$size) #make size numeric (rather than a chr) so can do means and such

#pull out just tagged fish
taggedfish <- allfish %>% filter(!is.na(tag_id))

########## Calculating values, constants, inputs from data
# Find max age (A) - maximum number of times a fish has been caught
tag_years <- taggedfish %>% 
  group_by(tag_id) %>%
  summarize(firstyear = min(year), lastyear = max(year), nyears = max(year) - min(year))

maxyears = max(tag_years$nyears)

# add site back in
tag_years_site <- left_join(tag_years, (taggedfish %>% select(tag_id, site)) %>% distinct(), by="tag_id") # 4 tags correspond to two sites - looked at below, emailed Michelle about checking data sheets

# # Note to self - if you just do distinct() above, it has 4 more rows - so looks like 4 tags that have multiple sites associated with them - should check into at some point...
# tagsatsites <- data.frame(table(tag_years_site$tag_id))
# tagsatsites %>% filter(Freq > 1)
# tag_years_site %>% filter(tag_id == "982000411818588")
# tag_years_site %>% filter(tag_id == "982000411818610")
# tag_years_site %>% filter(tag_id == "985153000401241")
# tag_years_site %>% filter(tag_id == "986112100172501")

# trying this with cap_id, since right now the max number of years recaught is the max possible for the # of years we've been tagging
cap_years <- allfish %>% 
  filter(!is.na(cap_id)) %>%
  group_by(cap_id) %>%
  summarize(firstyear = min(year), lastyear = max(year), nyears = max(year) - min(year), tag = tag_id[1])

maxyears_cap = max(cap_years$nyears)

# Do any fish in 2017 or 2018 have cap_ids?
# not when I try to find them like this...
test_caps <- allfish %>% 
  filter(!is.na(cap_id)) %>%
  filter(year %in% c(2017,2018))

# or like this
table((allfish %>% filter(year %in% c(2017,2018)))$cap_id)

# or like this
tags_in_cap_years <- cap_years %>% filter(!is.na(tag)) %>% select(cap_id, tag)
test_caps2 <- allfish %>%
  filter(cap_id %in% tags_in_cap_years$cap_id) 
table(test_caps2$year)

# maybe cap_ids haven't been populated yet for 2017 and 2018 (genetics not done for those years yet but could be fish sequenced and tagged pre-2016 and known to be recaptured via tag)?
# try to find that one old Visca fish to see - haven't quite gotten this working yet...
test_Visca <- allfish %>%
  filter(site == "Visca") %>%
  filter(!is.na(tag_id)) %>%
  filter(year == 2016)

test_Visca2 <- allfish %>%
  filter(anem_spp)

anems_Visca <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_spp) %>%
  collect() %>%
  filter(anem_table_id %in% (allfish %>% filter(site == "Visca"))$anem_table_id) 

########## Looking at eggs and recruits relationships


########## Assessing metrics


#################### Plots: ####################
##### Histograms of number of years fish are recaught
# all sites together
pdf(file = here("Plots/PersistenceMetrics", "NYearsRecaught_allsites.pdf"))
ggplot(data = tag_years, aes(nyears)) +
  geom_histogram(binwidth=1) +
  xlab("number of years recaught") + ggtitle("# years tagged fish recaught, all sites combined") +
  theme_bw()
dev.off()

# sites separated out
pdf(file = here("Plots/PersistenceMetrics", "NYearsRecaught_bysite.pdf"))
ggplot(data = tag_years_site, aes(nyears)) +
  geom_histogram(binwidth=1) +
  facet_grid(.~ site, labeller = label_wrap_gen(10)) + 
  xlab("number of years recaught") + ggtitle("# years tagged fish recaught by site") +
  theme_bw()
dev.off()


