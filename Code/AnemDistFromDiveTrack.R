# Script to find the distance from an anemone where a tagged fish was seen previously to the dive track for each sampling season
# Will be input into mark-recapture analysis to help with recapture estimates (ideally, to capture whether or not we even had a reasonable chance of catching the fish)
# Should think about how to best mesh this with AnemLocations.R and ClownfishMarkRecap.R - which ones should be loaded into others? Or should they all be part of one big script?

#################### Set-up: ####################
#Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(RMark)
library(lubridate)
#library(dbplyr)
library(ggplot2)
library(here)

# Load some data files while deciding how various scripts should interact/communicate
# Load output from AnemLocations.R (or source that file) - mean gps coordinates across multiple observations of anemones
load(file=here("Data",'AnemAllInfowLatLon2.RData')) #file with anems, after 2018 anems matched
load(file=here("Data", "AnemLatLonObsbyAnem.RData")) #file with lat-lon info appended

# Load encounter history data frame from ClownfishMarkRecap.R
load(file=here("Data", "encounters_all.RData"))

# Set a few things
sample_years = c(2015,2016,2017,2018)

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

#Go through table, match anem_id to entry in anem.Processed2, take lat/lon from there (should do this better, by using average lat/lons - partial work on this below)
addLatLons <- function(anem_vec, anemdf, latorlon) {
  out = rep(NA, length(anem_vec))
  
  for(i in 1:length(anem_vec)) {
    testid = anem_vec[i]
    
    if(!is.na(testid)) { #if the anem_id isn't NA
      matches = anemdf %>% filter(anem_id %in% testid) #see if there are matches
      if(length(matches$dive_table_id) != 0) {
        if(latorlon == "lat") {
          out[i] = matches$lat[1]
        } else if (latorlon == "lon") {
          out[i] = matches$lon[1]
        }
      } else if (length(matches$dive_table_id) == 0) {
        print(paste("No id matches for",testid, sep=" "))
      }
    }
  }
  return(out)
}

# Attempted function to get mean lat/lon of all anem observations rather than just the one from the year in consideration - mostly an issue w/2018 anems b/c not have anem_obs yet
# Could probably just add something to the other function (addLatLon) that pulls out anem_id_unq2 too, then could replace the anems in fish.Tagged with those and match up with anem.LatLon
# # Go through table, if there is an anem_id, match it to anem_obs or anem_id from the past (if 2018), so can link up with anem_id_unq2 in anem.LatLon
# updateAnemIds <- function(anem_vec, anem.Info, year) {
#   
#   out = rep(NA, length(anem_vec)) #output
#   
#   for(i in 1:length(anem_vec)) {
#     testid = anem_vec[i]
#     
#     if(year == 2018) {
#       
#     }
#     if(!is.na(testid)) { #if the anem_id in question is not NA
#       matches = anem.Info %>% filter(anem_id %in% testid) #see if there are any matches in anem_id
#       
#       if(length(matches$anem_table_id) != 0) { #if something comes up
#         if(!is.na(matches$anems_obs[1])) {
#           out[i] = paste("obs", matches$anem_obs[1], sep="") #make the id that matches the format in anem_id_unq2 if there is an anem_obs
#         } else {
#           out[i] = paste("id", testid, sep="") #if not an anem_obs, use the anem_id
#         }
#         
#       } else {
#         matches2 = anem.Info %>% filter(old_anem_id %in% testid) #see if any matches in old_anem_id
#       } if(length(matches2$anem_table_id) != 0) { #if something comes up
#         out[i] = paste("obs", matches2$anem_obs[2], sep="") #make the id that matches the format in anem_id_unq2
#       } else {
#         print(paste("Found no anem_id or old_anem_id matches with anem_obs for",testid))
#       }
#     } 
#     
#     
#     
#     matches = anem.Info %>% filter(anem_id %in% testid)
#     
#   }
# }
# for (i in 1:length())


#################### Running things: ####################
# Pull out info from data base
leyte <- read_db("Leyte")

fish.Info <- leyte %>% 
  tbl("clownfish") %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, cap_id, anem_table_id, recap, tag_id, color, size) %>%
  collect() %>%
  filter(!is.na(tag_id)) 

anem.Info <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, obs_time) %>%
  collect() #%>%
  #filter(anem_table_id %in% fish.Info$anem_table_id)

dive.Info <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, dive_type, date, site, gps) %>%
  collect() %>%
  mutate(year = as.integer(substring(date,1,4)))#%>%
  #filter(dive_table_id %in% allfish_anems$dive_table_id) %>%
  
gps.Info <- leyte %>%
  tbl("GPX") %>%
  select(lat, lon, time, unit) %>%
  collect(n = Inf) %>%
  mutate(obs_time = force_tz(ymd_hms(time), tzone = "UTC")) %>% #tell it that it is in UTC time zone
  mutate(month = month(obs_time), #and separate out useful components of the time (this and line above largely from Michelle's assign_db_gpx function)
         day = day(obs_time), 
         hour = hour(obs_time), 
         min = minute(obs_time), 
         sec = second(obs_time), 
         year = year(obs_time)) %>%
  separate(time, into = c("date", "time"), sep = " ") #pull out date separately as a chr string too

#merge fish, anem, and dive info into one data frame - list of each encounter of a tagged fish
fish.AllInfo <- left_join(fish.Info, (anem.Info %>% filter(anem_table_id %in% fish.Info$anem_table_id)), by="anem_table_id")
fish.AllInfo <- left_join(fish.AllInfo, (dive.Info %>% filter(dive_table_id %in% fish.AllInfo$dive_table_id)), by="dive_table_id")

# Pull out tagged fish - create a data frame that has tag, season tagged, seasons caught, anem_ids (first and others)
fish.Tagged <- encounters_all %>% select(tag_id, ch, site) %>%
  mutate(anem_2015 = rep(NA, length(tag_id)), anem_2016 = rep(NA, length(tag_id)), anem_2017 = rep(NA, length(tag_id)), anem_2018 = rep(NA, length(tag_id))) %>%
  mutate(year_tagged = rep(NA, length(tag_id)))

for(i in 1:length(fish.Tagged$tag_id)) {
  tag = fish.Tagged$tag_id[i]
  c2015 = as.integer(substring(fish.Tagged$ch[i],1,1)) #1/0 values for whether or not the fish was caught each sampling year
  c2016 = as.integer(substring(fish.Tagged$ch[i],2,2)) 
  c2017 = as.integer(substring(fish.Tagged$ch[i],3,3)) 
  c2018 = as.integer(substring(fish.Tagged$ch[i],4,4)) 
  year_tag = sample_years[min(which(c(c2015,c2016,c2017,c2018) > 0))] #gives year fish tagged
  fish.Tagged$year_tagged[i] = year_tag
  
  if(length((fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id) > 1) { print(paste("Multiple matching anems for fish",tag,"in year",year_tag))}
  
  if(year_tag == 2015){
    if(length((fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id) > 1) { print(paste("Multiple matching anems for fish",tag,"in year",year_tag))}
    fish.Tagged$anem_2015[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
    fish.Tagged$anem_2016[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2016))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2016))$anem_id, fish.Tagged$anem_2015[i])
    fish.Tagged$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id, fish.Tagged$anem_2016[i])
    fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
  } else if(year_tag == 2016){
    fish.Tagged$anem_2016[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
    fish.Tagged$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id, fish.Tagged$anem_2016[i])
    fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
  } else if(year_tag == 2017){
    fish.Tagged$anem_2017[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
    fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
  } else if(year_tag == 2018){
    fish.Tagged$anem_2018[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
  }
}

#Get about 50 print-outs of the message that multiple anems were found for a particular tag in a particular season. 
#Spot-checked a few and they had the same anem for each (usually 2nd capture was a recap dive) but should probably check more carefully at some point
#The above misses some anems (~14 in 2015) - there are 469 ch that start with 1 (like 1001, 1000), but only 455 anems in the anem_2015 vec - haven't looked into at all so don't know what's up

# Add lat/lon columns for each year
fish.Tagged <- fish.Tagged %>% 
  mutate(lat2015 = rep(NA, length(fish.Tagged$tag_id)), lon2015 = rep(NA, length(fish.Tagged$tag_id)), 
         lat2016 = rep(NA, length(fish.Tagged$tag_id)), lon2016 = rep(NA, length(fish.Tagged$tag_id)),
         lat2017 = rep(NA, length(fish.Tagged$tag_id)), lon2017 = rep(NA, length(fish.Tagged$tag_id)),
         lat2018 = rep(NA, length(fish.Tagged$tag_id)), lon2018 = rep(NA, length(fish.Tagged$tag_id)))

# Fill in lat/lon columns
fish.Tagged$lat2015 <- addLatLons(fish.Tagged$anem_2015, anem.Processed2, "lat")
fish.Tagged$lon2015 <- addLatLons(fish.Tagged$anem_2015, anem.Processed2, "lon")
fish.Tagged$lat2016 <- addLatLons(fish.Tagged$anem_2016, anem.Processed2, "lat")
fish.Tagged$lon2016 <- addLatLons(fish.Tagged$anem_2016, anem.Processed2, "lon")
fish.Tagged$lat2017 <- addLatLons(fish.Tagged$anem_2017, anem.Processed2, "lat")
fish.Tagged$lon2017 <- addLatLons(fish.Tagged$anem_2017, anem.Processed2, "lon")
fish.Tagged$lat2018 <- addLatLons(fish.Tagged$anem_2018, anem.Processed2, "lat")
fish.Tagged$lon2018 <- addLatLons(fish.Tagged$anem_2018, anem.Processed2, "lon")


# Go through table, if fish was caught in a year, find coordinates of anem from that year, then find distance from dives for that year (will be 0 usually, right?)
# If year was before a fish was tagged, put in NA
# If fish not caught in a year, find average coordinates of anem, then find distance from that to the tracks from sample year
# Not sure if this would work best as a function (where put in year? or maybe site?) or as series of dplyr commands...



##############################################################################
############# Feeder code below this point ###################################
tail_color_2015 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  filter(year == 2015) %>%
  group_by(tag_id) %>%
  summarize(tail_color_2015 = color[1])

#process anem.AllInfo (filter out anem_ids that are NA, add in unique anem_id, format date/time and switch to UTC time zone, add placeholder lat lon columns)
anem.Processed <- anem.AllInfo %>%
  filter(!is.na(anem_id) | anem_id != "-9999" | anem_id == "") %>% #filter out NAs, -9999s (if any left in there still...), and blanks; 10881 with NAs left in, 4977 after filtering (previously, before Michelle added 2018 and redid database to take out anem observations that were actually clownfish processing, had 9853 rows when NAs left in, 4056 when filtered out)
  filter(is.null(anem_spp) == FALSE) %>% #to remove "phantom" anem observations that were actually just clownfish processing, get rid of obs with anem_spp that is null...
  filter(is.na(anem_spp) == FALSE) %>% #or anem_spp that is NA
  mutate(anem_id = as.numeric(anem_id)) %>% #make anem_id numeric to make future joins/merges easier
  mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), paste("obs", anem_obs, sep=""))) %>% #add unique anem id to anem.AllInfo
  mutate(lat = as.numeric(rep(NA, length(anem_table_id)))) %>% #add in placeholder columns for lat and lon info
  mutate(lon = as.numeric(rep(NA, length(anem_table_id)))) %>%
  mutate(obs_time = force_tz(ymd_hms(str_c(date, obs_time, sep = " ")), tzone = "Asia/Manila")) %>% #tell it that it is currently in Asia/Manila time zone
  mutate(obs_time = with_tz(obs_time, tzone = "UTC")) %>% #convert to UTC so can compare with GPS data (this line and one above largely from Michelle's assign_db_gpx function)
  mutate(month = month(obs_time), #and separate out useful components of the time (this also from Michelle's assign_db_gpx function)
         day = day(obs_time), 
         hour = hour(obs_time), 
         min = minute(obs_time), 
         sec = second(obs_time), 
         year = year(obs_time)) %>%
  mutate(anem_id_unq2 = anem_id_unq) #add another anem_id_unq to use for matching up the anems seen in 2018 (since anem_obs hasn't been run for the new data yet)


  
  anem.Processed <- anem.AllInfo %>%
  filter(!is.na(anem_id) | anem_id != "-9999" | anem_id == "") %>% #filter out NAs, -9999s (if any left in there still...), and blanks; 10881 with NAs left in, 4977 after filtering (previously, before Michelle added 2018 and redid database to take out anem observations that were actually clownfish processing, had 9853 rows when NAs left in, 4056 when filtered out)
  filter(is.null(anem_spp) == FALSE) %>% #to remove "phantom" anem observations that were actually just clownfish processing, get rid of obs with anem_spp that is null...
  filter(is.na(anem_spp) == FALSE) %>% #or anem_spp that is NA
  mutate(anem_id = as.numeric(anem_id)) %>% #make anem_id numeric to make future joins/merges easier
  mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), paste("obs", anem_obs, sep=""))) %>% #add unique anem id to anem.AllInfo
  mutate(lat = as.numeric(rep(NA, length(anem_table_id)))) %>% #add in placeholder columns for lat and lon info
  mutate(lon = as.numeric(rep(NA, length(anem_table_id)))) %>%
  mutate(obs_time = force_tz(ymd_hms(str_c(date, obs_time, sep = " ")), tzone = "Asia/Manila")) %>% #tell it that it is currently in Asia/Manila time zone
  mutate(obs_time = with_tz(obs_time, tzone = "UTC")) %>% #convert to UTC so can compare with GPS data (this line and one above largely from Michelle's assign_db_gpx function)
  mutate(month = month(obs_time), #and separate out useful components of the time (this also from Michelle's assign_db_gpx function)
         day = day(obs_time), 
         hour = hour(obs_time), 
         min = minute(obs_time), 
         sec = second(obs_time), 
         year = year(obs_time)) %>%
  mutate(anem_id_unq2 = anem_id_unq) #add another anem_id_unq to use for matching up the anems seen in 2018 (since anem_obs hasn't been run for the new data yet)

# pull out just the year and put that in a separate column
allfish_dives$year <- as.integer(substring(allfish_dives$date,1,4))

#join together
allfish <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
allfish <- left_join(allfish, allfish_dives, by="dive_table_id")

allfish$size <- as.numeric(allfish$size) #make size numeric (rather than a chr) so can do means and such



# Pull out tagged fish - create a data frame that has tag, season tagged, seasons caught, anem_ids (first and others)

# Go through table, if fish was caught in a year, find coordinates of anem from that year, then find distance from dives for that year (will be 0 usually, right?)
# If year was before a fish was tagged, put in NA
# If fish not caught in a year, find average coordinates of anem, then find distance from that to the tracks from sample year
# Not sure if this would work best as a function (where put in year? or maybe site?) or as series of dplyr commands...


#pull out all lat/lon info so can feed into function above (rather than having to pull out from database each time)
alllatlons <- leyte %>%
  tbl("GPX") %>%
  select(lat, lon, time, unit) %>%
  collect(n = Inf) %>% 
  separate(time, into = c("date", "time"), sep = " ") 

# Distance to anemone where fish was previously caught
# If fish was caught in that season, use the coordinates of the anem from that season
# If not caught in that season or for some reason anem doesn't have coordinates, use mean of coordinates from previous sightings 
# Means I need to get that mean anem location code working again, go over it to make sure it makes sense....
# Should probably also try a run where mean is always used instead of the season-specific location when it is caught, to make sure that doesn't make a big difference

# One issue: tracks aren't continuous, b/c only take readings every 15 seconds...
# One option: just calculate the distance to every reading for the dives at that sight (just anem dives or all? could try both...) and then take the minimum

# Also need to figure out how to actually get it into the MARK input in the right way

# Seems like output for here should be data frame with tag ID, years as columns and tag numbers and distance in each year as input (NA for seasons before a fish was tagged)

#################### Plots: ####################
