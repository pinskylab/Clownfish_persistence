#Assess variance in anem gps location across visits
#File started: 11/30/17 (relevant code pulled from AnemVisitationTable.R)

#### To-dos:
# Check that lat/lons are pulled right (the last one on 1/14/15 for anem_obs 10 is different than what Michelle shows on the map)
# Check how I am calculating SD among measurements - is that actually the best way to do it? Would it be better to calculate average distance between points for the anem?

#################### Set-up: ####################
#Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
#library(dbplyr)
library(ggplot2)
#library(cowplot)
#library(fields)
library(here)

# #Load input files (code below creates these, just give the option to load here b/c takes several minutes to generate)
# load(file="AnemAllInfowLatLon.RData") #gives anem.AllInfo w/lat and lons averaged for each anem visit
# load(file="AnemLatLonObsbyAnem.RData") #gives anem.LatLon w/sd and var of lats and lons across vists to each anem 
# load(file="largeSDanems_100_dbdata.RData") #data from database of anems with SDlat between 100-200m
# load(file="largeSDanems_200_dbdata.RData") #data from database of anems with SDlat >200m

#Set constants
mperlat <- 111.111*1000 #meters per degree of latitude
distance_thresh <- 50 #distance threshold for zooming in on histogram of lat and lon standard devs plot

#Can disable diagnostics (or try updating RStudio?) if get tired of the many warnings about unitialized or unknown columns: https://stackoverflow.com/questions/39041115/fixing-a-multiple-warning-unknown-column

#################### Functions: ####################
# Functions and constants from my GitHub function/constant collection
script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Functions from Michelle's GitHub helpers script
#helper functions - do various tasks w/database (like assigning dates and site to fish and such)
#script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/helpers.R", ssl.verifypeer = FALSE)
#eval(parse(text = script))

#field_helpers (similar idea, in field repository) - this might be the newer version of the helpers collection?
script <- getURL("https://raw.githubusercontent.com/pinskylab/field/master/scripts/field_helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#function to find the lat and long for an anem_id (based off of Michelle's sample_latlon function)
anemid_latlong <- function(anem.table.id, anem.df, latlondata) { #anem.table.id is one anem_table_id value, anem.df is the anem.Processed data frame (so don't have to pull from db again here), latlondata is table of GPX data from database (rather than making the function call it each time); will need to think a bit more clearly about how to handle different locations read for different visits to the same anem_id (or different with same anem_obs); for now, just letting every row in anem.Info get a lat-long
  
  anem <- anem.df %>% filter(anem_table_id == anem.table.id) #get the relevant dive, time, site, etc. info for this anem_table_id
  
  # find the lat long for this anem observation
  latloninfo <- latlondata %>%
    filter(date %in% anem$date & unit == anem$gps) %>% #filter out just the GPS unit associated with this anem observation (added since previous time)
    filter(hour == anem$hour & min == anem$min) %>%
    mutate(lat = as.numeric(lat)) %>%
    mutate(lon = as.numeric(lon))
  
  #pull duplicates (so if sat in one place for enough time that multiple readings recorded there)
  #(there are more digits in the lats and lons than show up on the screen so sometimes things look like duplicates but aren't)
  dups_lat <- which(duplicated(latloninfo$lat)) #vector of positions of duplicate values 
  dups_lon <- which(duplicated(latloninfo$lon))
  
  #either take the mean of the lat/lon readings or the duplicated values, depending if there are duplicate points
  if(length(dups_lat) == 0) { #if all latitude points are different
    anem$lat <- round(mean(latloninfo$lat), digits = 5) #take the mean of the latitude values (digits = 5 b/c that is what Michelle had)
    anem$lon <- round(mean(latloninfo$lon), digits = 5) #take the mean of the longitude values
  }else{
    anem$lat <- latloninfo$lat[dups_lat[1]] #if are duplicates, take the value of the first duplicated point
    anem$lon <- latloninfo$lon[dups_lon[1]]
    print(paste("Dups in lat lons at anem_table_id", anem$anem_table_id, "on", anem$date, "with lat", anem$lat, sep = " ")) #just have this while trouble-shooting repeat entries in the database
  }
  
  return(anem)
  
}

# #function to see if rows with the same anem_obs also have the same site 
# compareObsandSite <- function(df) {
#   
# }

#################### Running things! ####################
leyte <- read_db("Leyte") 

#pull out GPS info from database, format + separate observation dates and times
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

#pull out dive info from database
dive.Info <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, date, dive_type, site, gps) %>%
  collect()

#pull out anem info
anem.Info <- leyte %>%
  tbl("anemones") %>%
  select(dive_table_id, anem_table_id, anem_id, anem_obs, old_anem_id, obs_time, anem_spp) %>%
  collect()

#merge dive and anem into one dataframe by dive_table_id (to assign date, year, site, dive type to each anem)
anem.AllInfo <- left_join(anem.Info, dive.Info, by="dive_table_id")

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
         year = year(obs_time))

#sometimes filtering issues so check to make sure various things gotten taken out
anemsppNA <- anem.Processed %>% filter(is.na(anem_spp)) #any nas in anem_spp?
if(length(anemsppNA$anem_table_id != 0)) {
  print("There are still NA anem_spp entries!")
}
anemsppNULL <- anem.Processed %>% filter(is.null(anem_spp)) #any nulls in anem_spp?
if(length(anemsppNULL$anem_table_id != 0)) {
  print("There are still NULL anem_spp entries!")
}

###### THE CODE BELOW PRODUCES THE RDATA FILES LOADED AT THE TOP (LOADED B/C TAKES SEVERAL MINUTES TO GENERATE)
#go through anem_table_ids and find the lat/lon for each anem observation (this is the step that takes a long time to do, prints row (out of 4977) so can track progress)
for (i in 1:length(anem.Processed$anem_table_id)) {
  outlatlon <- anemid_latlong(anem.Processed$anem_table_id[i], anem.Processed, gps.Info)
  anem.Processed$lat[i] <- outlatlon$lat
  anem.Processed$lon[i] <- outlatlon$lon
  print(i)
}

#save output, since takes several minutes to run the above for loop
save(anem.Processed, file=here("Data",'AnemAllInfowLatLon.RData'))

###### NOW BACK TO CODE THAT TAKES NO TIME TO RUN
#pull out list of anem_id_unqs and their sites with only one row per each
anem.IDUnqSites <- distinct(anem.Processed[c("anem_id_unq", "site")]) 

#filter out NAs in lat/lons, then find the mean, variance, and sd of lat observations and lon observations by anem_id_unq
#use Bessel's correction (https://en.wikipedia.org/wiki/Bessel%27s_correction) for the variance calculation (use (n-1) in demon instead of n, here mulitply var by n/(n-1)) b/c not know the true mean of the distribution and the sample size is small so otherwise variance biased low
anem.LatLon <- anem.Processed %>%
  filter(!is.na(lat) & !is.na(lon)) %>% #keep only non-NA lat and lon measurements
  group_by(anem_id_unq) %>% #group anems known to be the same through anem_id and anem_obs
  summarize(meanlat = mean(lat),
            varlat = var(lat)*(n()/(n()-1)),
            sdlat = sqrt(varlat),
            meanlon = mean(lon),
            varlon = var(lon)*(n()/(n()-1)),
            sdlon = sqrt(varlon),
            ngps = n()) #need to multiply variance by N/N-1 (same as multiplying by 1/N-1 when doing the sum) b/c not know the mean of the dist, estimating it, and small sample size (per Douglas discussion Thanksgiving weekend 2017)

#join sites back in
anem.LatLon <- left_join(anem.LatLon, anem.IDUnqSites, by = 'anem_id_unq') 

#convert sd from degrees of lat and lon to meters
anem.LatLon$sdlat_m <- anem.LatLon$sdlat*mperlat #convert latitudes
anem.LatLon$sdlon_m <- anem.LatLon$sdlon*mperlat*cos(anem.LatLon$meanlat*2*pi/360) #convert longitudes (distance of a degree depends on latitude)

#save data frame for ease in accessing in the future
save(anem.LatLon, file=here("Data", "AnemLatLonObsbyAnem.RData"))

##### NOW BACK TO CODE THAT IS NOT LOADED AT THE TOP
#filter out just the sds less than a threshold value to zoom in on histogram for plotting
anem.LatLon.Abbr <- anem.LatLon %>% filter(sdlat_m < distance_thresh & sdlon_m < distance_thresh)

#### THIS CODE PULLS OUT DATA FROM THE DATABASE FOR THE ANEMS WITH HIGH SDs (>100m), SAVED IN DATAFRAMES THAT CAN BE LOADED AT THE TOP
##### What are the anems that have large sdlat_m and sdlon_m? 
largeSDanems_200 <- anem.LatLon %>% filter(sdlat_m > 200) #filter out anems with sds > 200m (2 now)
largeSDanems_100 <- anem.LatLon %>% filter(sdlat_m > 100 & sdlat_m < 200) #filter out anems with sds between 100m and 200m (7 of these, was 20 before...)

##### Looking at the output and large SD anems in detail 
#look at them ones >200m individually to see if anything strange is going on
#obs24, obs30, obs307, obs330, obs423, obs513, obs523 (largeSDanems_100)
#obs1054, obs590

#go through these in detail
#id3225 - shows up at Tamakin Dacot and Haina on dives on 4/4/18
#obs24 - three observations, all at Palanas, in 2013, 2014, 2017 - coordinates are the same for 2013 and 2014 (anem_id 131), different for 2017 (anem_id 2628) - note says for 2017 dive, the dive_id was changed to 470 from 107 b/c dive was done across two sites (also Magbangon)
obs24 <- anem.Processed %>% filter(anem_obs == 24)
#obs30
obs30 <- anem.Processed %>% filter(anem_obs == 30)
#obs307 - two observations, one in 2017, one in 2015, at Wangag, anem_id 2064 noted for both
obs307 <- anem.Processed %>% filter(anem_obs == 307)
#obs330
obs330 <- anem.Processed %>% filter(anem_obs == 330)
#obs423
obs423 <- anem.Processed %>% filter(anem_obs == 423)
#obs513
obs513 <- anem.Processed %>% filter(anem_obs == 513)
#obs523
obs523 <- anem.Processed %>% filter(anem_obs == 523)
#obs1054
obs1054 <- anem.Processed %>% filter(anem_obs == 1054)
#obs590
obs590 <- anem.Processed %>% filter(anem_obs == 590)

#just checking stuff out... what about the anems that have a lot of observations? Do those look real or are phantom ones sneaking in?
table(anem.LatLon$ngps) #how many obs per anem are we getting? mostly 1, a couple with 7 (2), some with 6 (10)

#look at the ones with 7
anem.LatLon %>% filter(ngps == 7) #obs246 and obs97 - some have multiple obs within one dive (246 in 2017, both A), or seen in both 2015 seasons
anem.Processed %>% filter(anem_id_unq %in% (anem.LatLon %>% filter(ngps == 7))$anem_id_unq) %>% arrange(anem_obs, date)

#and the ones with 6
anem.LatLon %>% filter(ngps == 6)
anem.Processed %>% filter(anem_id_unq %in% (anem.LatLon %>% filter(ngps == 6))$anem_id_unq) %>% arrange(anem_obs, date)

#and the ones with 5
anem.LatLon %>% filter(ngps == 5)
anem.Processed %>% filter(anem_id_unq %in% (anem.LatLon %>% filter(ngps == 5))$anem_id_unq) %>% arrange(anem_obs, date)

#oldlist (from pre-2018 season when was working on this before)
# obs1014, obs1028, obs1034, obs1044, obs1046, obs135, obs231, obs306, obs327, obs334, obs379, obs516, obs533, obs582, obs663, obs672
# all at Magbangon and Wangag
# idtocheck <- 'obs672'
# anem.AllInfo %>% filter(anem_id_unq == idtocheck)
# 
# #pull the data from the database for the ones >200 and save
# #all of the anems with sdlat_m > 100m have anem_id_unq that start with obs so don't need to filter out for id (but might in future)
# largeSDanems_200_anemobstopull <- largeSDanems_200 %>% mutate(idtopull = substr(anem_id_unq, 4, length(anem_id_unq))) #separate out just the anem_obs number from "obsXXX"
# #largeSDanems_200_anemobstopull$idtopull <- as.numeric(largeSDanems_200_anemobstopull$idtopull) #convert from character to numeric
# 
# #now, pull info for those anems from the database to see if anything odd jumps out
# largeSDanems_200_dbdata_anems <- leyte %>% #pull anem info
#   tbl('anemones') %>%
#   select(anem_table_id, obs_time, dive_table_id, anem_id, anem_obs, old_anem_id) %>%
#   collect() %>%
#   filter(anem_obs %in% largeSDanems_200_anemobstopull$idtopull)
# 
# largeSDanems_200_dbdata_dive <- leyte %>% #pull dive info for those anem visits
#   tbl('diveinfo') %>%
#   select(dive_table_id, date, gps, site, dive_type) %>%
#   collect() %>%
#   filter(dive_table_id %in% largeSDanems_200_dbdata_anems$dive_table_id)
# 
# largeSDanems_200_dbdata <- left_join(largeSDanems_200_dbdata_anems, largeSDanems_200_dbdata_dive, by="dive_table_id") #join the two together
# largeSDanems_200_dbdata <- arrange(largeSDanems_200_dbdata, anem_obs) #arrange by anem_obs
# 
# #now pull the data from the database for the ones between 100 and 200m and save
# #all of the anems with sdlat_m > 100m have anem_id_unq that start with obs so don't need to filter out for id (but might in future)
# largeSDanems_100_anemobstopull <- largeSDanems_100 %>% mutate(idtopull = substr(anem_id_unq, 4, length(anem_id_unq))) #separate out just the anem_obs number from "obsXXX"
# #largeSDanems_100_anemobstopull$idtopull <- as.numeric(largeSDanems_100_anemobstopull$idtopull) #convert from character to numeric
# 
# #now, pull info for those anems from the database to see if anything odd jumps out
# largeSDanems_100_dbdata_anems <- leyte %>% #pull anem info
#   tbl('anemones') %>%
#   select(anem_table_id, obs_time, dive_table_id, anem_id, anem_obs, old_anem_id) %>%
#   collect() %>%
#   filter(anem_obs %in% largeSDanems_100_anemobstopull$idtopull)
# 
# largeSDanems_100_dbdata_dive <- leyte %>% #pull dive info for those anem visits
#   tbl('diveinfo') %>%
#   select(dive_table_id, date, gps, site, dive_type) %>%
#   collect() %>%
#   filter(dive_table_id %in% largeSDanems_100_dbdata_anems$dive_table_id)
# 
# largeSDanems_100_dbdata <- left_join(largeSDanems_100_dbdata_anems, largeSDanems_100_dbdata_dive, by="dive_table_id") #join the two together
# largeSDanems_100_dbdata <- arrange(largeSDanems_100_dbdata, anem_obs) #arrange by anem_obs
# 
# #save these two dataframes for future checking
# save(largeSDanems_100_dbdata, file="largeSDanems_100_dbdata.RData")
# save(largeSDanems_200_dbdata, file="largeSDanems_200_dbdata.RData")

##### Below was some code for checking to see if some anems were pulled from database incorrectly - not relevant now but don't want to delete yet (plot section is below)
# ##### Checking modified/synthesized/analyzed data against more raw data pulled straight from database - I think the whole issue here was that anem_id_unq was not always anem_id (so sometimes it looked like it was pulling the wrong thing if anem_obs was actually what it was using)
# # maybe something in here could get turned into a test? like checking that the time zones are handled appropriately (if that in fact turns out to be the problem)?
# #just pull out a couple of anem_ids to check
# # get anem info for just one anem_id
# ### There seems to be some confusion between whether a particular anem is coded by anem_id or anem_obs - 10 for anem_id is at Visca w/2 visits, 10 for anem_obs is at Wangag w/more visits
# ### Could some of those be getting combined to create anems with more visits than they actually have and across sites?
# ### Think more about how to check
# 
# # anem_id = 10: filtering from database, see two visits (5/22/13, 5/13/12), site is Visca; filtering in anem.LatLon or anem.AllInfo for anem_id_unq=10, see 7 visits, site is Wangag - has different anem_ids (2038, 890) and anem_obs is 10
# # also, anem_id = 10 has anem_obs = 297 (but don't see any other anem_ids associated with that anem_obs...)
# exanem = 10 # choose an anem_id to look at
# exanem <- c(2038, 890)
# anemidex <- leyte %>%
#   tbl("anemones") %>%
#   select(anem_table_id, obs_time, dive_table_id, anem_id, old_anem_id, anem_obs) %>%
#   collect() %>%
#   filter(anem_id %in% exanem)
# 
# anemidex <- leyte %>%
#   tbl("anemones") %>%
#   select(anem_table_id, obs_time, dive_table_id, anem_id, old_anem_id, anem_obs) %>%
#   collect() %>%
#   filter(anem_obs %in% exanem)
# 
# 
# # get the dive info for that anem_id
# anemiddiveinfoex <- leyte %>%
#   tbl("diveinfo") %>%
#   select(dive_table_id, date, gps, site) %>%
#   collect() %>%
#   filter(dive_table_id %in% anemidex$dive_table_id)
# 
# #and what does my anem.LatLon say?
# anemLatLonex <- anem.LatLon %>%
#   filter(anem_id_unq %in% exanem) 
# 
# #what does anem.AllInfo say?
# anemAllInfoex <- anem.AllInfo %>%
#   filter(anem_id_unq %in% exanem) #so the error is in anem.AllInfo - pulls visits from wrong site and too many
  
#################### Plots ####################
#scatter plot of all sds
pdf(file=here("Plots/AnemLocations", "AnemGPS_Estimates.pdf"))
ggplot(data=anem.LatLon, aes(sdlat_m, sdlon_m)) +
  geom_point(aes(color=site, size=ngps)) +
  #facet_grid(.~site_wrap, labeller=label_parsed) +
  xlab("latitude sd (m)") + ylab("longitude sd (m)") + ggtitle("Standard deviation of anem positions across visits") +
  theme_bw()
dev.off()

#scatter plot of just <50m sds 
pdf(file=here("Plots/AnemLocations", "AnemGPS_Estimates_Abbrv.pdf"))
ggplot(data=anem.LatLon.Abbr, aes(sdlat_m, sdlon_m)) +
  geom_point(aes(color=site, size=ngps)) +
  scale_x_continuous(limits = c(0, 50)) + scale_y_continuous(limits = c(0,50)) +
  #facet_grid(.~site_wrap, labeller=label_parsed) +
  xlab("latitude sd (m)") + ylab("longitude sd (m)") + ggtitle("Standard deviation of anem positions across visits for those <50m") +
  theme_bw()
dev.off()

#histogram of sdlat_m less than distance_thresh (here 50m)
pdf(file=here("Plots/AnemLocations","AnemGPS_Estimates_Hist.pdf"))
hist(anem.LatLon.Abbr$sdlat_m, breaks=50, xlab='Standard deviation of lat (m)', main=paste('Histogram of sd of lat (m) values < ', distance_thresh, sep=""))
dev.off()

#histogram of all sdlat_m
pdf(file=here("Plots/AnemLocations", "AnemGPS_Estimates_Hist_All.pdf"))
hist(anem.LatLon$sdlat_m, breaks=500, xlab='Standard deviation of lat (m)', main=paste('Histogram of sd of lat (m) values'))
dev.off()

#histogram of all ngps (number of times an anem visited)
pdf(file=here("Plots/AnemLocations","AnemVisits_Hist.pdf"))
hist(anem.LatLon$ngps, breaks=7, xlab='# visits with gps', main=paste('Histogram of times anems visited'))
dev.off()
