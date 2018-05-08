#Assess variance in anem gps location across visits
#File started: 11/30/17 (relevant code pulled from AnemVisitationTable.R)

#### To-dos:
# Check that lat/lons are pulled right (the last one on 1/14/15 for anem_obs 10 is different than what Michelle shows on the map)
# Could be that sometimes time - 8 hours is actually a different day -- should check for that possibility (not common but some)
# Check how I am calculating SD among measurements - is that actually the best way to do it? Would it be better to calculate average distance between points for the anem?

#################### Set-up: ####################
#Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(lubridate)
library(dbplyr)
library(ggplot2)
library(cowplot)
library(fields)
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
script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#function to find the lat and long for an anem_id (based off of Michelle's sample_latlon function)
anemid_latlong <- function(anem.table.id, latlondata) { #anem.table.id is one anem_table_id value, latlondata is table of GPX data from database (rather than making the function call it each time); will need to think a bit more clearly about how to handle different locations read for different visits to the same anem_id (or different with same anem_obs); for now, just letting every row in anem.Info get a lat-long
  
  #find the dive info and time for this anem observation
  dive <- leyte %>%
    tbl("anemones") %>%
    select(anem_table_id, obs_time, dive_table_id, anem_id) %>%
    collect() %>%
    filter(anem_table_id %in% anem.table.id)
  
  # find the date info and gps unit for this anem observation
  date <- leyte %>% 
    tbl("diveinfo") %>% 
    select(dive_table_id, date, gps, site) %>% 
    collect() %>% 
    filter(dive_table_id %in% dive$dive_table_id)
  
  #join with anem info, format obs time
  anem <- left_join(dive, date, by = "dive_table_id") %>% 
    separate(obs_time, into = c("hour", "minute", "second"), sep = ":") %>% #this line and the next directly from Michelle's code
    mutate(gpx_hour = as.numeric(hour) - 8)
  
  # find the lat long for this anem observation
  latloninfo <- latlondata %>%
    filter(date %in% anem$date) %>% 
    separate(time, into = c("hour", "minute", "second"), sep = ":") %>% 
    filter(as.numeric(hour) == anem$gpx_hour & as.numeric(minute) == anem$minute) 
  
  latloninfo$lat <- as.numeric(latloninfo$lat)
  latloninfo$lon <- as.numeric(latloninfo$lon)
  
  #often get multiple records for each anem_table_id (like if sit there for a while) - so multiple GPS points for same visit to an anemone, not differences across visits
  dups_lat <- which(duplicated(latloninfo$lat)) #vector of positions of duplicate values
  dups_lon <- which(duplicated(latloninfo$lon))
  
  #either take the mean of the lat/lon readings or the duplicated values, depending if there are duplicate points
  if(length(dups_lat) == 0) { #if all latitude points are different
    anem$lat <- round(mean(latloninfo$lat), digits = 5) #take the mean of the latitude values (digits = 5 b/c that is what Michelle had)
    anem$lon <- round(mean(latloninfo$lon), digits = 5) #take the mean of the longitude values
    #lat <- round(mean(latloninfo$lat), digits = 5) 
    #lon <- round(mean(latloninfo$lon), digits = 5)
  }else{
    anem$lat <- latloninfo$lat[dups_lat[1]] #if are duplicates, take the value of the first duplicated point
    anem$lon <- latloninfo$lon[dups_lon[1]]
    #lat <- latloninfo$lat[dups_lat[1]] #if are duplicates, take the value of the first duplicated point
    #lon <- latloninfo$lon[dups_lon[1]]
  }
  
  return(anem)
  
}

#function to see if rows with the same anem_obs also have the same site 
compareObsandSite <- function(df) {
  
}

#################### Running things! ####################
leyte <- read_db("Leyte") 

#pull out GPS, dive, anem info from database
alllatlons <- leyte %>%
  tbl("GPX") %>%
  select(lat, lon, time, unit) %>%
  collect(n = Inf) %>% 
  separate(time, into = c("date", "time"), sep = " ") 

dive.Info <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, date, dive_type, site) %>%
  collect() %>%
  mutate(year = as.integer(substring(date,1,4)))

anem.Info <- leyte %>%
  tbl("anemones") %>%
  select(dive_table_id, anem_table_id, anem_id, anem_obs, old_anem_id) %>%
  collect()

#merge into one dataframe by dive_table_id (to assign date, year, site, dive type to each anem), filter out anem_ids that are NA, add in unique anem id
anem.AllInfo <- left_join(anem.Info, dive.Info, by="dive_table_id") %>% 
  filter(!is.na(anem_id) | anem_id != "-9999" | anem_id == "") %>% #filter out NAs, -9999s (if any left in there still...), and blanks; 10881 with NAs left in, 4977 after filtering (previously, before Michelle added 2018 and redid database to take out anem observations that were actually clownfish processing, had 9853 rows when NAs left in, 4056 when filtered out)
  mutate(anem_id = as.numeric(anem_id)) %>% #make anem_id numeric to make future joins/merges easier
  mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), paste("obs", anem_obs, sep=""))) %>% #add unique anem id to anem.AllInfo
  mutate(lat = as.numeric(rep(NA, length(anem_table_id)))) %>% #add in placeholder columns for lat and lon info
  mutate(lon = as.numeric(rep(NA, length(anem_table_id))))






###### THE CODE BELOW PRODUCES THE RDATA FILES LOADED AT THE TOP (LOADED B/C TAKES SEVERAL MINUTES TO GENERATE)
###### THERE ARE A FEW LINES OF CODE AT THE BOTTOM OF THIS SECTION THAT FILTER OUT ANEMONES WITH LARGE SDs TO LOOK AT MORE CLOSELY

#find the lat/lon for each anem observation (could probably do this without a for loop but then the averaging across multiple gps records per observation got hairy when trying to do it all w/dplyr as vectors...)
anem.AllInfo$lat <- rep(NA, length(anem.AllInfo$anem_table_id))
anem.AllInfo$lon <- rep(NA, length(anem.AllInfo$anem_table_id))

#go through anem_table_ids and find the lat/lon for each anem observation (this is the step that takes a long time to do, prints row (out of 4044) so can track progress)
for (i in 1:length(anem.AllInfo$anem_table_id)) {
  outlatlon <- anemid_latlong(anem.AllInfo$anem_table_id[i], alllatlons)
  anem.AllInfo$lat[i] <- outlatlon$lat
  anem.AllInfo$lon[i] <- outlatlon$lon
  print(i)
}
#save output, since takes several minutes to run the above for loop
save(anem.AllInfo, file='AnemAllInfowLatLon.RData')

#pull out list of anem_id_unqs and their sites with only one row per each
anem.IDUnqSites <- distinct(anem.AllInfo[c("anem_id_unq", "site")]) 

#filter out NAs in lat/lons, then find the mean, variance, and sd of lat observations and lon observations by anem_id_unq
#use Bessel's correction (https://en.wikipedia.org/wiki/Bessel%27s_correction) for the variance calculation (use (n-1) in demon instead of n, here mulitply var by n/(n-1)) b/c not know the true mean of the distribution and the sample size is small so otherwise variance biased low
anem.LatLon <- anem.AllInfo %>% filter(!is.na(lat) & !is.na(lon)) %>% group_by(anem_id_unq) %>% summarize(meanlat = mean(lat), varlat = var(lat)*(n()/(n()-1)), sdlat = sqrt(varlat), meanlon = mean(lon), varlon = var(lon)*(n()/(n()-1)), sdlon = sqrt(varlon), ngps = n()) #need to multiply variance by N/N-1 (same as multiplying by 1/N-1 when doing the sum) b/c not know the mean of the dist, estimating it, and small sample size (per Douglas discussion Thanksgiving weekend 2017)
anem.LatLon <- left_join(anem.LatLon, anem.IDUnqSites, by = 'anem_id_unq') #add sites back in

#convert sd from degrees of lat and lon to meters
anem.LatLon$sdlat_m <- anem.LatLon$sdlat*mperlat #convert latitudes
anem.LatLon$sdlon_m <- anem.LatLon$sdlon*mperlat*cos(anem.LatLon$meanlat*2*pi/360) #convert longitudes (distance of a degree depends on latitude)

#filter out anem_obs = 7, which has an issue (shows up at two sites) 
anem.LatLon <- filter(anem.LatLon, anem_id_unq != "obs7")

#save data frame for ease in accessing in the future
save(anem.LatLon, file="AnemLatLonObsbyAnem.RData")

##### NOW BACK TO CODE THAT IS NOT LOADED AT THE TOP
#filter out just the sds less than a threshold value to zoom in on histogram for plotting
anem.LatLon.Abbr <- filter(anem.LatLon, sdlat_m < distance_thresh)

#### THIS CODE PULLS OUT DATA FROM THE DATABASE FOR THE ANEMS WITH HIGH SDs (>100m), SAVED IN DATAFRAMES THAT CAN BE LOADED AT THE TOP
##### What are the anems that have large sdlat_m and sdlon_m? 
largeSDanems_200 <- anem.LatLon %>% filter(sdlat_m > 200) #filter out anems with sds > 200m (16 of these)
largeSDanems_100 <- anem.LatLon %>% filter(sdlat_m > 100 & sdlat_m < 200) #filter out anems with sds between 100m and 200m (20 of these)

#look at them ones >200m individually to see if anything strange is going on
# obs1014, obs1028, obs1034, obs1044, obs1046, obs135, obs231, obs306, obs327, obs334, obs379, obs516, obs533, obs582, obs663, obs672
# all at Magbangon and Wangag
idtocheck <- 'obs672'
anem.AllInfo %>% filter(anem_id_unq == idtocheck)

#pull the data from the database for the ones >200 and save
#all of the anems with sdlat_m > 100m have anem_id_unq that start with obs so don't need to filter out for id (but might in future)
largeSDanems_200_anemobstopull <- largeSDanems_200 %>% mutate(idtopull = substr(anem_id_unq, 4, length(anem_id_unq))) #separate out just the anem_obs number from "obsXXX"
#largeSDanems_200_anemobstopull$idtopull <- as.numeric(largeSDanems_200_anemobstopull$idtopull) #convert from character to numeric

#now, pull info for those anems from the database to see if anything odd jumps out
largeSDanems_200_dbdata_anems <- leyte %>% #pull anem info
  tbl('anemones') %>%
  select(anem_table_id, obs_time, dive_table_id, anem_id, anem_obs, old_anem_id) %>%
  collect() %>%
  filter(anem_obs %in% largeSDanems_200_anemobstopull$idtopull)

largeSDanems_200_dbdata_dive <- leyte %>% #pull dive info for those anem visits
  tbl('diveinfo') %>%
  select(dive_table_id, date, gps, site, dive_type) %>%
  collect() %>%
  filter(dive_table_id %in% largeSDanems_200_dbdata_anems$dive_table_id)

largeSDanems_200_dbdata <- left_join(largeSDanems_200_dbdata_anems, largeSDanems_200_dbdata_dive, by="dive_table_id") #join the two together
largeSDanems_200_dbdata <- arrange(largeSDanems_200_dbdata, anem_obs) #arrange by anem_obs

#now pull the data from the database for the ones between 100 and 200m and save
#all of the anems with sdlat_m > 100m have anem_id_unq that start with obs so don't need to filter out for id (but might in future)
largeSDanems_100_anemobstopull <- largeSDanems_100 %>% mutate(idtopull = substr(anem_id_unq, 4, length(anem_id_unq))) #separate out just the anem_obs number from "obsXXX"
#largeSDanems_100_anemobstopull$idtopull <- as.numeric(largeSDanems_100_anemobstopull$idtopull) #convert from character to numeric

#now, pull info for those anems from the database to see if anything odd jumps out
largeSDanems_100_dbdata_anems <- leyte %>% #pull anem info
  tbl('anemones') %>%
  select(anem_table_id, obs_time, dive_table_id, anem_id, anem_obs, old_anem_id) %>%
  collect() %>%
  filter(anem_obs %in% largeSDanems_100_anemobstopull$idtopull)

largeSDanems_100_dbdata_dive <- leyte %>% #pull dive info for those anem visits
  tbl('diveinfo') %>%
  select(dive_table_id, date, gps, site, dive_type) %>%
  collect() %>%
  filter(dive_table_id %in% largeSDanems_100_dbdata_anems$dive_table_id)

largeSDanems_100_dbdata <- left_join(largeSDanems_100_dbdata_anems, largeSDanems_100_dbdata_dive, by="dive_table_id") #join the two together
largeSDanems_100_dbdata <- arrange(largeSDanems_100_dbdata, anem_obs) #arrange by anem_obs

#save these two dataframes for future checking
save(largeSDanems_100_dbdata, file="largeSDanems_100_dbdata.RData")
save(largeSDanems_200_dbdata, file="largeSDanems_200_dbdata.RData")

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
pdf(file="AnemGPS_Estimates.pdf")
ggplot(data=anem.LatLon, aes(sdlat_m, sdlon_m)) +
  geom_point(aes(color=site, size=ngps)) +
  #facet_grid(.~site_wrap, labeller=label_parsed) +
  xlab("latitude sd (m)") + ylab("longitude sd (m)") + ggtitle("Standard deviation of anem positions across visits") +
  theme_bw()
dev.off()

#histogram of sdlat_m less than distance_thresh (here 50m)
pdf(file="AnemGPS_Estimates_Hist.pdf")
hist(anem.LatLon.Abbr$sdlat_m, breaks=50, xlab='Standard deviation of lat (m)', main=paste('Histogram of sd of lat (m) values < ', distance_thresh, sep=""))
dev.off()

#histogram of all sdlat_m
pdf(file="AnemGPS_Estimates_Hist_All.pdf")
hist(anem.LatLon$sdlat_m, breaks=50, xlab='Standard deviation of lat (m)', main=paste('Histogram of sd of lat (m) values'))
dev.off()
