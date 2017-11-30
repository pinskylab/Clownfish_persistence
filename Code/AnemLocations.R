#Assess variance in anem gps location across visits
#File started: 11/30/17 (relevant code pulled from AnemVisitationTable.R)

rm(list=ls())

#################### Set-up: ####################
setwd("~/Box Sync/Rutgers postdoc/Clownfish_persistence/Code") #set working directory

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

#Load input files
load(file="AnemAllInfowLatLon.RData") #gives anem.AllInfo w/lat and lons averaged for each anem visit
load(file="AnemLatLonObsbyAnem.RData") #gives anem.LatLon w/sd and var of lats and lons across vists to each anem 

#Set constants
mperlat <- 111.111*1000 #meters per degree of latitude
distance_thresh <- 50 #distance threshold for zooming in on histogram of lat and lon standard devs plot

#Can disable diagnostics (or try updating RStudio?) if get tired of the many warnings about unitialized or unknown columns: https://stackoverflow.com/questions/39041115/fixing-a-multiple-warning-unknown-column

#################### Functions: ####################
#helper functions from Michelle's GitHub - do various tasks w/database (like assigning dates and site to fish and such)
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

#function to make vector of strings for column names for something done each year (like columns for sampling each year or minimum distance sampled to each anem each year, etc.)
makeYearlyColNames <- function(start.Year, end.Year, descriptor) { #start.Year is first year of sampling, end.Year is final year of sampling, descriptor is string to go before year in column name (like "min_dist_" if want columns like "min_dist_2012")
  
  if (start.Year > end.Year) {
    out <- NULL
  } 
  else {
    out <- as.vector(NA) #initalize output vector of column names
    year <- start.Year 
    
    for (i in 1:(end.Year - start.Year + 1)) {
      out[i] <- paste(descriptor, year, sep="") #create string like "min_dist_2012", where descriptor is something like "min_dist_" and year is the year sampled
      year <- year + 1
    }
  }
  return(out)
}

#################### Running things! ####################
leyte <- read_db("Leyte") 

#pull out all lat/lon info so can feed into function above (rather than having to pull out from database each time)
alllatlons <- leyte %>%
  tbl("GPX") %>%
  select(lat, lon, time, unit) %>%
  collect(n = Inf) %>% 
  separate(time, into = c("date", "time"), sep = " ") 

#the commented code below produces anem.AllInfo w/lat lon obs for each anem visit, currently commented out b/c can just load saved output file at the top
# #pull out relevant info from dive table and anemone table and merge so can match each anem_id to site and dates visited 
# dive.Info <- leyte %>% tbl("diveinfo") %>% select(date, dive_table_id, dive_type, site) %>% collect() #pull out list of dive dates, sites, and types
# dive.Info$year <- as.integer(substring(dive.Info$date, 1, 4)) #add a year column
# anem.Info <- leyte %>% tbl("anemones") %>% select(dive_table_id, anem_table_id, anem_id, anem_obs, old_anem_id) %>% collect() #pull out anem_id, dive_table_id, anem_table_id to match up with dive info
# 
# #merge into one dataframe by dive_table_id (to assign date, year, site, dive type to each anem), filter out anem_ids that are NA
# anem.AllInfo <- left_join(anem.Info, dive.Info, by="dive_table_id") %>% filter(!is.na(anem_id)) #9853 rows when NAs left in, 4056 when filtered out
# anem.AllInfo$anem_id <- as.numeric(anem.AllInfo$anem_id) #make anem_id numeric to make future joins/merges easier
# anem.AllInfo <- anem.AllInfo %>% filter(anem_id != "-9999") %>% collect() #remove -9999 anem_id
# anem.AllInfo <- anem.AllInfo %>% filter(anem_id != "") %>% collect() #remove blank anem_id
# 
# #grouping anem_ids by anem_obs values (so anemones we already know from tag surveys are actually the same one)
# anem.AllInfo <- anem.AllInfo %>% mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), anem_obs)) #add unique anem id to anem.AllInfo
# 
# #find the lat/lon for each anem observation (could probably do this without a for loop but then the averaging across multiple gps records per observation got hairy when trying to do it all w/dplyr as vectors...)
# anem.AllInfo$lat <- rep(NA, length(anem.AllInfo$anem_table_id))
# anem.AllInfo$lon <- rep(NA, length(anem.AllInfo$anem_table_id))
# 
# # this commented code produces the output file below that now I just load in to save time... creates on lat/lon observation per visit
# # #go through anem_table_ids and find the lat/lon for each anem observation 
# for (i in 1:length(anem.AllInfo$anem_table_id)) {
#   outlatlon <- anemid_latlong(anem.AllInfo$anem_table_id[i], alllatlons)
#   anem.AllInfo$lat[i] <- outlatlon$lat
#   anem.AllInfo$lon[i] <- outlatlon$lon
#   print(i)
# }
# # #save output, since takes several minutes to run the above for loop
# # #save(anem.AllInfo, file='AnemAllInfowLatLon.RData')

#pull out list of anem_id_unqs and their sites with only one row per each
anem.IDUnqSites <- distinct(anem.AllInfo[c("anem_id_unq", "site")]) 

#filter out NAs in lat/lons, then find the mean, variance, and sd of lat observations and lon observations by anem_id_unq
#use Bessel's correction (https://en.wikipedia.org/wiki/Bessel%27s_correction) for the variance calculation (use (n-1) in demon instead of n, here mulitply var by n/(n-1)) b/c not know the true mean of the distribution and the sample size is small so otherwise variance biased low
anem.LatLon <- anem.AllInfo %>% filter(!is.na(lat) & !is.na(lon)) %>% group_by(anem_id_unq) %>% summarize(meanlat = mean(lat), varlat = var(lat)*(n()/(n()-1)), sdlat = sqrt(varlat), meanlon = mean(lon), varlon = var(lon)*(n()/(n()-1)), sdlon = sqrt(varlon), ngps = n()) #need to multiply variance by N/N-1 (same as multiplying by 1/N-1 when doing the sum) b/c not know the mean of the dist, estimating it, and small sample size (per Douglas discussion Thanksgiving weekend 2017)
anem.LatLon <- left_join(anem.LatLon, anem.IDUnqSites, by = 'anem_id_unq') #add sites back in

#convert sd from degrees of lat and lon to meters
anem.LatLon$sdlat_m <- anem.LatLon$sdlat*mperlat #convert latitudes
anem.LatLon$sdlon_m <- anem.LatLon$sdlon*mperlat*cos(anem.LatLon$meanlat*2*pi/360) #convert longitudes (distance of a degree depends on latitude)

#filter out anem_id_unq = 7, which has an issue (shows up at two sites)
anem.LatLon <- filter(anem.LatLon, anem_id_unq != 7)
#save(anem.LatLon, file="AnemLatLonObsbyAnem.RData")

#filter out just the sds less than a threshold value to zoom in on histogram 
anem.LatLon.Abbr <- filter(anem.LatLon, sdlat_m < distance_thresh)

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
