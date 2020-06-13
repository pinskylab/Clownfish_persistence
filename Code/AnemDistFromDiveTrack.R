# Script to find the distance from an anemone where a tagged fish was seen previously to the dive track for each sampling season
# Takes encounter data produced in Clownfish_encounters.R, feeds into mark-recapture analysis in ClownfishMarkRecap.R

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))

# Load data frame with individual fish encounter histories and the one of marked fish with anemone + dive info (both from Clownfish_encounters.R)
load(file = here::here("Data/Script_outputs", "encounters_list.RData"))
load(file = here::here("Data/Script_outputs", "marked_fish.RData"))

# Load additional libraries
library(geosphere)
library(ggplot2)
library(cowplot)

#################### Functions: ####################
# Find the min distance from the anem where the fish was caught or last caught to a recorded GPS point in each year - from AnemDistFromDiveTrack.R
findDist <- function(sitevec, latvec, lonvec, capyearvec, gpxdf, divedf, sampleyear, site_visits){
  out <- rep(NA, length(sitevec))
  nsteps <- length(sitevec)
  
  for(i in 1:nsteps) {
    testsite = sitevec[i]
    testlat = as.numeric(latvec[i])
    testlon = as.numeric(lonvec[i])
    testcapyear = as.numeric(capyearvec[i])
    
    if(sampleyear >= testcapyear) {
      if(!is.na(testlat)) {
        visited <- site_visits$sampled[which(site_visits$site == testsite & site_visits$year == sampleyear)]  # pull out the sampled value from this year and site
        if(visited == 1) {  # if the sample site was visited that year, pull the dives from that site (for the spring months)
          relevantdives <- divedf %>%
            filter(site == testsite, as.numeric(year) == sampleyear, month %in% spring_months)
        } else {  # if not, pull dives from sites one position to the north and south
          site_N <- site_visits$site[which(site_visits$site == testsite & site_visits$year == sampleyear) - 1]  # site to the north
          site_S <- site_visits$site[which(site_visits$site == testsite & site_visits$year == sampleyear) + 1]  # site to the south
          
          relevantdives <- divedf %>%
            filter(site %in% c(site_N, site_S), as.numeric(year) == sampleyear, month %in% spring_months)
        }
        
        # Pull out gps points from the date(s) of those dives
        relevantgps <- gpxdf %>%
          filter(gps_date %in% relevantdives$dive_date)
        
        distvec <- rep(NA, length(relevantgps$lat))  # set a place to store distances
        
        for(j in 1:length(distvec)) {  # go through the list of coordinates and find the distance 
           #distvec2[j] = distGeo(c(testlon, testlat), c(relevantgps$lon[j], relevantgps$lat[j]))
          distvec[j] = distHaversine(c(testlon,testlat), c(as.numeric(relevantgps$lon[j]), as.numeric(relevantgps$lat[j])))
        }
        out[i] = min(distvec)  # select the minimum distance (closest a diver got to that anemone that year)
      }
    } else {
      out[i] = NA
    }
    print(paste(i, "out of", length(sitevec), "and", out[i]))
  } 
  return(out)
}

#################### Running things: ####################
# Save a new encounters data frame, add columns for anem_id, anem_obs, and anem_table_id for first capture anem 
encounters_dist <- encounters_list %>%
   mutate(capture_anem_id = NA,
          capture_anem_obs = NA,
          capture_anem_table_id = NA)

##### Find the anemone where each first was first caught and its coordinates
# Go through marked fish, find the anemone where they were first caught
for (i in 1:length(encounters_dist$fish_indiv)) {
  encounters_dist$capture_anem_id[i] = (marked_fish %>% filter(fish_indiv == encounters_dist$fish_indiv[i], year == encounters_dist$first_capture_year[i]))$anem_id[1] #record the anem_id where fish was first captured
  encounters_dist$capture_anem_obs[i] = (marked_fish %>% filter(fish_indiv == encounters_dist$fish_indiv[i], year == encounters_dist$first_capture_year[i]))$anem_obs[1]
  encounters_dist$capture_anem_table_id[i] = (marked_fish %>% filter(fish_indiv == encounters_dist$fish_indiv[i], year == encounters_dist$first_capture_year[i]))$anem_table_id[1]
}

# Go through and find lat lon for capture anems from anem_table_id (using anemid_latlong from Constants_database_common_functions)
capture_lat_vec <- rep(NA, length(encounters_dist$capture_anem_table_id))
capture_lon_vec <- rep(NA, length(encounters_dist$capture_anem_table_id))

for(i in 1:length(encounters_dist$capture_anem_table_id)) {
  anem_out <- anemid_latlong(encounters_dist$capture_anem_table_id[i], anems_Processed_all, gps_Info)
  capture_lat_vec[i] <- anem_out$lat
  capture_lon_vec[i] <- anem_out$lon
  print(i)
}

encounters_dist$capture_anem_lat <- capture_lat_vec
encounters_dist$capture_anem_lon <- capture_lon_vec

##### Find the distance from each anemone to the closest dive track in each sampling year
# replace the NA in site_visits with 0s
site_visits$sampled[is.na(site_visits$sampled)] <- 0  

# Calculate the distances from tracks in each year to the capture anem (NA for years prior to first capture year) - at some point, should check into why they're not all 0 or NA in 2012...
encounters_dist$dist2012 <- findDist(encounters_dist$site, encounters_dist$capture_anem_lat, encounters_dist$capture_anem_lon, encounters_dist$first_capture_year, gps_Info, dives_db_processed, 2012, site_visits)
encounters_dist$dist2013 <- findDist(encounters_dist$site, encounters_dist$capture_anem_lat, encounters_dist$capture_anem_lon, encounters_dist$first_capture_year, gps_Info, dives_db_processed, 2013, site_visits)
encounters_dist$dist2014 <- findDist(encounters_dist$site, encounters_dist$capture_anem_lat, encounters_dist$capture_anem_lon, encounters_dist$first_capture_year, gps_Info, dives_db_processed, 2014, site_visits)
encounters_dist$dist2015 <- findDist(encounters_dist$site, encounters_dist$capture_anem_lat, encounters_dist$capture_anem_lon, encounters_dist$first_capture_year, gps_Info, dives_db_processed, 2015, site_visits)
encounters_dist$dist2016 <- findDist(encounters_dist$site, encounters_dist$capture_anem_lat, encounters_dist$capture_anem_lon, encounters_dist$first_capture_year, gps_Info, dives_db_processed, 2016, site_visits)
encounters_dist$dist2017 <- findDist(encounters_dist$site, encounters_dist$capture_anem_lat, encounters_dist$capture_anem_lon, encounters_dist$first_capture_year, gps_Info, dives_db_processed, 2017, site_visits)
encounters_dist$dist2018 <- findDist(encounters_dist$site, encounters_dist$capture_anem_lat, encounters_dist$capture_anem_lon, encounters_dist$first_capture_year, gps_Info, dives_db_processed, 2018, site_visits)

# Save just the distances
min_survey_dist_to_anems <- encounters_dist %>% select(first_capture_year, capture_anem_id, capture_anem_obs, capture_anem_table_id, capture_anem_lat, capture_anem_lon,
                                                      dist2012, dist2013, dist2014, dist2015, dist2016, dist2017, dist2018)

##### Fill in the mean distance for pre-capture year anems (b/c cant have NAs in covariate columns)
mean2012 <- mean(encounters_dist$dist2012, na.rm=TRUE)  # do I need this?
mean2013 <- mean(encounters_dist$dist2013, na.rm=TRUE)
mean2014 <- mean(encounters_dist$dist2014, na.rm=TRUE)
mean2015 <- mean(encounters_dist$dist2015, na.rm=TRUE)
mean2016 <- mean(encounters_dist$dist2016, na.rm=TRUE)
mean2017 <- mean(encounters_dist$dist2017, na.rm=TRUE)
mean2018 <- mean(encounters_dist$dist2018, na.rm=TRUE)
mean_dist_allyears <- mean(c(encounters_dist$dist2012, encounters_dist$dist2013, encounters_dist$dist2014,
                             encounters_dist$dist2015, encounters_dist$dist2016, encounters_dist$dist2017,
                             encounters_dist$dist2018), na.rm=TRUE)

# Fill in the overall mean in one version and yearly mean in the other
encounters_dist_mean <- encounters_dist %>%
  mutate(dist2013 = replace(dist2013, is.na(dist2013), mean_dist_allyears),
         dist2014 = replace(dist2014, is.na(dist2014), mean_dist_allyears),
         dist2015 = replace(dist2015, is.na(dist2015), mean_dist_allyears),
         dist2016 = replace(dist2016, is.na(dist2016), mean_dist_allyears),
         dist2017 = replace(dist2017, is.na(dist2017), mean_dist_allyears),
         dist2018 = replace(dist2018, is.na(dist2018), mean_dist_allyears))

encounters_dist_mean_by_year <- encounters_dist %>%
  mutate(dist2013 = replace(dist2013, is.na(dist2013), mean2013),
         dist2014 = replace(dist2014, is.na(dist2014), mean2014),
         dist2015 = replace(dist2015, is.na(dist2015), mean2015),
         dist2016 = replace(dist2016, is.na(dist2016), mean2016),
         dist2017 = replace(dist2017, is.na(dist2017), mean2017),
         dist2018 = replace(dist2018, is.na(dist2018), mean2018))

#################### Plots: ####################

##### Histograms of distances from tracks to capture anems for tagged fish by year
d_2012_plot <- ggplot(data = encounters_dist, aes(dist2012)) +
  geom_histogram(aes(dist2012), bins=30) +
  xlab("distance (m)") + ggtitle("2012") +
  theme_bw()
    
d_2013_plot <- ggplot(data = encounters_dist, aes(dist2013)) +
  geom_histogram(aes(dist2013), bins=50) +
  xlab("distance (m)") + ggtitle("2013") +
  theme_bw()

d_2014_plot <- ggplot(data = encounters_dist, aes(dist2014)) +
  geom_histogram(aes(dist2014), bins=50) +
  xlab("distance (m)") + ggtitle("2014") +
  theme_bw()

d_2015_plot <- ggplot(data = encounters_dist, aes(dist2015)) +
  geom_histogram(aes(dist2015), bins=50) +
  xlab("distance (m)") + ggtitle("2015") +
  theme_bw()

d_2016_plot <- ggplot(data = encounters_dist, aes(dist2016)) +
  geom_histogram(aes(dist2016), bins=50) +
  xlab("distance (m)") + ggtitle("2016") +
  theme_bw()

d_2017_plot <- ggplot(data = encounters_dist, aes(dist2017)) +
  geom_histogram(aes(dist2017), bins=50) +
  xlab("distance (m)") + ggtitle("2017") +
  theme_bw()

d_2018_plot <- ggplot(data = encounters_dist, aes(dist2018)) +
  geom_histogram(aes(dist2018), bins=50) +
  xlab("distance (m)") + ggtitle("2018") +
  theme_bw()

pdf(file = here("Plots/DistFromTrackToAnems", "Distance_to_marked_fish_capture_anems_allyears.pdf"))
plot_grid(d_2012_plot, d_2013_plot, d_2014_plot, d_2015_plot, d_2016_plot, d_2017_plot, d_2018_plot, 
          labels = c("a","b","c","d","e","f","g"), ncol=2)
dev.off()

#################### Saving output: ####################
save(min_survey_dist_to_anems, file = here::here("Data/Script_outputs", "min_survey_dist_to_anems.RData"))  # data frame with anem_table_ids for capture anems plus their lat and lon and distance to sampling track in each year
save(encounters_dist, file = here::here("Data/Script_outputs", "encounters_dist.RData"))  # encounter histories with distances to dive tracks in each year
save(encounters_dist_mean, file = here::here("Data/Script_outputs", "encounters_dist_mean.RData"))  # encounter histories with overall mean dist filled in for NAs (pre-capture years, plus about 5 missing coords)
save(encounters_dist_mean_by_year, file = here::here("Data/Script_outputs", "encounters_dist_mean_by_year.RData"))  # encounter histories with mean dist by year filled in for NAs


# 
# 
# ##### IN THE PROCESS OF EDITING THIS FUNCTION, REALIZED anem_Processed doesn't have lat/lon, will need to run function in Constants_database_common_functions but updating it first
# 
# addLatLons <- function(encounters_df, anemdf) {
#   nentries = length(encounters_df$capture_anem_id)  # number of anems to go through, variable to use as short-hand length
#   out_lat = rep(NA, nentries)  # latitude output vector
#   out_lon = rep(NA, nentries)  # longitude output vector
#   
#   for(i in 1:nentries) {
#     test_obs = encounters_df$capture_anem_obs[i]  # anem_obs of the anem in question
#     test_table_id = encounters_df$capture_anem_table_id[i]  # anem_table_id of the anem in question
#     test_id = encounters_df$capture_anem_id[i]  # anem_id of the anem in question (if it has one)
#     
#     # try to get coordinates from the anem_obs first
#     if(!is.na(anem_obs)) {  # if the anem_obs isn't NA
#       matches = anemdf %>% filter(anem_obs == test_obs)  # see if there are any matches in the anem info data table
#       if(length(matches$dive_table_id != 0)) {  # if there are matches
#         out_lat[i] = mean(matches$lat, rm.na = TRUE)  # put the mean of the lat coordinate
#       }
#       
#     }
#     if(!is.na(test_id)) { #if the anem_id isn't NA
#       matches = anemdf %>% filter(anem_id == test_id) #see if there are any matches
#       if(length(matches$dive_table_id != 0)) { #if there are matches
#         if(latorlon == "lat") {
#           out[i] = matches$lat[1] #put the lat coord in if lat
#         } else if (latorlon == "lon") {
#           out[i] = matches$lon[1] #put the lon coord in if lon
#         }
#       }
#     }
#     
#     # check to see if got a coordinate from the anem_id, if not try the anem_table_id
#     if(is.na(out[i])) { #if the coordinate is still NA
#       if(!is.na(test_table_id)) { #and there is an anem_table_id
#         matches = anemdf %>% filter(anem_table_id == test_table_id) #see if there are any matches
#         if(length(matches$dive_table_id != 0)) { #if there are matches
#           if(latorlon == "lat") {
#             out[i] = matches$lat[1] #put the lat coord in if lat
#           } else if (latorlon == "lon") {
#             out[i] = matches$lon[1] #put the lon coord in if lon
#           }
#         }
#       }
#     }
#     
#     # check if either of those got a coordinate
#     if(is.na(out[i])) {
#       print(paste("No coordinates for anem_id",test_id,'or anem_table_id',test_table_id, sep=" "))
#     }
#   }
#   return(out)
# }
# 
# encounters_all$dist_2012 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2012, site_visits)
# 
# # Find the min distance from the anem where the fish was caught or last caught to a recorded GPS point in each year - from AnemDistFromDiveTrack.R
# findDist <- function(sitevec, latvec, lonvec, capyearvec, gpxdf, divedf, sampleyear, site_visits){
#   out <- rep(NA, length(sitevec))
#   nsteps <- length(sitevec)
#   
#   for(i in 1:nsteps) {
#     testsite = sitevec[i]
#     testlat = as.numeric(latvec[i])
#     testlon = as.numeric(lonvec[i])
#     testcapyear = as.numeric(capyearvec[i])
#     
#     if(sampleyear >= testcapyear) {
#       if(!is.na(testlat)) {
#         visited <- site_visits$sampled[which(site_visits$site == testsite & site_visits$year == sampleyear)] # pull out the sampled value from this year and site
#         if(visited == 1) { # if the sample site was visited that year, pull the dives from that site
#           relevantdives <- divedf %>%
#             filter(site == testsite, as.numeric(year) == sampleyear)
#         } else { #if not, pull dives from sites 1 to the north and south
#           site_N <- site_visits$site[which(site_visits$site == testsite & site_visits$year == sampleyear) - 1] #site to the north
#           site_S <- site_visits$site[which(site_visits$site == testsite & site_visits$year == sampleyear) + 1] #site to the south
#           
#           relevantdives <- divedf %>%
#             filter(site %in% c(site_N, site_S), as.numeric(year) == sampleyear)
#         }
#         
#         #filter out gps points from those dives
#         # relevantgps <- gpxdf %>%
#         #   filter(gps_date %in% relevantdives$date & gps_year == sampleyear)
#         relevantgps <- gpxdf %>%
#           filter(gps_date %in% relevantdives$dive_date)
#         
#         distvec <- rep(NA, length(relevantgps$lat)) #set a place to store distances
#         #distvec2 <- rep(NA, length(relevantgps$lat))
#         
#         for(j in 1:length(distvec)) { #go through the list of coordinates
#           #distvec2[j] = distGeo(c(testlon, testlat), c(relevantgps$lon[j], relevantgps$lat[j]))
#           distvec[j] = distHaversine(c(testlon,testlat), c(as.numeric(relevantgps$lon[j]), as.numeric(relevantgps$lat[j])))
#           #print(paste(j, "of", length(distvec)))
#         }
#         out[i] = min(distvec)
#         #print(min(distvec2))
#       }
#     } else {
#       out[i] = NA
#     }
#     print(paste(i, "out of", length(sitevec), "and", out[i]))
#   } 
#   return(out)
# }
# 
# 
# ##### Find distance from anemone in previous years
# encounters_all <- encounters_all %>%
#   mutate(capture_anem_id = NA,
#          capture_anem_table_id = NA,
#          capture_anem_lat = NA,
#          capture_anem_lon = NA,
#          dist_2012 = NA,
#          dist_2013 = NA,
#          dist_2014 = NA,
#          dist_2015 = NA,
#          dist_2016 = NA,
#          dist_2017 = NA,
#          dist_2018 = NA)
# 
# 
# ######################## OLD CODE ######################### Previous version of script 
# 
# # Go through list of tagged fish, find distance from anem where first caught for all years after capture (originally, updated anem each year to where it had been caught (code in AmenDistFromDiveTrack) but doesn't that bias fish that have been recaught to have lower distances?)
# for(i in 1:length(encounters_all$fish_indiv)) {
#   
#   # #Check for multiple anems recorded in capture year, print message if find some - this just records mulitple observations, not actually different anem_ids... amend to fix that...
#   # if(length((allfish_mark %>% filter(fish_id == encounters_all$fish_id[i], year == encounters_all$first_capture_year[i]))$anem_id) > 1) { 
#   #   mult_anems_test <- (allfish_mark %>% filter(fish_id == encounters_all$fish_id[i], year == encounters_all$first_capture_year[i]))$anem_id
#   #   for(l in 1:length(mult_anems_test))
#   # 
#   #   print(paste("Multiple matching anems for fish",encounters_all$fish_id[i],"in year",encounters_all$first_capture_year[i]))
#   # }
#   
#   # Record capture anem (id and anem_table_id) in first capture year and find lat/lon coordinates
#   encounters_all$capture_anem_id[i] = (marked_fish %>% filter(fish_indiv == encounters_all$fish_indiv[i], year == encounters_all$first_capture_year[i]))$anem_id[1] #record the anem where fish was first captured
#   encounters_all$capture_anem_table_id[i] = (marked_fish %>% filter(fish_indiv == encounters_all$fish_indiv[i], year == encounters_all$first_capture_year[i]))$anem_table_id[1]
# }
# 
# # Now fill in the lat/lon coordinates for the capture anem
# #encounters_all$capture_anem_lat <- addLatLons(encounters_all, anems_Processed, "lat")
# encounters_all$capture_anem_lat <- addLatLons(encounters_all, anem.Processed_full, "lat") #now 5 anems without coordinates (which seems much more reasonable) and same ones for both lat and lon: anem_table_id 363, 364, 365, 366, 8094, none with anem_id
# encounters_all$capture_anem_lon <- addLatLons(encounters_all, anem.Processed_full, "lon")
# # talked to Michelle about these ones - probably from when the battery died at one point (2012) or the gps unit was submerged (2016)
# 
# # replace the NA in site_visits with 0s
# site_visits$sampled[is.na(site_visits$sampled)] <- 0  
# 
# ## NEED TO RUN THESE OVERNIGHT - TOO SLOW TO DO DURING THE DAY!
# # Calculate the distances from tracks in each year to the capture anem (NA for years prior to first capture year) - at some point, should check into why they're not all 0 or NA in 2012...
# encounters_all$dist_2012 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2012, site_visits)
# encounters_all$dist_2013 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2013, site_visits)
# encounters_all$dist_2014 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2014, site_visits)
# encounters_all$dist_2015 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2015, site_visits)
# encounters_all$dist_2016 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2016, site_visits)
# encounters_all$dist_2017 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2017, site_visits)
# encounters_all$dist_2018 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2018, site_visits)
# 
# addLatLons <- function(encounters_df, anemdf, latorlon) {
#   nentries = length(encounters_df$capture_anem_id) #number of anems to go through, variable to use as short-hand length
#   out = rep(NA, nentries) #output vector
#   
#   for(i in 1:nentries) {
#     test_table_id = encounters_df$capture_anem_table_id[i] #anem_table_id of the anem in question
#     test_id = encounters_df$capture_anem_id[i] #anem_id of the anem in question (if it has one)
#     # 
#     # # see if can get a lat/lon from the anem_id
#     # if(is.na(test_id)) { #if the anem_id is NA
#     #   if(is.na(test_table_id)) { #and the anem_table_id is NA
#     #     print('Both anem_id',)
#     #   }
#     # }
#     # 
#     
#     if(!is.na(test_id)) { #if the anem_id isn't NA
#       matches = anemdf %>% filter(anem_id == test_id) #see if there are any matches
#       if(length(matches$dive_table_id != 0)) { #if there are matches
#         if(latorlon == "lat") {
#           out[i] = matches$lat[1] #put the lat coord in if lat
#         } else if (latorlon == "lon") {
#           out[i] = matches$lon[1] #put the lon coord in if lon
#         }
#       }
#     }
#     
#     # check to see if got a coordinate from the anem_id, if not try the anem_table_id
#     if(is.na(out[i])) { #if the coordinate is still NA
#       if(!is.na(test_table_id)) { #and there is an anem_table_id
#         matches = anemdf %>% filter(anem_table_id == test_table_id) #see if there are any matches
#         if(length(matches$dive_table_id != 0)) { #if there are matches
#           if(latorlon == "lat") {
#             out[i] = matches$lat[1] #put the lat coord in if lat
#           } else if (latorlon == "lon") {
#             out[i] = matches$lon[1] #put the lon coord in if lon
#           }
#         }
#       }
#     }
#     
#     # check if either of those got a coordinate
#     if(is.na(out[i])) {
#       print(paste("No coordinates for anem_id",test_id,'or anem_table_id',test_table_id, sep=" "))
#     }
#   }
#   return(out)
# }
# 
# encounters_all$dist_2012 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2012, site_visits)
# 
# # Find the min distance from the anem where the fish was caught or last caught to a recorded GPS point in each year - from AnemDistFromDiveTrack.R
# findDist <- function(sitevec, latvec, lonvec, capyearvec, gpxdf, divedf, sampleyear, site_visits){
#   out <- rep(NA, length(sitevec))
#   nsteps <- length(sitevec)
#   
#   for(i in 1:nsteps) {
#     testsite = sitevec[i]
#     testlat = as.numeric(latvec[i])
#     testlon = as.numeric(lonvec[i])
#     testcapyear = as.numeric(capyearvec[i])
#     
#     if(sampleyear >= testcapyear) {
#       if(!is.na(testlat)) {
#         visited <- site_visits$sampled[which(site_visits$site == testsite & site_visits$year == sampleyear)] # pull out the sampled value from this year and site
#         if(visited == 1) { # if the sample site was visited that year, pull the dives from that site
#           relevantdives <- divedf %>%
#             filter(site == testsite, as.numeric(year) == sampleyear)
#         } else { #if not, pull dives from sites 1 to the north and south
#           site_N <- site_visits$site[which(site_visits$site == testsite & site_visits$year == sampleyear) - 1] #site to the north
#           site_S <- site_visits$site[which(site_visits$site == testsite & site_visits$year == sampleyear) + 1] #site to the south
#           
#           relevantdives <- divedf %>%
#             filter(site %in% c(site_N, site_S), as.numeric(year) == sampleyear)
#         }
#         
#         #filter out gps points from those dives
#         # relevantgps <- gpxdf %>%
#         #   filter(gps_date %in% relevantdives$date & gps_year == sampleyear)
#         relevantgps <- gpxdf %>%
#           filter(gps_date %in% relevantdives$dive_date)
#         
#         distvec <- rep(NA, length(relevantgps$lat)) #set a place to store distances
#         #distvec2 <- rep(NA, length(relevantgps$lat))
#         
#         for(j in 1:length(distvec)) { #go through the list of coordinates
#           #distvec2[j] = distGeo(c(testlon, testlat), c(relevantgps$lon[j], relevantgps$lat[j]))
#           distvec[j] = distHaversine(c(testlon,testlat), c(as.numeric(relevantgps$lon[j]), as.numeric(relevantgps$lat[j])))
#           #print(paste(j, "of", length(distvec)))
#         }
#         out[i] = min(distvec)
#         #print(min(distvec2))
#       }
#     } else {
#       out[i] = NA
#     }
#     print(paste(i, "out of", length(sitevec), "and", out[i]))
#   } 
#   return(out)
# }
# 
# #test <- findDist(encounters_all$site[1:10], encounters_all$capture_anem_lat[1:10], encounters_all$capture_anem_lon[1:10], encounters_all$first_capture_year[1:10], gps_Info, dives_db, 2012, site_visits)
# #save(encounters_all, file=here('Data','encounters_all.RData'))
# saveRDS(encounters_all, file=here('Data','encounters_all_with_dist.RData'))
# 
# # save just distances
# min_survey_dist_to_anems <- encounters_all %>% select(first_capture_year, capture_anem_id, capture_anem_table_id, capture_anem_lat, capture_anem_lon,
#                                                       dist_2012, dist_2013, dist_2014, dist_2015, dist_2016, dist_2017, dist_2018)
# 
# # Will be input into mark-recapture analysis to help with recapture estimates (ideally, to capture whether or not we even had a reasonable chance of catching the fish)
# # Should think about how to best mesh this with AnemLocations.R and ClownfishMarkRecap.R - which ones should be loaded into others? Or should they all be part of one big script?
# 
# 
# #################### Set-up: ####################
# #Load relevant libraries
# library(RCurl) #allows running R scripts from GitHub
# library(RMySQL) #might need to load this to connect to the database?
# library(dplyr)
# library(tidyr)
# #library(RMark)
# library(lubridate)
# library(geosphere)
# #library(dbplyr)
# library(ggplot2)
# library(cowplot)
# library(here)
# 
# # Load some data files while deciding how various scripts should interact/communicate
# # Load output from AnemLocations.R (or source that file) - mean gps coordinates across multiple observations of anemones
# load(file=here("Data",'AnemAllInfowLatLon2.RData')) #file with anems, after 2018 anems matched
# load(file=here("Data", "AnemLatLonObsbyAnem.RData")) #file with lat-lon info appended
# 
# # Load encounter history data frame from ClownfishMarkRecap.R
# load(file=here("Data", "encounters_all.RData"))
# 
# # Set a few things
# sample_years = c(2012,2013,2014,2015,2016,2017,2018)
# large_dist = 500 #threshold large distance for checking into anemsß
# 
# #################### Functions: ####################
# # Functions and constants from my GitHub function/constant collection
# # script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQAXExRItr2tnepmkRr7NRt4hylZrks5aciBtwA%3D%3D", ssl.verifypeer = FALSE)
# # eval(parse(text = script))
# script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
# eval(parse(text = script))
# 
# # Functions from Michelle's GitHub helpers script
# #field_helpers (similar idea to helpers, in field repository) - this might be the newer version of the helpers collection?
# script <- getURL("https://raw.githubusercontent.com/pinskylab/field/master/scripts/field_helpers.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
# 
# #Go through table, match anem_id to entry in anem.Processed2, take lat/lon from there (should do this better, by using average lat/lons - partial work on this below)
# addLatLons <- function(anem_vec, anemdf, latorlon) {
#   out = rep(NA, length(anem_vec))
#   
#   for(i in 1:length(anem_vec)) {
#     testid = anem_vec[i]
#     
#     if(!is.na(testid)) { #if the anem_id isn't NA
#       matches = anemdf %>% filter(anem_id %in% testid) #see if there are matches
#       if(length(matches$dive_table_id) != 0) {
#         if(latorlon == "lat") {
#           out[i] = matches$lat[1]
#         } else if (latorlon == "lon") {
#           out[i] = matches$lon[1]
#         }
#       } else if (length(matches$dive_table_id) == 0) {
#         print(paste("No id matches for",testid, sep=" "))
#       }
#     }
#   }
#   return(out)
# }
# 
# # Find the min distance from the anem where the fish was caught or last caught to a recorded GPS point in each year
# findDist <- function(sitevec, latvec, lonvec, gpxdf, divedf, sampleyear){
#   out <- rep(NA, length(sitevec))
#   nsteps <- length(sitevec)
#   
#   for(i in 1:nsteps) {
#     testsite = sitevec[i]
#     testlat = as.numeric(latvec[i])
#     testlon = as.numeric(lonvec[i])
#     
#     if(!is.na(testlat)) {
#       #filter out dives at the right site in the right year
#       relevantdives <- divedf %>% 
#         filter(site == testsite, year == sampleyear)
#       
#       #filter out gps points from those dives
#       relevantgps <- gpxdf %>%
#         filter(date %in% relevantdives$date, year == sampleyear)
#       
#       distvec <- rep(NA, length(relevantgps$lat)) #set a place to store distances
#       
#       for(j in 1:length(distvec)) { #go through the list of coordinates
#         distvec[j] = distHaversine(c(testlon,testlat), c(as.numeric(relevantgps$lon[j]), as.numeric(relevantgps$lat[j])))
#         #print(paste(j, "of", length(distvec)))
#       }
#       out[i] = min(distvec)
#     }
#     print(paste(i, "and", out[i]))
#   } 
#   return(out)
# }
# 
# # Attempted function to get mean lat/lon of all anem observations rather than just the one from the year in consideration - mostly an issue w/2018 anems b/c not have anem_obs yet
# # Could probably just add something to the other function (addLatLon) that pulls out anem_id_unq2 too, then could replace the anems in fish.Tagged with those and match up with anem.LatLon
# # # Go through table, if there is an anem_id, match it to anem_obs or anem_id from the past (if 2018), so can link up with anem_id_unq2 in anem.LatLon
# # updateAnemIds <- function(anem_vec, anem.Info, year) {
# #   
# #   out = rep(NA, length(anem_vec)) #output
# #   
# #   for(i in 1:length(anem_vec)) {
# #     testid = anem_vec[i]
# #     
# #     if(year == 2018) {
# #       
# #     }
# #     if(!is.na(testid)) { #if the anem_id in question is not NA
# #       matches = anem.Info %>% filter(anem_id %in% testid) #see if there are any matches in anem_id
# #       
# #       if(length(matches$anem_table_id) != 0) { #if something comes up
# #         if(!is.na(matches$anems_obs[1])) {
# #           out[i] = paste("obs", matches$anem_obs[1], sep="") #make the id that matches the format in anem_id_unq2 if there is an anem_obs
# #         } else {
# #           out[i] = paste("id", testid, sep="") #if not an anem_obs, use the anem_id
# #         }
# #         
# #       } else {
# #         matches2 = anem.Info %>% filter(old_anem_id %in% testid) #see if any matches in old_anem_id
# #       } if(length(matches2$anem_table_id) != 0) { #if something comes up
# #         out[i] = paste("obs", matches2$anem_obs[2], sep="") #make the id that matches the format in anem_id_unq2
# #       } else {
# #         print(paste("Found no anem_id or old_anem_id matches with anem_obs for",testid))
# #       }
# #     } 
# #     
# #     
# #     
# #     matches = anem.Info %>% filter(anem_id %in% testid)
# #     
# #   }
# # }
# # for (i in 1:length())
# 
# 
# #################### Running things: ####################
# # Pull out info from data base
# leyte <- read_db("Leyte")
# 
# fish.Info <- leyte %>% 
#   tbl("clownfish") %>%
#   select(fish_table_id, anem_table_id, fish_spp, sample_id, gen_id, anem_table_id, recap, tag_id, color, size) %>%
#   collect() %>%
#   #filter(!is.na(tag_id) | !is.na(sample_id) | !is.na(cap_id)) 
#   filter(!is.na(tag_id) | !is.na(sample_id) | !is.na(gen_id))
# 
# anem.Info <- leyte %>%
#   tbl("anemones") %>%
#   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, obs_time) %>%
#   collect() #%>%
#   #filter(anem_table_id %in% fish.Info$anem_table_id)
# 
# dive.Info <- leyte %>%
#   tbl("diveinfo") %>%
#   select(dive_table_id, dive_type, date, site, gps) %>%
#   collect() %>%
#   mutate(year = as.integer(substring(date,1,4)))#%>%
#   #filter(dive_table_id %in% allfish_anems$dive_table_id) %>%
#   
# gps.Info <- leyte %>%
#   tbl("GPX") %>%
#   select(lat, lon, time, unit) %>%
#   collect(n = Inf) %>%
#   mutate(obs_time = force_tz(ymd_hms(time), tzone = "UTC")) %>% #tell it that it is in UTC time zone
#   mutate(month = month(obs_time), #and separate out useful components of the time (this and line above largely from Michelle's assign_db_gpx function)
#          day = day(obs_time), 
#          hour = hour(obs_time), 
#          min = minute(obs_time), 
#          sec = second(obs_time), 
#          year = year(obs_time)) %>%
#   separate(time, into = c("date", "time"), sep = " ") #pull out date separately as a chr string too
# 
# #merge fish, anem, and dive info into one data frame - list of each encounter of a tagged fish
# fish.AllInfo <- left_join(fish.Info, (anem.Info %>% filter(anem_table_id %in% fish.Info$anem_table_id)), by="anem_table_id")
# fish.AllInfo <- left_join(fish.AllInfo, (dive.Info %>% filter(dive_table_id %in% fish.AllInfo$dive_table_id)), by="dive_table_id")
# 
# # Pull out tagged fish - create a data frame that has tag, season tagged, seasons caught, anem_ids (first and others)
# fish.Tagged <- encounters_all %>% select(tag_id, ch, site) %>%
#   mutate(anem_2015 = rep(NA, length(tag_id)), anem_2016 = rep(NA, length(tag_id)), anem_2017 = rep(NA, length(tag_id)), anem_2018 = rep(NA, length(tag_id))) %>%
#   mutate(year_tagged = rep(NA, length(tag_id)))
# 
# # Go through the list of tagged fish, if it could have been caught in a year, assign it the anem where it was caught or most recently caught
# for(i in 1:length(fish.Tagged$tag_id)) {
#   tag = fish.Tagged$tag_id[i]
#   c2015 = as.integer(substring(fish.Tagged$ch[i],1,1)) #1/0 values for whether or not the fish was caught each sampling year
#   c2016 = as.integer(substring(fish.Tagged$ch[i],2,2)) 
#   c2017 = as.integer(substring(fish.Tagged$ch[i],3,3)) 
#   c2018 = as.integer(substring(fish.Tagged$ch[i],4,4)) 
#   year_tag = sample_years[min(which(c(c2015,c2016,c2017,c2018) > 0))] #gives year fish tagged
#   fish.Tagged$year_tagged[i] = year_tag
#   
#   if(length((fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id) > 1) { print(paste("Multiple matching anems for fish",tag,"in year",year_tag))}
#   
#   if(year_tag == 2015){
#     if(length((fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id) > 1) { print(paste("Multiple matching anems for fish",tag,"in year",year_tag))}
#     fish.Tagged$anem_2015[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#     fish.Tagged$anem_2016[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2016))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2016))$anem_id, fish.Tagged$anem_2015[i])
#     fish.Tagged$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id, fish.Tagged$anem_2016[i])
#     fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
#   } else if(year_tag == 2016){
#     fish.Tagged$anem_2016[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#     fish.Tagged$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id, fish.Tagged$anem_2016[i])
#     fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
#   } else if(year_tag == 2017){
#     fish.Tagged$anem_2017[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#     fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
#   } else if(year_tag == 2018){
#     fish.Tagged$anem_2018[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#   }
# }
# 
# #Get about 50 print-outs of the message that multiple anems were found for a particular tag in a particular season. 
# #Spot-checked a few and they had the same anem for each (usually 2nd capture was a recap dive) but should probably check more carefully at some point
# #The above misses some anems (~14 in 2015) - there are 469 ch that start with 1 (like 1001, 1000), but only 455 anems in the anem_2015 vec - haven't looked into at all so don't know what's up
# 
# # Add lat/lon columns for each year
# fish.Tagged <- fish.Tagged %>% 
#   mutate(lat2015 = rep(NA, length(fish.Tagged$tag_id)), lon2015 = rep(NA, length(fish.Tagged$tag_id)), 
#          lat2016 = rep(NA, length(fish.Tagged$tag_id)), lon2016 = rep(NA, length(fish.Tagged$tag_id)),
#          lat2017 = rep(NA, length(fish.Tagged$tag_id)), lon2017 = rep(NA, length(fish.Tagged$tag_id)),
#          lat2018 = rep(NA, length(fish.Tagged$tag_id)), lon2018 = rep(NA, length(fish.Tagged$tag_id)))
# 
# # Fill in lat/lon columns - right now, just fills in the lat/lon for that anem observation, not the mean observation where old_anem_id and anem_obs info is taken into account (which would be in anem.LatLon)
# fish.Tagged$lat2015 <- addLatLons(fish.Tagged$anem_2015, anem.Processed2, "lat")
# fish.Tagged$lon2015 <- addLatLons(fish.Tagged$anem_2015, anem.Processed2, "lon")
# fish.Tagged$lat2016 <- addLatLons(fish.Tagged$anem_2016, anem.Processed2, "lat")
# fish.Tagged$lon2016 <- addLatLons(fish.Tagged$anem_2016, anem.Processed2, "lon")
# fish.Tagged$lat2017 <- addLatLons(fish.Tagged$anem_2017, anem.Processed2, "lat")
# fish.Tagged$lon2017 <- addLatLons(fish.Tagged$anem_2017, anem.Processed2, "lon")
# fish.Tagged$lat2018 <- addLatLons(fish.Tagged$anem_2018, anem.Processed2, "lat")
# fish.Tagged$lon2018 <- addLatLons(fish.Tagged$anem_2018, anem.Processed2, "lon")
# 
# # Find min distance (in m) from anem to a recorded gps point from that year (for now, no sorting by dive type or anything like that, just site and year)
# dist2015 <- findDist(fish.Tagged$site, fish.Tagged$lat2015, fish.Tagged$lon2015, gps.Info, dive.Info, 2015)
# dist2016 <- findDist(fish.Tagged$site, fish.Tagged$lat2016, fish.Tagged$lon2016, gps.Info, dive.Info, 2016)
# dist2017 <- findDist(fish.Tagged$site, fish.Tagged$lat2017, fish.Tagged$lon2017, gps.Info, dive.Info, 2017)
# dist2018 <- findDist(fish.Tagged$site, fish.Tagged$lat2018, fish.Tagged$lon2018, gps.Info, dive.Info, 2018)
# 
# # Add distance columns to fish.Tagged
# fish.Tagged <- fish.Tagged %>%
#   mutate(dist_2015 = dist2015, dist_2016 = dist2016, dist_2017 = dist2017, dist_2018 = dist2018)
# 
# ##### Check distances that seem quite large
# # Find the ones that are big outliers to check into a bit more
# far_2015 <- fish.Tagged %>% filter(dist_2015 >= large_dist) #should't be any of these (and there aren't!)
# far_2016 <- fish.Tagged %>% filter(dist_2016 >= large_dist) #none of these either
# far_2017 <- fish.Tagged %>% filter(dist_2017 >= large_dist) #none of these either
# far_2018 <- fish.Tagged %>% filter(dist_2018 >= large_dist) #one of these - this is the fish that moved from Wangag to Hicgop South!
# 
# # Filter out just that fish, in case want to do some runs without it
# fish.Tagged_nomovement <- fish.Tagged %>% filter(tag_id != 982000411818588)
# 
# 
# #################### Plots: ####################
# # Histogram of distances from tracks to anems with tagged fish
# d_2015 <- ggplot(data = fish.Tagged, aes(dist2015)) +
#   geom_histogram(aes(dist2015), bins=30) +
#   xlab("distance (m) in 2015") + ggtitle("2015 track distance to anems") +
#   theme_bw()
# d_2016 <- ggplot(data = fish.Tagged, aes(dist2016)) +
#   geom_histogram(aes(dist2016), bins=30) +
#   xlab("distance (m) in 2016") + ggtitle("2016 track distance to anems") +
#   theme_bw()
# d_2017 <- ggplot(data = fish.Tagged, aes(dist2017)) +
#   geom_histogram(aes(dist2017), bins=100) +
#   xlab("distance (m) in 2017") + ggtitle("2017 track distance to anems") +
#   theme_bw()
# d_2018 <- ggplot(data = fish.Tagged, aes(dist2018)) +
#   geom_histogram(aes(dist2018), bins=200) +
#   xlab("distance (m) in 2018") + ggtitle("2018 track distance to anems") +
#   theme_bw()
# 
# pdf(file = here("Plots/DistFromTrackToAnems", "Distance_to_tagged_fish_anems_allyears.pdf"))
# plot_grid(d_2015, d_2016, d_2017, d_2018, labels = "AUTO")
# dev.off()
# 
# # # Can't quite get this to run but seems like it might be close...
# # ggplot(data = fish.Tagged, aes(dist_2015,dist_2016,dist_2017,dist_2018)) +
# #   geom_histogram(aes(dist_2015), bins = 15) +
# #   geom_histogram(aes(dist_2016), bins = 30) +
# #   geom_histogram(aes(dist_2017), bins = 50) +
# #   geom_histogram(aes(dist_2017), bins = 200) +
# #   facet_grid(.~variable) +
# #   ggtitle("Distance (m) from anem where fish was previously caught") +
# #   theme_bw()
# 
# # Histograms of distances from tracks to anems with tagged fish, each year separately
# hist(fish.Tagged$dist_2015, xlab="distance (m) from track to fish anem", main="2015 distances")
# hist(fish.Tagged$dist_2016, xlab="distance (m) from track to fish anem", main="2016 distances")
# hist(fish.Tagged$dist_2017, breaks=100, xlab="distance (m) from track to fish anem", main="2017 distances")
# hist(fish.Tagged$dist_2018, xlab="distance (m) from track to fish anem", main="2018 distances")
# 
# # Histogram of distances from tracks to anems with tagged fish, without the fish that moved from Wangag (2017) to Hicgop South (2018)
# d_2015_2 <- ggplot(data = fish.Tagged_nomovement, aes(dist_2015)) +
#   geom_histogram(aes(dist_2015), bins=50) +
#   xlab("distance (m) in 2015") + ggtitle("2015 track distance to anems") +
#   theme_bw()
# d_2016_2 <- ggplot(data = fish.Tagged_nomovement, aes(dist_2016)) +
#   geom_histogram(bins=50) +
#   xlab("distance (m) in 2016") + ggtitle("2016 track distance to anems") +
#   theme_bw()
# d_2017_2 <- ggplot(data = fish.Tagged_nomovement, aes(dist_2017)) +
#   geom_histogram(aes(dist_2017), bins=50) +
#   xlab("distance (m) in 2017") + ggtitle("2017 track distance to anems \nno Wangag-HS fish") +
#   theme_bw()
# d_2018_2 <- ggplot(data = fish.Tagged_nomovement, aes(dist_2018)) +
#   geom_histogram(aes(dist_2018), bins=50) +
#   xlab("distance (m) in 2018") + ggtitle("2018 track distance to anems \nno Wangag-HS fish") +
#   theme_bw()
# 
# pdf(file = here("Plots/DistFromTrackToAnems", "Distance_to_tagged_fish_anems_allyears_withoutmovingfish.pdf"))
# plot_grid(d_2015_2, d_2016_2, d_2017_2, d_2018_2, labels = "AUTO")
# dev.off()
# 
# 
# #################### Saving output: ####################
# save(dist2015, file=here("Data", "dist2015.RData"))
# save(dist2016, file=here("Data", "dist2016.RData"))
# save(dist2017, file=here("Data", "dist2017.RData"))
# save(dist2018, file=here("Data", "dist2018.RData"))
# save(fish.Tagged, file=here("Data", "fish_Tagged.RData"))
# save(fish.Tagged_nomovement, file=here("Data", "fish_tagged_nomovement.RData"))
# 
# 
# ##############################################################################
# ############# Feeder code below this point ###################################
# # Go through table, if fish was caught in a year, find coordinates of anem from that year, then find distance from dives for that year (will be 0 usually, right?)
# # If year was before a fish was tagged, put in NA
# # If fish not caught in a year, find average coordinates of anem, then find distance from that to the tracks from sample year
# # Not sure if this would work best as a function (where put in year? or maybe site?) or as series of dplyr commands...
# 
# # tail_color_2015 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
# #   filter(year == 2015) %>%
# #   group_by(tag_id) %>%
# #   summarize(tail_color_2015 = color[1])
# # 
# # #process anem.AllInfo (filter out anem_ids that are NA, add in unique anem_id, format date/time and switch to UTC time zone, add placeholder lat lon columns)
# # anem.Processed <- anem.AllInfo %>%
# #   filter(!is.na(anem_id) | anem_id != "-9999" | anem_id == "") %>% #filter out NAs, -9999s (if any left in there still...), and blanks; 10881 with NAs left in, 4977 after filtering (previously, before Michelle added 2018 and redid database to take out anem observations that were actually clownfish processing, had 9853 rows when NAs left in, 4056 when filtered out)
# #   filter(is.null(anem_spp) == FALSE) %>% #to remove "phantom" anem observations that were actually just clownfish processing, get rid of obs with anem_spp that is null...
# #   filter(is.na(anem_spp) == FALSE) %>% #or anem_spp that is NA
# #   mutate(anem_id = as.numeric(anem_id)) %>% #make anem_id numeric to make future joins/merges easier
# #   mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), paste("obs", anem_obs, sep=""))) %>% #add unique anem id to anem.AllInfo
# #   mutate(lat = as.numeric(rep(NA, length(anem_table_id)))) %>% #add in placeholder columns for lat and lon info
# #   mutate(lon = as.numeric(rep(NA, length(anem_table_id)))) %>%
# #   mutate(obs_time = force_tz(ymd_hms(str_c(date, obs_time, sep = " ")), tzone = "Asia/Manila")) %>% #tell it that it is currently in Asia/Manila time zone
# #   mutate(obs_time = with_tz(obs_time, tzone = "UTC")) %>% #convert to UTC so can compare with GPS data (this line and one above largely from Michelle's assign_db_gpx function)
# #   mutate(month = month(obs_time), #and separate out useful components of the time (this also from Michelle's assign_db_gpx function)
# #          day = day(obs_time), 
# #          hour = hour(obs_time), 
# #          min = minute(obs_time), 
# #          sec = second(obs_time), 
# #          year = year(obs_time)) %>%
# #   mutate(anem_id_unq2 = anem_id_unq) #add another anem_id_unq to use for matching up the anems seen in 2018 (since anem_obs hasn't been run for the new data yet)
# # 
# # 
# #   
# #   anem.Processed <- anem.AllInfo %>%
# #   filter(!is.na(anem_id) | anem_id != "-9999" | anem_id == "") %>% #filter out NAs, -9999s (if any left in there still...), and blanks; 10881 with NAs left in, 4977 after filtering (previously, before Michelle added 2018 and redid database to take out anem observations that were actually clownfish processing, had 9853 rows when NAs left in, 4056 when filtered out)
# #   filter(is.null(anem_spp) == FALSE) %>% #to remove "phantom" anem observations that were actually just clownfish processing, get rid of obs with anem_spp that is null...
# #   filter(is.na(anem_spp) == FALSE) %>% #or anem_spp that is NA
# #   mutate(anem_id = as.numeric(anem_id)) %>% #make anem_id numeric to make future joins/merges easier
# #   mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), paste("obs", anem_obs, sep=""))) %>% #add unique anem id to anem.AllInfo
# #   mutate(lat = as.numeric(rep(NA, length(anem_table_id)))) %>% #add in placeholder columns for lat and lon info
# #   mutate(lon = as.numeric(rep(NA, length(anem_table_id)))) %>%
# #   mutate(obs_time = force_tz(ymd_hms(str_c(date, obs_time, sep = " ")), tzone = "Asia/Manila")) %>% #tell it that it is currently in Asia/Manila time zone
# #   mutate(obs_time = with_tz(obs_time, tzone = "UTC")) %>% #convert to UTC so can compare with GPS data (this line and one above largely from Michelle's assign_db_gpx function)
# #   mutate(month = month(obs_time), #and separate out useful components of the time (this also from Michelle's assign_db_gpx function)
# #          day = day(obs_time), 
# #          hour = hour(obs_time), 
# #          min = minute(obs_time), 
# #          sec = second(obs_time), 
# #          year = year(obs_time)) %>%
# #   mutate(anem_id_unq2 = anem_id_unq) #add another anem_id_unq to use for matching up the anems seen in 2018 (since anem_obs hasn't been run for the new data yet)
# # 
# # # pull out just the year and put that in a separate column
# # allfish_dives$year <- as.integer(substring(allfish_dives$date,1,4))
# # 
# # #join together
# # allfish <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
# # allfish <- left_join(allfish, allfish_dives, by="dive_table_id")
# # 
# # allfish$size <- as.numeric(allfish$size) #make size numeric (rather than a chr) so can do means and such
# # 
# # 
# # 
# # # Pull out tagged fish - create a data frame that has tag, season tagged, seasons caught, anem_ids (first and others)
# # 
# # # Go through table, if fish was caught in a year, find coordinates of anem from that year, then find distance from dives for that year (will be 0 usually, right?)
# # # If year was before a fish was tagged, put in NA
# # # If fish not caught in a year, find average coordinates of anem, then find distance from that to the tracks from sample year
# # # Not sure if this would work best as a function (where put in year? or maybe site?) or as series of dplyr commands...
# # 
# # 
# # #pull out all lat/lon info so can feed into function above (rather than having to pull out from database each time)
# # alllatlons <- leyte %>%
# #   tbl("GPX") %>%
# #   select(lat, lon, time, unit) %>%
# #   collect(n = Inf) %>% 
# #   separate(time, into = c("date", "time"), sep = " ") 
# 
# # Distance to anemone where fish was previously caught
# # If fish was caught in that season, use the coordinates of the anem from that season
# # If not caught in that season or for some reason anem doesn't have coordinates, use mean of coordinates from previous sightings 
# # Means I need to get that mean anem location code working again, go over it to make sure it makes sense....
# # Should probably also try a run where mean is always used instead of the season-specific location when it is caught, to make sure that doesn't make a big difference
# 
# # One issue: tracks aren't continuous, b/c only take readings every 15 seconds...
# # One option: just calculate the distance to every reading for the dives at that sight (just anem dives or all? could try both...) and then take the minimum
# 
# # Also need to figure out how to actually get it into the MARK input in the right way
# 
# # Seems like output for here should be data frame with tag ID, years as columns and tag numbers and distance in each year as input (NA for seasons before a fish was tagged)
# 
# # pdf(file=here("Plots/AnemLocations","AnemVisits_Hist.pdf"))
# # hist(anem.LatLon$ngps, breaks=7, xlab='# visits with gps', main=paste('Histogram of times anems visited'))
# # dev.off()
# # 
# # # Comparing C and A dives (plotted by year), adding in D,E,M for 2012
# # pdf(file = here("Plots", "DistanceFromFishAnen_allsites_byyear.pdf"))
# # ggplot(data = fish.Tagged, aes(lat2015)) +
# #   geom_histogram(data = fish.Tagged, alpha=0.5, binwidth=1) +
# #   facet_grid(.~year)
# # 
# # ggplot(data = (fishInfo %>% filter(dive_type %in% c("A","C","D","E","M"))), aes(size_num, fill = dive_type)) +
# #   geom_histogram(data = (fishInfo %>% filter(dive_type == "C")), alpha = 0.5, binwidth = 1) +
# #   geom_histogram(data = (fishInfo %>% filter(dive_type == "A")), alpha = 0.5, binwidth = 1) +
# #   geom_histogram(data = (fishInfo %>% filter(dive_type == "D")), alpha = 0.5, binwidth = 1) +
# #   geom_histogram(data = (fishInfo %>% filter(dive_type == "E")), alpha = 0.5, binwidth = 1) +
# #   geom_histogram(data = (fishInfo %>% filter(dive_type == "M")), alpha = 0.5, binwidth = 1) +
# #   facet_grid(.~ year) +
# #   xlab("size (cm)") + ylab("# fish") + ggtitle("Size histograms for fish from A,C,D,E,M dives (all sites combined)") +
# #   theme_bw()
# # dev.off()
# # 
# # 
