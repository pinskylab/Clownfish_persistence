# Script to find the distance from an anemone where a tagged fish was caught to the dive track for each sampling season
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


