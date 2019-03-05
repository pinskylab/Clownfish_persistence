#Useful functions and constants for clownfish work - copied from Common_constants_and_functions.R in Clownfish_data_analysis 10/15/18
#Database pulls, constants, functions I use throughout the scripts in this repository

#################### Install packages ####################
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(geosphere)
library(varhandle)
library(reshape)
library(here)

#################### Constants, indices, etc. ####################
##### Vector and indices to reference sites easily in plots and filtering and such (can call things like site_vec[Cabatoan] or site_vec[1]))
site_vec <- c("Cabatoan", "Caridad Cemetery", "Caridad Proper", "Elementary School", "Gabas", "Haina", "Hicgop South",
              "N. Magbangon", "S. Magbangon", "Palanas", "Poroc Rose", "Poroc San Flower", "San Agustin", "Sitio Baybayon", "Sitio Lonas", 
              "Sitio Tugas", "Tamakin Dacot", "Visca", "Wangag")
Cabatoan <- 1
CaridadCemetery <- 2
CaridadProper <- 3
ElementarySchool <- 4
Gabas <- 5
Haina <- 6
HicgopSouth <- 7
N_Magbangon <- 8
S_Magbangon <- 9
Palanas <- 10
PorocRose <- 11
PorocSanFlower <- 12
SanAgustin <- 13
SitioBaybayon <- 14
SitioLonas <- 15
SitioTugas <- 16
TamakinDacot <- 17
Visca <- 18
Wangag <- 19

#vector of sites from north to south
# this one has Sitio Hicgop in in, which as far as I can tell doesn't have any fish associated to it in the db
# site_vec_NS <- c('Palanas', 'Wangag', 'N. Magbangon', 'S. Magbangon' , 'Cabatoan',
#                  'Caridad Cemetery', 'Caridad Proper', 'Sitio Hicgop', 'Hicgop South',
#                  'Sitio Tugas', 'Elementary School', 'Sitio Lonas', 'San Agustin',
#                  'Poroc San Flower', 'Poroc Rose', 'Visca', 'Gabas', 'Tamakin Dacot',
#                  'Haina', 'Sitio Baybayon')

site_vec_NS <- c('Palanas', 'Wangag', 'N. Magbangon', 'S. Magbangon' , 'Cabatoan',
                 'Caridad Cemetery', 'Caridad Proper', 'Hicgop South',
                 'Sitio Tugas', 'Elementary School', 'Sitio Lonas', 'San Agustin',
                 'Poroc San Flower', 'Poroc Rose', 'Visca', 'Gabas', 'Tamakin Dacot',
                 'Haina', 'Sitio Baybayon')

# data frame with site names and order alphabetically and geographically (N-S)
site_vec_order <- data.frame(site_name = site_vec)
site_vec_order$alpha_order <- seq(1:length(site_vec))
site_vec_order$geo_order <- c(5, 6, 7, 10, 16, 18, 8, 3, 4, 1, 14, 13, 12, 19, 11, 9, 17, 15, 2)
  
#tagging and fin-clipping thresholds
min_tag_size <- 6.0 #minimum size for tagging   
min_clip_size <- 3.5 #minimum size for fin-clip

#size thresholds for determining stage (just made up based on gut for now) - update based on data - Michelle has a boxplot somewhere?
min_breeding_F_size <- 6
min_breeding_M_size <- 6
breeding_F_YR_cutoff <- 9 #for now, saying if a fish is greater than 9cm but marked YR, probably a female
female_size_cutoff <- 10 #for now, saying if a fish is >10cm but we don't know anything about the color, probably a female

#first anemone tag in 2018
tag1_2018 <- 2938

#first metal anemone tag number (started using in May 2015)
first_metal_tag <- 2001 

#number of years sampled
years_sampled <- c(2012, 2013, 2014, 2015, 2016, 2017, 2018)

#years with tagging
tag_sample_years = c(2015,2016,2017,2018)

#winter months
winter_months <- c(1,2) #to pull out winter 2015 surveys - check that they didn't go into March too
spring_months <- c(3,4,5,6,7,8) #to pull out non-winter 2015 surveys

#recaptured fish known to be caught at two or more sites (been checked for typos)
multiple_site_recaps <- data.frame(tag_id = c('982000411818588', '982000411818610', '985153000401241'),
                                   gen_id = c(NA, NA, NA))

#anems at N and S ends of sites (visually from QGIS, in particular mid anems totally eyeballed)
Cabatoan_N <- 198 # 3195, 951 other options, 2502
Cabatoan_mid <- 2241
Cabatoan_S <- 904 #2250, 185 other options
CaridadCemetery_N <- 695
CaridadCemetery_mid <- 2683 #2256
CaridadCemetery_S <- 292
CaridadCemetery_mid <- 289
CaridadProper_N <- 291
CaridadProper_S <- 290
CaridadProper_mid <- NA
ElementarySchool_N <- 2303
ElementarySchool_mid <- 1306
ElementarySchool_S <- 702
Gabas_N <- 2667 #2271
Gabas_mid <- 1288
Gabas_S <- 1340
Haina_E <- 2138
Haina_mid <- 441 #2144, 2934
Haina_W <- 3080
HicgopSouth_N <- 300
HicgopSouth_mid <- 2128
HicgopSouth_S <- 697 #2296
Magbangon_N <- 2079
Magbangon_mid <- 680
Magbangon_S <- 213
Magbangon_N_N <- 2079
Magbangon_N_mid <- #2257 is mid of northern-most chunk; 475 or 2952 N of mid-chunk
  #Magbangon_N_S <- 680 #680 is S of hull, 1114 is another option; 212 is S end of northern-most chunk, 1391 is another option
  Magbangon_N_S <- 1114 #680 doesn't seem to have a lon value?
Magbangon_S_N <- 1113 #209 is another option
Magbangon_S_mid <- 2437
Magbangon_S_S <- 213 #214, 215 other options 
Palanas_N <- 2001 #1030 also quite close
Palanas_mid <- 2632 #totally eyeballed...
Palanas_S <- 876 #426 also a good end point
PorocRose_N <- 24 #2650
PorocRose_mid <- 2310
PorocRose_S <- 724
PorocSanFlower_N <- 902 #2315 another option
PorocSanFlower_mid <- 2646
PorocSanFlower_S <- 377 #2319
SanAgustin_N <- 2662
SanAgustin_mid <- 2660 #711 is mid of hull
SanAgustin_S <- 705 #outside of hull, 2129 is at bottom of hull 
SitioBaybayon_N <- 1302
SitioBaybayon_mid <- 538
SitioBaybayon_S <- 2747 #outside hull (maybe KC sites?), #805, 2148 w/in hull
SitioLonas_N <- 48
SitioLonas_S <- 48
SitioLonas_mid <- 48
SitioTugas_N <- 54
SitioTugas_S <- 59
SitioTugas_mid <- NA
TamakinDacot_N <- 2554
TamakinDacot_mid <- 2861
TamakinDacot_S <- 2270 #2147
Visca_N <- 4
Visca_mid <- 8 #776
Visca_S <- 2314
Wangag_N <- 2985
Wangag_mid <- 2734
Wangag_S <- 2063 #1034 also a good end point

# Put north, south, mid anems at each site into a dataframe
site_edge_anems <- data.frame(site = site_vec_NS, site_geo_order = c(5, 6, 7, 10, 16, 18, 8, 3, 4, 1, 14, 13, 12, 19, 11, 9, 17, 15, 2),
                              north_anem = c(Palanas_N, Wangag_N, Magbangon_N_N, Magbangon_S_N, Cabatoan_N, CaridadCemetery_N, CaridadProper_N,
                                             HicgopSouth_N, SitioTugas_N, ElementarySchool_N, SitioLonas_N, SanAgustin_N, PorocSanFlower_N,
                                             PorocRose_N, Visca_N, Gabas_N, TamakinDacot_N, Haina_W, SitioBaybayon_N),
                              south_anem = c(Palanas_S, Wangag_S, Magbangon_N_S, Magbangon_S_S, Cabatoan_S, CaridadCemetery_S, CaridadProper_S,
                                             HicgopSouth_S, SitioTugas_S, ElementarySchool_S, SitioLonas_S, SanAgustin_S, PorocSanFlower_S,
                                             PorocRose_S, Visca_S, Gabas_S, TamakinDacot_S, Haina_E, SitioBaybayon_S),
                              mid_anem = c(Palanas_mid, Wangag_mid, Magbangon_N_mid, Magbangon_S_mid, Cabatoan_mid, CaridadCemetery_mid, CaridadProper_mid,
                                HicgopSouth_mid, SitioTugas_mid, ElementarySchool_mid, SitioLonas_mid, SanAgustin_mid, PorocSanFlower_mid,
                                PorocRose_mid, Visca_mid, Gabas_mid, TamakinDacot_mid, Haina_mid, SitioBaybayon_mid))

# Prob of catching a fish by site, from KC script: https://github.com/katcatalano/parentage/blob/master/notebooks/proportion_sampled_allison.ipynb
prob_r <- c(0.5555556, 0.2647059, 0.8888889, 0.6666667, 0.2000000, #2016 recapture dives
            0.8333333, 0.4666667, 0.2000000, 0.8333333, 1.0000000, #2017 recapture dives
            0.3333333, 0.5789474, 0.6250000, 0.4090909) #2018 recapture dives


#################### Functions: ####################
# Functions from Michelle's GitHub helpers script
#field_helpers (similar idea to helpers, in field repository) - this might be the newer version of the helpers collection?
script <- getURL("https://raw.githubusercontent.com/pinskylab/field/master/scripts/field_helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

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

# Finds the real parameter estimate from the logit estimate
logit_recip <- function(logitval) {
  recip = (exp(logitval))/(1 + exp(logitval))
  return(recip)
}

# Function to attach 2018 anems to anems previously seen (since right now, anem_obs hasn't been extended to all of them), from AnemLocations.R script originally, this version from TotalAnemsAtSite.R script
attach2018anems <- function(anemdf) {
  
  #make anem_id numeric, create anem_id_unq
  anemdf <- anemdf %>%
    mutate(anem_id = as.numeric(anem_id)) %>% #make anem_id numeric to make future joins/merges easier
    mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), paste("obs", anem_obs, sep="")))  #add unique anem id so can track anems easier (has obs or id in front of whichever number it is reporting: eg. obs192)
  
  #filter out anems that might be repeat visits (2018, tag_id less than the new tags for 2018 or new tag and old tag present) and don't already have an anem_obs
  anems2018 <- anemdf %>% #646
    filter(year == 2018) %>%
    filter(anem_id < tag1_2018 | (anem_id >= tag1_2018 & !is.na(old_anem_id))) %>% #filter out the ones that could could have other anem_ids (sighting of existing anem_id tag or new tag with old_anem_id filled in), checked this (below) and think it is filtering out correctly...
    filter(is.na(anem_obs)) #some anem_obs already filled in so take those out... 
    
  #other 2018 anems that aren't candidates for repeat obs (so can add back in later)
  anems2018_2 <- anemdf %>% #270
    filter(year == 2018) %>%
    filter(anem_id >= tag1_2018 & is.na(old_anem_id))
  
  #other 2018 anems that already have an anem_obs (so can add back in later) - checked (once...) that anems2018, anems2018_2, and anems2018_3 covers all 2018 anems
  anems2018_3 <- anemdf %>%
    filter(year == 2018) %>%
    filter(!is.na(anem_obs))
  
  #filter out anems that are candidates for revisits (just anything from before 2018...), and that will get added back into final df later
  otheranems <- anemdf %>% #4068
    filter(year != 2018)
  
  #the filtering above covers all anems except anem_table_id 10473 with anem_id 2535, which has no year associated with it
  #test <- anemdf %>% filter(!year %in% c(2012,2013,2014,2015,2016,2017,2018)) #this is how I found that anem
  
  #go through anems2018 that might be revisits and see if there are anems from previous years that match
  for (i in 1:length(anems2018$anem_id)) {
    
    testid <- anems2018$anem_id[i] #pull out anem_id to compare
    testoldid <- anems2018$old_anem_id[i] #pull out old_anem_id
    
    matchanem <- filter(otheranems, anem_id == testid)  #does the anem_id match an anem_id record from the past?
    matcholdanem <- filter(otheranems, anem_id == testoldid) #does the old_anem_id match an anem_id from the past?
    
    # does the anem_id match an anem_id from the past? 
    if (length(matchanem$anem_id) > 0) { #if the anem_id matches an old one
      # if so, does the site match?
      if (matchanem$site[1] == anems2018$site[i]) { #make sure the site matches
        anems2018$anem_id_unq[i] = matchanem$anem_id_unq[1] #if there are multiple records from the past that match, take the first one
      } else {
        print(paste("Site does not match for anem_id", testid)) #print message if site doesn't match
      }
      # if not, does an old_anem_id match an anem_id from the past?
    } else if (length(matcholdanem$anem_id) > 0) { #if the old_anem_id matches one from the past
      # if so, does the site match?
      if (matcholdanem$site[1] == anems2018$site[i]) { #check site
        anems2018$anem_id_unq[i] = matcholdanem$anem_id_unq[1]
      } else {
        print(paste("Site does not match for old_anem_id", testoldid, "and anem_id", testid))
      }
    } else {
      anems2018$anem_id_unq[i] = anems2018$anem_id_unq[i]
      print(paste("No past anem matches found for testid", testid, "and testoldid", testoldid))
    }
  }
  out <- rbind(anems2018, anems2018_2, anems2018_3, otheranems)
  
  if(length(out$anem_table_id) == (length(anemdf$anem_table_id)-1)) { #1 anem (anem_table_id 10473, anem_id 2535 has no year listed - investigate further later) (maybe one of those 3 anems I sighted that used to get lost in the filtering b/c not have a year associated with them (NA)?)
    print("All anems accounted for.")
  } else {
    print("Some anems missing or double-counted.")
    print(length(out$anem_table_id))
    print(length(anemdf$anem_table_id))
  }
  return(out)
}

# Function to find the lat and long for an anem_id (based off of Michelle's sample_latlon function) - copied from AnemLocations.R
anemid_latlong <- function(anem.table.id, anem.df, latlondata) { #anem.table.id is one anem_table_id value, anem.df is the anem.Processed data frame (so don't have to pull from db again here), latlondata is table of GPX data from database (rather than making the function call it each time); will need to think a bit more clearly about how to handle different locations read for different visits to the same anem_id (or different with same anem_obs); for now, just letting every row in anem.Info get a lat-long
  
  #this is what causes the multiple entries - pulls multiple rows for a few anems (81) that have multiple entries for the same anem_table_id in the database
  anem <- anem.df %>% 
    filter(anem_table_id == anem.table.id) %>% #get the relevant dive, time, site, etc. info for this anem_table_id
    distinct(anem_table_id, .keep_all = TRUE) #added this in to get remove multiple entries that exist for some 2018 anem_table_ids
  
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

# Function with myconfig file info in it, for some reason new version of R/RStudio can't find the database...
read_db <- function(x) {
  db <- src_mysql(dbname = x,
                  port = 3306,
                  create = F,
                  host = "amphiprion.deenr.rutgers.edu",
                  user = "allisond",
                  password = "fish!NM?717")
  return(db)
}
#################### Pull out database info ####################
leyte <- read_db("Leyte") 

anem_db <- leyte %>% tbl("anemones") %>% collect()
fish_db <- leyte %>% tbl("clownfish") %>% collect()
fish_seen_db <- leyte %>% tbl('clown_sightings') %>% collect()
dives_db <- leyte %>% tbl("diveinfo") %>% collect()
gps_db <- leyte %>% tbl("GPX") %>% collect()

##### Pull all APCL caught or otherwise in the clownfish table
allfish_fish <- leyte %>% 
  tbl("clownfish") %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, gen_id, anem_table_id, recap, tag_id, color, size, fish_obs_time, fish_notes) %>%
  collect() %>%
  filter(fish_spp == 'APCL') %>%
  mutate(size = as.numeric(size))  # make the size numeric (rather than chr) so can do means and such

# and their corresponding anemones
allfish_anems <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_notes) %>%
  collect() %>%
  filter(anem_table_id %in% allfish_fish$anem_table_id)

# and the corresponding dive info
allfish_dives <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, dive_type, date, site, gps, dive_notes) %>%
  collect() %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id) 

# pull out just the year and put that in a separate column
allfish_dives$year <- as.integer(substring(allfish_dives$date,1,4))
 
# join together
allfish_caught <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
allfish_caught <- left_join(allfish_caught, allfish_dives, by="dive_table_id")

##### Create list of sites and the years they were sampled
# Summarize sites sampled each year
dives_processed <- dives_db %>%
  mutate(year = as.integer(substring(date,1,4))) %>%
  group_by(year, site) %>% 
  summarize(ndives = n()) %>%
  mutate(sampled = ifelse(ndives > 0, 1, 0)) # 1 if sampled (all 1 here b/c no entries for sites not sampled a particular year...)

site_visits <- data.frame(site = rep(site_vec_NS, length(years_sampled)))
site_visits$year <- c(rep(2012, length(site_vec_NS)), rep(2013, length(site_vec_NS)), rep(2014, length(site_vec_NS)), rep(2015, length(site_vec_NS)),
                      rep(2016, length(site_vec_NS)), rep(2017, length(site_vec_NS)), rep(2018, length(site_vec_NS)))
site_visits <- left_join(site_visits, dives_processed %>% select(year, site, sampled), by=c('year','site'))  # 1 if sampled in a particular year, NA if not

#################### Save files ####################
save(allfish_caught, file = here::here("Data", "allfish_caught.RData"))  # all caught APCL



# ####### Attempts at joining together APCL seen and caught on the same anems from the two tables, for Hannah,... see Katrina's email about how she did this for Adriana
# 
# # joining dives + anems to all fish from clownfish database for Hannah
# allfish_fishdb <- fish_db %>%
#   select(fish_table_id, anem_table_id, fish_spp, sample_id, gen_id, anem_table_id, recap, tag_id, color, size) 
#   
# allfish_caught_anems <- anem_db %>%
#   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_spp, anem_dia, anem_obs_time, notes) %>%
#   filter(anem_table_id %in% allfish_fishdb$anem_table_id)
# 
# allfish_caught_dives <- dives_db %>%
#   select(dive_table_id, dive_type, date, site, gps) %>%
#   filter(dive_table_id %in% allfish_caught_anems$dive_table_id)
# 
# allfish_caught <- left_join(allfish_fishdb, allfish_caught_anems, by='anem_table_id')
# allfish_caught <- left_join(allfish_caught, allfish_caught_dives, by='dive_table_id') %>%
#   mutate(year = substring(date,1,4))
# 
# # joining dives + anems to all fish from clown_sighted database for Hannah
# allfish_seendb <- fish_seen_db %>%
#   select(est_table_id, anem_table_id, fish_spp, size, notes, collector) 
# 
# allfish_seen_anems <- anem_db %>%
#   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_obs_time) %>%
#   filter(anem_table_id %in% allfish_seendb$anem_table_id)
# 
# allfish_seen_dives <- dives_db %>%
#   select(dive_table_id, dive_type, date, site, gps) %>%
#   filter(dive_table_id %in% allfish_seen_anems$dive_table_id)
# 
# allfish_seen <- left_join(allfish_seendb, allfish_seen_anems, by='anem_table_id')
# allfish_seen <- left_join(allfish_seen, allfish_seen_dives, by='dive_table_id') %>%
#   mutate(year = substring(date,1,4))
# 
# 
# # Try joining the two by anem_table_id - just presence/absence of fish spp
# allfish_seen_fish_spp <- allfish_seen %>%
#   dplyr::distinct(anem_table_id, fish_spp, .keep_all=TRUE)
# 
# allfish_caught_fish_spp <- allfish_caught %>%
#   distinct(anem_table_id, fish_spp, .keep_all=TRUE)
# 
# allfish_fish_spp_PA <- rbind(allfish_seen_fish_spp %>% select(anem_table_id, fish_spp, year, site),
#                             allfish_caught_fish_spp %>% select(anem_table_id, fish_spp, year, site))
#       
# 
# # Keeping anem_ids in the mix (b/c otherwise, how can track anems through time?) - doesn't filter out any differently...
# allfish_seen_fish_spp_anem_id <- allfish_seen %>%
#   dplyr::distinct(anem_table_id, fish_spp, anem_id, .keep_all=TRUE)
# 
# allfish_caught_fish_spp_anem_id <- allfish_caught %>%
#   distinct(anem_table_id, fish_spp, anem_id, .keep_all=TRUE)
# 
# # Are there any overlapping anem_table_ids between the two dfs? (shouldn't be, right?)
# test3 <- allfish_caught_fish_spp %>% filter(anem_table_id %in% allfish_seen_fish_spp)
# 
# # This data frame 
# allfish_fish_spp_PA_distinct <- allfish_fish_spp_PA %>%
#   distinct(anem_table_id, fish_spp, year, site)
# 
# allfish_fish_spp_PA %>%
#   filter(!anem_table_id %in% allfish_fish_spp_PA_distinct$anem_table_id)
# 
# !where_case_travelled_1 %in%
# 
# test2 <- as.data.frame(table(allfish_fish_spp_PA$anem_table_id))
# 
# ###################
# 
#   
# allAPCL_fish <- leyte %>% 
#   tbl("clownfish") %>%
#   select(fish_table_id, anem_table_id, fish_spp, sample_id, gen_id, anem_table_id, recap, tag_id, color, size) %>%
#   collect() %>%
#   filter(fish_spp == "APCL")
# 
# allAPCL_anems <- leyte %>%
#   tbl("anemones") %>%
#   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
#   collect() %>%
#   filter(anem_table_id %in% allAPCL_fish$anem_table_id)
# 
# allAPCL_dives <- leyte %>%
#   tbl("diveinfo") %>%
#   select(dive_table_id, dive_type, date, site, gps) %>%
#   collect() %>%
#   filter(dive_table_id %in% allAPCL_anems$dive_table_id)
# 
# # pull out just the year and put that in a separate column
# allAPCL_dives$year <- as.integer(substring(allAPCL_dives$date,1,4))
# allAPCL_dives$month <- as.integer(substring(allAPCL_dives$date,6,7))
# 
# # join together
# allAPCL <- left_join(allAPCL_fish, allAPCL_anems, by="anem_table_id")
# allAPCL <- left_join(allAPCL, allAPCL_dives, by="dive_table_id")
# 
# allAPCL$size <- as.numeric(allAPCL$size) #make size numeric (rather than a chr) so can do means and such
# 
# #pull out just tagged fish
# taggedfish <- allAPCL %>% filter(!is.na(tag_id))
# 
# # all anems that had APCL at some point (have a tag)
# anems_tagged <- anem_db %>%
#   filter(!is.na(anem_id) | !is.na(anem_obs)) %>%
#   mutate(anem_id_unq = if_else(is.na(anem_obs), paste('id',anem_id,sep=''), paste('obs',anem_obs,sep='')))
# 

# #################### Save files ####################
# save(leyte, file = here::here("Data/Database_backups", "leyte.RData"))
# save(anem_db, file = here::here("Data/Database_backups", "anem_db.RData"))
# save(fish_db, file = here::here("Data/Database_backups", "fish_db.RData"))
# save(fish_seen_db, file = here::here('Data/Database_backups', 'fish_seen_db.RData'))
# save(dives_db, file = here::here("Data/Database_backups", "dives_db.RData"))
# save(gps_db, file = here::here("Data/Database_backups", "gps_db.RData"))
# 
# save(allAPCL, file=here::here('Data', 'allAPCL.RData'))
# save(taggedfish, file=here::here('Data', 'taggedfish.RData'))
# 
# 
# # Data for Hannah
# save(anems_tagged, file=here::here('Data','anems_tagged.RData'))
# save(allfish_caught, file=here::here('Data','allfish_caught.RData'))
# save(allfish_seen, file=here::here('Data','allfish_seen.RData'))
# 
# load(file = here("Data", "anem_db.RData"))
# load(file = here("Data", "fish_db.RData"))
# load(file = here("Data", "dives_db.RData"))
# load(file = here("Data", "gps_db.RData"))
