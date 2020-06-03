# Functions and constants used throughout clownfish persistence scripts 

#################### Install packages: ####################
library(RCurl) #allows running R scripts from GitHub - now reading dispersal kernel tables from KC parentage repository
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(here)

#################### Load data from database and other analyses: ####################
### Option 1: pull from database + GitHub
# source(here::here("Code", "Pull_data_from_database.R"))

### Option 2: load previously-pulled saved files
# Database files
load(file = here::here("Data/Data_from_database", "anem_db.RData"))
load(file = here::here("Data/Data_from_database", "fish_db.RData"))
load(file = here::here("Data/Data_from_database", "fish_seen_db.RData"))
load(file = here::here("Data/Data_from_database", "dives_db.RData"))
load(file = here::here("Data/Data_from_database", "gps_db.RData"))

# Files from other analyses (parentage, matching up fish, etc.)
load(file = here::here("Data/From_other_analyses", "fish-obs.RData"))  # fish obs table linking together individual fish marked and recaught by both genetic and tag methods (https://github.com/pinskylab/genomics/blob/master/data/fish-obs.csv)
load(file = here::here("Data/From_other_analyses", "site_centers.RData"))  # site centers (eyeballed in QGIS by KC)
load(file = here::here("Data/From_other_analyses", "kernel_summary.RData"))  # kernel fits (https://github.com/katcatalano/parentage/blob/master/kernel_fitting/1340_loci/final_results/tables/kernel_fitting_summary.csv)
load(file = here::here("Data/From_other_analyses", "all_parents.RData"))  # all parents in parentage analysis (https://github.com/katcatalano/parentage/blob/master/colony2/20190523_1340loci/input/all_parents_corrected.txt)
load(file = here::here("Data/From_other_analyses", "all_offspring.RData"))  # all offspring in parentage analysis (https://github.com/katcatalano/parentage/blob/master/colony2/20190523_1340loci/input/all_offspring_corrected.txt)
load(file = here::here("Data/From_other_analyses", "k_theta_allyear_95CI_values.RData"))  # 95% k and theta fits (https://github.com/katcatalano/parentage/blob/master/kernel_fitting/1340_loci/final_results/likelihood_profiles_grid_search/profile95CI_AllYears.csv)
load(file = here::here("Data/From_other_analyses", "length_count8llEA.RData"))  # size-fecundity model (from Adam Yawdoszyn work)

load(file = here::here("Data/From_other_analyses", "old_parent_file.RData"))  # loading old parent file to see if can just keep using it

#################### Set constants and parameters: ####################

##### Sampling information
years_parentage <- c(2012, 2013, 2014, 2015, 2016, 2017, 2018)  # years included in parentage analyses
years_sampled <- c(2012, 2013, 2014, 2015, 2016, 2017, 2018)  # years sampled
tag_sample_years = c(2015,2016,2017,2018)  # years with PIT tagging

winter_months <- c(1,2)  # to pull out winter 2015 surveys
spring_months <- c(3,4,5,6,7,8)  # to pull out non-winter 2015 surveys

# dive types where clownfish were sampled
clown_sample_dive_types <- c("0","A","C","D","E","F","M")  # everything but R 

##### Fish information (all sizes in cm)
min_tag_size <- 6  # minimum size for tagging   
min_clip_size <- 3.5  # minimum size for fin-clip
min_size = 0  # minimum fish size 
max_size = 15  # largest fish we've seen

# #size thresholds for determining stage - CHECK IF THESE ARE ACTUALLY STILL USED!!
# min_breeding_F_size <- 6  # size thresholds for determining stage (made up based on gut for now, should update based on data - Michelle boxplot?)
# min_breeding_M_size <- 6  # size thresholds for determining stage (made up based on gut for now, should update based on data - Michelle boxplot?)
# breeding_F_YR_cutoff <- 9  # for now, saying if a fish is greater than 9cm but marked YR, probably a female
# female_size_cutoff <- 10  # for now, saying if a fish is >10cm but we don't know anything about the color, probably a female

##### Anemone information
first_metal_tag <- 2001  # first metal anemone tag number (started using metal tags in May 2015)
#tag1_2018 <- 2938  # first anemone tag in 2018

##### Other values
clutches_per_year_mean = 11.9  # from Holtswarth et al. 2017 paper
n_runs = 1000  # number of runs for uncertainty analyses
max_larval_nav_buffer = 1000 # in m, maximum amount of buffer added to site edges to account for larval navigation
northern_site_pos = 1  # when sites are ordered N-S, number of northern-most site (Palanas)
southern_site_pos = 19  # when sites are ordered N-S, number of southern-most site (Sitio Baybayon)

##### Site information
# Vector of sites, in alphabetical order, for referencing sites easily in plots and filtering and such (does not include Sitio Hicgop, which has no fish)
site_vec <- c("Cabatoan", "Caridad Cemetery", "Caridad Proper", "Elementary School", "Gabas", "Haina", "Hicgop South",
              "N. Magbangon", "Palanas", "Poroc Rose", "Poroc San Flower", "S. Magbangon", "San Agustin", "Sitio Baybayon", "Sitio Lonas",
              "Sitio Tugas", "Tamakin Dacot", "Visca", "Wangag")

# Set sites to include for total possible sampling area (all site areas times all years sampled) - excluding Caridad Proper, Sitio Lonas, Sitio Tugas from this because they disappeared partway through (check with Michelle on that) (should I exclude the Sitio Lonas match too?)
sites_for_total_areas <- c('Cabatoan', 'Caridad Cemetery', 'Elementary School', 'Gabas', 
                           'Haina', 'Hicgop South', 'N. Magbangon', 'Palanas', 'Poroc Rose',
                           'Poroc San Flower', 'San Agustin', 'Sitio Baybayon', 'Tamakin Dacot',
                           'Visca', 'Wangag', 'S. Magbangon')  

# Vector of sites from north to south (does not include Sitio Hicgop, which has no fish)
site_vec_NS <- c('Palanas', 'Wangag', 'N. Magbangon', 'S. Magbangon' , 'Cabatoan',
                 'Caridad Cemetery', 'Caridad Proper', 'Hicgop South',
                 'Sitio Tugas', 'Elementary School', 'Sitio Lonas', 'San Agustin',
                 'Poroc San Flower', 'Poroc Rose', 'Visca', 'Gabas', 'Tamakin Dacot',
                 'Haina', 'Sitio Baybayon')

# Vector of sites in alphabetical order when no spaces (for MARK, otherwise issues running some of the models with site)
no_space_sites_alpha <- c("Cabatoan", "Caridad Cemetery", "Caridad Proper", "Elementary School", "Gabas", "Haina", "Hicgop South",
                          "N. Magbangon", "Palanas", "Poroc Rose", "Poroc San Flower", "San Agustin", "Sitio Baybayon", 
                          "Sitio Lonas", "Sitio Tugas", "S. Magbangon", "Tamakin Dacot", "Visca", "Wangag")

# No-space sites, without ones we didn't revisit
no_space_sites_revisited = c("Cabatoan", "Caridad Cemetery", "Elementary School", "Gabas", "Haina",
                             "Hicgop South", "N. Magbangon", "Palanas", "Poroc Rose", "Poroc San Flower", 
                             "San Agustin", "Sitio Baybayon", "S. Magbangon", "Tamakin Dacot", "Visca", "Wangag")

# Sites not revisited (only sampled in 1 year)
sites_not_revisited = c("Caridad Proper", "Sitio Lonas", "Sitio Tugas")

##### Anemones at edges of sites (identified visually from QGIS, in particular mid anems)
Cabatoan_N <- 198  # 3195, 951 other options, 2502
Cabatoan_mid <- 2241
Cabatoan_S <- 904  # 2250, 185 other options
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
Magbangon_N_mid <- 2257  # is mid of northern-most chunk; 475 or 2952 N of mid-chunk 
Magbangon_N_S <- 1114  # 680 doesn't seem to have a lon value? #680 is S of hull, 1114 is another option; 212 is S end of northern-most chunk, 1391 is another option
Magbangon_S_N <- 1113  # 209 is another option
Magbangon_S_mid <- 2437
Magbangon_S_S <- 213  # 214, 215 other options 
Palanas_N <- 2001  # 1030 also quite close
Palanas_mid <- 2632  # totally eyeballed...
Palanas_S <- 876  # 426 also a good end point
PorocRose_N <- 24  # 2650
PorocRose_mid <- 2310
PorocRose_S <- 724
PorocSanFlower_N <- 902  # 2315 another option
PorocSanFlower_mid <- 2646
PorocSanFlower_S <- 377  # 2319
SanAgustin_N <- 2662
SanAgustin_mid <- 2660  # 711 is mid of hull
SanAgustin_S <- 705  # outside of hull, 2129 is at bottom of hull 
SitioBaybayon_N <- 1302
SitioBaybayon_mid <- 538
SitioBaybayon_S <- 2747  # outside hull (maybe KC sites?), #805, 2148 w/in hull
SitioLonas_N <- 48
SitioLonas_S <- 48
SitioLonas_mid <- 48
SitioTugas_N <- 54
SitioTugas_S <- 59
SitioTugas_mid <- NA
TamakinDacot_N <- 2554
TamakinDacot_mid <- 2861
TamakinDacot_S <- 2270  # 2147
Visca_N <- 4
Visca_mid <- 8  # 776
Visca_S <- 2314
Wangag_N <- 2985
Wangag_mid <- 2734
Wangag_S <- 2063  # 1034 also a good end point

# Put north, south, mid anems at each site into a dataframe
site_edge_anems <- data.frame(site = site_vec_NS, site_geo_order = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
                              anem_id = c(Palanas_N, Wangag_N, Magbangon_N_N, Magbangon_S_N, Cabatoan_N, CaridadCemetery_N, CaridadProper_N,
                                          HicgopSouth_N, SitioTugas_N, ElementarySchool_N, SitioLonas_N, SanAgustin_N, PorocSanFlower_N,
                                          PorocRose_N, Visca_N, Gabas_N, TamakinDacot_N, Haina_W, SitioBaybayon_N, Palanas_S, Wangag_S, Magbangon_N_S, Magbangon_S_S, Cabatoan_S, CaridadCemetery_S, CaridadProper_S,
                                          HicgopSouth_S, SitioTugas_S, ElementarySchool_S, SitioLonas_S, SanAgustin_S, PorocSanFlower_S,
                                          PorocRose_S, Visca_S, Gabas_S, TamakinDacot_S, Haina_E, SitioBaybayon_S, Palanas_mid, Wangag_mid, Magbangon_N_mid, Magbangon_S_mid, Cabatoan_mid, CaridadCemetery_mid, CaridadProper_mid,
                                          HicgopSouth_mid, SitioTugas_mid, ElementarySchool_mid, SitioLonas_mid, SanAgustin_mid, PorocSanFlower_mid,
                                          PorocRose_mid, Visca_mid, Gabas_mid, TamakinDacot_mid, Haina_mid, SitioBaybayon_mid),
                              anem_loc = c(rep('north',length(site_vec_NS)), rep('south',length(site_vec_NS)), rep('mid', length(site_vec_NS))))

#################### Functions: ####################

# Function to make vector of strings for column names for something done each year (like columns for sampling each year or minimum distance sampled to each anem each year, etc.)
makeYearlyColNames <- function(start.Year, end.Year, descriptor) {  # start.Year is first year of sampling, end.Year is final year of sampling, descriptor is string to go before year in column name (like "min_dist_" if want columns like "min_dist_2012")
  
  if (start.Year > end.Year) {
    out <- NULL
  } 
  else {
    out <- as.vector(NA)  # initalize output vector of column names
    year <- start.Year 
    
    for (i in 1:(end.Year - start.Year + 1)) {
      out[i] <- paste(descriptor, year, sep="")  # create string like "min_dist_2012", where descriptor is something like "min_dist_" and year is the year sampled
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


# Function to find the lat and long for an anem_id (based off of Michelle's sample_latlon function) 
anemid_latlong <- function(anem.table.id, anem.df, latlondata) {  # anem.no is an anem_table_id or anem_obs or anem_id; anem.no.type is which of those it is, anem.df is the anem.Processed data frame, latlondata is table of GPX data from database 

  # this is what causes the multiple entries - pulls multiple rows for a few anems (81) that have multiple entries for the same anem_table_id in the database
  anem <- anem.df %>%
    filter(anem_table_id == anem.table.id) %>%  # get the relevant dive, time, site, etc. info for this anem_table_id
    distinct(anem_table_id, .keep_all = TRUE)  # added this in to get remove multiple entries that exist for some 2018 anem_table_ids

  latloninfo <- latlondata %>%
    filter(as.character(gps_date) %in% anem$date & unit == anem$gps) %>%  # filter out just the GPS unit associated with this anem observation and on the right date
    filter(gps_hour == anem$anem_hour & gps_min == anem$anem_min) %>%
    mutate(lat = as.numeric(lat)) %>%
    mutate(lon = as.numeric(lon))
  
  anem$lat <- round(mean(latloninfo$lat, na.rm = TRUE), digits = 5)  # take the mean of the latitude values (digits = 5 b/c that is what Michelle had)
  anem$lon <- round(mean(latloninfo$lon, na.rm = TRUE), digits = 5)  # take the mean of the longitude values
  
  return(anem)
}


#################### Edit + clean data from other analyses (if necessary): ####################
 
##### Prob of catching a fish by site from the recapture dives, from KC script: https://github.com/katcatalano/parentage/blob/master/notebooks/proportion_sampled_allison.ipynb
prob_r <- c(0.5555556, 0.2647059, 0.8888889, 0.6666667, 0.2000000,  # 2016 recapture dives
            0.8333333, 0.4666667, 0.2000000, 0.8333333, 1.0000000,  # 2017 recapture dives
            0.3333333, 0.5789474, 0.6250000, 0.4090909) #2018 recapture dives

# Raw egg counts from photos (from egg_data2018f.csv in Adam's repository)
egg_counts_AY_data <- c(479, 590, 586, 305, 679, 683, 387, 720, 427, 688, 169, 655, 414, 352, 1102, 265, 1886, 904,
                        851, 160, 648, 766, 1060, 670, 351, 557)  

##### Site centers - edit so site names match the ones I use (add spaces in N. Magbangon and S. Magbangon)
site_centers$site[which(site_centers$site == "N.Magbangon")] = "N. Magbangon"
site_centers$site[which(site_centers$site == "S.Magbangon")] = "S. Magbangon"

#################### Pull out relevant values: ####################

# Pull out number of parents, number of offspring genotyped, number of offspring matched
n_offspring_genotyped <- (kernel_summary %>% filter(Year == "2012-2018"))$NumOffsSampled  # number of potential offspring genotyped (whether or not they were matched to parents)
n_offspring_matched <- (kernel_summary %>% filter(Year == "2012-2018"))$NumParentageMatches  # number of offspring matched to parents
n_parents_genotyped <- length((all_parents %>% distinct(fish_indiv))$fish_indiv)  # number of potential parents genotyped (whether or not they were matched to offspring)

# Pull out kernel fits (from 2D fitting where both k and theta fit simultaneously)
theta_allyears <- (kernel_summary %>% filter(Year == "2012-2018"))$best_theta
k_allyears <- (kernel_summary %>% filter(Year == "2012-2018"))$best_k

# Mean eggs per clutch from photo counts
mean_eggs_per_clutch_from_counts <- mean(egg_counts_AY_data)

#################### Process data ####################

##### Data frame with site names and order both alphabetically and geographically (N-S)
site_vec_order <- data.frame(site_name = site_vec, stringsAsFactors = FALSE)
site_vec_order$alpha_order <- seq(1:length(site_vec))
site_vec_order$geo_order <- c(5, 6, 7, 10, 16, 18, 8, 3, 1, 14, 13, 4, 12, 19, 11, 9, 17, 15, 2)

##### Match up other relevant info (site, date, fish_indiv, etc.) to fish in the clownfish table
# Pull out year and month into a separate column in dives_db
dives_db_processed <- dives_db %>%
  mutate(year = as.integer(substring(date,1,4))) %>%
  mutate(month = as.integer(substring(date,6,7))) %>%
  mutate(dive_date = date(date))

# Pull all APCL caught or otherwise in the clownfish table
allfish_fish <- fish_db %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, anem_table_id, recap, tag_id, color, sex, size, fish_obs_time, fish_notes) %>%
  filter(fish_spp == 'APCL') %>%
  mutate(size = as.numeric(size))  # make the size numeric (rather than chr) so can do means and such

# and their corresponding anemones
allfish_anems <- anem_db %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_notes) %>%
  filter(anem_table_id %in% allfish_fish$anem_table_id)

# and the corresponding dive info
allfish_dives <- dives_db_processed %>%
  select(dive_table_id, dive_type, date, year, month, site, gps, dive_notes) %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id) 

# join together
allfish_caught <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
allfish_caught <- left_join(allfish_caught, allfish_dives, by="dive_table_id")

# add in the gen_ids and fish_indiv (now in a separate gen_id table) - gen_id only comes in the time the fish was sequenced, not at all captures
allfish_caught <- left_join(allfish_caught, (fish_obs %>% select(fish_table_id, gen_id, fish_indiv)), by = "fish_table_id")

# remove intermediate dataframes, just to keep things tidy
rm(allfish_fish, allfish_anems, allfish_dives) 

##### Create list of sites and the years they were sampled 
# Summarize sites sampled each year
dives_processed <- dives_db_processed %>%
  group_by(year, site) %>% 
  summarize(ndives = n()) %>%
  mutate(sampled = ifelse(ndives > 0, 1, 0))  # 1 if sampled (all 1 here b/c no entries for sites not sampled a particular year...)

# Summarize sites sampled for clownfish each year (so only include spring dives, not the 2015 winter anemone survey dives)
dives_processed_clownfish <- dives_db_processed %>%
  filter(month %in% spring_months) %>%  # only include dives in the spring months, 
  group_by(year, site) %>% 
  summarize(ndives = n()) %>%
  mutate(sampled = ifelse(ndives > 0, 1, 0))  # 1 if sampled (all 1 here b/c no entries for sites not sampled a particular year...)

# Put sampling indicator into dataframe with sites
site_visits <- data.frame(site = rep(site_vec_NS, length(years_sampled)))
site_visits$year <- c(rep(2012, length(site_vec_NS)), rep(2013, length(site_vec_NS)), rep(2014, length(site_vec_NS)), rep(2015, length(site_vec_NS)),
                      rep(2016, length(site_vec_NS)), rep(2017, length(site_vec_NS)), rep(2018, length(site_vec_NS)))
site_visits <- left_join(site_visits, dives_processed_clownfish %>% select(year, site, sampled), by=c('year','site'))  # 1 if sampled in a particular year, NA if not

##### Process anems to merge with dive info, get times in same time zone as gps
# Merge with dive info
anems_Processed_dive <- left_join(anem_db, dives_db_processed, by="dive_table_id")

# Take out anems with no anem_id or species, convert time zones so matches with GPS 
anems_Processed <- anems_Processed_dive %>%
  filter(!is.na(anem_id) | anem_id != "-9999" | anem_id == "") %>%  # filter out NAs, -9999s (if any left in there still), and blanks
  filter(is.null(anem_spp) == FALSE) %>%  # to remove "phantom" anem observations that were actually just clownfish processing, get rid of obs with anem_spp that is null...
  filter(is.na(anem_spp) == FALSE) %>%  # or anem_spp that is NA
  mutate(anem_id = as.numeric(anem_id)) %>%  # make anem_id numeric to make future joins/merges easier
  mutate(anem_id_unq = anem_obs) %>%  # now can just use anem_obs, keeping this in there so don't have to change rest of code
  mutate(obs_time = force_tz(ymd_hms(str_c(date, anem_obs_time, sep = " ")), tzone = "Asia/Manila")) %>%  # tell it that it is currently in Asia/Manila time zone
  mutate(obs_time = with_tz(obs_time, tzone = "UTC")) %>%  # convert to UTC so can compare with GPS data (this line and one above largely from Michelle's assign_db_gpx function)
  mutate(anem_month = month(obs_time),  # and separate out useful components of the time (this also from Michelle's assign_db_gpx function)
         anem_day = day(obs_time), 
         anem_hour = hour(obs_time), 
         anem_min = minute(obs_time), 
         anem_sec = second(obs_time))
# When run the above, get a warning that two failed to parse - assume those are the two that don't have anem_obs_times

# Do a version for all anems (even those without anem_ids)
anems_Processed_all <- anems_Processed_dive %>%
  mutate(obs_time = force_tz(ymd_hms(str_c(date, anem_obs_time, sep = " ")), tzone = "Asia/Manila")) %>%  # tell it that it is currently in Asia/Manila time zone
  mutate(obs_time = with_tz(obs_time, tzone = "UTC")) %>%  # convert to UTC so can compare with GPS data (this line and one above largely from Michelle's assign_db_gpx function)
  mutate(anem_month = month(obs_time),  # and separate out useful components of the time (this also from Michelle's assign_db_gpx function)
         anem_day = day(obs_time), 
         anem_hour = hour(obs_time), 
         anem_min = minute(obs_time), 
         anem_sec = second(obs_time))
# When run the above, get a warning that three failed to parse - assume those are the two that don't have anem_obs_times plus some other one

##### Process gps data so date and time is on there and more accessible for comparison
gps_Info <- gps_db %>%
  mutate(gps_date = date(time),
         gps_day = day(time),
         gps_hour = hour(time),
         gps_min = minute(time),
         gps_sec = second(time),
         gps_year = year(time),
         gps_month = month(time))

# Attach lat/lons to anems - takes a long time to run...
# anems_Processed_latlon <- anems_Processed %>%
#   mutate(lat = NA,
#          lon = NA)
# 
# for(i in 1:length(anems_Processed_latlon$anem_table_id)) {
#   output = anemid_latlong(anems_Processed_latlon$anem_table_id[i], anems_Processed, gps_Info)
#   anems_Processed_latlon$lat[i] = output$lat
#   anems_Processed_latlon$lon[i] = output$lon
# }

##### Add sites to fish in the parents file (new file)
# Make fish_indiv a character so can mesh with allfish_caught
all_parents_edited <- all_parents %>%
  mutate(fish_indiv = as.character(fish_indiv))

# Then join the two so parents have site assigned to them
all_parents_by_site <- left_join(all_parents_edited %>% dplyr::rename(fish_indiv_parent = fish_indiv, gen_id_parent = gen_id), 
                                 allfish_caught %>% select(fish_indiv, site, gen_id, sample_id), by = "sample_id") %>%
  distinct(fish_indiv, .keep_all = TRUE)  # this doesn't actually remove any for now, think about whether it should be included

# Now do it with old version of parents file too so can compare
# Make fish_indiv a character so can mesh with allfish_caught
all_parents_edited_old <- old_parent_file %>%
  mutate(fish_indiv = as.character(fish_indiv))

# Then join the two so parents have site assigned to them
all_parents_by_site_old <- left_join(all_parents_edited_old %>% dplyr::rename(fish_indiv_parent = fish_indiv, gen_id_parent = gen_id), 
                                 allfish_caught %>% select(fish_indiv, site, gen_id, sample_id), by = "sample_id") %>%
  distinct(fish_indiv, .keep_all = TRUE)  # this doesn't actually remove any for now, think about whether it should be included

# Compare number of parents by site
old_parents_by_site <- data.frame(table(all_parents_by_site_old$site))
new_parents_by_site <- data.frame(table(all_parents_by_site$site))

#################### Save files processed in this script ####################
save(site_vec_order, file = here::here("Data/Script_outputs", "site_vec_order.RData"))
save(allfish_caught, file = here::here("Data/Script_outputs", "allfish_caught.RData"))  # all caught APCL 
save(dives_processed, file = here::here("Data/Script_outputs", "dives_processed.RData"))  # number of dives at each site in each year (for determining if a site was sampled)
save(dives_processed_clownfish, file = here::here("Data/Script_outputs", "dives_processed_clownfish.RData"))  # number of dives during clownfish sampling time in each year (not necessarily clownfish dives specifically)
save(site_visits, file = here::here("Data/Script_outputs", "site_visits.RData"))  # whether or not each site was sampled in each year
save(anems_Processed, file = here::here("Data/Script_outputs", "anems_Processed.RData"))  # anems with ids, times processed and such
save(anems_Processed_all, file = here::here("Data/Script_outputs", "anems_Processed_all.RData"))  # times processed for all anems, even without ids
save(all_parents_by_site, file = here::here("Data/Script_outputs", "all_parents_by_site.RData"))  # number of parents by site
save(site_edge_anems, file = here::here("Data/Script_outputs", "site_edge_anems.RData"))  # anem ids at edges and middle of each site (eyeballed)

