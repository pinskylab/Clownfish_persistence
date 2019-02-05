# Script to calculate persistence metrics
# For now, using fake data + parameters where necessary to get structure right - will note clearly, format might shift as real stuff comes in
# 
# #################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))

# #Load relevant libraries
# library(RCurl) #allows running R scripts from GitHub
# library(RMySQL) #might need to load this to connect to the database?
# library(dplyr)
# library(tidyr)
# library(lubridate)
# library(stringr)
# library(geosphere)
# library(varhandle)
# library(ggplot2)
# library(reshape)
# library(here)

# Load some data files while deciding how various scripts should interact/communicate
# Load output from AnemLocations.R (or source that file) - mean gps coordinates across multiple observations of anemones
#load(file=here("Data",'AnemAllInfowLatLon2.RData')) #file with anems, after 2018 anems matched
#load(file=here("Data", "AnemLatLonObsbyAnem.RData")) #file with lat-lon info appended

# Load encounter history data frame from ClownfishMarkRecap.R
#load(file=here("Data", "encounters_all.RData"))

# # Load mig_est data
# # #load(file=here("Data", "MigEstPooled_out_justimmigration.txt"))
# # # Old version - MigEst pooled matrix but w/out ghost populations
# # migEst_conn <- read.table(file=here("Data", "MigEstPooled_out_justimmigration.txt"), 
# #            header = TRUE, sep = "")
# migEst_sites <- c("Cabatoan", "Caridad Cemetery", "Caridad Proper", "Elementary School", "Gabas",
#                   "Haina", "Hicgop South", "Magbangon", "Palanas", "Poroc Rose", "Poroc San Flower",
#                   "San Agustin", "Sitio Baybayon", "Sitio Lonas", "Sitio Tugas", "Tamakin Dacot",
#                   "Visca", "Wangag")
# migEst_conn <- read.table(file=here::here('Data','allyears_average_migest_justmeantable.txt'), skip=2, header=TRUE) #just the mean all-years table from allyears_average_migest.txt, run on 9/25/18
# # While using this migEst output that has Magbangon as one site (not broken into N and S), split it by assuming the arrivals at each are the same and half of the dispersers from Mag are from N and half from S
# # Duplicate the Magbangon row in the migEst output to have one for N. Mag and one for S. Mag (eventually, will get a run with them as separate sites), and remove the column with pop number
# migEst_conn <- migEst_conn[,-1] #remove Pop number column
# Mag_pos = 8 #position of the Magbangon row and column (once pop column is removed)
# n_migEstrow <- dim(migEst_conn)[1] #number of rows total
# n_migEstcol <- dim(migEst_conn)[2]
# # duplicate the row showing composition of arrivers to Mag so same for N and S
# migEst_1 <- migEst_conn[1:Mag_pos,] #get Mag row once here (as well as those before it)
# migEst_2 <- migEst_conn[Mag_pos:n_migEstrow,] #get Mag row once here (as well as those after it)
# migEst_conn <- rbind(migEst_1, migEst_2)
# # duplicate the column showing outgoers from Mag and halve the values (so same total percentage from Mag, just half are from N and half are from S)
# migEst_a <- migEst_conn[,1:Mag_pos] #get Mag column once here (plus all the ones to the left of it)
# migEst_a[,Mag_pos] = migEst_a[,Mag_pos]/2 # divide values in half
# migEst_b <- migEst_conn[,Mag_pos:n_migEstcol] #get Mag column again here, plus all the columns to the right of it
# migEst_b[,1] <- migEst_b[,1]/2
# migEst_conn <- cbind(migEst_a, migEst_b)
# 
# rm(migEst_1, migEst_2, migEst_a, migEst_b) #remove in-between data frames to clear up clutter...

# Load in site_areas
load(file=here::here('Data','site_areas.RData'))

# Load file with proportion habitat sampled estimates (from TotalAnemsAtSite.R)
load(file=here("Data",'anem_sampling_table.RData')) #file with anems, after 2018 anems matched

# Load in new MigEst output (with N. Magbangon and S. Magbangon separated) - all years combined and confidence intervals
migEst_conn <- read.table(file=here::here('Data', '20181029_connmat_allyears.txt'), header=TRUE) #read in all years combined values
migEst_conn <- migEst_conn[,-1] #remove Pop number column
migEst_conn_CI <- read.table(file=here::here('Data', '20181029_confid_allyears.txt'), header=TRUE) #read in confidence intervals
migEst_conn_UCI <- migEst_conn_CI[, grepl("CI_UJ_", names(migEst_conn_CI))] #pull out just upper confidence estimates
migEst_conn_LCI <- migEst_conn_CI[, grepl('CI_LJ_', names(migEst_conn_CI))] #pull out just lower confidence estimates

# Load all parentage matches (as of Nov 2018, N and S Mag separated but before 2016, 2017, 2018 genotypes are in)
parentage_moms <- read.csv(file=here('Data','20181017colony_migest_mums_allyears.csv'), stringsAsFactors = FALSE)
parentage_dads <- read.csv(file=here('Data','20181017colony_migest_dads_allyears.csv'), stringsAsFactors = FALSE)
parentage_trios <- read.csv(file=here('Data','20181017colony_migest_trios_allyears.csv'), stringsAsFactors = FALSE)

# # old files before N and S Mag were separated and before 2016, 2017, 2018 genotypes are in)
# parentage_moms <- read.csv(file=here('Data','20180713colony_migest_mums_allyears.csv'), stringsAsFactors = FALSE)
# parentage_dads <- read.csv(file=here('Data','20180713colony_migest_dads_allyears.csv'), stringsAsFactors = FALSE)
# parentage_trios <- read.csv(file=here('Data','20180713colony_migest_trios_allyears.csv'), stringsAsFactors = FALSE)

# Load recruit estimates
load(file=here("Data", "recruits_info.RData"))
load(file=here('Data', 'breedingF_info.RData'))
#load(file=here("Data", "females_recruits_summary.RData"))
#load(file=here("Data", "metapop_level.RData"))

# Load egg-recruit linear model estimates
#load(file=here("Data", "metapop_lm.RData")) #estimate at metapop level, pretty low R2 and not significant
#load(file=here("Data", "egg_recruit_lm.RData")) #estimate using each site and year as a point, R2 of about 0.53, significant
load(file=here('Data', 'egg_recruit_lm_mod1.RData'))  # estimate using each site and year as a point, scaled by m2
load(file=here('Data', 'egg_recruit_lm_mod6.RData'))  # same as egg_recruit_lm_mod1 buth with intercept forced to be 0
load(file=here('Data', 'metapop_lm_mod5.RData'))

# Load LEP estimates
load(file=here("Data","LEP_out.RData")) #LEP estimates from LEP_estimate.R

# Load proportion habitat sampled info
load(file=here('Data','anem_sampling_table.RData'))

# Dispersal kernel parameters from Katrina (estimates as of 12/15/18 from KC paper draft) 
k_2012 = -2.67
theta_2012 = 3
k_2013 = -3.27
theta_2013 = 3
k_2014 = -2.38
theta_2014 = 2
k_2015 = -2.73
theta_2015 = 2
k_allyears = -1.36  # with 2012-2015 data
theta_allyears = 0.5  # with 2012-2015 data

# #estimates as of 11/01/2018, in 20181101_allkernels.pdf 
# k_allyears = 0.07 #2012-2015 data
# theta_allyears = 0.5 #2012-2015 data
# k_2012 =-2.35
# theta_2012 = 0.5
# k_2013 = -1.83
# theta_2013 = 0.5
# k_2014 = -2.23
# theta_2014 = 2
# k_2015 = -2.7
# theta_2015 = 2

# PRELIM OR FAKE FOR NOW THAT WILL GET REPLACED BY REAL DATA
eggs_per_clutch = 1763 #from LEP_calc_WillWhite.R
clutch_per_year = 11.9 #from LEP_calc_WillWhite.R
eggs_intercept = 100
eggs_slope = 107 #from Adam Y Aresty poster, #eggs/cm in female size

# LEP
LEP_Will <- 1780 #eggs/recruit
LEP_Will_eggs_per_clutch <- LEP_out$LEP_W #my size-survival relationship, Will's eggs per clutch
LEP <- LEP_out$LEP #AY's prelim egg data, my size-survival relationship

# Process saved gps database info
gps.Info <- gps_db %>%  
  select(lat, lon, time, unit) %>%
  mutate(obs_time = force_tz(ymd_hms(time), tzone = "UTC")) %>% #tell it that it is in UTC time zone
  mutate(month = month(obs_time), #and separate out useful components of the time (this and line above largely from Michelle's assign_db_gpx function)
         day = day(obs_time),
         hour = hour(obs_time),
         min = minute(obs_time),
         sec = second(obs_time),
         year = year(obs_time)) %>%
  separate(time, into = c("date", "time"), sep = " ") #pull out date separately as a chr string too

# # prelim dispersal kernel parameters from Katrina (from email sent 7/5/18)
# k_2012 = 0.1
# theta_2012 = 5
# k_2013 = 0.025
# theta_2013 = 2 #or 1
# k_2014 = 0.1
# theta_2014 = 2
# k_2015 = 0.4
# theta_2015 = 1

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

# Eggs by female size function (from Adam Y's work), could change to have default intercept or slope in there once have those better estimated
eggsBySize <- function(intercept, slope, fish_size) {
  eggs <- intercept + slope*fish_size
  return(eggs)
}

# # Calculate LEP - consider updating later to include some of the calculating of survivals and such in here instead of as inputs? Probably better as separate function, though
# findLEP <- function(l_vec, f_vec) { #l_vec: vector of survivals-to-age, f_vec: vector of fecundities-at-age. Need to make sure first age of reproduction and max are are included in bounds, otherwise can define starting age as would like.
#   out <- sum(l_vec*f_vec)
#   return(out)
# }
# 
# # Calculate LEP via an IPM (like Will did in LEP_calc_WillWhite.R)

# Calculate self-persistence - one version
findSP_v1 <- function(LEP, home_recruits, egg_prod) { #LEP: lifetime egg production (see above), home_recruits: number o
  SP <- LEP * (home_recruits/egg_prod)
  return(SP)
}

# Another version of calculating self-persistence
findSP_v2 <- function(LEP, R_E_slope, prob_disp_home) { #LEP: lifetime egg production (see above), R_E_slope: slope converting eggs to recruits next year, prob_disp_home: probability of a recruit dispersing home
  SP <- LEP * R_E_slope * prob_disp_home
  return(SP)
}

# Calculate dispersal kernels (would need to be changed if theta values changed)
disp_allyears <- function(d) {  # theta = 0.5, equation for p(d) in eqn. 6c in Bode et al. 2018
  z = exp(k_allyears)
  disp = (z/2)*exp(-(z*d)^(theta_allyears))
  return(disp)
}

disp_2012 <- function(d) {  # theta = 3, equation for p(d) in eqn. 6b in Bode et al. 2018
  z = exp(k_2012)
  disp = (3*z)/(gamma(1/theta_2012))*exp(-(z*d)^theta_2012)
  return(disp)
}

disp_2013 <- function(d) {  # theta = 3, equation for p(d) in eqn. 6b in Bode et al. 2018
  z = exp(k_2013)
  disp = (3*z)/(gamma(1/theta_2013))*exp(-(z*d)^theta_2013)
  return(disp)
}

disp_2014 <- function(d) {  # theta = 2, equation for p(d) in eqn. 6a in Bode et al. 2018
  z = exp(k_2014)
  disp = ((2*z)/gamma(1/theta_2014))*exp(-(z*d)^theta_2014)
  return(disp)
} 

disp_2015 <- function(d) {  # theta = 2, equation for p(d) in eqn. 6a in Bode et al. 2018
  z = exp(k_2015)
  disp = ((2*z)/gamma(1/theta_2015))*exp(-(z*d)^theta_2015)
  return(disp)
} 

# # Calculate dispersal kernel (equation from Katrina email 7/5/18, based on Bode et al. 2017 paper)
# dispKernel <- function(k, theta, d) {
#   disp <- k*exp(-(k*d)^theta)
#   return(disp)
# }
# 
# # Dispersal kernel 2013
# dispKernel2013 <- function(d) {
#   k_2013*exp(-(k_2013*d)^theta_2013)
# }
# 
# # Dispersal kernel 2014
# dispKernel2014 <- function(d) {
#   k_2014*exp(-(k_2014*d)^theta_2014)
# }
# 
# # Dispersal kernel 2015
# dispKernel2015 <- function(d) {
#   k_2015*exp(-(k_2015*d)^theta_2015)
# }
# 
# # Dispersal kernel overall
# dispKernelallyears <- function(d) {
#   k_allyears*exp(-(k_allyears*d)^theta_allyears)
# }
# 
# 
# # Integrate dispersal kernel  (for theta = 1)
# integrateDK_theta1 <- function(d1,d2,k) {
#   int_d2 <- -exp(-k*d2)
#   int_d1 <- -exp(-k*d1)
#   int_out <- int_d2-int_d1
#   return(int_out)
# }

# Find lat lon for an anem (very similar function to anemid_latlong in AnemLocations.R - should probably combine the two and put in my common constants and code script...)
anemid_latlong_2 <- function(anem.table.id, anemdf, latlondata) { #anem.table.id is one anem_table_id value, anem.df is the anem.Processed data frame (so don't have to pull from db again here), latlondata is table of GPX data from database (rather than making the function call it each time); will need to think a bit more clearly about how to handle different locations read for different visits to the same anem_id (or different with same anem_obs); for now, just letting every row in anem.Info get a lat-long
  
  anem <- anemdf %>% filter(anem_table_id == anem.table.id) #get the relevant dive, time, site, etc. info for this anem_table_id
  
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

# for findSP_v2 - inputs are LEP, R_E_slope, prob_disp_home - find recruits from eggs
findRfromE <- function(m,eggs,b) {
  recruits = m*eggs + b
  return(recruits)
}

#################### Running things: ####################
# ########## Pull info from database
# #leyte <- read_db("Leyte")
# leyte <- read_db("Leyte")
# 
# allfish_fish <- leyte %>% 
#   tbl("clownfish") %>%
#   select(fish_table_id, anem_table_id, fish_spp, sample_id, gen_id, anem_table_id, recap, tag_id, color, size) %>%
#   collect() 
# 
# allfish_anems <- leyte %>%
#   tbl("anemones") %>%
#   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
#   collect() %>%
#   filter(anem_table_id %in% allfish_fish$anem_table_id)
# 
# allfish_dives <- leyte %>%
#   tbl("diveinfo") %>%
#   select(dive_table_id, dive_type, date, site, gps) %>%
#   collect() %>%
#   filter(dive_table_id %in% allfish_anems$dive_table_id)
# 
# # pull out just the year and put that in a separate column
# allfish_dives$year <- as.integer(substring(allfish_dives$date,1,4))
# 
# #join together
# allfish <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
# allfish <- left_join(allfish, allfish_dives, by="dive_table_id")
# 
# allfish$size <- as.numeric(allfish$size) #make size numeric (rather than a chr) so can do means and such
# 
# #pull GPS info
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

# #pull out just tagged fish
# taggedfish <- allfish %>% filter(!is.na(tag_id))

########## Calculating values, constants, inputs from data
# ##### Find max age (A) - maximum number of times a fish has been caught
# tag_years <- taggedfish %>% 
#   group_by(tag_id) %>%
#   summarize(firstyear = min(year), lastyear = max(year), nyears = max(year) - min(year))
# 
# maxyears = max(tag_years$nyears)
# 
# # add site back in
# tag_years_site <- left_join(tag_years, (taggedfish %>% select(tag_id, site)) %>% distinct(), by="tag_id") # 4 tags correspond to two sites - looked at below, emailed Michelle about checking data sheets

# # Note to self - if you just do distinct() above, it has 4 more rows - so looks like 4 tags that have multiple sites associated with them - should check into at some point...
# tagsatsites <- data.frame(table(tag_years_site$tag_id))
# tagsatsites %>% filter(Freq > 1)
# tag_years_site %>% filter(tag_id == "982000411818588")
# tag_years_site %>% filter(tag_id == "982000411818610")
# tag_years_site %>% filter(tag_id == "985153000401241")
# tag_years_site %>% filter(tag_id == "986112100172501")

# At some point, go back and figure this out - check to see if can find that fish at Visca and if it even is the same fish we catch each time on that STME
# # trying this with cap_id, since right now the max number of years recaught is the max possible for the # of years we've been tagging
# cap_years <- allfish %>% 
#   filter(!is.na(cap_id)) %>%
#   group_by(cap_id) %>%
#   summarize(firstyear = min(year), lastyear = max(year), nyears = max(year) - min(year), tag = tag_id[1])
# 
# maxyears_cap = max(cap_years$nyears)
# 
# # Do any fish in 2017 or 2018 have cap_ids?
# # not when I try to find them like this...
# test_caps <- allfish %>% 
#   filter(!is.na(cap_id)) %>%
#   filter(year %in% c(2017,2018))
# 
# # or like this
# table((allfish %>% filter(year %in% c(2017,2018)))$cap_id)
# 
# # or like this
# tags_in_cap_years <- cap_years %>% filter(!is.na(tag)) %>% select(cap_id, tag)
# test_caps2 <- allfish %>%
#   filter(cap_id %in% tags_in_cap_years$cap_id) 
# table(test_caps2$year)
# 
# # maybe cap_ids haven't been populated yet for 2017 and 2018 (genetics not done for those years yet but could be fish sequenced and tagged pre-2016 and known to be recaptured via tag)?
# # try to find that one old Visca fish to see - haven't quite gotten this working yet...
# test_Visca <- allfish %>%
#   filter(site == "Visca") %>%
#   filter(!is.na(tag_id)) %>%
#   filter(year == 2016)
# 
# test_Visca2 <- allfish %>%
#   filter(anem_spp)
# 
# anems_Visca <- leyte %>%
#   tbl("anemones") %>%
#   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_spp) %>%
#   collect() %>%
#   filter(anem_table_id %in% (allfish %>% filter(site == "Visca"))$anem_table_id) 


##### Find widths of sites (for SP calculation, used to calculate estimate home recruits from dispersal kernel)
# make a data frame of sites, boundary anems, lat lons, widths
site_width_info <- data.frame(site = site_vec) #initialize dataframe

# anem_id of northern-most anem at site (or eastern-most in case of Haina, which runs E-W)
site_width_info$N_anem <- c(Cabatoan_N, CaridadCemetery_N, CaridadProper_N, ElementarySchool_N, Gabas_N, Haina_E, HicgopSouth_N,
                            Magbangon_N_N, Magbangon_S_N, Palanas_N, PorocRose_N, PorocSanFlower_N, SanAgustin_N, SitioBaybayon_N, SitioLonas_N,
                            SitioTugas_N, TamakinDacot_N, Visca_N, Wangag_N)

# anem_id of southern-most anem at site (or western-most in case of Haina, which runs E-W)
site_width_info$S_anem <- c(Cabatoan_S, CaridadCemetery_S, CaridadProper_S, ElementarySchool_S, Gabas_S, Haina_W, HicgopSouth_S, Magbangon_N_S,
                            Magbangon_S_S, Palanas_S, PorocRose_S, PorocSanFlower_S, SanAgustin_S, SitioBaybayon_S, SitioLonas_S,
                            SitioTugas_S, TamakinDacot_S, Visca_S, Wangag_S)

# anem_id of mid anem at site (or western-most in case of Haina, which runs E-W) #except CardiadProper, Sitio Lonas, Sitio Tugas, used N
site_width_info$M_anem <- c(Cabatoan_mid, CaridadCemetery_mid, CaridadProper_N, ElementarySchool_mid, Gabas_mid, Haina_mid, HicgopSouth_mid,
                            Magbangon_N_mid, Magbangon_S_mid, Palanas_mid, PorocRose_mid, PorocSanFlower_mid, SanAgustin_mid, SitioBaybayon_mid, SitioLonas_N,
                            SitioTugas_N, TamakinDacot_mid, Visca_mid, Wangag_mid)

# pull out anem_table_ids for those anems so can use anemid_latlong2 function to find lat lon for boundary anems
site_width_anems <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_obs_time) %>%
  collect() %>%
  filter(anem_id %in% c(site_width_info$N_anem, site_width_info$S_anem, site_width_info$M_anem)) 

# pull out the dive info associated with those anems
site_width_dives <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, date, site, gps) %>%
  collect() %>%
  filter(dive_table_id %in% site_width_anems$dive_table_id)

# join into one data frame to put into MS function (anem_latlong), get the anem obs_time in the same time zone (UTC) as the GPX data, create placeholder lat and lon columns
site_width_anemdives <- left_join(site_width_anems, site_width_dives, by="dive_table_id") %>%
  mutate(obs_time = force_tz(ymd_hms(str_c(date, anem_obs_time, sep = " ")), tzone = "Asia/Manila")) %>% #tell it that it is currently in Asia/Manila time zone
  mutate(obs_time = with_tz(obs_time, tzone = "UTC")) %>% #convert to UTC so can compare with GPS data (this line and one above largely from Michelle's assign_db_gpx function)
  mutate(month = month(obs_time), #and separate out useful components of the time (this also from Michelle's assign_db_gpx function)
         day = day(obs_time), 
         hour = hour(obs_time), 
         min = minute(obs_time), 
         sec = second(obs_time), 
         year = year(obs_time)) %>%
  mutate(lat = as.numeric(NA), lon = as.numeric(NA)) #put in a placeholder column for lat and lon

# find lat lon for the boundary anems -- THIS DOESN'T WORK FOR THE MID ANEMS
for(i in 1:length(site_width_anemdives$anem_table_id)) {
  out_anemLL <- anemid_latlong_2(site_width_anemdives$anem_table_id[i], site_width_anemdives, gps.Info)  # SHOULD COMBINE/DECIDE BETWEEN THE TWO ANEM_LATLON FUNCTIONS AT SOME POINT
  #out_anemLL <- anemid_latlong(site_width_anemdives$anem_table_id, site_width_anemdives, gps.Info)
  site_width_anemdives$lat[i] = out_anemLL$lat[1] #figure out why this is sometimes the wrong replacement length... (put [1] in to try to solve that issue but not sure why it's needed)
  site_width_anemdives$lon[i] = out_anemLL$lon[1]
}

# just pull out one row for each anem_id 
site_width_anemdives_short <- site_width_anemdives %>%
  distinct(anem_id, .keep_all = TRUE)

# looks like anem_id is a chr now? switch to numeric
site_width_anemdives_short <- site_width_anemdives_short %>%
  mutate(anem_id = as.numeric(anem_id))

# put lat and lons into info table (requires some column renaming so needs to be run in order)
# lat/lon for N_anem
site_width_info <- left_join(site_width_info %>% dplyr::rename(anem_id = N_anem), site_width_anemdives_short %>% select(anem_id, lat, lon), by="anem_id") #join in the lat lons for the N anem
site_width_info <- site_width_info %>% dplyr::rename(N_anem = anem_id, N_lat = lat, N_lon = lon) #rename the anem_id, lat, lon columns to be N-specific
# lat/lon for S_anem
site_width_info <- left_join(site_width_info %>% dplyr::rename(anem_id = S_anem), site_width_anemdives_short %>% select(anem_id, lat, lon), by="anem_id") #join in the lat lons for the S anem
site_width_info <- site_width_info %>% dplyr::rename(S_anem = anem_id, S_lat = lat, S_lon = lon) #rename the anem_id, lat, lon columns to be S-specific 
# lat/lon for M_anem
site_width_info <- left_join(site_width_info %>% dplyr::rename(anem_id = M_anem), site_width_anemdives_short %>% select(anem_id, lat, lon), by="anem_id") #join in the lat lons for the mid anem
site_width_info <- site_width_info %>% dplyr::rename(M_anem = anem_id, M_lat = lat, M_lon = lon) #rename the anem_id, lat, lon columns to be M-specific 

#and calculate the distance!
for(i in 1:length(site_width_info$site)) {
  site_width_info$width_m[i] = distHaversine(c(site_width_info$N_lon[i], site_width_info$N_lat[i]), c(site_width_info$S_lon[i], site_width_info$S_lat[i]))
  site_width_info$width_km[i] = site_width_info$width_m[i]/1000 
}

# Find distances among sites
# Set up data frame
site_dist_info <- data.frame(org_site = c(rep(site_vec[1],length(site_vec)), rep(site_vec[2],length(site_vec)),
                                          rep(site_vec[3],length(site_vec)), rep(site_vec[4],length(site_vec)),
                                          rep(site_vec[5],length(site_vec)), rep(site_vec[6],length(site_vec)),
                                          rep(site_vec[7],length(site_vec)), rep(site_vec[8],length(site_vec)),
                                          rep(site_vec[9],length(site_vec)), rep(site_vec[10],length(site_vec)),
                                          rep(site_vec[11],length(site_vec)), rep(site_vec[12],length(site_vec)),
                                          rep(site_vec[13],length(site_vec)), rep(site_vec[14],length(site_vec)),
                                          rep(site_vec[15],length(site_vec)), rep(site_vec[16],length(site_vec)),
                                          rep(site_vec[17],length(site_vec)), rep(site_vec[18],length(site_vec)), 
                                          rep(site_vec[19],length(site_vec))), stringsAsFactors = FALSE)
site_dist_info$dest_site <- rep(site_vec, length(site_vec))
site_dist_info$dist_N_to_S_m <- rep(NA, length(site_dist_info$org_site))
site_dist_info$dist_N_to_S_km <- rep(NA, length(site_dist_info$org_site))
site_dist_info$dist_S_to_N_m <- rep(NA, length(site_dist_info$org_site))
site_dist_info$dist_S_to_N_km <- rep(NA, length(site_dist_info$org_site))
site_dist_info$dist_avg <- rep(NA, length(site_dist_info$org_site))
site_dist_info$dest_width <- rep(NA, length(site_dist_info$org_site))
site_dist_info$d1_km <- rep(NA, length(site_dist_info$org_site))
site_dist_info$d2_km <- rep(NA, length(site_dist_info$org_site))

# and find site distances! (for now, doing N to N, should do mid to mid...)
for(i in 1:length(site_dist_info$org_site)) {
  site_org = site_dist_info$org_site[i]
  N_lat_org = (site_width_info %>% filter(site == site_org))$N_lat 
  N_lon_org = (site_width_info %>% filter(site == site_org))$N_lon
  S_lat_org = (site_width_info %>% filter(site == site_org))$S_lat 
  S_lon_org = (site_width_info %>% filter(site == site_org))$S_lon
  site_dest = site_dist_info$dest_site[i]
  N_lat_dest = (site_width_info %>% filter(site == site_dest))$N_lat
  N_lon_dest = (site_width_info %>% filter(site == site_dest))$N_lon
  S_lat_dest = (site_width_info %>% filter(site == site_dest))$S_lat
  S_lon_dest = (site_width_info %>% filter(site == site_dest))$S_lon
  site_dist_info$dist_N_to_S_m[i] = distHaversine(c(N_lon_org, N_lat_org), c(S_lon_dest, S_lat_dest))
  site_dist_info$dist_N_to_S_km[i] = site_dist_info$dist_N_to_S_m[i]/1000 
  site_dist_info$dist_S_to_N_m[i] = distHaversine(c(N_lon_dest, N_lat_dest), c(S_lon_org, S_lat_org))
  site_dist_info$dist_S_to_N_km[i] = site_dist_info$dist_S_to_N_m[i]/1000 
  
  site_dist_info$dist_avg[i] = (site_dist_info$dist_N_to_S_km[i] + site_dist_info$dist_S_to_N_km[i])/2
  site_dist_info$dest_width[i] <- (site_width_info %>% filter(site == site_dest))$width_km
  
  site_dist_info$d1_km[i] <- site_dist_info$dist_avg[i] - site_dist_info$dest_width[i]/2
  site_dist_info$d2_km[i] <- site_dist_info$dist_avg[i] + site_dist_info$dest_width[i]/2
  # if(site_dist_info$dist_N_to_S_km[i] >= site_dist_info$dist_S_to_N_km[i]){
  #   site_dist_info$d1_km[i] <- site_dist_info$dist_S_to_N_km[i]
  #   site_dist_info$d2_km[i] <- site_dist_info$dist_N_to_S_km[i]
  # } else {
  #   site_dist_info$d1_km[i] <- site_dist_info$dist_N_to_S_km[i]
  #   site_dist_info$d2_km[i] <- site_dist_info$dist_S_to_N_km[i]
  # }
}

# For sites going to self, put N-S width in
for (i in 1:length(site_dist_info$org_site)) {
  site_org = site_dist_info$org_site[i]
  site_dest = site_dist_info$dest_site[i]
  if (site_org == site_dest) {
    N_lat_org = (site_width_info %>% filter(site == site_org))$N_lat 
    N_lon_org = (site_width_info %>% filter(site == site_org))$N_lon
    S_lat_org = (site_width_info %>% filter(site == site_org))$S_lat 
    S_lon_org = (site_width_info %>% filter(site == site_org))$S_lon
    site_dist_info$dist_N_to_S_m[i] = distHaversine(c(N_lon_org, N_lat_org), c(S_lon_org, S_lat_org))
    site_dist_info$dist_N_to_S_km[i] = site_dist_info$dist_N_to_S_m[i]/1000 
    site_dist_info$dist_S_to_N_m[i] = distHaversine(c(N_lon_org, N_lat_org), c(S_lon_org, S_lat_org))
    site_dist_info$dist_S_to_N_km[i] = site_dist_info$dist_S_to_N_m[i]/1000 
    site_dist_info$d1_km[i] = 0
    site_dist_info$d2_km[i] = site_dist_info$dist_N_to_S_km[i]
  }
}

########## Looking at eggs and recruits relationships
# Egg-recruit slope with LEP - the models don't seem to be loading correctly... give different values when use summary() here vs. EggRecruitRelationship.R - check into this!!
# intercept_mod1 <- 1.856e+01 #need to figure out how to extract this from egg_recruits_est_mod1  # these commented values are what show up when use summary() on RData file loaded in this script...
# slope_mod1 <- 5.430e-05 #need to figure out how to extract this from egg_recruits_est_mod1

# these values written in from egg_recruits_est_metalTA_m2_mod1 when summary() called in EggRecruitsRelationship.R - Adjusted R2 of 0.1369, both intercept (0.05) and slope (0.01) significant
intercept_mod1 <- 1.203e-03
slope_mod1 <- 3.117e-05

# these values written in from egg_recruits_est_metalTA_0intercept_m2_mod6 when summary() called in EggRecruitsRelationship.R - Adjusted R2 of 0.4682, slope (0.0001) significant
intercept_mod6 <- 0
slope_mod6 <- 4.711e-05

# These were taken from loaded-in file, also probably wrong, though haven't checked...
# intercept_meta_mod1 <- -2.657e+02
# slope_meta_mod1 <- 6.624e-05

recruits_per_LEP <- findRfromE(slope_mod1, LEP, intercept_mod1)  # 18! Means persistence is possible! THIS SEEMS SUPER WRONG... pretty high... and now down to 0.007, very low...
recruits_per_LEP_mod6 <- findRfromE(slope_mod6, LEP, intercept_mod6)  # 0.009...

recruits_per_LEP_Will_eggs <- findRfromE(slope_mod1, LEP_Will_eggs_per_clutch, intercept_mod1)  # using Will eggs-per-clutch estimate
recruits_per_LEP_Will_LEP <- findRfromE(slope_mod1, LEP_Will, intercept_mod1)  # using Will LEP estimate
#recruits_per_LEP_meta <- findRfromE(slope_meta_mod1, LEP, intercept_meta_mod1)

plot(0:2000, intercept_mod1 + slope_mod1*(0:2000))

########## Pull out summary # recruits and # eggs by site
# find summary N recruits and N eggs by site, averaged across years sampled (also take mid, high, low?) - lots of NAs in recruits_info in prob_hab sampled... should check into that (maybe b/c 0 metal-tagged anems?)
# recruits info - something is wrong with a the N. and S. Magbangon proportion habitat sampled.....
demog_info_recruits <- recruits_info %>%
  group_by(site) %>%
  summarize(mean_est_R = mean(totalR_est_metalTA, na.rm=TRUE),
            high_est_R = max(totalR_est_metalTA, na.rm=TRUE),
            low_est_R = min(totalR_est_metalTA, na.rm=TRUE),
            mean_raw_R = mean(Nrecruits, na.rm=TRUE),
            high_raw_R = max(Nrecruits, na.rm=TRUE),
            low_raw_R = min(Nrecruits, na.rm=TRUE))

demog_info_eggs <- breedingF_info %>%
  group_by(site) %>%
  summarize(mean_est_eggs = mean(est_eggs_metalTA, na.rm=TRUE),
            high_est_eggs = max(est_eggs_metalTA, na.rm=TRUE),
            low_est_eggs = min(est_eggs_metalTA, na.rm=TRUE),
            mean_raw_F = mean(NbreedingF_combo, na.rm=TRUE),
            high_raw_F = max(NbreedingF_combo, na.rm=TRUE),
            low_raw_F = min(NbreedingF_combo, na.rm=TRUE))

########## Converting migEst "migration rates" (proportion of settlers at dest that came from each of the other sites) into pij values
migEst_pijmat <- site_dist_info %>% select(org_site, dest_site) 
migEst_pijmat <- left_join(migEst_pijmat, site_vec_order, by=c("org_site" = "site_name")) #coerces factor into character vector...
migEst_pijmat <- migEst_pijmat %>%
  dplyr::rename(org_alpha_order = alpha_order, org_geo_order = geo_order) %>% #rename order columns so clear that they are for the origin site
  mutate(mig_rate = rep(NA, length(org_site)), prop_disp = rep(NA, length(org_site)),  #columns to put in MigEst migration rates and the conversions to pijs
         mig_rate_LCI = rep(NA, length(org_site)), mig_rate_UCI = rep(NA, length(org_site)), #columns for upper and lower confidence intervals of migEst output
         prop_disp_LCI = rep(NA, length(org_site)), prop_disp_UCI = rep(NA, length(org_site)), #columns for upper and lower confidence intervals of converted pijs
         dest_alpha_order = rep(1:19, 19), #destination site number (by alpha order)
         dest_geo_order = rep(site_vec_order$geo_order, 19)) #destination site number by geographical order 

# Put the migEst estimates in the new data frame - cycle through the org and dest sites 
for(i in 1:19) { #19 sites right now
  org_num = i
  for(j in 1:19) {
    dest_num = j
    
    # pull out migEst estimates
    migRate = migEst_conn[j,i] # migration rate from MigEst
    migRate_LCI = migEst_conn_LCI[j,i] # lower confidence interval estimate of migration rate from MigEst
    migRate_UCI = migEst_conn_UCI[j,i] # upper confidence interval estimate of migration rate from MigEst
    
    # put them in the right row for origin and destination
    migEst_pijmat$mig_rate[which(migEst_pijmat$org_alpha_order == org_num & migEst_pijmat$dest_alpha_order == dest_num)] = migRate #migration rate
    migEst_pijmat$mig_rate_LCI[which(migEst_pijmat$org_alpha_order == org_num & migEst_pijmat$dest_alpha_order == dest_num)] = migRate_LCI #lower bound of confidence interval for migration rate
    migEst_pijmat$mig_rate_UCI[which(migEst_pijmat$org_alpha_order == org_num & migEst_pijmat$dest_alpha_order == dest_num)] = migRate_UCI #upper bound of confidence interval for migration rate
  }
}

# Convert migEst output to proportion of recruits from site i settling at site j
for(i in 1:length(migEst_pijmat$org_site)){
  org_site = migEst_pijmat$org_site[i]
  dest_site = migEst_pijmat$dest_site[i]
  mig_rate = migEst_pijmat$mig_rate[i]
  mig_rate_LCI = migEst_pijmat$mig_rate_LCI[i]
  mig_rate_UCI = migEst_pijmat$mig_rate_UCI[i]
  recruits_to_dest = (demog_info_recruits %>% filter(site == dest_site))$mean_est_R
  eggs_from_org = (demog_info_eggs %>% filter(site == org_site))$mean_est_eggs
  recruits_from_org = findRfromE(slope_mod1, eggs_from_org, intercept_mod1)

  migEst_pijmat$prop_disp[i] = (mig_rate * recruits_to_dest)/recruits_from_org
  migEst_pijmat$prop_disp_LCI[i] = (mig_rate_LCI * recruits_to_dest)/recruits_from_org
  migEst_pijmat$prop_disp_UCI[i] = (mig_rate_UCI * recruits_to_dest)/recruits_from_org
}

##### ADD IN HIGH AND LOW ESTIMATES OF PROB_DISP WITH HIGH/LOW est R and est eggs - all combos?

##### Scale up parentage numbers - want number of recruits arriving home to patch i - find by taking number of parentage matches going home and scale up by percentage hab sampled
# combine parentage files (mums, dads, trios) - first rename columns so they match across the files, add a column for match type, then rbind
parentage_dads <- parentage_dads %>% 
  dplyr::rename(parent_site = par2_site, nmatches = n_dad, offspring_site = offs_site) %>% 
  mutate(match_type = rep('dad', dim(parentage_dads)[1]))
parentage_moms <- parentage_moms %>% 
  dplyr::rename(parent_site = par1_site, nmatches = n_mum, offspring_site = offs_site) %>%
  mutate(match_type = rep('mom', dim(parentage_moms)[1]))
parentage_trios <- parentage_trios %>% 
  dplyr::rename(parent_site = par1_site, nmatches = n_trios, offspring_site = offs_site) %>%
  mutate(match_type = rep('trio', dim(parentage_trios)[1]))

parentage_matches_raw <- rbind(parentage_dads, parentage_moms, parentage_trios)

# select out just the self matches, group by site and year, and sum across match types
parentage_matches_self <- parentage_matches_raw %>% 
  filter(offspring_site == parent_site) %>% 
  group_by(year, offspring_site) %>%
  summarize(n_matches = sum(nmatches))

# join with proportion habitat sampled
parentage_matches_self <- left_join(parentage_matches_self, anems_table %>% select(prop_hab_sampled_metal_TA, year, site), by=c('year' = 'year', 'offspring_site' = 'site'))

# join with egg production numbers
parentage_matches_self <- left_join(parentage_matches_self, breedingF_info %>% select(site, year, est_eggs_metalTA), by = c('year' = 'year', 'offspring_site' = 'site'))

# # Don't need to do this anymore b/c have new parentage files now with N and S Mag separated
# # WHILE PARENTAGE MATCHES ARE JUST MAG, MAKE PROP HAB SAMPLED THE AVERAGE OF N MAG AND S MAG 
# prop_hab_Mag2013 <- mean((anems_table %>% filter(year == 2013, site %in% c('N. Magbangon', 'S. Magbangon')))$prop_hab_sampled_metal_TA)
# prop_hab_Mag2014 <- mean((anems_table %>% filter(year == 2014, site %in% c('N. Magbangon', 'S. Magbangon')))$prop_hab_sampled_metal_TA)
# prop_hab_Mag2015 <- mean((anems_table %>% filter(year == 2015, site %in% c('N. Magbangon', 'S. Magbangon')))$prop_hab_sampled_metal_TA)
# 
# parentage_matches_self$prop_hab_sampled_metal_TA[2] = prop_hab_Mag2013
# parentage_matches_self$prop_hab_sampled_metal_TA[7] = prop_hab_Mag2014
# parentage_matches_self$prop_hab_sampled_metal_TA[12] = prop_hab_Mag2015

# scale up raw matches by proportion habitat sampled
parentage_matches_self <- parentage_matches_self %>%
  mutate(nrecruits_scaled = n_matches/prop_hab_sampled_metal_TA,
         nrecruits_rounded = round(n_matches/prop_hab_sampled_metal_TA)) 

########## Assess survival from egg-recruit in a new way - use total parentage matches/N parents genotypes
# Total number of offspring identified
n_offspring_parentage <- sum(parentage_matches_raw$nmatches)  # all offspring identified via parentage
n_parents_parentage_df <- allfish_caught %>%
  filter(!is.na(gen_id)) %>%
  filter(size >= min_breeding_M_size | color == 'YP' | color == 'O') %>%
  distinct(gen_id, .keep_all = TRUE) 
n_parents_parentage <- length(n_parents_parentage_df$gen_id)

# Sum up total site area (all site areas times all years sampled) - total possible sampling area
sites_for_total_areas <- c('Cabatoan', 'Caridad Cemetery', 'Elementary School', 'Gabas', 
                           'Haina', 'Hicgop South', 'N. Magbangon', 'Palanas', 'Poroc Rose',
                           'Poroc San Flower', 'San Agustin', 'Sitio Baybayon', 'Tamakin Dacot',
                           'Visca', 'Wangag', 'S. Magbangon')  # excluding Caridad Proper, Sitio Lonas, Sitio Tugas from this (should I exclude the Sitio Lonas match too?)
site_areas_modified <- site_areas %>% filter(site %in% sites_for_total_areas)

total_area_all_years <- sum(site_areas_modified$kmsq_area)*length(years_sampled)

# Find proportion of habitat sampled overall - sum of area sampled in each
sampled_area_each_year <- left_join(site_areas_modified, anems_table %>% select(site, year, prop_hab_sampled_metal_TA), by = 'site')
sampled_area_each_year <- sampled_area_each_year %>%
  mutate(area_sampled = prop_hab_sampled_metal_TA*kmsq_area)
total_area_sampled <- sum(sampled_area_each_year$area_sampled)

total_prop_hab_sampled <- total_area_sampled/total_area_all_years

# Estimate survival from adult-recruit
prop_F_M <- 0.5  # saying 50% of the "adults" we clip are males that won't make it to females -- reasonable? could check this. But LEP takes that into account, right?
tagged_offspring <- n_parents_parentage*LEP
recruited_offspring <- n_offspring_parentage/total_prop_hab_sampled
surv_egg_recruit <- recruited_offspring/tagged_offspring
  
#How does that compare to what we had from the estimated egg-recruit relationship?
surv_RfromE <- findRfromE(slope_mod1, 1, intercept_mod1)

###### Re-Run METRICS! BUT NEED TO RE-RUN LEP FIRST!

########## Assessing metrics

##### Do pops look persistence just based on the number of recruits? (estimating this by calculating SP using all the recruits coming in)
SP_allRecruits <- left_join(demog_info_recruits, demog_info_eggs, by = 'site') %>%
  mutate(SP_mean_est = LEP*(mean_est_R/mean_est_eggs),
         SP_high_est = LEP*(high_est_R/mean_est_eggs))

SP_allRecruits_2 <- left_join(demog_info_recruits, demog_info_eggs, by = 'site') %>%
  mutate(SP_mean_est = 5000*(mean_est_R/mean_est_eggs),
         SP_high_est = 5000*(high_est_R/mean_est_eggs))

##### Self-persistence, using migEst estimates
SP_migEst <- migEst_pijmat %>%
  filter(org_site == dest_site) %>% #just the self-self info
  mutate(prop_disp = if_else(prop_disp > 1, 1, prop_disp), #if converted prop_disp is > 1, put at 1
         prop_disp_LCI = if_else(prop_disp_LCI > 1, 1, prop_disp_LCI), #if converted prop_disp_LCI is > 1, put at 1
         prop_disp_UCI = if_else(prop_disp_UCI > 1, 1, prop_disp_UCI), #if converted prop_disp_LCI is > 1, put at 1
         SP = prop_disp*recruits_per_LEP,
         SP_mrLCI = prop_disp_LCI*recruits_per_LEP,
         SP_mrUCI = prop_disp_UCI*recruits_per_LEP) 

# Separate out migEst point, upper, and lower CI intervals so can rbind all the SP estimates together
SP_migEst_est <- SP_migEst %>%
  select(org_site, org_alpha_order, org_geo_order, SP) %>%
  dplyr::rename(site = org_site, alpha_order = org_alpha_order, geo_order = org_geo_order) %>%
  mutate(est_method = 'migEst_pointEst',
         data = 'migEst',
         year = 'combined')

SP_migEst_LCI <- SP_migEst %>%
  select(org_site, org_alpha_order, org_geo_order, SP_mrLCI) %>%
  dplyr::rename(site = org_site, alpha_order = org_alpha_order, geo_order = org_geo_order, SP = SP_mrLCI) %>%
  mutate(est_method = 'migEst_LCI',
         data = 'migEst',
         year = 'combined')

SP_migEst_UCI <- SP_migEst %>%
  select(org_site, org_alpha_order, org_geo_order, SP_mrUCI) %>%
  dplyr::rename(site = org_site, alpha_order = org_alpha_order, geo_order = org_geo_order, SP = SP_mrUCI) %>%
  mutate(est_method = 'migEst_UCI',
         data = 'migEst',
         year = 'combined')

##### Self-persistence, using raw and scaled up parentage matches
SP_parentage <- parentage_matches_self %>%
  mutate(SP_raw = LEP*(n_matches/est_eggs_metalTA),
         SP_scaled = LEP*(nrecruits_rounded/est_eggs_metalTA))

# Join with site_vec_order to get geo order column
SP_parentage <- left_join(SP_parentage, site_vec_order, by = c("offspring_site" = "site_name"))

# Separate out raw and scaled so can join all SP methods together
SP_parentage_raw <- SP_parentage %>%
  select(year, offspring_site, SP_raw, alpha_order, geo_order) %>%
  dplyr::rename(site = offspring_site, SP = SP_raw) %>%
  mutate(est_method = 'parentage_raw',
         data = 'parentage')

SP_parentage_scaled <- SP_parentage %>%
  select(year, offspring_site, SP_scaled, alpha_order, geo_order) %>%
  dplyr::rename(site = offspring_site, SP = SP_scaled) %>%
  mutate(est_method = 'parentage_scaled',
         data = 'parentage')

##### Self-persistence, using dispersal kernels
SP_kernels <- site_dist_info %>%
  filter(org_site == dest_site) %>%
  select(org_site, dest_width) %>%
  mutate(pdisp_allyears = as.numeric(NA),
         pdisp_2012 = as.numeric(NA),
         pdisp_2013 = as.numeric(NA),
         pdisp_2014 = as.numeric(NA),
         pdisp_2015 = as.numeric(NA))

# Fill in prob of dispersing (tried using mutate instead but integrate doesn't work well within dplyr)
for(i in 1:length(SP_kernels$org_site)) {
  SP_kernels$pdisp_allyears[i] = integrate(disp_allyears, 0, SP_kernels$dest_width[i])$value
  SP_kernels$pdisp_2012[i] = integrate(disp_2012, 0, SP_kernels$dest_width[i])$value
  SP_kernels$pdisp_2013[i] = integrate(disp_2013, 0, SP_kernels$dest_width[i])$value
  SP_kernels$pdisp_2014[i] = integrate(disp_2014, 0, SP_kernels$dest_width[i])$value
  SP_kernels$pdisp_2015[i] = integrate(disp_2015, 0, SP_kernels$dest_width[i])$value
}

# Calculate SP
SP_kernels <- SP_kernels %>%
  mutate(SP_allyears = pdisp_allyears*recruits_per_LEP,
         SP_2012 = pdisp_2012*recruits_per_LEP,
         SP_2013 = pdisp_2013*recruits_per_LEP,
         SP_2014 = pdisp_2014*recruits_per_LEP,
         SP_2015 = pdisp_2015*recruits_per_LEP)

# Join with site_vec_info
SP_kernels <- left_join(SP_kernels, site_vec_order, by = c('org_site' = 'site_name'))

# Separate out SP by year so can bind together with other SP estimates
SP_kernels_allyears <- SP_kernels %>%
  select(org_site, SP_allyears, alpha_order, geo_order) %>%
  dplyr::rename(site = org_site, SP = SP_allyears) %>%
  mutate(year = 'combined',
         est_method = 'all years kernel',
         data = 'kernel')
  
SP_kernels_2012 <- SP_kernels %>%
  select(org_site, SP_2012, alpha_order, geo_order) %>%
  dplyr::rename(site = org_site, SP = SP_2012) %>%
  mutate(year = '2012',
         est_method = '2012 kernel',
         data = 'kernel')

SP_kernels_2013 <- SP_kernels %>%
  select(org_site, SP_2013, alpha_order, geo_order) %>%
  dplyr::rename(site = org_site, SP = SP_2013) %>%
  mutate(year = '2013',
         est_method = '2013 kernel',
         data = 'kernel')

SP_kernels_2014 <- SP_kernels %>%
  select(org_site, SP_2014, alpha_order, geo_order) %>%
  dplyr::rename(site = org_site, SP = SP_2014) %>%
  mutate(year = '2014',
         est_method = '2014 kernel',
         data = 'kernel')

SP_kernels_2015 <- SP_kernels %>%
  select(org_site, SP_2015, alpha_order, geo_order) %>%
  dplyr::rename(site = org_site, SP = SP_2015) %>%
  mutate(year = '2015',
         est_method = '2015 kernel',
         data = 'kernel')

##### Bind self-persistence estimates together so can easily plot to compare
SP_all <- rbind(SP_migEst_est, SP_migEst_LCI, SP_migEst_UCI, 
                as.data.frame(SP_parentage_raw), as.data.frame(SP_parentage_scaled),
                SP_kernels_allyears, SP_kernels_2012, SP_kernels_2013,
                SP_kernels_2014, SP_kernels_2015)

##### Connectivity matrix, using dispersal kernels
# Set up matrices
conn_mat_2012K <- matrix(NA, nrow = length(site_vec_order$site_name), ncol = length(site_vec_order$site_name))
conn_mat_2013K <- matrix(NA, nrow = length(site_vec_order$site_name), ncol = length(site_vec_order$site_name))
conn_mat_2014K <- matrix(NA, nrow = length(site_vec_order$site_name), ncol = length(site_vec_order$site_name))
conn_mat_2015K <- matrix(NA, nrow = length(site_vec_order$site_name), ncol = length(site_vec_order$site_name))
conn_mat_allyearsK <- matrix(NA, nrow = length(site_vec_order$site_name), ncol = length(site_vec_order$site_name))

site_dist_info <- site_dist_info %>%
  mutate(prob_disp_2012 = as.numeric(NA),
         prob_disp_2013 = as.numeric(NA),
         prob_disp_2014 = as.numeric(NA),
         prob_disp_2015 = as.numeric(NA),
         prob_disp_allyears = as.numeric(NA))

# Find prob of dispersing within the width of the site - just using integrate function here, should swap out for analytical integration...
for(i in 1:length(site_dist_info$org_site)){ 
  site_dist_info$prob_disp_2012[i] = integrate(disp_2012, site_dist_info$d1_km[i], site_dist_info$d2_km[i])$value 
  site_dist_info$prob_disp_2013[i] = integrate(disp_2013, site_dist_info$d1_km[i], site_dist_info$d2_km[i])$value 
  site_dist_info$prob_disp_2014[i] = integrate(disp_2014, site_dist_info$d1_km[i], site_dist_info$d2_km[i])$value 
  site_dist_info$prob_disp_2015[i] = integrate(disp_2015, site_dist_info$d1_km[i], site_dist_info$d2_km[i])$value 
  site_dist_info$prob_disp_allyears[i] = integrate(disp_allyears, site_dist_info$d1_km[i], site_dist_info$d2_km[i])$value 
}

# Add in origin site alpha and geo order
site_dist_info <- left_join(site_dist_info, site_vec_order, by = c('org_site' = 'site_name'))
site_dist_info <- site_dist_info %>%
  dplyr::rename(org_alpha_order = alpha_order, org_geo_order = geo_order)

# Add in destination site alpha and geo order
site_dist_info <- left_join(site_dist_info, site_vec_order, by = c('dest_site' = 'site_name'))
site_dist_info <- site_dist_info %>%
  dplyr::rename(dest_alpha_order = alpha_order, dest_geo_order = geo_order)

# Multiply by recruits_per_LEP to get realized connectivity
site_dist_info <- site_dist_info %>%
  mutate(realizedC_2012 = prob_disp_2012*recruits_per_LEP,
         realizedC_2013 = prob_disp_2013*recruits_per_LEP,
         realizedC_2014 = prob_disp_2014*recruits_per_LEP,
         realizedC_2015 = prob_disp_2015*recruits_per_LEP,
         realizedC_allyears = prob_disp_allyears*recruits_per_LEP)

# Create connectivity and realized connectivity matrices as matrices rather than data frames so can find eigenvalues
rCmat_2012 = matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))
rCmat_2013 = matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))
rCmat_2014 = matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))
rCmat_2015 = matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))
rCmat_allyears = matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))
rCmat_allyears_test =  matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))

Cmat_allyears_dK_justconn = matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))
Cmat_2012_dK_justconn = matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))
Cmat_2013_dK_justconn = matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))
Cmat_2014_dK_justconn = matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))
Cmat_2015_dK_justconn = matrix(NA, length(site_vec_order$site_name), length(site_vec_order$site_name))

for(i in 1:length(site_dist_info$org_site)) {
  column = site_dist_info$org_geo_order[i]  # column is origin 
  row = site_dist_info$dest_geo_order[i]  # row is destination
  rCmat_2012[row, column] = site_dist_info$realizedC_2012[i]
  rCmat_2013[row, column] = site_dist_info$realizedC_2013[i]
  rCmat_2014[row, column] = site_dist_info$realizedC_2014[i]
  rCmat_2015[row, column] = site_dist_info$realizedC_2015[i]
  rCmat_allyears[row, column] = site_dist_info$realizedC_allyears[i]
  
  Cmat_allyears_dK_justconn[row, column] = site_dist_info$prob_disp_allyears[i]
  Cmat_2012_dK_justconn[row, column] = site_dist_info$prob_disp_2012[i]
  Cmat_2013_dK_justconn[row, column] = site_dist_info$prob_disp_2013[i]
  Cmat_2014_dK_justconn[row, column] = site_dist_info$prob_disp_2014[i]
  Cmat_2015_dK_justconn[row, column] = site_dist_info$prob_disp_2015[i]
  
  rCmat_allyears_test[column, row] = site_dist_info$realizedC_allyears[i]  # threw this in just to make sure it didn't seem to matter which order the dest vs. org were (row vs. column) - doesn't seem to, has same eigenvalues as the one with the rows and columns switched
}

# Find eigenvalues
eig_2012_dK = eigen(rCmat_2012)
eig_2013_dK = eigen(rCmat_2013)
eig_2014_dK = eigen(rCmat_2014)
eig_2015_dK = eigen(rCmat_2015)
eig_allyears_dK = eigen(rCmat_allyears)
eig_allyears_dK_test = eigen(rCmat_allyears_test)

eig_allyears_dK_conn = eigen(Cmat_allyears_justconn)
eig_2012_dK_conn = eigen(Cmat_2012_dK_justconn)
eig_2013_dK_conn = eigen(Cmat_2013_dK_justconn)
eig_2014_dK_conn = eigen(Cmat_2014_dK_justconn)
eig_2015_dK_conn = eigen(Cmat_2015_dK_justconn)

##### Connectivity matrix, using migEst
# Set up matrices
rCmat_allyearsME <- matrix(NA, nrow = length(site_vec_order$site_name), ncol = length(site_vec_order$site_name))
rCmat_allyearsME_lowCI <- matrix(NA, nrow = length(site_vec_order$site_name), ncol = length(site_vec_order$site_name))
rCmat_allyearsME_highCI <- matrix(NA, nrow = length(site_vec_order$site_name), ncol = length(site_vec_order$site_name))
Cmat_allyearsME_justconn <- matrix(NA, nrow = length(site_vec_order$site_name), ncol = length(site_vec_order$site_name))

# Make realized connectivity from mig Est prop disp (and make any > 1 = 1)
migEst_pijmat <- migEst_pijmat %>%
  mutate(prop_disp_edit = if_else(prop_disp > 1, 1, prop_disp),
         prop_disp_LCI_edit = if_else(prop_disp_LCI > 1, 1, prop_disp_LCI),
         prop_disp_UCI_edit = if_else(prop_disp_UCI > 1, 1, prop_disp_UCI)) %>%
  mutate(r_prop_disp = prop_disp_edit*recruits_per_LEP,
         r_prop_disp_LCI = prop_disp_LCI_edit*recruits_per_LEP,
         r_prop_disp_UCI = prop_disp_UCI_edit*recruits_per_LEP)
  
# Fill in the matrices with the realized connetivity values
for(i in 1:length(migEst_pijmat$org_site)) {
  column = migEst_pijmat$org_geo_order[i]  # column is origin 
  row = migEst_pijmat$dest_geo_order[i]  # row is destination
  
  rCmat_allyearsME[row, column] = migEst_pijmat$r_prop_disp[i]
  rCmat_allyearsME_lowCI[row, column] = migEst_pijmat$r_prop_disp_LCI[i]
  rCmat_allyearsME_highCI[row, column] = migEst_pijmat$r_prop_disp_UCI[i]
  Cmat_allyearsME_justconn = migEst_pijmat$prop_disp_edit[i]
}

# Find eigenvalues
eig_allyears_ME = eigen(rCmat_allyearsME)
eig_allyears_ME_LCI = eigen(rCmat_allyearsME_lowCI)
eig_allyears_ME_UCI = eigen(rCmat_allyearsME_highCI)
eig_allyears_ME_conn = eigen(Cmat_allyearsME_justconn)

##### Pull out first eigenvalues so can plot
NP_eigs <- data.frame(eigen_val = c(eig_2012_dK$values[1], 
                                    eig_2013_dK$values[1],
                                    eig_2014_dK$values[1],
                                    eig_2015_dK$values[1],
                                    eig_allyears_dK$values[1],
                                    eig_allyears_dK_conn$values[1],
                                    eig_2012_dK_conn$values[1],
                                    eig_2013_dK_conn$values[1],
                                    eig_2014_dK_conn$values[1],
                                    eig_2015_dK_conn$values[1],
                                    Re(eig_allyears_ME$values[1]),
                                    Re(eig_allyears_ME_LCI$values[1]),
                                    Re(eig_allyears_ME_UCI$values[1]),
                                    Re(eig_allyears_ME_conn$values[1])),
                      origin = c('2012 kernel', 
                                 '2013 kernel',
                                 '2014 kernel',
                                 '2015 kernel',
                                 'all years kernel',
                                 'all years kernel conn',
                                 '2012 kernel conn',
                                 '2013 kernel conn',
                                 '2014 kernel conn',
                                 '2015 kernel conn',
                                 'all years MigEst',
                                 'all years MigEst LCI',
                                 'all years MigEst UCI', 
                                 'all years MigEst just conn'),
                      type = c('realized','realized','realized','realized','realized',
                               'connectivity','connectivity','connectivity','connectivity','connectivity',
                               'realized','realized','realized','connectivity'))
                      #type = c('kernel', 'kernel', 'kernel', 'kernel', 'kernel', 'kernel', 'kernel', 
                      #         'MigEst', 'MigEst', 'MigEst', 'MigEst'))

##### Calculating an example dispersal kernel to plot
distance_dp <- data.frame(distance = seq(0,50,0.5))
distance_dp$prob_2014 <- dispKernel(k_2014, theta_2014, distance_dp$distance)
distance_dp$prob_2013 <- dispKernel(k_2013, 1, distance_dp$distance)
distance_dp$prob_2015 <- dispKernel(k_2015, theta_2015, distance_dp$distance)


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
pdf(file = here("Plots/PersistenceMetrics", "NYearsRecaught_bysite.pdf"), width=11,height=5)
ggplot(data = tag_years_site, aes(nyears)) +
  geom_histogram(binwidth=1) +
  facet_grid(.~ site, labeller = label_wrap_gen(8)) + 
  xlab("number of years recaught") + ggtitle("# years tagged fish recaught by site") +
  theme_bw()
dev.off()

##### Realized connectivity matrix -- all years dispersal kernel
pdf(file=here("Plots/PersistenceMetrics", "Realized_connectivity_allyears_kernel.pdf"),width=7,height=5,useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=site_dist_info, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=realizedC_allyears)) +
  scale_fill_gradient(high='black',low="white") +
  #scale_fill_gradient(high='lightgray',low='black') +
  labs(x='origin',y='destination') +
  theme_bw() +
  #theme(text = element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### Connectivity matrix, all years dispersal kernel
pdf(file=here("Plots/PersistenceMetrics", "Connectivity_allyears_kernel.pdf"),width=7,height=5,useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=site_dist_info, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_allyears)) +
  scale_fill_gradient(high='black',low="white") +
  #scale_fill_gradient(high='lightgray',low='black') +
  labs(x='origin',y='destination') + ggtitle('All years dispersal kernel connectivity') +
  theme_bw() +
  #theme(text = element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()


##### Proportion of recruits from i arriving at j - estimates converted from Mig Est output - some prop_disps are bigger than one... think through whether that's a quirk (and can just cap them at 1) or a structural error...
pdf(file = here::here('Plots/PersistenceMetrics', 'MigEst_propdisp_mat.pdf'),width=7,height=5,useDingbats=F)
ggplot(data = migEst_pijmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order), fill=prop_disp)) +
  geom_tile() +
  scale_fill_gradient(high='black', low='white') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab("origin") + ylab("destination") + ggtitle("Prop. of recruits from origin settling at dest, converted from migEst")
dev.off()

# Lower confidence interval proportion of recruits from i arriving at j - estimates converted from Mig Est output LCI
pdf(file = here::here('Plots/PersistenceMetrics', 'MigEst_propdisp_mat_LCI.pdf'))
ggplot(data = migEst_pijmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order), fill=prop_disp_LCI)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab("origin site") + ylab("destination site") + ggtitle("LCI prop recruits from origin settling at dest, from migEst LCI") 
dev.off()

# Upper confidence interval proportion of recruits from i arriving at j - estimates converted from Mig Est output UCI
pdf(file = here::here('Plots/PersistenceMetrics', 'MigEst_propdisp_mat_UCI.pdf'))
ggplot(data = migEst_pijmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order), fill=prop_disp_UCI)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab("origin site") + ylab("destination site") + ggtitle("UCI prop recruits from origin settling at dest, from migEst UCI") 
dev.off()

# Bar plot of migEst-converted pijs with upper and lower confidence intervals
pdf(file = here::here('Plots/PersistenceMetrics', 'MigEst_propdist_by_site_withCI.pdf'))
ggplot(data = (migEst_pijmat %>% filter(org_site == dest_site)), aes(x=reorder(org_site, org_geo_order), y=prop_disp, ymin=prop_disp_LCI, ymax=prop_disp_UCI)) +
  geom_errorbar() +
  geom_point() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab("site") + ylab("prop disp with CI") + ggtitle("Prop disp with CI converted from MigEst") 
dev.off()

# Bar plot of migEst mig rates with upper and lower confidence intervals
pdf(file = here::here('Plots/PersistenceMetrics', 'MigEst_migrate_by_site.pdf'))
ggplot(data = (migEst_pijmat %>% filter(org_site == dest_site)), aes(x=reorder(org_site, org_geo_order), y=mig_rate, ymin=mig_rate_LCI, ymax=mig_rate_UCI)) +
  geom_errorbar() +
  geom_point() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab("site") + ylab("mig rate with CI") + ggtitle("Mig rate with CI from MigEst") 
dev.off()

##### SP metric - all methods
pdf(file = here('Plots/PersistenceMetrics', 'SP_allmethods.pdf'))
ggplot(data = SP_all, aes(x = reorder(site, geo_order), y=SP, color = est_method)) +
  geom_bar(aes(fill = data), position = 'dodge', stat='identity') +
  #geom_col() +
  geom_hline(yintercept = 1) +
  xlab("site") + ylab("self-persistence") + ggtitle("SP by site and method") +
  theme_bw() +
  #theme(text =  element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### SP metric - all methods - zoomed
pdf(file = here('Plots/PersistenceMetrics', 'SP_allmethods_zoomed.pdf'))
ggplot(data = SP_all, aes(x = reorder(site, geo_order), y=SP, color = est_method)) +
  geom_bar(aes(fill = data), position = 'dodge', stat='identity') +
  #geom_col() +
  #geom_hline(yintercept = 1) +
  xlab("site") + ylab("self-persistence") + ggtitle("SP by site and method") +
  theme_bw() +
  #theme(text =  element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### SP metric - all methods - this one looks weird...
pdf(file = here('Plots/PersistenceMetrics', 'SP_allmethods_site.pdf'))
ggplot(data = SP_all, aes(x = reorder(site, geo_order), y=SP, color = est_method)) +
  geom_bar(aes(fill = data), position = 'dodge', stat='identity') +
  #geom_col() +
  geom_hline(yintercept = 1) +
  facet_wrap(~site) +
  xlab("site") + ylab("self-persistence") + ggtitle("SP by site and method") +
  theme_bw() +
  #theme(text =  element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### SP metric - all methods
pdf(file = here('Plots/PersistenceMetrics', 'SP_allmethods_method.pdf'))
ggplot(data = SP_all, aes(x = reorder(site, geo_order), y=SP, color = est_method)) +
  geom_bar(aes(fill = data), position = 'dodge', stat='identity') +
  #geom_col() +
  geom_hline(yintercept = 1) +
  facet_wrap(~ est_method) +
  xlab("site") + ylab("self-persistence") + ggtitle("SP by site and method") +
  theme_bw() +
  #theme(text =  element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### SP metric - all methods
pdf(file = here('Plots/PersistenceMetrics', 'SP_allmethods_method_zoomed.pdf'))
ggplot(data = SP_all, aes(x = reorder(site, geo_order), y=SP, fill = est_method)) +
  geom_bar(position = 'dodge', stat='identity') +
  #geom_col() +
  #geom_hline(yintercept = 1) +
  facet_wrap(~ est_method) +
  xlab("site") + ylab("self-persistence") + ggtitle("SP by site and method") +
  theme_bw() +
  #theme(text =  element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()


##### SP metric - migEst  - these all look way high... something is wrong... and now look super low once I changed the input values from the egg-recruit relationship...
pdf(file = here("Plots/PersistenceMetrics", "SP_by_site_migEst.pdf"))
ggplot(data = SP_migEst, aes(x=reorder(org_site, org_geo_order), y=SP)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 1) +
  xlab("site") + ylab("self-persistence") + ggtitle("SP by site, migEst") +
  theme_bw() +
  theme(text =  element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### NP metrics
# NP metric (first eigenvalue) from dispersal kernels and MigEst reverse estimates - realized connectivity
pdf(file = here::here('Plots/PersistenceMetrics', 'NP_realizedconnectivity.pdf'))
ggplot(data = NP_eigs %>% filter(type == 'realized'), aes(x= origin, y = eigen_val, fill = origin)) +
  geom_bar(stat='identity') +
  xlab('estimate type') + ylab('first eigenvalue') + ggtitle('Eigenvalue of realized connectivity matrix') +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()
  
# NP metric (first eigenvalue) from dispersal kernels and MigEst reverse estimates - just connectivity
pdf(file = here::here('Plots/PersistenceMetrics', 'Eig_connmats.pdf'))
ggplot(data = NP_eigs %>% filter(type == 'connectivity'), aes(x= origin, y = eigen_val, fill = origin)) +
  geom_bar(stat='identity') +
  xlab('estimate origin') + ylab('first eigenvalue') + ggtitle('Eigenvalue of connectivity matrix') +
  #geom_hline(yintercept = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()


# ##### SP metric (FAKE DATA RIGHT NOW)
# pdf(file = here("Plots/PersistenceMetrics", "SP_by_site.pdf"))
# ggplot(data = SP_site, aes(x=site, y=SP_avg)) +
#   geom_bar(stat="identity") +
#   geom_hline(yintercept = 1) +
#   xlab("site") + ylab("self-persistence") + ggtitle("FAKE DATA SP by site plot") +
#   theme_bw() +
#   theme(text =  element_text(size=20)) +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# dev.off()

# ##### SP metric (prelim data for ESA!)
# pdf(file = here("Plots/PersistenceMetrics", "SP_by_site_ESAdraft.pdf"))
# ggplot(data = SP_site, aes(x=site, y=SP_avg)) +
#   geom_bar(stat="identity") +
#   geom_hline(yintercept = 1) +
#   xlab("site") + ylab("self-persistence") +
#   theme_bw() +
#   theme(text =  element_text(size=20)) +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# dev.off()

##### Prob of self-dispersing (from dispersal kernel) (prelim data for ESA!)
pdf(file = here("Plots/PersistenceMetrics", "Disp_home_by_site_ESAdraft.pdf"))
ggplot(data = SP_site, aes(x=site, y=prob_disperse)) +
  geom_bar(stat="identity") +
  xlab("site") + ylab("prob of dispersing home") +
  theme_bw() +
  theme(text =  element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### Example dispersal kernel (for MPE poster)
pdf(file = here("Plots/PersistenceMetrics", "Ex_disp_kernel.pdf"))
ggplot(data = distance_dp, aes(x=distance, y=prob)) +
  geom_line(size = 5) +
  xlab("distance") + ylab("probability of recruiting") +
  theme_bw() +
  theme(text = element_text(size=40)) +
  theme(axis.text.x = element_text(size=30)) + theme(axis.text.y = element_text(size=30))
dev.off()

##### Example dispersal kernel (for ESA talk) 
pdf(file = here("Plots/PersistenceMetrics", "Ex_disp_kernel_ESA_1.pdf"))
ggplot(data = distance_dp, aes(x=distance, y=prob_2013)) +
  geom_line(size = 5, color = "black") +
  xlab("distance") + ylab("probability of recruiting") +
  theme_bw() +
  theme(text = element_text(size=40)) +
  theme(axis.text.x = element_text(size=30)) + theme(axis.text.y = element_text(size=30))
dev.off()

##### Example dispersal kernel (for ESA talk) - multiple kernels
pdf(file = here("Plots/PersistenceMetrics", "Ex_disp_kernel_ESA_3.pdf"))
ggplot(data = distance_dp, aes(x=distance)) +
  geom_line(aes(x=distance, y=prob_2013), size = 3, color = "black") +
  geom_line(aes(x=distance, y=prob_2014), size = 3, color = "blue") +
  geom_line(aes(x=distance, y=prob_2015), size = 3, color = "red") +
  xlab("distance") + ylab("probability of recruiting") +
  theme_bw() +
  theme(text = element_text(size=40)) +
  theme(axis.text.x = element_text(size=30)) + theme(axis.text.y = element_text(size=30))
dev.off()

# Realized connectivity matrix
site_dist_info_plot <- site_dist_info %>%
  dplyr::rename(p2013 = prob_disperse_2013, p2014 = prob_disperse_2014_oldway, p2015 = prob_disperse_2015)

# 2013 kernel, analytical integration
pdf(file=here("Plots/PersistenceMetrics", "Realized_connectivity_2013kernel.pdf"),useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=site_dist_info_plot, aes(x=org_site,y=dest_site)) +
  geom_tile(aes(fill=realizedC_2013)) +
  scale_fill_gradient(high='black',low="white") +
  #scale_fill_gradient(high='lightgray',low='black') +
  labs(x='origin',y='destination') +
  theme_bw() +
  theme(text = element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

# Realized connectivity matrix
# 2013 kernel, non-analytical integration
pdf(file=here("Plots/PersistenceMetrics", "Realized_connectivity_2013kernel_nonan.pdf"),width=7,height=5,useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=site_dist_info, aes(x=org_site,y=dest_site)) +
  geom_tile(aes(fill=realizedC_2013_oldway)) +
  scale_fill_gradient(high='black',low="white") +
  #scale_fill_gradient(high='lightgray',low='black') +
  labs(x='origin',y='destination') +
  theme_bw() +
  theme(text = element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

# 2014 kernel, non-analytical integration
pdf(file=here("Plots/PersistenceMetrics", "Realized_connectivity_2014kernel_nonan.pdf"),width=7,height=5,useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=site_dist_info, aes(x=org_site,y=dest_site)) +
  geom_tile(aes(fill=realizedC_2014)) +
  scale_fill_gradient(high='black',low="white") +
  #scale_fill_gradient(high='lightgray',low='black') +
  labs(x='origin',y='destination') +
  theme_bw() +
  theme(text = element_text(size=15)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(legend.title=element_text(size=0))
dev.off()

pdf(file=here("Plots/PersistenceMetrics", "Realized_connectivity_2014kernel_nonan_orderedsites.pdf"),width=7,height=5,useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=site_dist_info, aes(x=origin,y=destination)) +
  geom_tile(aes(fill=realizedC_2014)) +
  scale_fill_gradient(high='black',low="white") +
  #scale_fill_gradient(high='lightgray',low='black') +
  labs(x='origin',y='destination') +
  theme_bw() +
  theme(text = element_text(size=15)) +
  #theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(legend.title=element_text(size=0))
dev.off()

# 2015 kernel, analytical integration
pdf(file=here("Plots/PersistenceMetrics", "Realized_connectivity_2015kernel.pdf"),width=7,height=5,useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=site_dist_info, aes(x=org_site,y=dest_site)) +
  geom_tile(aes(fill=realizedC_2015)) +
  scale_fill_gradient(high='black',low="white") +
  #scale_fill_gradient(high='lightgray',low='black') +
  labs(x='origin',y='destination') +
  theme_bw() +
  theme(text = element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

# Just prob dispersing
# 2014 kernel, non-analytical integration
pdf(file=here("Plots/PersistenceMetrics", "Connectivity_2014kernel_nonan.pdf"),useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=site_dist_info_plot, aes(x=org_site,y=dest_site)) +
  geom_tile(aes(fill=p2014)) +
  scale_fill_gradient(high='black',low="white") +
  #scale_fill_gradient(high='lightgray',low='black') +
  labs(x='origin',y='destination') +
  theme_bw() +
  theme(text = element_text(size=15)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

# 2013 kernel, analytical integration
pdf(file=here("Plots/PersistenceMetrics", "Connectivity_2013kernel_2.pdf"),useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=site_dist_info_plot, aes(x=org_site,y=dest_site)) +
  geom_tile(aes(fill=p2013)) +
  scale_fill_gradient(high='black',low="white") +
  #scale_fill_gradient(high='lightgray',low='black') +
  labs(x='origin',y='destination') +
  theme_bw() +
  theme(text = element_text(size=15)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

# 2015 kernel, analytical integration
pdf(file=here("Plots/PersistenceMetrics", "Connectivity_2015kernel_2.pdf"),width=7,height=5,useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=site_dist_info, aes(x=org_site,y=dest_site)) +
  geom_tile(aes(fill=prob_disperse_2015)) +
  scale_fill_gradient(high='black',low="white") +
  #scale_fill_gradient(high='lightgray',low='black') +
  labs(x='origin',y='destination') +
  theme_bw() +
  theme(text = element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

#################### Saving things: ####################
save(allfish, file=here("Data", "allfish.RData"))
save(allfish_fish, file=here("Data", "allfish_fish.RData"))
save(allfish_anems, file=here("Data", "allfish_anems.RData"))
save(allfish_dives, file=here("Data", "allfish_dives.RData"))
save(gps.Info, file=here("Data", "gps.Info.RData"))
save(site_dist_info, file=here("Data", "site_dist_info.RData"))

# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}
# ----- END plot function ----- #

#################### Old code: ####################


# #### STOPPED EDITING HERE FOR NP STUFF!
# site_dist_info <- site_dist_info %>%
#   mutate(realizedC_2014 = prob_disperse_2014_oldway*recruits_per_LEP) %>%
#   mutate(realizedC_2013 = prob_disperse_2013*recruits_per_LEP) %>%
#   mutate(realizedC_2015 = prob_disperse_2015*recruits_per_LEP) %>%
#   mutate(realizedC_2013_oldway = prob_disperse_2013_oldway*recruits_per_LEP) 

# # Fill in with prob_disp
# for(i in 1:length(site_vec_order$site_name)) {
#   row_val <- site_vec_order$geo_order[i]
#   for(j in 1:length(site_vec_order$site_name)) {
#     col_val <- site_vec_order$geo_order[j]
#     
#     conn_mat
#   }
#   
# }



# ##### Connectivity matrix, using dispersal kernels
# 
# # # Find prob of dispersing within the width of the site
# # for(i in 1:length(site_width_info$site)){ #this doesn't work
# #   site_width_info$prob_disperse_2013[i] = integrateDK_theta1(0, site_width_info$width_km[i], k_2013) #2013 kernel using theta=1
# #   site_width_info$prob_disperse_2015[i] = integrateDK_theta1(0, site_width_info$width_km[i], k_2015) #2015 kernel
# #   site_width_info$prob_disperse_2014_oldway[i] = integrate(dispKernel2014, 0, site_width_info$width_km[i])$value #this is almost certainly a terrible way to integrate - make this better! Could probably even do it analytically!
# #   site_width_info$prob_disperse_2013_oldway[i] = integrate(dispKernel2013, 0, site_width_info$width_km[i])$value
# # }
# 
# test_1 <- seq(0,50,by=0.5)
# test_2 <- dispKernelallyears(test_1)
# plot(test_1, test_2)
# 
# # # Find prob of dispersing within the width of the site - just using integrate function here, should swap out for analytical integration...
# # for(i in 1:length(site_width_info$site)){ #this doesn't work
# #   site_width_info$prob_disperse_2013[i] = integrate(dispKernel2012, 0, site_width_info$width_km[i])$value 
# #   site_width_info$prob_disperse_2015[i] = integrate(dispKernel2013, 0, site_width_info$width_km[i])$value 
# #   site_width_info$prob_disperse_2014[i] = integrate(dispKernel2014, 0, site_width_info$width_km[i])$value #this is almost certainly a terrible way to integrate - make this better! Could probably even do it analytically!
# #   site_width_info$prob_disperse_2013[i] = integrate(dispKernel2015, 0, site_width_info$width_km[i])$value
# # }
# # 
# # # Find prob of dispersing among sites
# # for(i in 1:length(site_dist_info$org_site)){ 
# #   site_dist_info$prob_disperse_2013[i] = integrateDK_theta1(site_dist_info$d1_km[i], site_dist_info$d2_km[i], k_2013) #2013 kernel using theta=1
# #   site_dist_info$prob_disperse_2015[i] = integrateDK_theta1(site_dist_info$d1_km[i], site_dist_info$d2_km[i], k_2015) #2015 kernel
# #   site_dist_info$prob_disperse_2014_oldway[i] = integrate(dispKernel2014, site_dist_info$d1_km[i], site_dist_info$d2_km[i])$value #this is almost certainly a terrible way to integrate - make this better! Could probably even do it analytically!
# #   site_dist_info$prob_disperse_2013_oldway[i] = integrate(dispKernel2013, site_dist_info$d1_km[i], site_dist_info$d2_km[i])$value
# # }
# 
# # Add in index values
# site_dist_info$origin <- c(rep(4, length(site_vec)), rep(5, length(site_vec)), rep(6, length(site_vec)), #Cabatoan, Caridad Cemetery, Caridad Proper
#                            rep(9, length(site_vec)), rep(15, length(site_vec)), rep(17, length(site_vec)), #Elementary School, Gabas, Haina
#                            rep(7, length(site_vec)), rep(3, length(site_vec)), rep(1, length(site_vec)), #Hicgop South, Magbangon, Palanas
#                            rep(13, length(site_vec)), rep(12, length(site_vec)), rep(11, length(site_vec)), #Poroc Rose, Poroc San Flower, San Agustin
#                            rep(18, length(site_vec)), rep(10, length(site_vec)), rep(8, length(site_vec)), #Sitio Baybayon, Sitio Lonas, Sitio Tugas
#                            rep(16, length(site_vec)), rep(14, length(site_vec)), rep(2, length(site_vec))) #Tamakin Dacot, Visca, Wangag
# 
# site_dist_info$destination <- rep(c(4,5,6,9,15,17,7,3,1,13,12,11,18,10,8,16,14,2), length(site_vec))
# # list of sites N-S: 
# # 1) Palanas, 2) Wangag, 3) Magbangon, 4) Cabatoan, 5) Caridad Cemetery, 6) Caridad Proper
# # 7) Hicgop South, 8) Sitio Tugas, 9) Elementary School 10) Sitio Lonas
# # 11) San Agustin, 12) Poroc San Flower, 13) Poroc Rose, 14) Visca, 15) Gabas  1288, 2273 2271 anems
# # 16) Tamakin Dacot, 17) Haina, 18) Sitio Baybayon
# 
# # Turn site_dist_info into a matrix
# disp_mat_2014 <- matrix(NA, 18, 18)
# realized_C_2014 <- matrix(NA, 18, 18)
# realized_C_2013 <- matrix(NA, 18, 18)
# realized_C_2015 <- matrix(NA, 18, 18)
# for(i in 1:length(site_dist_info$org_site)) {
#   column = site_dist_info$origin[i]
#   row = site_dist_info$destination[i]
#   disp_mat_2014[row, column] = site_dist_info$prob_disperse_2014_oldway[i]
#   realized_C_2014[row, column] = site_dist_info$realizedC_2014[i]
#   realized_C_2015[row, column] = site_dist_info$realizedC_2015[i]
#   realized_C_2013[row, column] = site_dist_info$realizedC_2013[i]
# }
# 
# # Find eigenvalues
# eig_2013 <- eigen(realized_C_2013)
# eig_2014 <- eigen(realized_C_2014)
# eig_2015 <- eigen(realized_C_2015)
# 
# 
# 
# # Try to melt this into a matrix
# #disp_mat_2014 <- melt(site_dist_info %>% select(org_site, dest_site, prob_disperse_2014_oldway), id.vars = c("org_site", "dest_site"))
# 
# #SHOULD TRY ASSIGNING EACH A ROW AND COLUMN VALUE (CAN BE N-S), THEN MAKE THOSE THE INDICES
# 
# 
# # # Old way (for practice presentation)
# # for(i in 1:length(site_width_info$site)){ #this doesn't work
# #   site_width_info$prob_disperse[i] = integrate(dispKernel2014, 0, site_width_info$width_km[i])$value #this is almost certainly a terrible way to integrate - make this better! Could probably even do it analytically!
# #   site_width_info$prob_disperse_abserror[i] = integrate(dispKernel2014, 0, site_width_info$width_km[i])$abs.error
# # }
# # site_width_info$prob_disperse <- unlist(site_width_info$prob_disperse) #makes prob_disperse a list for some reason, make it not anymore
# 
# #save(site_width_info, file=here("Data", "site_width_info.RData"))
# 
# ## Calculating the metric!
# SP_site <- site_width_info %>% select(site, width_km, prob_disperse, prob_disperse_abserror) %>%
#   mutate(SP_avg = prob_disperse*recruits_per_LEP)
# 
# # ## Try network persistence! - uh-oh, migEst gives self-recruitment...
# # realizedC <- migEst_conn[,2:19]*recruits_per_LEP #multiply connectivity matrix by recruits_from_LEP
# # realizedCmat <- as.matrix(realizedC)
# # realizedC_mat_eig <- eigen(realizedC)
# 
# ## Calculate realized connectivity matrix in plotting way (dataframe)
# site_dist_info <- site_dist_info %>%
#   mutate(realizedC_2014 = prob_disperse_2014_oldway*recruits_per_LEP) %>%
#   mutate(realizedC_2013 = prob_disperse_2013*recruits_per_LEP) %>%
#   mutate(realizedC_2015 = prob_disperse_2015*recruits_per_LEP) %>%
#   mutate(realizedC_2013_oldway = prob_disperse_2013_oldway*recruits_per_LEP) 
#   
# 
# 
# plot(realizedCmat)
# myImagePlot(realizedCmat, xlabels, ylabels, zlim, title=c("my title")) 

# ## FAKE DATA FOR NOW, JUST TO SEE
# # Average number of recruits per site (ALL FAKE DATA!!!!!)
# recruits_to_site <- data.frame(site = site_vec)
# recruits_to_site$avg_recruit_pop <- c(45,15,15,20,20,55,45,200,350,10,10,15,450,5,5,75,30,400)
# recruits_to_site$avg_recruits_from_site <- recruits_to_site$avg_recruit_pop*10
# 
# # Average number of recruits per site (trying out real data)
# recruits_info_input <- recruits_info %>% 
#   group_by(site) %>%
#   summarise(avg_annual_recruits_est = mean(totalR_est), avg_annual_recruits_raw = mean(Nrecruits))
#   
# recruits_to_site <- data.frame(site = site_vec)
# recruits_to_site$avg_recruit_pop <- c(45,15,15,20,20,55,45,200,350,10,10,15,450,5,5,75,30,400)
# recruits_to_site$avg_recruits_from_site <- recruits_to_site$avg_recruit_pop*10
# 
# # Average egg output a year per site (ALL FAKE DATA!!!!!!)
# egg_output <- data.frame(site = site_vec)
# egg_output$breedingF <- c(30,10,10,5,10,40,30,100,200,5,5,5,300,0,0,50,15,200)
# egg_output$eggs_per_year <- egg_output$breedingF*eggs_per_clutch*clutch_per_year
# 
# ## Finding number of recruits returning home for each site (numerator of LR fraction) 
# # Find prob of dispersing within the width of the site
# for(i in 1:length(site_width_info$site)){ #this doesn't work
#   site_width_info$prob_disperse[i] = integrate(dispKernel2014, 0, site_width_info$width_km[i])$value #this is almost certainly a terrible way to integrate - make this better! Could probably even do it analytically!
#   site_width_info$prob_disperse_abserror[i] = integrate(dispKernel2014, 0, site_width_info$width_km[i])$abs.error
# }
# site_width_info$prob_disperse <- unlist(site_width_info$prob_disperse) #makes prob_disperse a list for some reason, make it not anymore
# 
# save(site_width_info, file=here("Data", "site_width_info.RData"))
# 
# ## Calculating the metric!
# SP_site <- data.frame(site = site_vec)
# SP_site <- left_join(SP_site, egg_output, by="site")
# SP_site <- left_join(SP_site, recruits_to_site, by="site")
# SP_site <- left_join(SP_site, site_width_info %>% select(site, prob_disperse, prob_disperse_abserror), by="site")
# SP_site$home_recruits <- SP_site$avg_recruits_from_site*SP_site$prob_disperse
# SP_site$home_recuits_high <- SP_site$avg_recruits_from_site*(SP_site$prob_disperse + SP_site$prob_disperse_abserror)
# SP_site$home_recuits_low <- SP_site$avg_recruits_from_site*(SP_site$prob_disperse - SP_site$prob_disperse_abserror)
# SP_site$SP_avg <- findSP(rep(LEP_Will, length(SP_site$site)), SP_site$home_recruits, SP_site$eggs_per_year)
# SP_site$SP_avgU <- findSP(rep(LEP_Will, length(SP_site$site)), SP_site$home_recuits_high, SP_site$eggs_per_year)
# SP_site$SP_avgL <- findSP(rep(LEP_Will, length(SP_site$site)), SP_site$home_recuits_low, SP_site$eggs_per_year)
# SP_site$site <- unfactor(SP_site$site)
# SP_site$SP_avg[14:15] <- c(0,0) #make sites w/0 or NaN SP have 0
