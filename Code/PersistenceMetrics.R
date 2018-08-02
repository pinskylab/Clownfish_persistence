# Script to calculate persistence metrics
# For now, using fake data + parameters where necessary to get structure right - will note clearly, format might shift as real stuff comes in

#################### Set-up: ####################
#Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(geosphere)
library(varhandle)
library(ggplot2)
library(here)

# Load some data files while deciding how various scripts should interact/communicate
# Load output from AnemLocations.R (or source that file) - mean gps coordinates across multiple observations of anemones
#load(file=here("Data",'AnemAllInfowLatLon2.RData')) #file with anems, after 2018 anems matched
#load(file=here("Data", "AnemLatLonObsbyAnem.RData")) #file with lat-lon info appended

# Load encounter history data frame from ClownfishMarkRecap.R
#load(file=here("Data", "encounters_all.RData"))

# Load mig_est data
load(file=here("Data", "MigEstPooled_out_justimmigration.txt"))

migEst_conn <- read.table(file=here("Data", "MigEstPooled_out_justimmigration.txt"), 
           header = TRUE, sep = "")

migEst_sites <- c("Cabatoan", "Caridad Cemetery", "Caridad Proper", "Elementary School", "Gabas",
                  "Haina", "Hicgop South", "Magbangon", "Palanas", "Poroc Rose", "Poroc San Flower",
                  "San Agustin", "Sitio Baybayon", "Sitio Lonas", "Sitio Tugas", "Tamakin Dacot",
                  "Visca", "Wangag")
migEst_site <- data.frame(site = migEst_sites, pop = seq(from=1,to=18,by=1), stringsAsFactors = FALSE)

# Load recruit estimates
load(file=here("Data", "recruits_info.RData"))
load(file=here("Data", "females_recruits_summary.RData"))
load(file=here("Data", "metapop_level.RData"))

# Load egg-recruit linear model estimates
load(file=here("Data", "metapop_lm.RData")) #estimate at metapop level, pretty low R2 and not significant
load(file=here("Data", "egg_recruit_lm.RData")) #estimate using each site and year as a point, R2 of about 0.53, significant

# Load LEP estimates
load(file=here("Data","LEP_out.RData")) #LEP estimates from LEP_estimate.R

# Set a few things
sample_years = c(2015,2016,2017,2018)

# PRELIM OR FAKE FOR NOW THAT WILL GET REPLACED BY REAL DATA
eggs_per_clutch = 1763 #from LEP_calc_WillWhite.R
clutch_per_year = 11.9 #from LEP_calc_WillWhite.R
eggs_intercept = 100
eggs_slope = 107 #from Adam Y Aresty poster, #eggs/cm in female size

# prelim dispersal kernel parameters from Katrina (from email sent 7/5/18)
k_2012 = 0.1
theta_2012 = 5
k_2013 = 0.025
theta_2013 = 2 #or 1
k_2014 = 0.1
theta_2014 = 2
k_2015 = 0.4
theta_2015 = 1

# LEP
LEP_Will <- 1780 #eggs/recruit
LEP_Will_eggs_per_clutch <- LEP_out$LEP_W #my size-survival relationship, Will's eggs per clutch
LEP <- LEP_out$LEP #AY's prelim egg data, my size-survival relationship

# anems at N and S ends of sites (visually from QGIS, in particular mid anems totally eyeballed)
Cabatoan_N <- 198 # 3195, 951 other options, 2502
Cabatoan_mid <- 2241
Cabatoan_S <- 904 #2250, 185 other options
CaridadCemetery_N <- 695
CaridadCemetery_mid <- 2683 #2256
CaridadCemetery_S <- 292
CaridadCemetery_mid <- 289
CaridadProper_N <- 291
CaridadProper_S <- 290
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
Magbangon_N_S <- 680 #680 is S of hull, 1114 is another option; 212 is S end of northern-most chunk, 1391 is another option
Magbagon_S_N <- 1113 #209 is another option
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
SitioBaybaon_mid <- 538
SitioBaybayon_S <- 2747 #outside hull (maybe KC sites?), #805, 2148 w/in hull
SitioLonas_N <- 48
SitioLonas_S <- 48
SitioTugas_N <- 54
SitioTugas_S <- 59
TamakinDacot_N <- 2554
TamakinDacot_mid <- 2861
TamakinDacot_S <- 2270 #2147
Visca_N <- 4
Visca_mid <- 8 #776
Visca_S <- 2314
Wangag_N <- 2985
Wangag_mid <- 2734
Wangag_S <- 2063 #1034 also a good end point


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
findLEP <- function(l_vec, f_vec) { #l_vec: vector of survivals-to-age, f_vec: vector of fecundities-at-age. Need to make sure first age of reproduction and max are are included in bounds, otherwise can define starting age as would like.
  out <- sum(l_vec*f_vec)
  return(out)
}

# Calculate LEP via an IPM (like Will did in LEP_calc_WillWhite.R)

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

# Calculate dispersal kernel (equation from Katrina email 7/5/18, based on Bode et al. 2017 paper)
dispKernel <- function(k, theta, d) {
  disp <- k*exp(-(k*d)^theta)
  return(disp)
}

# Dispersal kernel 2014
dispKernel2014 <- function(d) {
  k_2014*exp(-(k_2014*d)^theta_2014)
}

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

#pull GPS info
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

#pull out just tagged fish
taggedfish <- allfish %>% filter(!is.na(tag_id))

########## Calculating values, constants, inputs from data
##### Find max age (A) - maximum number of times a fish has been caught
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
                            Magbangon_N, Palanas_N, PorocRose_N, PorocSanFlower_N, SanAgustin_N, SitioBaybayon_N, SitioLonas_N,
                            SitioTugas_N, TamakinDacot_N, Visca_N, Wangag_N)

# anem_id of southern-most anem at site (or western-most in case of Haina, which runs E-W)
site_width_info$S_anem <- c(Cabatoan_S, CaridadCemetery_S, CaridadProper_S, ElementarySchool_S, Gabas_S, Haina_W, HicgopSouth_S,
                            Magbangon_S, Palanas_S, PorocRose_S, PorocSanFlower_S, SanAgustin_S, SitioBaybayon_S, SitioLonas_S,
                            SitioTugas_S, TamakinDacot_S, Visca_S, Wangag_S)

# pull out anem_table_ids for those anems so can use anemid_latlong2 function to find lat lon for boundary anems
site_width_anems <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, obs_time) %>%
  collect() %>%
  filter(anem_id %in% c(site_width_info$N_anem, site_width_info$S_anem)) 

# pull out the dive info associated with those anems
site_width_dives <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, date, site, gps) %>%
  collect() %>%
  filter(dive_table_id %in% site_width_anems$dive_table_id)

# join into one data frame to put into MS function (anem_latlong), get the anem obs_time in the same time zone (UTC) as the GPX data, create placeholder lat and lon columns
site_width_anemdives <- left_join(site_width_anems, site_width_dives, by="dive_table_id") %>%
  mutate(obs_time = force_tz(ymd_hms(str_c(date, obs_time, sep = " ")), tzone = "Asia/Manila")) %>% #tell it that it is currently in Asia/Manila time zone
  mutate(obs_time = with_tz(obs_time, tzone = "UTC")) %>% #convert to UTC so can compare with GPS data (this line and one above largely from Michelle's assign_db_gpx function)
  mutate(month = month(obs_time), #and separate out useful components of the time (this also from Michelle's assign_db_gpx function)
         day = day(obs_time), 
         hour = hour(obs_time), 
         min = minute(obs_time), 
         sec = second(obs_time), 
         year = year(obs_time)) %>%
  mutate(lat = as.numeric(NA), lon = as.numeric(NA)) #put in a placeholder column for lat and lon

# find lat lon for the boundary anems
for(i in 1:length(site_width_anemdives$anem_table_id)) {
  out_anemLL <- anemid_latlong_2(site_width_anemdives$anem_table_id[i], site_width_anemdives, gps.Info)
  site_width_anemdives$lat[i] = out_anemLL$lat
  site_width_anemdives$lon[i] = out_anemLL$lon
}

# just pull out one row for each anem_id 
site_width_anemdives_short <- site_width_anemdives %>%
  distinct(anem_id, .keep_all = TRUE)

# put lat and lons into info table (requires some column renaming so needs to be run in order)
# lat/lon for N_anem
site_width_info <- left_join(site_width_info %>% rename(anem_id = N_anem), site_width_anemdives_short %>% select(anem_id, lat, lon), by="anem_id") #join in the lat lons for the N anem
site_width_info <- site_width_info %>% rename(N_anem = anem_id, N_lat = lat, N_lon = lon) #rename the anem_id, lat, lon columns to be N-specific
# lat/lon for S_anem
site_width_info <- left_join(site_width_info %>% rename(anem_id = S_anem), site_width_anemdives_short %>% select(anem_id, lat, lon), by="anem_id") #join in the lat lons for the S anem
site_width_info <- site_width_info %>% rename(S_anem = anem_id, S_lat = lat, S_lon = lon) #rename the anem_id, lat, lon columns to be S-specific 

#and calculate the distance!
for(i in 1:length(site_width_info$site)) {
  site_width_info$width_m[i] = distHaversine(c(site_width_info$N_lon[i], site_width_info$N_lat[i]), c(site_width_info$S_lon[i], site_width_info$S_lat[i]))
  site_width_info$width_km[i] = site_width_info$width_m[i]/1000 
}

########## Looking at eggs and recruits relationships


########## Assessing metrics

##### Self-persistence
# for findSP_v2 - inputs are LEP, R_E_slope, prob_disp_home
findRfromE <- function(m,eggs,b) {
  recruits = m*eggs + b
  return(recruits)
}

# Egg-recruit slope with LEP
intercept_mod1 <- 1.856e+01 #need to figure out how to extract this from egg_recruits_est_mod1
slope_mod1 <- 5.430e-05 #need to figure out how to extract this from egg_recruits_est_mod1
recruits_per_LEP <- findRfromE(slope_mod1, LEP, intercept_mod1) #18! Means persistence is possible!

# Find prob of dispersing within the width of the site
for(i in 1:length(site_width_info$site)){ #this doesn't work
  site_width_info$prob_disperse[i] = integrate(dispKernel2014, 0, site_width_info$width_km[i])$value #this is almost certainly a terrible way to integrate - make this better! Could probably even do it analytically!
  site_width_info$prob_disperse_abserror[i] = integrate(dispKernel2014, 0, site_width_info$width_km[i])$abs.error
}
site_width_info$prob_disperse <- unlist(site_width_info$prob_disperse) #makes prob_disperse a list for some reason, make it not anymore

#save(site_width_info, file=here("Data", "site_width_info.RData"))

## Calculating the metric!
SP_site <- site_width_info %>% select(site, width_km, prob_disperse, prob_disperse_abserror) %>%
  mutate(SP_avg = prob_disperse*recruits_per_LEP)

## Try network persistence!
realizedC <- migEst_conn[,2:19]*recruits_per_LEP #multiply connectivity matrix by recruits_from_LEP
realizedCmat <- as.matrix(realizedC)
realizedC_mat_eig <- eigen(realizedC)

pdf(file=here("Plots/PersistenceMetrics", "Realized_connectivity_pooled_migEst.pdf"),width=7,height=5,useDingbats=F) #looks like the lines don't intercept the y-axis together...
ggplot(data=realizedC, aes(x=Spread,y=Magnitude)) +
  geom_tile(aes(fill=Abundance)) +
  #scale_fill_gradient(high='lightgray',low="steelblue") +
  scale_fill_gradient(high='lightgray',low='black') +
  labs(x='Spread (patch lengths) ',y=expression(paste('Magnitude (',nu,')'))) +
  #labs(title='25% habitat',x='Spread (patch lengths) ',y=expression(paste('Magnitude (',nu,')'))) +
  theme_bw() +
  theme(text = element_text(size=20)) 
dev.off()

plot(realizedCmat)
myImagePlot(realizedCmat, xlabels, ylabels, zlim, title=c("my title")) 

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

##### Calculating an example dispersal kernel to plot
distance_dp <- data.frame(distance = seq(0,50,0.5))
distance_dp$prob <- dispKernel(k_2014, theta_2014, distance_dp$distance)



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

##### SP metric (FAKE DATA RIGHT NOW)
pdf(file = here("Plots/PersistenceMetrics", "SP_by_site.pdf"))
ggplot(data = SP_site, aes(x=site, y=SP_avg)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 1) +
  xlab("site") + ylab("self-persistence") + ggtitle("FAKE DATA SP by site plot") +
  theme_bw() +
  theme(text =  element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### SP metric (prelim data for ESA!)
pdf(file = here("Plots/PersistenceMetrics", "SP_by_site_ESAdraft.pdf"))
ggplot(data = SP_site, aes(x=site, y=SP_avg)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = 1) +
  xlab("site") + ylab("self-persistence") +
  theme_bw() +
  theme(text =  element_text(size=20)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

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
