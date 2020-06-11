# Find proportion of anemones occupied by APCL vs. "empty"

#################### Set-up: ####################
# Load packages
#library(stringr)
#library(ggplot2)

# Pull data, functions, constants 
source(here::here('Code', 'Constants_database_common_functions.R'))

# Could have an alternate section that runs the relevant constants/etc. and loads the saved data files... think about organization a bit more...

# Set vectors of months, years, and dive types
all_months <- c(1,2,3,4,5,6,7,8,9,10,11,12)  # for use in anemone occupancy function
#all_years <- c(2012, 2013, 2014, 2015, 2016, 2017, 2018)
anem_occ_dives <- c("A","C","D","E","F","M","R")  # use all dive types for now, except R (think through this more)

# This is in Constants_database_common_functions - should probably just use this? But what are 0 dives?
clown_sample_dive_types <- c("0","A","C","D","E","F","M")  # everything but R 

#################### Functions: ####################

# Find the number and proportion of anemones in each site and year occupied by APCL vs. occupied by other clownfish spp. vs. empty (modified from findCompOccupancyByYear in Metapop_exploration.R in Clownfish_metapop repository)
anemOcupancyByYear <- function(year_val, month_vals, divetype_vals, allanemsdf) {  # inputs: year_val - years to consider, months_val - months to consider
  
  #filter out just the relevant year and months
  allanemsdf_filt <- allanemsdf %>% filter(year == year_val) %>%  # filter by year
    filter(month %in% month_vals) %>%  # filter by months
    filter(dive_type %in% divetype_vals)  # filter by dive types
      
  #limit to one obs per anemone_table_id (think this through a bit more - ideally this gets one row per fish spp per anem; are there some anems showing up more than that? like if they were visited multiple times per season?) - actually, this seems a bit low... only 652 anem/spp combos?
  #allanemsdf_distinct <- allanemsdf_filt %>% distinct(anem_table_id, dive_table_id, anem_id, anem_obs, old_anem_id, anem_spp, obs_time, fish_spp.x, fish_spp.y, dive_type, date, site, gps, year, month)
  allanemsdf_distinct <- allanemsdf_filt %>% distinct(anem_table_id, fish_spp_clownfish, fish_spp_sighted, .keep_all = TRUE)
  
  # #filter for anem_table_ids that have multiple rows (multiple spp on that anem)
  # allanemsdf_mult <- as.data.frame(table(allanemsdf_distinct$anem_table_id)) #make a table of the anem_table_id frequencies
  # allanemsdf_mult2 <- allanemsdf_mult %>% filter(Freq >= 2) #filter out just those that appear at least twice
  # 
  # #check that it's just fish_spp that is the different and not anything else...
  # sppcheck <- allanemsdf_distinct %>% filter(anem_table_id %in% allanemsdf_mult2$Var1) #pull out just the entries with the anem_table_ids that appear at least twice
  # sppcheck_distinct <- sppcheck %>% distinct(anem_table_id, dive_table_id, anem_id, anem_obs, old_anem_id, anem_spp, obs_time, dive_type, date, site, gps, year, month) #pull out distinct rows for all columns except fish_spp (b/c we assume that is the one differing between the two rows for each anem_table_id)
  # n_nonfishsppdiffs <- length(allanemsdf_mult2$Var1) - length(sppcheck_distinct$anem_table_id) #should be 0 if fish_spp is the only thing changing
  # if (n_nonfishsppdiffs != 0) { #print a message if the number of obs not accounted for by fish_spp isn't 0 (so something other than fish_spp is different)
  #   print("Something other than fish_spp is causing multiple rows per anem_table_id.")
  # } else {
  #   print("Fish spp accounted for multiple anem_table_id observations.")
  # }
  # 
  # #check that anems aren't being included multiple times (multiple obs per season, say) 
  # anem_ids_mult <- as.data.frame(table(sppcheck_distinct$anem_id))
  # anem_ids_mult2 <- anem_ids_mult %>% filter(Freq >= 2) #filter out just the anem_ids observed more than once
  # n_multanemids <- length(anem_ids_mult2$Var1) #find out how many anems are observed more than once
  # print(n_multanemids) 
  # #FIGURE OUT HOW TO HANDLE MULTIPLE ANEM SIGHTINGS!
  # #test2 <- allfishdf_distinct %>% filter(anem_id == 553)
  # 
  # #find the total number of anems sampled that year and season (since some have multiple clownfish spp on them)
  # n_2visits <- allanemsdf_mult2 %>% filter(Freq == 2)
  # n_3visits <- allanemsdf_mult2 %>% filter(Freq == 3)
  # n_4visits <- allanemsdf_mult2 %>% filter(Freq == 4)
  # 
  # totalanems <- length(allanemsdf_distinct$anem_table_id) - length(n_2visits$Var1) - 2*length(n_3visits$Var1) - 3*length(n_4visits$Var1)
  totalanems <- length(allanemsdf_distinct$anem_table_id)
  
  #find frequency and percentage of each spp on anems 
  anems_fishspp_combined <- allanemsdf_distinct %>%
    mutate(fish_spp = case_when(is.na(fish_spp_clownfish) & is.na(fish_spp_sighted) ~ 'UNOC',  # no fish
                                is.na(fish_spp_clownfish) & !is.na(fish_spp_sighted) ~ fish_spp_sighted,  # no fish caught but fish seen
                                !is.na(fish_spp_clownfish) & is.na(fish_spp_sighted) ~ fish_spp_clownfish,  # no fish seen but fish caught
                                !is.na(fish_spp_clownfish) & !is.na(fish_spp_sighted) & fish_spp_clownfish == fish_spp_sighted ~ fish_spp_clownfish,  # fish caught and seen but both same species
                                !is.na(fish_spp_clownfish) & !is.na(fish_spp_sighted) & fish_spp_clownfish != fish_spp_sighted ~ paste(fish_spp_clownfish, fish_spp_sighted, sep="+")))  # fish caught and seen but different species
  
  anems_Cspp <- as.data.frame(table(anems_fishspp_combined$fish_spp, useNA = "ifany"), stringsAsFactors = FALSE)
  anems_Cspp$Perc <- anems_Cspp$Freq/totalanems
  anems_Cspp <- anems_Cspp %>%
    mutate(year = year_val)
  
  # anems_Cspp_by_site <- anems_fishspp_combined %>%
  #   group_by(site) %>%
  #   mutate(totalsiteanems = n()) %>%
  #   ungroup() %>%
  #   group_by(site, fish_spp) %>%
  #   summarize(nAnems = n(),
  #             percAnems = n()/totalsiteanems,
  #             year = year_val) 
  
  anems_Cspp_by_site <- anems_fishspp_combined %>%
    group_by(site, fish_spp) %>%
    summarize(nAnems = n(),
              year = year_val)
  
  out <- list(anems_Cspp = anems_Cspp, anems_Cspp_by_site = anems_Cspp_by_site)
  
  return(out)
  
}

#################### Running things: ####################

##### Data set up

#### Testing 4 combos
# 1) all_anems_with_fish - original one in script, uses anem_db and allfish_caught
# 2) all_anems_with_fish_2 - uses anem_db and fish_db
# 3) all_anems_with_fish_3 - uses anems_Processed_all and allfish_caught
# 4) all_anems_with_fish_4 - uses anems_Processed_all and fish_db

# Join together all anems with fish_caught -- anems_db is all anems and allfish_caught is all APCL in "clownfish table", both created in Constants_database_common_functions
all_anems_with_fish <- left_join(anem_db, allfish_caught %>% select(anem_table_id, fish_table_id, fish_spp, color, sex, size, fish_obs_time, fish_notes), by = "anem_table_id")
all_anems_with_fish_2 <- left_join(anem_db, fish_db %>% select(anem_table_id, fish_table_id, fish_spp, color, sex, size, fish_obs_time, fish_notes), by = "anem_table_id")
all_anems_with_fish_3 <- left_join(anems_Processed_all, allfish_caught %>% select(anem_table_id, fish_table_id, fish_spp, color, sex, size, fish_obs_time, fish_notes), by = "anem_table_id")
all_anems_with_fish_4 <- left_join(anems_Processed_all, fish_db %>% select(anem_table_id, fish_table_id, fish_spp, color, sex, size, fish_obs_time, fish_notes), by = "anem_table_id")


# Rename fish_spp and size columns so clear which fish data table it came from
all_anems_with_fish <- all_anems_with_fish %>%
  dplyr::rename(fish_spp_clownfish = fish_spp, size_clownfish = size)
all_anems_with_fish_2 <- all_anems_with_fish_2 %>%
  dplyr::rename(fish_spp_clownfish = fish_spp, size_clownfish = size)
all_anems_with_fish_3 <- all_anems_with_fish_3 %>%
  dplyr::rename(fish_spp_clownfish = fish_spp, size_clownfish = size)
all_anems_with_fish_4 <- all_anems_with_fish_4 %>%
  dplyr::rename(fish_spp_clownfish = fish_spp, size_clownfish = size)

# Now join in fish seen
all_anems_with_fish <- left_join(all_anems_with_fish, fish_seen_db %>% select(fish_spp, size, anem_table_id, est_table_id), by = "anem_table_id")
all_anems_with_fish_2 <- left_join(all_anems_with_fish_2, fish_seen_db %>% select(fish_spp, size, anem_table_id, est_table_id), by = "anem_table_id")
all_anems_with_fish_3<- left_join(all_anems_with_fish_3, fish_seen_db %>% select(fish_spp, size, anem_table_id, est_table_id), by = "anem_table_id")
all_anems_with_fish_4 <- left_join(all_anems_with_fish_4, fish_seen_db %>% select(fish_spp, size, anem_table_id, est_table_id), by = "anem_table_id")

# Rename fish_spp and size columns so clear came from clown_sightings data table
all_anems_with_fish <- all_anems_with_fish %>%
  dplyr::rename(fish_spp_sighted = fish_spp, size_sighted = size)
all_anems_with_fish_2 <- all_anems_with_fish_2 %>%
  dplyr::rename(fish_spp_sighted = fish_spp, size_sighted = size)
all_anems_with_fish_3 <- all_anems_with_fish_3 %>%
  dplyr::rename(fish_spp_sighted = fish_spp, size_sighted = size)
all_anems_with_fish_4 <- all_anems_with_fish_4 %>%
  dplyr::rename(fish_spp_sighted = fish_spp, size_sighted = size)

# Join in dive info too
all_anems_with_fish <- left_join(all_anems_with_fish, dives_db_processed, by = "dive_table_id")
all_anems_with_fish_2 <- left_join(all_anems_with_fish_2, dives_db_processed, by = "dive_table_id")

##### Do occupancy calcs 
# Test whether it matters if I use 0 dives or not - not at all for 2015, only very slightly for 2017
anemOcc_2015_W_a <- anemOcupancyByYear(2015, winter_months, anem_occ_dives, all_anems_with_fish)  # why nothing but APCL and UNOC? Had the other species in the earlier version of this (in Assessing_anemone_occupancy)
anemOcc_2017_a <- anemOcupancyByYear(2017, all_months, anem_occ_dives, all_anems_with_fish)

anemOcc_2015_W_b <- anemOcupancyByYear(2015, winter_months, clown_sample_dive_types, all_anems_with_fish)  # why nothing but APCL and UNOC? Had the other species in the earlier version of this (in Assessing_anemone_occupancy)
anemOcc_2017_b <- anemOcupancyByYear(2017, all_months, clown_sample_dive_types, all_anems_with_fish)

anemOcc_2015_W_c <- anemOcupancyByYear(2015, winter_months, anem_occ_dives, all_anems_with_fish_2)  # why nothing but APCL and UNOC? Had the other species in the earlier version of this (in Assessing_anemone_occupancy)
anemOcc_2017_c <- anemOcupancyByYear(2017, all_months, anem_occ_dives, all_anems_with_fish_2)

anemOcc_2015_W_d <- anemOcupancyByYear(2015, winter_months, clown_sample_dive_types, all_anems_with_fish_2)  # why nothing but APCL and UNOC? Had the other species in the earlier version of this (in Assessing_anemone_occupancy)
anemOcc_2017_d <- anemOcupancyByYear(2017, all_months, clown_sample_dive_types, all_anems_with_fish_2)

# Test whether using anem_db vs. anems_Processed_all matters and fish_db vs. allfish_caught
### Doesn't matter whether use anems_Processed_all or anems_db, using fish_db rather than allfish_caught does restore the missing non-APCL clownfish in 2015. 
# The number of APCL anems stays the same in 2015, regardless of anem or fish df used, but total number increases by 70 from 3635 when use allfish_caught to 3705 when use fish_db
# Not sure exactly why? Seems like there should have just been a shuffling from UNOC to other fish spp?
## the 2017 numbers are unchanged, doesn't matter if use allfish_caught vs. fish_db or anems_Processed_all vs. anem_db
anemOcc_2015_W <- anemOcupancyByYear(2015, winter_months, anem_occ_dives, all_anems_with_fish)  # why nothing but APCL and UNOC? Had the other species in the earlier version of this (in Assessing_anemone_occupancy)
anemOcc_2017 <- anemOcupancyByYear(2017, all_months, anem_occ_dives, all_anems_with_fish)
anemOcc_2015_W_2 <- anemOcupancyByYear(2015, winter_months, anem_occ_dives, all_anems_with_fish_2)  # why nothing but APCL and UNOC? Had the other species in the earlier version of this (in Assessing_anemone_occupancy)
anemOcc_2017_2 <- anemOcupancyByYear(2017, all_months, anem_occ_dives, all_anems_with_fish_2)
anemOcc_2015_W_3 <- anemOcupancyByYear(2015, winter_months, anem_occ_dives, all_anems_with_fish_3)  # why nothing but APCL and UNOC? Had the other species in the earlier version of this (in Assessing_anemone_occupancy)
anemOcc_2017_3 <- anemOcupancyByYear(2017, all_months, anem_occ_dives, all_anems_with_fish_3)
anemOcc_2015_W_4 <- anemOcupancyByYear(2015, winter_months, anem_occ_dives, all_anems_with_fish_4)  # why nothing but APCL and UNOC? Had the other species in the earlier version of this (in Assessing_anemone_occupancy)
anemOcc_2017_4 <- anemOcupancyByYear(2017, all_months, anem_occ_dives, all_anems_with_fish_4)

# Why are there 50 more anems in 2015 when use fish_db instead of allfish_caught?
test_fish_db_anems <- fish_db %>% filter(anem_table_id %in% (anems_Processed_all %>% filter(year == 2015))$anem_table_id) 
test_fish_db_anems %>% distinct(anem_table_id)  # 1847 anem_table_ids
test_allfish_caught_anems <- allfish_caught %>% filter(anem_table_id %in% (anems_Processed_all %>% filter(year == 2015))$anem_table_id)
test_allfish_caught_anems %>% distinct(anem_table_id)  # 1162 anem_table_ids
anem_diffs <- setdiff(test_fish_db_anems$anem_table_id, test_allfish_caught_anems$anem_table_id)  # 685 differences

test_allfish_caught_fish <- all_anems_with_fish %>% filter(anem_table_id %in% anem_diffs) %>% select(anem_table_id, anem_id, fish_table_id, fish_spp_clownfish, est_table_id, fish_spp_sighted, dive_type, date, site, year)
test_fish_db_fish <- all_anems_with_fish_2 %>% filter(anem_table_id %in% anem_diffs) %>% select(anem_table_id, anem_id, fish_table_id, fish_spp_clownfish, est_table_id, fish_spp_sighted, dive_type, date, site, year)

table(test_allfish_caught_fish$fish_spp_clownfish)  # no fish recorded on these anems from the fish_db
table(test_allfish_caught_fish$fish_spp_sighted)  # no records

table(test_fish_db_fish$fish_spp_clownfish)  # all but two have fish from fish_db (remaining two are NA)
table(test_fish_db_fish$fish_spp_sighted)  # no records

# Convert into available habitat just considering UNOC as available and considering all non-APCL anems available
non_APCL_clown <- c("APFR","APOC","APPE","APSE","PRBI")
non_APCL_clown_and_UNOC <- c(non_APCL_clown, "UNOC")
APCL_clown <- c("APCL", "APCL+APFR", "APCL+APPE", "APCL+APSA")

anemOcc_2017_nonAPCLclownperc = sum((anemOcc_2017$anems_Cspp %>% filter(Var1 %in% non_APCL_clown))$Perc)
anemOcc_2017_nonAPCLperc = sum((anemOcc_2017$anems_Cspp %>% filter(Var1 %in% non_APCL_clown_and_UNOC))$Perc)
anemOcc_2017_APCLperc = sum((anemOcc_2017$anems_Cspp %>% filter(Var1 %in% APCL_clown))$Perc)

anemOcc_2017_nonAPCLclownperc_2 = sum((anemOcc_2017_2$anems_Cspp %>% filter(Var1 %in% non_APCL_clown))$Perc)
anemOcc_2017_nonAPCLperc_2 = sum((anemOcc_2017_2$anems_Cspp %>% filter(Var1 %in% non_APCL_clown_and_UNOC))$Perc)
anemOcc_2017_APCLperc_2 = sum((anemOcc_2017_2$anems_Cspp %>% filter(Var1 %in% APCL_clown))$Perc)

anemOcc_2017_nonAPCLclownperc_3 = sum((anemOcc_2017_3$anems_Cspp %>% filter(Var1 %in% non_APCL_clown))$Perc)
anemOcc_2017_nonAPCLperc_3 = sum((anemOcc_2017_3$anems_Cspp %>% filter(Var1 %in% non_APCL_clown_and_UNOC))$Perc)
anemOcc_2017_APCLperc_3 = sum((anemOcc_2017_3$anems_Cspp %>% filter(Var1 %in% APCL_clown))$Perc)

anemOcc_2017_nonAPCLclownperc_4 = sum((anemOcc_2017_4$anems_Cspp %>% filter(Var1 %in% non_APCL_clown))$Perc)
anemOcc_2017_nonAPCLperc_4 = sum((anemOcc_2017_4$anems_Cspp %>% filter(Var1 %in% non_APCL_clown_and_UNOC))$Perc)
anemOcc_2017_APCLperc_4 = sum((anemOcc_2017_4$anems_Cspp %>% filter(Var1 %in% APCL_clown))$Perc)

# Hmm, these aren't quite right - would need to calculate perc nonAPCLclown for 2015 now too
anems_APCL_and_not_byyear <- data.frame(year = c(2015, 2017), 
                                     perc_APCL = c(anemOcc_2015_W$anems_Cspp$Perc[1], anemOcc_2017_APCLperc), 
                                     perc_non_APCL_clown = c(0, anemOcc_2017_nonAPCLclownperc),
                                     perc_UNOC = c(as.numeric((anemOcc_2015_W$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc), as.numeric((anemOcc_2017$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc)),
                                     perc_non_APCL_clown_or_UNOC = c(0, anemOcc_2017_nonAPCLperc))

anems_APCL_and_not_byyear_2 <- data.frame(year = c(2015, 2017), 
                                        perc_APCL = c(anemOcc_2015_W_2$anems_Cspp$Perc[1], anemOcc_2017_APCLperc_2), 
                                        perc_non_APCL_clown = c(0, anemOcc_2017_nonAPCLclownperc_2),
                                        perc_UNOC = c(as.numeric((anemOcc_2015_W_2$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc), as.numeric((anemOcc_2017_2$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc)),
                                        perc_non_APCL_clown_or_UNOC = c(0, anemOcc_2017_nonAPCLperc_2))

anems_APCL_and_not_byyear <- data.frame(year = c(2015, 2017), 
                                        perc_APCL = c(anemOcc_2015_W$anems_Cspp$Perc[1], anemOcc_2017_APCLperc), 
                                        perc_non_APCL_clown = c(0, anemOcc_2017_nonAPCLclownperc),
                                        perc_UNOC = c(as.numeric((anemOcc_2015_W$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc), as.numeric((anemOcc_2017$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc)),
                                        perc_non_APCL_clown_or_UNOC = c(0, anemOcc_2017_nonAPCLperc))

anems_APCL_and_not_byyear <- data.frame(year = c(2015, 2017), 
                                        perc_APCL = c(anemOcc_2015_W$anems_Cspp$Perc[1], anemOcc_2017_APCLperc), 
                                        perc_non_APCL_clown = c(0, anemOcc_2017_nonAPCLclownperc),
                                        perc_UNOC = c(as.numeric((anemOcc_2015_W$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc), as.numeric((anemOcc_2017$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc)),
                                        perc_non_APCL_clown_or_UNOC = c(0, anemOcc_2017_nonAPCLperc))

# Try with 0 dives for 2017 to compare
anemOcc_2017_nonAPCLclownperc_2 = sum((anemOcc_2017_2$anems_Cspp %>% filter(Var1 %in% non_APCL_clown))$Perc)
anemOcc_2017_nonAPCLperc_2 = sum((anemOcc_2017_2$anems_Cspp %>% filter(Var1 %in% non_APCL_clown_and_UNOC))$Perc)
anemOcc_2017_APCLperc_2 = sum((anemOcc_2017_2$anems_Cspp %>% filter(Var1 %in% APCL_clown))$Perc)

anems_APCL_and_not_byyear_2 <- data.frame(year = c(2015, 2017), 
                                        perc_APCL = c(anemOcc_2015_W$anems_Cspp$Perc[1], anemOcc_2017_APCLperc_2), 
                                        perc_non_APCL_clown = c(0, anemOcc_2017_nonAPCLclownperc_2),
                                        perc_UNOC = c(as.numeric((anemOcc_2015_W$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc), as.numeric((anemOcc_2017_2$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc)),
                                        perc_non_APCL_clown_or_UNOC = c(0, anemOcc_2017_nonAPCLperc_2))


# Where are the other clownfish spp in 2015?
test_for_2015 <- left_join(anems_Processed_all, fish_seen_db, by = "anem_table_id") %>%
  dplyr::rename(fish_spp_seen_db = fish_spp)
test_for_2015 <- left_join(test_for_2015, fish_db, by = "anem_table_id") %>%
  dplyr::rename(fish_spp_fish_db = fish_spp)
# there aren't any other clownfish spp besides APCL for 2015 in fish_seen_db but there are in fish_db - does this mean I need to update things?

# currently averaging between 2015W and 2017 values (but 2015W ones have some issues b/c not seeing any spp other than APCL in the data but did before...)
anems_APCL_and_not <- data.frame(perc_hab = c("APCL", "UNOC"), 
                                 value = c(mean(anems_APCL_and_not_byyear$perc_APCL),
                                           mean(anems_APCL_and_not_byyear$perc_UNOC)))

#################### Save output: ####################
saveRDS(anems_APCL_and_not, here::here("Data/Script_outputs", "anems_APCL_and_not.RData"))
