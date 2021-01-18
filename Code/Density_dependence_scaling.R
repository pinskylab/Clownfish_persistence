# Find proportion of anemones occupied by APCL vs. "empty"

#################### Set-up: ####################

# Pull data, functions, constants 
source(here::here('Code', 'Constants_database_common_functions.R'))

# Set vectors of months, years, and dive types
all_months <- c(1,2,3,4,5,6,7,8,9,10,11,12)  # for use in anemone occupancy function
anem_occ_dives <- clown_sample_dive_types  # use all dive types except for R (clown_sample_dive_types is in Constants_database_common_functions.R)

#################### Functions: ####################

# Find the number and proportion of anemones in each site and year occupied by APCL vs. occupied by other clownfish spp. vs. empty 
anemOcupancyByYear <- function(year_val, month_vals, divetype_vals, allanemsdf) {  # inputs: year_val - years to consider, months_val - months to consider, divetype_vals - dive types to consider, allanemsdf - df of anem sightings and fish on them (created below)
  
  # filter out just the relevant year and months
  allanemsdf_filt <- allanemsdf %>% filter(year == year_val) %>%  # filter by year
    filter(month %in% month_vals) %>%  # filter by months
    filter(dive_type %in% divetype_vals)  # filter by dive types
      
  # try to get one row per anem_table_id 
  # first, get distinct anem_table_id, fish_spp_clownfish, fish_spp_sighted combos
  allanemsdf_distinct <- allanemsdf_filt %>% distinct(anem_table_id, fish_spp_clownfish, fish_spp_sighted, .keep_all = TRUE) %>%  # this limits down some (gets rid of multiple obs because one per measured fish on the anem) but doesn't work if multiple species were seen (most anems appear once, some twice)
    group_by(anem_table_id) %>%
    mutate(nrow = n()) %>%
    ungroup()
  
  # for ones that only appear once after that, find different combos of species on anems
  anems_fish_spp_combined_1 <- allanemsdf_distinct %>% 
    filter(nrow == 1) %>%
    mutate(fish_spp = case_when(is.na(fish_spp_clownfish) & is.na(fish_spp_sighted) ~ 'UNOC',  # no fish
                                is.na(fish_spp_clownfish) & !is.na(fish_spp_sighted) ~ fish_spp_sighted,  # no fish caught but fish seen
                                !is.na(fish_spp_clownfish) & is.na(fish_spp_sighted) ~ fish_spp_clownfish,  # no fish seen but fish caught
                                !is.na(fish_spp_clownfish) & !is.na(fish_spp_sighted) & fish_spp_clownfish == fish_spp_sighted ~ fish_spp_clownfish,  # fish caught and seen but both same species
                                !is.na(fish_spp_clownfish) & !is.na(fish_spp_sighted) & fish_spp_clownfish != fish_spp_sighted ~ paste(fish_spp_clownfish, fish_spp_sighted, sep="+")))  # fish caught and seen but different species
  
  # for ones that appear multiple times, easier to sort by year
  # in 2015, all fish_spp_sighted values are NA
  if(year_val == 2015) {
    anems_fish_spp_combined_2 <- allanemsdf_distinct %>%
      filter(nrow == 2) %>%
      group_by(anem_table_id) %>%
      mutate(fish_spp = case_when(is.na(fish_spp_clownfish[1]) & is.na(fish_spp_clownfish[2]) ~ "UNOC",  # if both are NA, UNOC - this shouldn't happen
                                  !is.na(fish_spp_clownfish[1]) & is.na(fish_spp_clownfish[2]) ~ fish_spp_clownfish[1],
                                  is.na(fish_spp_clownfish[1]) & !is.na(fish_spp_clownfish[2]) ~ fish_spp_clownfish[2],
                                  !is.na(fish_spp_clownfish[1]) & !is.na(fish_spp_clownfish[2]) & fish_spp_clownfish[1] != fish_spp_clownfish[2] ~ paste(fish_spp_clownfish[1], fish_spp_clownfish[2], sep="+"),
                                  !is.na(fish_spp_clownfish[1]) & !is.na(fish_spp_clownfish[2]) & fish_spp_clownfish[1] == fish_spp_clownfish[2] ~ fish_spp_clownfish[1])) %>%  # this one shouldn't happen either 
      ungroup() %>%
      distinct(anem_table_id, .keep_all = TRUE)
    
    anems_fish_spp_combined_multiple <- anems_fish_spp_combined_2
  }
  
  # in 2017, sometimes APCL appears both in fish_spp_clownfish and fish_spp_sighted and there are some anem_table_ids that appear 3 times, fish_spp_clownfish is only ever APCL or NA, fish_spp_sighted is never NA
  if(year_val == 2017) { 
    # for anem_table_ids that appear 2 times
    anems_fish_spp_combined_2 <- allanemsdf_distinct %>% 
      filter(nrow == 2) %>%
      group_by(anem_table_id) %>%
      mutate(fish_spp = case_when(is.na(fish_spp_clownfish[1]) & is.na(fish_spp_clownfish[2]) & is.na(fish_spp_sighted[1]) & is.na(fish_spp_sighted[2]) ~ "UNOC",  # no fish, shouldn't happen
                                  is.na(fish_spp_clownfish[1]) & is.na(fish_spp_clownfish[2]) & !is.na(fish_spp_sighted[1]) & !is.na(fish_spp_sighted[2]) & fish_spp_sighted[1] != fish_spp_sighted[2] ~ paste(fish_spp_sighted[1],fish_spp_sighted[2],sep="+"),  # both clown caughts are NA, clown sighteds are different
                                  (fish_spp_clownfish[1] == "APCL" | fish_spp_clownfish[2] == "APCL") & fish_spp_sighted[1] != "APCL" & fish_spp_sighted[2] != "APCL" & fish_spp_sighted[1] != fish_spp_sighted[2] ~ paste("APCL",fish_spp_sighted[1],fish_spp_sighted[2],sep="+"),  # one of the clown caughts is APCL, neither clown sighted is but both are different
                                  (fish_spp_clownfish[1] == "APCL" | fish_spp_clownfish[2] == "APCL") & fish_spp_sighted[1] != "APCL" & fish_spp_sighted[2] != "APCL" & fish_spp_sighted[1] == fish_spp_sighted[2] ~ paste("APCL",fish_spp_sighted[1], sep = "+"),  # one of the clown caughts is APCL, neither clown sighted is but they are the same
                                  (fish_spp_clownfish[1] == "APCL" | fish_spp_clownfish[2] == "APCL") & (fish_spp_sighted[1] == "APCL" | fish_spp_sighted[2] == "APCL") & fish_spp_sighted[1] != fish_spp_sighted[2] ~ paste(fish_spp_sighted[1],fish_spp_sighted[2],sep="+"),  # one of the clown caughts is APCL and one of the clown sighted is but not both
                                  (fish_spp_clownfish[1] == "APCL" | fish_spp_clownfish[2] == "APCL") & fish_spp_sighted[1] == "APCL" & fish_spp_sighted[2] == "APCL" ~ "APCL")) %>%
      ungroup() %>%
      distinct(anem_table_id, .keep_all = TRUE)
                                  
    # for anem_table_ids that appear 3 times (only one case so not doing a thorough search of possibilities in case_when)
    anems_fish_spp_combined_3 <- allanemsdf_distinct %>%
      filter(nrow == 3) %>%
      group_by(anem_table_id) %>%
      mutate(fish_spp = case_when((fish_spp_sighted[1] == "APCL" | fish_spp_sighted[2] == "APCL" | fish_spp_sighted[3] == "APCL") ~ paste(fish_spp_sighted[1],fish_spp_sighted[2],fish_spp_sighted[3],sep="+"))) %>%
      ungroup() %>%
      distinct(anem_table_id, .keep_all = TRUE)
    
    anems_fish_spp_combined_multiple = rbind(anems_fish_spp_combined_2, anems_fish_spp_combined_3)
                                 
    }
  
   # combine those that appear once with those that appear multiple times
  anems_fish_spp_combined <- rbind(anems_fish_spp_combined_1, anems_fish_spp_combined_multiple)
  
  # find the total number of anems
  totalanems <- length(anems_fish_spp_combined$anem_table_id)
  
  # match up equivalent combination categories
  anems_fish_spp_combined <- anems_fish_spp_combined %>%
    mutate(fish_spp_tidied = case_when(fish_spp == "APCL" ~ "APCL",
                                       fish_spp == "APFR" ~ "APFR",
                                       fish_spp == "APOC" ~ "APOC",
                                       fish_spp == "APPO" ~ "APPO",
                                       fish_spp == "APPE" ~ "APPE",
                                       fish_spp == "APSE" ~ "APSA",  # assuming this one is a typo...
                                       fish_spp == "APSA" ~ "APSA",
                                       fish_spp == "PRBI" ~ "PRBI",
                                       fish_spp == "UNOC" ~ "UNOC",
                                       (fish_spp == "APCL+APFR" | fish_spp == "APFR+APCL") ~ "APCL+APFR",
                                       (fish_spp == "APCL+APPE" | fish_spp == "APPE+APCL") ~ "APCL+APPE",
                                       (fish_spp == "APCL+APSA" | fish_spp == "APSA+APCL") ~ "APCL+APSA",
                                       (fish_spp == "APCL+APPO" | fish_spp == "APPO+APCL") ~ "APCL+APPO",
                                       (fish_spp == "APPE+APFR" | fish_spp == "APFR+APPE") ~ "APFR+APPE",
                                       fish_spp == "APCL+APSA+APPE" ~ "APCL+APSA+APPE"))
    
  #find frequency and percentage of each spp on anems 
  anems_Cspp <- as.data.frame(table(anems_fish_spp_combined$fish_spp_tidied, useNA = "ifany"), stringsAsFactors = FALSE)
  anems_Cspp$Perc <- anems_Cspp$Freq/totalanems
  anems_Cspp <- anems_Cspp %>%
    mutate(year = year_val)
  
  anems_Cspp_by_site <- anems_fish_spp_combined %>%
    group_by(site, fish_spp_tidied) %>%
    summarize(nAnems = n(),
              year = year_val)
  
  out <- list(anems_Cspp = anems_Cspp, anems_Cspp_by_site = anems_Cspp_by_site, totalanems = totalanems)
  
  return(out)
  
}

#################### Running things: ####################

##### Data set up
# Join together all anems with fish in fish_db
all_anems_with_fish <- left_join(anem_db, fish_db %>% select(anem_table_id, fish_table_id, fish_spp, color, sex, size, fish_obs_time, fish_notes), by = "anem_table_id")

# Rename fish_spp and size columns so clear which fish data table it came from
all_anems_with_fish <- all_anems_with_fish %>%
  dplyr::rename(fish_spp_clownfish = fish_spp, size_clownfish = size)

# Now join in fish seen
all_anems_with_fish <- left_join(all_anems_with_fish, fish_seen_db %>% select(fish_spp, size, anem_table_id, est_table_id), by = "anem_table_id")

# Rename fish_spp and size columns so clear came from clown_sightings data table
all_anems_with_fish <- all_anems_with_fish %>%
  dplyr::rename(fish_spp_sighted = fish_spp, size_sighted = size)

# Join in dive info too
all_anems_with_fish <- left_join(all_anems_with_fish, dives_db_processed, by = "dive_table_id")

##### Do occupancy calcs 
anemOcc_2015_W <- anemOcupancyByYear(2015, winter_months, anem_occ_dives, all_anems_with_fish)
anemOcc_2017 <- anemOcupancyByYear(2017, all_months, anem_occ_dives, all_anems_with_fish)

# For revisions, find perc UNOC, perc occupied by non-APCL clownfish, and perc occupied by APCL clownfish in 2015 and 2017
APCL_list <- c("APCL","APCL+APFR","APCL+APPE","APCL+APSA")
other_clown_list <- c("APFR","APFR+APPE","APOC","APPE","APPO","APSA","PRBI")

perc_APCL_2015 <- anemOcc_2015_W$anems_Cspp %>% filter(Var1 %in% APCL_list) %>% summarize(sumAPCL = sum(Perc))
perc_otherclown_2015 <- anemOcc_2015_W$anems_Cspp %>% filter(Var1 %in% other_clown_list) %>% summarize(sumOther = sum(Perc))
perc_UNOC_2015 <- anemOcc_2015_W$anems_Cspp %>% filter(Var1 == "UNOC") %>% summarize(sumUNOC = sum(Perc))

perc_APCL_2017 <- anemOcc_2017$anems_Cspp %>% filter(Var1 %in% APCL_list) %>% summarize(sumAPCL = sum(Perc))
perc_otherclown_2017 <- anemOcc_2017$anems_Cspp %>% filter(Var1 %in% other_clown_list) %>% summarize(sumOther = sum(Perc))
perc_UNOC_2017 <- anemOcc_2017$anems_Cspp %>% filter(Var1 == "UNOC") %>% summarize(sumUNOC = sum(Perc))

# Convert into available habitat just considering UNOC as available and considering all non-APCL anems available
non_APCL_clown <- c("APFR","APOC","APPE","APSE","PRBI","APSA","APPO","APFR+APPE")
non_APCL_clown_and_UNOC <- c(non_APCL_clown, "UNOC")
APCL_clown <- c("APCL", "APCL+APFR", "APCL+APPE", "APCL+APSA", "APCL+APPO", "APCL+APSA+APPE")

# Find total perc and number of anems occupied by different catgories
# Perc
anemOcc_2015_W_nonAPCLclownperc = sum((anemOcc_2015_W$anems_Cspp %>% filter(Var1 %in% non_APCL_clown))$Perc)  # anems occupied by clownfish that aren't APCL
anemOcc_2015_W_nonAPCLperc = sum((anemOcc_2015_W$anems_Cspp %>% filter(Var1 %in% non_APCL_clown_and_UNOC))$Perc)  # anems not occupied by APCL
anemOcc_2015_W_APCLperc = sum((anemOcc_2015_W$anems_Cspp %>% filter(Var1 %in% APCL_clown))$Perc)  # anems occupied by APCL (even if other clownfish present too)

anemOcc_2017_nonAPCLclownperc = sum((anemOcc_2017$anems_Cspp %>% filter(Var1 %in% non_APCL_clown))$Perc)  # anems occupied by clownfish that aren't APCL
anemOcc_2017_nonAPCLperc = sum((anemOcc_2017$anems_Cspp %>% filter(Var1 %in% non_APCL_clown_and_UNOC))$Perc)  # anems not occupied by APCL
anemOcc_2017_APCLperc = sum((anemOcc_2017$anems_Cspp %>% filter(Var1 %in% APCL_clown))$Perc)  # anems occupied by APCL (even if other clownfish present too)

# Number anems
anemOcc_2015_W_nonAPCLclown_anems = sum((anemOcc_2015_W$anems_Cspp %>% filter(Var1 %in% non_APCL_clown))$Freq)  # anems occupied by clownfish that aren't APCL
anemOcc_2015_W_nonAPCL_anems = sum((anemOcc_2015_W$anems_Cspp %>% filter(Var1 %in% non_APCL_clown_and_UNOC))$Freq)  # anems not occupied by APCL
anemOcc_2015_W_APCL_anems = sum((anemOcc_2015_W$anems_Cspp %>% filter(Var1 %in% APCL_clown))$Freq)  # anems occupied by APCL (even if other clownfish present too)

anemOcc_2017_nonAPCLclown_anems = sum((anemOcc_2017$anems_Cspp %>% filter(Var1 %in% non_APCL_clown))$Freq)  # anems occupied by clownfish that aren't APCL
anemOcc_2017_nonAPCL_anems = sum((anemOcc_2017$anems_Cspp %>% filter(Var1 %in% non_APCL_clown_and_UNOC))$Freq)  # anems not occupied by APCL
anemOcc_2017_APCL_anems = sum((anemOcc_2017$anems_Cspp %>% filter(Var1 %in% APCL_clown))$Freq)  # anems occupied by APCL (even if other clownfish present too)

# Just consider anems occupied by APCL (somehow) or UNOC as the "total habitat"
anems_APCL_and_not_by_year <- data.frame(year = c(2015, 2017),
                                         anems_APCL = c(anemOcc_2015_W_APCL_anems, anemOcc_2017_APCL_anems),
                                         anems_other_clown = c(anemOcc_2015_W_nonAPCLclown_anems, anemOcc_2017_nonAPCLclown_anems),
                                         anems_UNOC = c(as.numeric((anemOcc_2015_W$anems_Cspp %>% filter(Var1 == "UNOC"))$Freq), as.numeric((anemOcc_2017$anems_Cspp %>% filter(Var1 == "UNOC"))$Freq)),
                                         anems_non_APCL = c(anemOcc_2015_W_nonAPCL_anems, anemOcc_2017_nonAPCL_anems),
                                         perc_APCL = c(anemOcc_2015_W_APCLperc, anemOcc_2017_APCLperc),
                                         perc_other_clown = c(anemOcc_2015_W_nonAPCLclownperc, anemOcc_2017_nonAPCLclownperc),
                                         perc_UNOC = c(as.numeric((anemOcc_2015_W$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc), as.numeric((anemOcc_2017$anems_Cspp %>% filter(Var1 == "UNOC"))$Perc)),
                                         perc_non_APCL = c(anemOcc_2015_W_nonAPCLperc, anemOcc_2017_nonAPCLperc),
                                         total_anems = c(anemOcc_2015_W$totalanems, anemOcc_2017$totalanems), stringsAsFactors = FALSE) %>%
  mutate(total_APCL_and_UNOC_anems = anems_APCL + anems_UNOC,
         total_APCL_and_UNOC_perc = perc_APCL + perc_UNOC)

# currently averaging between 2015W and 2017 values 
anems_APCL_and_not <- data.frame(perc_hab = c("APCL", "UNOC"), 
                                 value = c(mean(anems_APCL_and_not_by_year$perc_APCL),
                                           mean(anems_APCL_and_not_by_year$perc_UNOC)))

#################### Save output: ####################
save(anems_APCL_and_not, file=here::here("Data/Script_outputs", "anems_APCL_and_not.RData"))
save(anems_APCL_and_not_by_year, file=here::here("Data/Script_outputs", "anems_APCL_and_not_by_year.RData"))
