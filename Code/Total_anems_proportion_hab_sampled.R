# Calculate total number of anemones at each site and proportion of habitat sampled each year based on anemone numbers

#### LOOKS LIKE PLOTS DON'T QUITE MATCH PREVIOUS ESTIMATES - WHY NOT?

##### What if instead of doing total area across years, do total anemones across years (like we do for within-year prop hab sampled)

#################### Set-up: ####################
# Load packages
library(stringr)
library(ggplot2)

# Pull data, functions, constants 
source(here::here('Code', 'Constants_database_common_functions.R'))

# Set sites to include for total possible sampling area (all site areas times all years sampled) - excluding Caridad Proper, Sitio Lonas, Sitio Tugas from this because they disappeared partway through (check with Michelle on that) (should I exclude the Sitio Lonas match too?)
sites_for_total_areas <- c('Cabatoan', 'Caridad Cemetery', 'Elementary School', 'Gabas', 
                           'Haina', 'Hicgop South', 'N. Magbangon', 'Palanas', 'Poroc Rose',
                           'Poroc San Flower', 'San Agustin', 'Sitio Baybayon', 'Tamakin Dacot',
                           'Visca', 'Wangag', 'S. Magbangon')  

#################### Functions: ####################
# Find the number of tagged anemones visited by each site, each year 
pull_anems_by_year <- function(anemsdf, allanems, divesdf, year_i, site_i, survey_months, dive_types) {  # anemsdf is processed anems (anems_Processed), allanems are anems associated with APCL - regardless if tagged, divesdf is all dives, year_i is year of interest, site_i is site, survey_months is list of months when clownfish were sampled
  
  # Pull out dives relevant to year and site
  dives <- divesdf %>%
    filter(site == site_i, year == year_i, month %in% survey_months) %>%  # pull out dives of the right site, year, and season (relevant for 2015)
    filter(dive_type %in% dive_types)  # pull out dives of the right type
  
  # Pull out tagged anems
  tagged_anems_from_dives <- anemsdf %>%
    filter(dive_table_id %in% dives$dive_table_id) %>%  # just the anems from those dives
    distinct(anem_id_unq, .keep_all = TRUE) %>%  # only keep one row per anem (since some seen multiple times a season)
    summarize(anems_sampled = n())  # count the total number of anems
  
 # Pull out sampled anems
  sampled_anems_from_dives <- allanems %>%
    filter(site == site_i, year == year_i, month %in% survey_months) %>%
    summarize(anems_sampled =  n())
  
  # Anemones weren't tagged in 2012, so check if the year is 2012 to know what kind of anem observations to count in 'best estimate' of anems sampled that year
  if(year_i != 2012) {  
    
    # If not 2012, pull out tagged anems seen on those dives
    anems_visited = as.integer(tagged_anems_from_dives)
  } else if (year_i == 2012) {
    
    # If it is 2012, use anems sampled without specifying that they are tagged (do I need to worry about dive types?)
    anems_visited = as.integer(sampled_anems_from_dives)
  }

  # Put anem counts together in a list (tagged_anems for tagged ones, sampled_anems for sampled, visited_anems for best est of those visited)
  out <- list(tagged_anems = as.integer(tagged_anems_from_dives), sampled_anems = as.integer(sampled_anems_from_dives), visited_anems = anems_visited)
  
  return(out) 
}

#################### Running things: ####################
##### Pull out anems associated with fish from clownfish table
# Pull out fish
APCL_caught <- fish_db %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, gen_id, anem_table_id, recap, tag_id, color, size) %>%
  filter(fish_spp == "APCL")

# and their associated anems 
APCL_anems <- anem_db %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_obs_time, anem_id, anem_obs, old_anem_id) %>%
  filter(anem_table_id %in% APCL_caught$anem_table_id)

# and their associated dives
APCL_dives <- dives_db %>%
  select(dive_table_id, dive_type, date, site, gps, year, month) %>%
  filter(dive_table_id %in% APCL_anems$dive_table_id)

# Join together anem and dive info
all_APCL_anems <- left_join(APCL_anems, APCL_dives, by="dive_table_id")  # this is about 250 fewer than anems_Processed... not sure why...

# Remove intermediate data frames for neatness
rm(APCL_caught, APCL_anems, APCL_dives)

##### Find total habitat at, as indicated by anemones, via four methods
# Create data frames to store the output for various methods of assessing total anemones
methods = c("metal tags", "all tags", "seen twice", "2015 survey")

total_anems_by_site_metal <- data.frame(site = site_vec_order$site_name) %>%
  mutate(method = "metal tags",
         n_total_anems = NA)
total_anems_by_site_alltags <- data.frame(site = site_vec_order$site_name) %>%
  mutate(method = "all tags",
         n_total_anems = NA)
total_anems_by_site_seen2x <- data.frame(site = site_vec_order$site_name) %>%
  mutate(method = "seen twice",
         n_total_anems = NA)
total_anems_by_site_2015W <- data.frame(site = site_vec_order$site_name) %>%
  mutate(method = "2015 survey",
         n_total_anems = NA)

# Go through sites, pull out total number of anems by the different methods
for(i in 1:length(site_vec_order$site_name)) {
  
  # anemones tagged with metal tags (less likely to get overgrown/lost)
  metal_anems_at_site <- anems_Processed %>%
    filter(site == site_vec_order$site_name[i]) %>%  # just pull out anemones at the relevant site
    filter(anem_id >= first_metal_tag) %>%  # just those with metal tags
    distinct(anem_id_unq, .keep_all = TRUE) %>%  # keep only one row for each (if seen multiple times across years)
    summarize(total_anems = n())
  
  total_anems_by_site_metal$n_total_anems[i] = as.integer(metal_anems_at_site)  # as integer, rather than data frame, so doesn't turn into a list
  
  # all tagged anemones
  all_tagged_anems_at_site <- anems_Processed %>%
    filter(site == site_vec_order$site_name[i]) %>%  # just pull out anemones at the relevant site
    distinct(anem_id_unq, .keep_all = TRUE) %>%  # keep only one row for each (should all have anem_ids b/c that was already filtered for in anems_Processed)
    summarize(total_anems = n())
  
  total_anems_by_site_alltags$n_total_anems[i] = as.integer(all_tagged_anems_at_site)
  
  # anemones seen twice
  anems_seen_twice <- anems_Processed %>%
    filter(site == site_vec_order$site_name[i]) %>%  # just pull out anemones at the relevant site
    distinct(anem_id_unq, year, .keep_all = TRUE) %>%  # pull out just one observation of each tagged anem in a particular year
    group_by(anem_id_unq) %>%  # group by anem
    summarize(n_years = n()) %>%  # count how many years they were observed
    filter(n_years >= 2) %>%  # filter out just those seen two or more years
    summarize(total_anems = n())  
  
  total_anems_by_site_seen2x$n_total_anems[i] = as.integer(anems_seen_twice)

  # anemones seen on the 2015 winter anemone survey
  anems_seen_2015W <- anems_Processed %>%
    filter(site == site_vec_order$site_name[i]) %>%  # just pull out anemones at the relevant site
    filter(year == 2015, anem_month %in% winter_months) %>%  # just pull out anemones seen during the 2015 winter field season anemone survey
    distinct(anem_id_unq, .keep_all = TRUE) %>%  # filter out one observation per anemone
    summarize(total_anems = n())
  
  total_anems_by_site_2015W$n_total_anems[i] = as.integer(anems_seen_2015W)
  
}

# Bind together the four methods into one data frame
total_anems_by_site <- rbind(total_anems_by_site_metal, total_anems_by_site_alltags, total_anems_by_site_seen2x, total_anems_by_site_2015W)

# Remove the individual method data frames, for neatness
rm(total_anems_by_site_metal, total_anems_by_site_alltags, total_anems_by_site_seen2x, total_anems_by_site_2015W)

##### Find the amount of the site visited, as estimated by number of anems visited each year in each site
# Set up data frames for each year
anems_visited_2012 <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
  mutate(year = 2012,
         n_tagged_anems = NA,
         n_sampled_anems = NA,
         n_anems = NA)

anems_visited_2013 <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
  mutate(year = 2013,
         n_tagged_anems = NA,
         n_sampled_anems = NA,
         n_anems = NA)

anems_visited_2014 <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
  mutate(year = 2014,
         n_tagged_anems = NA,
         n_sampled_anems = NA,
         n_anems = NA)

anems_visited_2015 <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
  mutate(year = 2015,
         n_tagged_anems = NA,
         n_sampled_anems = NA,
         n_anems = NA)

anems_visited_2016 <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
  mutate(year = 2016,
         n_tagged_anems = NA,
         n_sampled_anems = NA,
         n_anems = NA)

anems_visited_2017 <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
  mutate(year = 2017,
         n_tagged_anems = NA,
         n_sampled_anems = NA,
         n_anems = NA)

anems_visited_2018 <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
  mutate(year = 2018,
         n_tagged_anems = NA,
         n_sampled_anems = NA,
         n_anems = NA)

# Run through sites each year to find number of anemones visited each year at each site
for(i in 1:length(site_vec_order$site_name)) {
  # 2012
  anems_visited_2012$n_tagged_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2012, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$tagged_anems
  anems_visited_2012$n_sampled_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2012, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$sampled_anems
  anems_visited_2012$n_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2012, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$visited_anems

  # 2013
  anems_visited_2013$n_tagged_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2013, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$tagged_anems
  anems_visited_2013$n_sampled_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2013, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$sampled_anems
  anems_visited_2013$n_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2013, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$visited_anems

  # 2014
  anems_visited_2014$n_tagged_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2014, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$tagged_anems
  anems_visited_2014$n_sampled_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2014, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$sampled_anems
  anems_visited_2014$n_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2014, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$visited_anems
  
  # 2015
  anems_visited_2015$n_tagged_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2015, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$tagged_anems
  anems_visited_2015$n_sampled_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2015, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$sampled_anems
  anems_visited_2015$n_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2015, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$visited_anems
  
  # 2016
  anems_visited_2016$n_tagged_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2016, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$tagged_anems
  anems_visited_2016$n_sampled_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2016, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$sampled_anems
  anems_visited_2016$n_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2016, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$visited_anems
  
  # 2017
  anems_visited_2017$n_tagged_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2017, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$tagged_anems
  anems_visited_2017$n_sampled_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2017, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$sampled_anems
  anems_visited_2017$n_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2017, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$visited_anems
  
  # 2018
  anems_visited_2018$n_tagged_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2018, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$tagged_anems
  anems_visited_2018$n_sampled_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2018, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$sampled_anems
  anems_visited_2018$n_anems[i] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db, 2018, site_vec_order$site_name[i], spring_months, clown_sample_dive_types))$visited_anems
}

# Bind the years together
anems_visited_by_year <- rbind(anems_visited_2012, anems_visited_2013, anems_visited_2014, anems_visited_2015,
                     anems_visited_2016, anems_visited_2017, anems_visited_2018)

# Remove intermediate data frames for neatness
rm(anems_visited_2012, anems_visited_2013, anems_visited_2014, anems_visited_2015, anems_visited_2016, anems_visited_2017, anems_visited_2018)

# Add in total anems per site
anems_visited_by_year <- left_join(anems_visited_by_year, total_anems_by_site, by = "site")

##### Calculate proportion habitat sampled and create a tidied version where NaN, Inf, and values >1 are edited
anems_visited_by_year <- anems_visited_by_year %>%
  mutate(prop_hab_sampled = n_anems/n_total_anems) %>%  # for all methods of determining total number of tags, calculate prop_hab_sampled
  mutate(n_anems_tidied = case_when(n_anems <= n_total_anems ~ n_anems,
                                   n_anems > n_total_anems ~ n_total_anems)) %>%  # to avoid prop_hab_sampled > 1, restrict n_anems in each year and site to the total at that site
  mutate(prop_hab_sampled_tidied = case_when(n_anems_tidied/n_total_anems <= 1 ~ n_anems_tidied/n_total_anems,  
                                              n_anems_tidied/n_total_anems == Inf | is.na(n_anems_tidied/n_total_anems) ~ 0))  # make a tidier version of prop_hab_sampled (no Inf, NaN, > 1 values), using n_anems_tidied

# # How I was doing this tidying before - moved to tidying anems so it carries over into the cumulative counts (anems over the total at a site don't count in the total sampled at all sites across years)
# mutate(prop_hab_sampled_tidied = case_when(prop_hab_sampled <= 1 ~ prop_hab_sampled,  # if it's a real number and less than one, stick with original prop_hab_sampled
#                                            prop_hab_sampled == Inf | is.na(prop_hab_sampled) ~ 0,  # if it's infinite or otherwise not a number (b/c total anems is 0), use 0
#                                            prop_hab_sampled > 1 ~ 1)) %>%  # if it's bigger than 1, round down to 1 

##### Calculate the total amount of habitat visited over time (used in estimating egg-recruit survival)
# Sum up total site area (all site areas times all years sampled) - total possible sampling area 
site_areas_modified <- site_areas %>% filter(site %in% sites_for_total_areas)  # Pull just the area from those sites

# Find amount of habitat sampled overall - sum of area sampled in each year
sampled_area_each_year <- left_join(site_areas_modified, anems_visited_by_year, by = 'site') %>%  # join with total area
  mutate(area_sampled = prop_hab_sampled_tidied*kmsq_area)  # for each method, site, and year, find area sampled in kmsq, using tidied-up prop hab sampled

time_frames <- c("2012", "2012-2013", "2012-2014", "2012-2015", "2012-2016", "2012-2017", "2012-2018")

# Set up a data frame for collecting total area sampled across different time frames
time_frame_list <- rep(time_frames[1], length(methods))
for (i in 2:length(time_frames)) {
  time_frame_list <- c(time_frame_list, rep(time_frames[i], length(methods)))
}

total_area_sampled_through_time <- data.frame(method = rep(methods, length(time_frames)),
                                              time_frame = time_frame_list, stringsAsFactors = FALSE) %>%
  mutate(total_possible_sample_area_km2 = NA,
         total_area_sampled_km2 = NA,
         total_prop_hab_sampled_area = NA,
         end_year = str_sub(time_frame, -4, -1),
         total_possible_sample_anems = NA,  # total number of anems that would have been possible to sample in the time frame (sum of total at all sites X number of years of sampling)
         total_anems_sampled = NA,  # sum of total number of anems sampled at all sites in that time period
         total_anems_sampled_tidied = NA,  # sum of total number of anems sampled at all sites in that time period but not including anems that exceed the total at a site (so if 15 sampled recorded at a site with 10 estimated total, this number sums using 10)
         total_prop_hab_sampled_anems = NA,
         total_prop_hab_sampled_anems_tidied = NA)  # doesn't include anems sampled at a site that exceed the total anems at that site estimated by each method

# total_area_sampled <- data.frame(method = methods) %>%
#   mutate(total_area_sampled = NA,
#          total_prop_hab_sampled = NA,
#          total_area_sampled_2012to2015 = NA,
#          total_prop_hab_sampled_2012to2015 = NA)

# Find total area that could have been sampled, adding a sample year each time (for parentage)
for (i in 1:length(total_area_sampled_through_time$method)) {
  end_year <- as.integer(total_area_sampled_through_time$end_year[i])
  years_to_include <- seq(2012, end_year, by=1)
  
  # Find possible area to sample, area sampled, and prop_hab sampled using area
  total_area_sampled_through_time$total_possible_sample_area_km2[i] = sum(site_areas_modified$kmsq_area*(length(years_to_include)))
  total_area_sampled_through_time$total_area_sampled_km2[i] = sum((sampled_area_each_year %>% filter(year %in% years_to_include & method == total_area_sampled_through_time$method[i]))$area_sampled)
  #total_area_sampled_through_time$total_area_sampled[i] = sum((sampled_area_each_year %>% filter(method == total_area_sampled_through_time$methods[i], year %in% years_to_include))$area_sampled)
  total_area_sampled_through_time$total_prop_hab_sampled_area[i] = total_area_sampled_through_time$total_area_sampled_km2[i]/total_area_sampled_through_time$total_possible_sample_area_km2[i]

  # Find possible area to sample, area sampled, and prob_hab sampled using anems and proportion of anems, rather than converting to area (like do for within-year proportion hab sampled)
  total_area_sampled_through_time$total_possible_sample_anems[i] = sum((sampled_area_each_year %>% filter(year %in% years_to_include & method == total_area_sampled_through_time$method[i]))$n_total_anems)
  total_area_sampled_through_time$total_anems_sampled[i] = sum((sampled_area_each_year %>% filter(year %in% years_to_include & method == total_area_sampled_through_time$method[i]))$n_anems)
  total_area_sampled_through_time$total_anems_sampled_tidied[i] = sum((sampled_area_each_year %>% filter(year %in% years_to_include & method == total_area_sampled_through_time$method[i]))$n_anems_tidied)
  total_area_sampled_through_time$total_prop_hab_sampled_anems[i] = total_area_sampled_through_time$total_anems_sampled[i]/total_area_sampled_through_time$total_possible_sample_anems[i]
  total_area_sampled_through_time$total_prop_hab_sampled_anems_tidied[i] = total_area_sampled_through_time$total_anems_sampled_tidied[i]/total_area_sampled_through_time$total_possible_sample_anems[i]
}

# Make into a data frame for easier comparison plotting
total_sampling_across_years <- data.frame(total_anems_method = rep(total_area_sampled_through_time$method, 3),
                                          time_frame = rep(total_area_sampled_through_time$time_frame, 3),
                                          total_prop_hab_sampled = c(total_area_sampled_through_time$total_prop_hab_sampled_area, total_area_sampled_through_time$total_prop_hab_sampled_anems, total_area_sampled_through_time$total_prop_hab_sampled_anems_tidied),
                                          total_area_method = c(rep("area", length(total_area_sampled_through_time$method)), rep("anems", length(total_area_sampled_through_time$method)), rep("anems tidied", length(total_area_sampled_through_time$method))))

# total_area_all_years <- sum(site_areas_modified$kmsq_area)*length(years_sampled)
# 
# total_area_2012to2015 <- sum(site_areas_modified$kmsq_area)*length(years_parentage)  # until all genetic data is in parentage, only act as if we sampled in 2012-2015 (for egg-recruit survival estimate)
# 
# # Using each method of determining number of anemones at a site, go through and find total area sampled in each year
# # (Is this the best way of doing it? Total area sampled doesn't change within a year and should be some way of determing that like from tracks or something... but total area of site would...))
# for(i in 1:length(methods)) {
#   # for all sampling years (2012-2018)
#   total_area_sampled$total_area_sampled[i] = sum((sampled_area_each_year %>% filter(method == methods[i]))$area_sampled)
#   total_area_sampled$total_prop_hab_sampled[i] = total_area_sampled$total_area_sampled[i]/total_area_all_years
#   
#   # for early sampling years (2012-2015), currently the ones in parentage analysis and kernel estimates
#   total_area_sampled$total_area_sampled_2012to2015[i] = sum((sampled_area_each_year %>% filter(method == methods[i]) %>% filter(year %in% years_parentage))$area_sampled)
#   total_area_sampled$total_prop_hab_sampled_2012to2015[i] = total_area_sampled$total_area_sampled_2012to2015[i]/total_area_2012to2015
# }

# Remove intermediate data frames for neatness
rm(site_areas_modified, sampled_area_each_year)

#################### Plots: #################### (already have versions of these in this folder from the other script but now can compare)
# Look at cumulative proportion total habitat sampled, using metal tags as total anemones
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "Cumulative_prop_hab_sampled_method_comp.pdf"))
ggplot(data = total_sampling_across_years %>% filter(total_anems_method == "metal tags"), aes(x = time_frame, y = total_prop_hab_sampled, fill = total_area_method)) +
  geom_bar(position = "dodge", stat = "identity") +
  ggtitle("Cumulative proportion habitat sampled (metal tags)") + xlab("Sampling years") + ylab("Proportion sampled") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
dev.off()
  
# Look at cumulative proportion total habitat sampled, using metal tags as total anemones
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "Cumulative_prop_hab_sampled_method_comp_all_anems_methods.pdf"))
ggplot(data = total_sampling_across_years, aes(x = time_frame, y = total_prop_hab_sampled, fill = total_area_method)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~total_anems_method) +
  ggtitle("Cumulative proportion habitat sampled") + xlab("Sampling years") + ylab("Proportion sampled") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
dev.off()

# Look at proportion sampled, using metal tags as total anemones
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "Metal_TA_prop_hab_sampled.pdf"))
ggplot(data = anems_visited_by_year %>% filter(method=="metal tags"), aes(x=site, y=prop_hab_sampled)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  facet_wrap(~year) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of proportion habitat sampled, metal tags as TA")
dev.off()

# Compare different ways of estimating total anems at site
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "Method_comparison_total_anems_at_site.pdf"))
ggplot(data = anems_visited_by_year, aes(site, n_total_anems)) +
  geom_bar(aes(fill = method), position = "dodge", stat = "identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  #theme(legend.position = "bottom") +
  ggtitle("Estimates of total anemones at site")
dev.off()

# Compare different estimates of proportion habitat sampled
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "Method_comparison_prop_hab_sampled_across_sites_and_years.pdf"))
ggplot(data = anems_visited_by_year, aes(x=site, y=prop_hab_sampled, fill=method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  facet_wrap(~year) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of proportion habitat sampled")
dev.off()

# Compare different estimates of proportion habitat sampled, but with edited values (to get rid of NaN, Inf, and > 1)
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "Method_comparison_edited_prop_hab_sampled_across_sites_and_years.pdf"))
ggplot(data = anems_visited_by_year, aes(x=site, y=prop_hab_sampled_tidied, fill=method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  scale_y_continuous(limits = c(0, 1.5)) +
  facet_wrap(~year) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of cleaned-up proportion habitat sampled")
dev.off()

# And looking at each year individually
# 2012
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "2012_prop_hab_sampled.pdf"))
ggplot(data = anems_visited_by_year %>% filter(year == 2012), aes(x=site, y=prop_hab_sampled, fill=method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of proportion habitat sampled - 2012")
dev.off()

# 2013
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "2013_prop_hab_sampled.pdf"))
ggplot(data = anems_visited_by_year %>% filter(year == 2013), aes(x=site, y=prop_hab_sampled, fill=method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of proportion habitat sampled - 2013")
dev.off()

# 2014
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "2014_prop_hab_sampled.pdf"))
ggplot(data = anems_visited_by_year %>% filter(year == 2014), aes(x=site, y=prop_hab_sampled, fill=method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of proportion habitat sampled - 2014")
dev.off()

# 2015
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "2015_prop_hab_sampled.pdf"))
ggplot(data = anems_visited_by_year %>% filter(year == 2015), aes(x=site, y=prop_hab_sampled, fill=method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of proportion habitat sampled - 2015")
dev.off()

# 2016
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "2016_prop_hab_sampled.pdf"))
ggplot(data = anems_visited_by_year %>% filter(year == 2016), aes(x=site, y=prop_hab_sampled, fill=method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of proportion habitat sampled - 2016")
dev.off()

# 2017
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "2017_prop_hab_sampled.pdf"))
ggplot(data = anems_visited_by_year %>% filter(year == 2017), aes(x=site, y=prop_hab_sampled, fill=method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of proportion habitat sampled - 2017")
dev.off()

# 2018
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "2018_prop_hab_sampled.pdf"))
ggplot(data = anems_visited_by_year %>% filter(year == 2018), aes(x=site, y=prop_hab_sampled, fill=method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of proportion habitat sampled - 2018")
dev.off()

# Look at overall proportion sampled, using metal tags as total anemones, by time frame, area vs. anems method
pdf(file = here::here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_through_time_area_vs_anems.pdf"))
ggplot(data = total_sampling_across_years %>% filter(total_anems_method=="metal tags"), aes(x=time_frame, y=total_prop_hab_sampled, fill=total_area_method)) +
  geom_bar(position = "dodge", stat = "identity") +
  #geom_bar(aes(x=time_frame, y=total_prop_hab_sampled_anems), position = "dodge", stat = "identity", fill = "dark blue") +
  geom_hline(yintercept = 1) +
  xlab("sampling years") + ylab("total prop habitat sampled") +
  #facet_wrap(~year) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(text = element_text(size=12)) +
  ggtitle("Estimates of proportion habitat sampled, anems vs area method")
dev.off()

#################### Saving output: ####################
save(anems_visited_by_year, file=here::here("Data/Script_outputs", "anems_visited_by_year.RData"))  # file with total number of anems and prop hab sampled by method
#save(sampled_area_each_year, file=here::here("Data/Script_outputs", "sampled_area_each_year.RData"))
#save(total_area_sampled, file=here::here("Data", "total_area_sampled.RData"))  # file with total area sampled through time by method
save(total_area_sampled_through_time, file=here::here("Data/Script_outputs", "total_area_sampled_through_time.RData"))
save(total_sampling_across_years, file=here::here("Data/Script_outputs", "total_sampling_across_years.RData"))  # summary of prop hab sampled by different total area and total anem methods

# # # Load helpful functions from other scripts
# # script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R", ssl.verifypeer = FALSE)
# # eval(parse(text = script))
# # # script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
# # # eval(parse(text = script))
# # 
# # # Functions from Michelle's GitHub helpers script
# # #helper functions - do various tasks w/database (like assigning dates and site to fish and such)
# # script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/helpers.R", ssl.verifypeer = FALSE)
# # eval(parse(text = script))
# # 
# # # function to attach 2018 anems to anems previously seen (since right now, anem_obs hasn't been extended to all of them), from AnemLocations.R script
# # attach2018anems <- function(anemdf) {
# #   
# #   #filter out anems that might be repeat visits (2018, tag_id less than the new tags for 2018 or new tag and old tag present) and don't already have an anem_obs
# #   anems2018 <- anemdf %>% #646
# #     filter(year == 2018) %>%
# #     filter(anem_id < tag1_2018 | (anem_id >= tag1_2018 & !is.na(old_anem_id))) %>% #filter out the ones that could could have other anem_ids (sighting of existing anem_id tag or new tag with old_anem_id filled in), checked this (below) and think it is filtering out correctly...
# #     filter(is.na(anem_obs)) #some anem_obs already filled in so take those out...
# #   
# #   #other 2018 anems that aren't candidates for repeat obs (so can add back in later)
# #   anems2018_2 <- anemdf %>% #270
# #     filter(year == 2018) %>%
# #     filter(anem_id >= tag1_2018 & is.na(old_anem_id))
# #   
# #   #other 2018 anems that already have an anem_obs (so can add back in later) - checked (once...) that anems2018, anems2018_2, and anems2018_3 covers all 2018 anems
# #   anems2018_3 <- anemdf %>%
# #     filter(year == 2018) %>%
# #     filter(!is.na(anem_obs))
# #     
# #   #filter out anems that are candidates for revisits (just anything from before 2018...), and that will get added back into final df later
# #   otheranems <- anemdf %>% #4068
# #     filter(year != 2018)
# #   
# #   #the filtering above covers all anems except anem_table_id 10473 with anem_id 2535, which has no year associated with it
# #   #test <- anemdf %>% filter(!year %in% c(2012,2013,2014,2015,2016,2017,2018)) #this is how I found that anem
# #   
# #   #go through anems2018 that might be revisits and see if there are anems from previous years that match
# #   for (i in 1:length(anems2018$anem_id)) {
# #     
# #     testid <- anems2018$anem_id[i] #pull out anem_id to compare
# #     testoldid <- anems2018$old_anem_id[i] #pull out old_anem_id
# #     
# #     matchanem <- filter(otheranems, anem_id == testid)  #does the anem_id match an anem_id record from the past?
# #     matcholdanem <- filter(otheranems, anem_id == testoldid) #does the old_anem_id match an anem_id from the past?
# #     
# #     # does the anem_id match an anem_id from the past? 
# #     if (length(matchanem$anem_id) > 0) { #if the anem_id matches an old one
# #       # if so, does the site match?
# #       if (matchanem$site[1] == anems2018$site[i]) { #make sure the site matches
# #          anems2018$anem_id_unq2[i] = matchanem$anem_id_unq[1] #if there are multiple records from the past that match, take the first one
# #       } else {
# #         print(paste("Site does not match for anem_id", testid)) #print message if site doesn't match
# #       }
# #       # if not, does an old_anem_id match an anem_id from the past?
# #     } else if (length(matcholdanem$anem_id) > 0) { #if the old_anem_id matches one from the past
# #       # if so, does the site match?
# #       if (matcholdanem$site[1] == anems2018$site[i]) { #check site
# #         anems2018$anem_id_unq2[i] = matcholdanem$anem_id_unq[1]
# #       } else {
# #         print(paste("Site does not match for old_anem_id", testoldid, "and anem_id", testid))
# #       }
# #     } else {
# #       anems2018$anem_id_unq2[i] = anems2018$anem_id_unq[i]
# #       print(paste("No past anem matches found for testid", testid, "and testoldid", testoldid))
# #     }
# #   }
# #   out <- rbind(anems2018, anems2018_2, anems2018_3, otheranems)
# #   
# #   if(length(out$anem_table_id) == (length(anemdf$anem_table_id)-1)) { #1 anem (anem_table_id 10473, anem_id 2535 has no year listed - investigate further later) (maybe one of those 3 anems I sighted that used to get lost in the filtering b/c not have a year associated with them (NA)?)
# #     print("All anems accounted for.")
# #   } else {
# #     print("Some anems missing or double-counted.")
# #     print(length(out$anem_table_id))
# #     print(length(anemdf$anem_table_id))
# #   }
# #   return(out)
# # }
# 
# # Find the number of tagged anemones visited by each site, each year 
# pull_tagged_anems_by_year <- function(anemsdf, year_i) {
#   
#   # pull all APCL seen
#   allfish_fish <- leyte %>% 
#     tbl("clownfish") %>%
#     select(fish_table_id, anem_table_id, fish_spp, sample_id, gen_id, anem_table_id, recap, tag_id, color, size) %>%
#     collect() %>%
#     filter(fish_spp == "APCL")
#   
#   # and anems associated with that
#   allfish_anems <- leyte %>%
#     tbl("anemones") %>%
#     select(anem_table_id, dive_table_id, anem_obs, obs_time, anem_id, anem_obs, old_anem_id) %>%
#     collect() %>%
#     filter(anem_table_id %in% allfish_fish$anem_table_id)
#   
#   allfish_dives <- leyte %>%
#     tbl("diveinfo") %>%
#     select(dive_table_id, dive_type, date, site, gps) %>%
#     collect() %>%
#     filter(dive_table_id %in% allfish_anems$dive_table_id)
#   
#   all_anems <- left_join(allfish_anems, allfish_dives, by="dive_table_id") 
#   all_anems <- all_anems %>% mutate(year = as.integer(substring(all_anems$date,1,4))) # this line not in the code KC sent me
#   
#   # for each site in the passed dataframe
#   for(i in 1:length(anemsdf$site)) {
#     
#     # pull out all the dives at that site in the selected year that were clownfish dives - none listed for 2012
#     # changed to pull out all dives except those that are R - should check and see how this might affect things...
#     dives <- alldives %>%
#       mutate(month = as.integer(substring(alldives$date,6,7))) %>%
#       filter(site == anemsdf$site[i]) %>%
#       filter(year == year_i) %>%
#       filter(month %in% month_list) %>% #filter out 2015 winter dives (just anem survey, no clownfish caught)
#       #filter(dive_type == "C")
#       filter(dive_type != "R") 
#     
#     # test different ways of calculating anems visited 
#     # pull out anems seem on those dives, using tagged anems (this doesn't find any from 2012 - even though there are some tagged at Visca - dive_table_id = 20)
#     anems_from_dives <- allanems_2018obsfixed %>%
#       filter(dive_table_id %in% dives$dive_table_id) %>% #just the anems on those dives
#       distinct(anem_id_unq2, .keep_all = TRUE) #only keep one row of each anem_id_unq2 (which is obs if had it, id if not, updated for 2018 anems) (could probably just use anem_id here too)
#     
#     anems_total <- sum(count(anems_from_dives, anem_id_unq2)$n)
#     
#     anemsdf$n_tagged_anems[i] = anems_total
#     
#     # now just try anems sampled without specifying that they are tagged (this is how KC was doing it, built off of her emailed code, in Testing_KC_anems_sampled.R)
#     n_anemones_sampled <- all_anems %>%
#       filter(site == anemsdf$site[i]) %>%
#       filter(year == year_i) %>%
#       summarise(n_anem_sampled = n())
#     
#     anemsdf$n_sampled_anems[i] = n_anemones_sampled$n_anem_sampled
#           #anems_total <- as.integer(anems_from_dives %>% summarise(n())) #count up the number, save it as a number rather than a dataframe
#     
#     # Put together "best estimate" of anemones visited at a site in year - tagged anemones seen (n_anems_tagged) for all years except 2012, where use other method (n_anems_seen) because anems were only tagged at Visca in 2012
#     if(year_i == 2012) {
#       anemsdf$n_anems[i] = anemsdf$n_sampled_anems[i]
#     }
#     else {
#       anemsdf$n_anems[i] = anemsdf$n_tagged_anems[i]
#     }
#   }
#   return(anemsdf)
# }
# 
# #################### Running things: ####################
# ##### Pull out info from database
# leyte <- read_db("Leyte")
# 
# # filter out all anems that have an anem_id (which means they were an APCL anem at some point, right?)
# allanems <- leyte %>%
#   tbl("anemones") %>%
#   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
#   collect() %>%
#   filter(!is.na(anem_id))
# 
# allanems <- anem_db %>%
#   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
#   filter(!is.na(anem_id))
# 
# # filter out all dives of any type
# alldives <- leyte %>%
#   tbl("diveinfo") %>%
#   select(dive_table_id, dive_type, date, site, gps) %>%
#   collect() %>%
#   mutate(year = as.integer(substring(date,1,4)))
# 
# # join dives in with anems so can get years associated with each anem (for adding anem_obs to 2018 anems) (some of this code originated in AnemLocations.R)
# allanemswithdives <- left_join(allanems, alldives, by="dive_table_id") %>%
#   #mutate(year = as.integer(substring(date,1,4))) %>%
#   mutate(anem_id = as.numeric(anem_id)) %>% #make anem_id numeric to make future joings/merges easier
#   mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), paste("obs", anem_obs, sep=""))) %>% #add unique anem id so can track anems easier (has obs or id in front of whichever number it is reporting: eg. obs192)
#   mutate(anem_id_unq2 = anem_id_unq) #add another anem_id_unq to use for matching up the anems seen in 2018 (since anem_obs hasn't been run for the new data yet
# 
# # match up anems from 2018 with those seen in previous years - anem_id_unq2 now the column that should identify anems, with same anems having the same string
# allanems_2018obsfixed <- attach2018anems(allanemswithdives) %>% #this will give warning messages and print out a couple of errors - just means one anem is off but not a big deal for now
#   mutate(month = as.integer(substring(date,6,7))) 
# 
# ##### Decide the "total" number of anemones at a site (even though, in reality it changes some across years), using a few methods
# total_anems_by_site <- data.frame(site = site_vec) #initialize a data frame with list of sites
# total_anems_by_site$total_anems <- rep(NA, length(total_anems_by_site$site)) #add column for total anems (Strategy 1)
# total_anems_by_site$total_anems_metal_tags <- rep(NA, length(total_anems_by_site$site)) #column for total anems (Strategy 2: just metal tags)
# total_anems_by_site$total_anems_seen2x <- rep(NA, length(total_anems_by_site$site)) #column for total anems (Strategy 3: just metal tags but only with those seen at least twice)
# total_anems_by_site$total_anems_2015W <- rep(NA, length(total_anems_by_site$site)) #column for total anems (Strategy 4: just 2015 anem survey tags)
# total_anems_by_site$low_est <- rep(NA, length(total_anems_by_site$site)) 
# total_anems_by_site$high_est <- rep(NA, length(total_anems_by_site$site))
# total_anems_by_site$mean_est <- rep(NA, length(total_anems_by_site$site))
# 
# ### Strategy 1: Pull all tagged anems from all dives for each site, then weed out duplicates
# 
# for(i in 1:length(site_vec)) { #go through the sites
#   
#   # pull out all dives at that site
#   dives <- alldives %>%
#     filter(site == site_vec[i])
#   
#   # pull out all anems seen on those dives
#   anems_from_dives <- allanems_2018obsfixed %>%
#     filter(dive_table_id %in% dives$dive_table_id) %>% #just the anems on those dives
#     distinct(anem_id_unq2, .keep_all = TRUE) #only keep one row of each anem_id_unq2 (which is obs if had it, id if not, updated for 2018 anems)
# 
#   anems_total <- as.integer(anems_from_dives %>% summarise(n())) #count up the number, save it as a number rather than a dataframe
#   
#   total_anems_by_site$site[i] = site_vec[i]
#   total_anems_by_site$total_anems[i] = anems_total 
# }
# 
# ### Strategy 2: Pull all metal tags from all dives for each site
# # # Find first metal anem tag (placed in May 2015, according to Malin notes in data doc)
# # metal_tags_2015 <- allanems_2018obsfixed %>%
# #   mutate(month = as.integer(substring(date,6,7))) %>%
# #   filter(year == 2015) %>% 
# #   filter(month == 5) %>%
# #   arrange(anem_id) # based on looking through these, think the first metal tag is 2001 (but there are a lot of other tags noted before and after that one that aren't in the 2000s and don't have old_anem_ids - emailed Michelle/Malin to clarify what tagging that season looked like)
# # 
# # #first_metal_tag <- min(metal_tags_2015$anem_id) #this didn't work b/c there are still pre-2000s tags scene (or placed?) in May 2015
# 
# for(i in 1:length(site_vec)) { #go through the sites
#   metal_anems_at_site <- allanems_2018obsfixed %>%
#     filter(site == site_vec[i]) %>% #filter out just dives from that site
#     filter(anem_id >= first_metal_tag) %>% #filter out just anems with metal tags at that site 
#     distinct(anem_id_unq2, .keep_all = TRUE) #only keep one row of each (if seen multiple times across years)
# 
#   anems_total <- as.integer(metal_anems_at_site %>% summarize(n())) #count up how many there are)
#   
#   total_anems_by_site$total_anems_metal_tags[i] = anems_total
#     
# } #also tried this in a similar for loop way to Strategy 1, where you pull out the relevant dives, then pull out the anems associated with those dives and got the same answer
# 
# ### Strategy 3: Looking at anems with any type of tags that have been seen for at least two years
# 
# for(i in 1:length(site_vec)) { #go through sites
#   tagged_anems_at_site <- allanems_2018obsfixed %>%
#     filter(site == site_vec[i]) %>% #filter out just anems from that site
#     distinct(anem_id_unq2, year, .keep_all = TRUE) %>% #pull out just one obs of each anem in a particular year
#     group_by(anem_id_unq2) %>% #group by anem
#     summarize(n_years = n())
#   
#   anems_total <- as.integer(length((tagged_anems_at_site %>% filter(n_years >= 2))$n_years))
#   
#   total_anems_by_site$total_anems_seen2x[i] = anems_total
# }
# 
# ### Strategy 4: Anems tagged in the anem survey in winter 2015 - but should supplement for a few sites with other years, like Haina and Sitio Baybayon, right?
# 
# for(i in 1:length(site_vec)) {
#   anems_2015W <- allanems_2018obsfixed %>%
#     filter(site == site_vec[i]) %>% 
#     filter(year == 2015, month %in% c(1,2,3)) %>%
#     distinct(anem_id_unq2, .keep_all = TRUE) 
#     
#   anems_total <- as.integer(anems_2015W %>% summarize(n())) 
#   
#   total_anems_by_site$total_anems_2015W[i] <- anems_total
# }
# 
# ### Compare estimates
# for(i in 1:length(total_anems_by_site$site)) {
#   total_anems_by_site$low_est[i] = min(c(total_anems_by_site$total_anems[i], total_anems_by_site$total_anems_metal_tags[i], total_anems_by_site$total_anems_seen2x[i], total_anems_by_site$total_anems_2015W[i]))
#   total_anems_by_site$high_est[i] = max(c(total_anems_by_site$total_anems[i], total_anems_by_site$total_anems_metal_tags[i], total_anems_by_site$total_anems_seen2x[i], total_anems_by_site$total_anems_2015W[i]))
#   total_anems_by_site$mean_est[i] = round(mean(c(total_anems_by_site$total_anems[i], total_anems_by_site$total_anems_metal_tags[i], total_anems_by_site$total_anems_seen2x[i], total_anems_by_site$total_anems_2015W[i]), na.rm = TRUE))
# }
# 
# ##### Find the number of tagged anems visited each year in each site
# #anems_visited_by_site_and_year <- data.frame(site = rep(site_vec, 7)) #initialize a data frame with list of sites for each year
# #anems_visited_by_site_and_year$year <- c(rep(2012, length(site_vec)), rep(2013, length(site_vec)), rep(2014, length(site_vec)), rep(2015, length(site_vec)), rep(2016, length(site_vec)), rep(2017, length(site_vec)), rep(2018, length(site_vec)))
# 
# # Trying by year instead...
# anems_visited_2012 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
# anems_visited_2012$year <- rep(2012, length(site_vec))
# anems_visited_2012$n_tagged_anems <- rep(NA, length(site_vec))
# anems_visited_2012$n_sampled_anems <- rep(NA, length(site_vec))
# anems_visited_2012$n_anems <- rep(NA, length(site_vec))
# anems_visited_2012 <- pull_tagged_anems_by_year(anems_visited_2012, 2012)
# 
# anems_visited_2013 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
# anems_visited_2013$year <- rep(2013, length(site_vec))
# anems_visited_2013$n_tagged_anems <- rep(NA, length(site_vec))
# anems_visited_2013$n_sampled_anems <- rep(NA, length(site_vec))
# anems_visited_2013$n_anems <- rep(NA, length(site_vec))
# anems_visited_2013 <- pull_tagged_anems_by_year(anems_visited_2013, 2013)
# 
# anems_visited_2014 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
# anems_visited_2014$year <- rep(2014, length(site_vec))
# anems_visited_2014$n_tagged_anems <- rep(NA, length(site_vec))
# anems_visited_2014$n_sampled_anems <- rep(NA, length(site_vec))
# anems_visited_2014$n_anems <- rep(NA, length(site_vec))
# anems_visited_2014 <- pull_tagged_anems_by_year(anems_visited_2014, 2014)
# 
# anems_visited_2015 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
# anems_visited_2015$year <- rep(2015, length(site_vec))
# anems_visited_2015$n_tagged_anems <- rep(NA, length(site_vec))
# anems_visited_2015$n_sampled_anems <- rep(NA, length(site_vec))
# anems_visited_2015$n_anems <- rep(NA, length(site_vec))
# anems_visited_2015 <- pull_tagged_anems_by_year(anems_visited_2015, 2015)
# 
# anems_visited_2016 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
# anems_visited_2016$year <- rep(2016, length(site_vec))
# anems_visited_2016$n_tagged_anems <- rep(NA, length(site_vec))
# anems_visited_2016$n_sampled_anems <- rep(NA, length(site_vec))
# anems_visited_2016$n_anems <- rep(NA, length(site_vec))
# anems_visited_2016 <- pull_tagged_anems_by_year(anems_visited_2016, 2016)
# 
# anems_visited_2017 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
# anems_visited_2017$year <- rep(2017, length(site_vec))
# anems_visited_2017$n_tagged_anems <- rep(NA, length(site_vec))
# anems_visited_2017$n_sampled_anems <- rep(NA, length(site_vec))
# anems_visited_2017$n_anems <- rep(NA, length(site_vec))
# anems_visited_2017 <- pull_tagged_anems_by_year(anems_visited_2017, 2017)
# 
# anems_visited_2018 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
# anems_visited_2018$year <- rep(2018, length(site_vec))
# anems_visited_2018$n_tagged_anems <- rep(NA, length(site_vec))
# anems_visited_2018$n_sampled_anems <- rep(NA, length(site_vec))
# anems_visited_2018$n_anems <- rep(NA, length(site_vec))
# anems_visited_2018 <- pull_tagged_anems_by_year(anems_visited_2018, 2018)
# 
# # and put all the years together
# anems_table <- rbind(anems_visited_2012, anems_visited_2013, anems_visited_2014, anems_visited_2015,
#                      anems_visited_2016, anems_visited_2017, anems_visited_2018)
# 
# # add back in total anems per site
# anems_table <- left_join(anems_table, total_anems_by_site, by = "site")
# 
# # and calculate proportion habitat sampled
# anems_table$prop_hab_sampled = anems_table$n_anems/anems_table$total_anems  #this is what it was before (but num might have changed in a couple cases)
# anems_table$prop_hab_sampled_high_TA = anems_table$n_anems/anems_table$high_est
# anems_table$prop_hab_sampled_low_TA = anems_table$n_anems/anems_table$low_est
# anems_table$prop_hab_sampled_mid_TA = anems_table$n_anems/anems_table$mean_est
# anems_table$prop_hab_sampled_metal_TA = anems_table$n_anems/anems_table$total_anems_metal_tags
# 
# ### Re-format data frame for total anem estimates estimates to make plotting easier
# total_anem_methods <- c("all tags", "metal tags", "seen 2x", "2015 W") # list of methods used to estimate total anems at a site
# total_anems_for_plotting <- data_frame(site = rep(site_vec, length(total_anem_methods)))
# total_anems_for_plotting$total_anems <- rep(NA, length(total_anems_for_plotting$site))
# total_anems_for_plotting$method <- rep(NA, length(total_anems_for_plotting$site))
# 
# for(i in 1:length(total_anem_methods)) {
#   start_index = (i-1)*length(site_vec)+1
#   end_index = i*length(site_vec)
#   total_anems_for_plotting$method[start_index:end_index] <- rep(total_anem_methods[i], length(site_vec)) #fill in the method
#   
#   if(i == 1) {
#     total_anems_for_plotting$total_anems[start_index:end_index] = anems_table$total_anems[1:length(site_vec)]
#   } else if (i == 2) {
#     total_anems_for_plotting$total_anems[start_index:end_index] = anems_table$total_anems_metal_tags[1:length(site_vec)]
#   } else if (i == 3) {
#     total_anems_for_plotting$total_anems[start_index:end_index] = anems_table$total_anems_seen2x[1:length(site_vec)]
#   } else if (i == 4) {
#     total_anems_for_plotting$total_anems[start_index:end_index] = anems_table$total_anems_2015W[1:length(site_vec)]
#   }
# }
# 
# ### Re-format data frame for proportion habitat sampled estimates estimates to make plotting easier
# # Pull out and reorganize anems_table by year and prop_hab_sampled method
# prop_hab_for_plotting <- data.frame()
# for(i in 1:length(years_sampled)) {
#   
#   # high TA
#   anems_table_highTA <- anems_table %>% 
#     filter(year == years_sampled[i]) %>% #pull out the right year
#     select(site, year, prop_hab_sampled_high_TA) %>% #and the relevant prop hab sampled method
#     rename(prop_hab_sampled = prop_hab_sampled_high_TA) %>% #rename
#     mutate(method = "high TA") #add method column
#   
#   # low TA
#   anems_table_lowTA <- anems_table %>% 
#     filter(year == years_sampled[i]) %>% #pull out the right year
#     select(site, year, prop_hab_sampled_low_TA) %>% #and the relevant prop hab sampled method
#     rename(prop_hab_sampled = prop_hab_sampled_low_TA) %>% #rename
#     mutate(method = "low TA") #add method column
#   
#   # mid TA
#   anems_table_meanTA <- anems_table %>% 
#     filter(year == years_sampled[i]) %>% #pull out the right year
#     select(site, year, prop_hab_sampled_mid_TA) %>% #and the relevant prop hab sampled method
#     rename(prop_hab_sampled = prop_hab_sampled_mid_TA) %>% #rename
#     mutate(method = "mean TA") #add method column
#   
#   # metal TA
#   anems_table_metalTA <- anems_table %>% 
#     filter(year == years_sampled[i]) %>% #pull out the right year
#     select(site, year, prop_hab_sampled_metal_TA) %>% #and the relevant prop hab sampled method
#     rename(prop_hab_sampled = prop_hab_sampled_metal_TA) %>% #rename
#     mutate(method = "metal TA") #add method column
#   
#   prop_hab_for_plotting <- rbind(prop_hab_for_plotting, anems_table_highTA, anems_table_lowTA, anems_table_meanTA, anems_table_metalTA)
# }
# 
# # Get rid of Inf for sites that have no total anems the way we are counting and make NA (but some years should make 1 - figure out which sites and years that applies to)
# is.na(prop_hab_for_plotting$prop_hab_sampled) <- sapply(prop_hab_for_plotting$prop_hab_sampled, is.infinite)
# 
# ########## KC - if you nuance how to handle Inf (since sometimes should be 0, other times 1), add it here too! #################
# 
# # #################### Plots: ####################
# # Look at proportion sampled, using metal tags as total anemones
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_metalTA.pdf"))
# ggplot(data = anems_table, aes(x=site, y=prop_hab_sampled_metal_TA)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   facet_wrap(~year) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled, metal tags as TA")
# dev.off()
# 
# 
# # Compare different ways of estimating total anems at site
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Total_anems_at_site_methods_comp.pdf"))
# ggplot(data = total_anems_for_plotting, aes(site, total_anems)) +
#   geom_bar(aes(fill = method), position = "dodge", stat = "identity") +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   #theme(legend.position = "bottom") +
#   ggtitle("Estimates of total anemones at site")
# dev.off()
# 
# # Compare different estimates of proportion habitat sampled
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_methods_across_sites_and_years.pdf"))
# ggplot(data = prop_hab_for_plotting, aes(x=site, y=prop_hab_sampled, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   facet_wrap(~year) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled")
# dev.off()
# 
# # Compare different estimates of proportion habitat sampled, but zooming in a bit
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_methods_across_sites_and_years_zoomed.pdf"))
# ggplot(data = prop_hab_for_plotting, aes(x=site, y=prop_hab_sampled, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   scale_y_continuous(limits = c(0, 1.5)) +
#   facet_wrap(~year) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled - zoomed")
# dev.off()
# 
# # And looking at each year individually
# # 2012
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_methods_2012.pdf"))
# ggplot(data = prop_hab_for_plotting %>% filter(year == 2012), aes(x=site, y=prop_hab_sampled, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled - 2012")
# dev.off()
# 
# # 2013
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_methods_2013.pdf"))
# ggplot(data = prop_hab_for_plotting %>% filter(year == 2013), aes(x=site, y=prop_hab_sampled, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled - 2013")
# dev.off()
# 
# # 2014
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_methods_2014.pdf"))
# ggplot(data = prop_hab_for_plotting %>% filter(year == 2014), aes(x=site, y=prop_hab_sampled, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled - 2014")
# dev.off()
# 
# # 2015
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_methods_2015.pdf"))
# ggplot(data = prop_hab_for_plotting %>% filter(year == 2015), aes(x=site, y=prop_hab_sampled, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled - 2015")
# dev.off()
# 
# # 2016
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_methods_2016.pdf"))
# ggplot(data = prop_hab_for_plotting %>% filter(year == 2016), aes(x=site, y=prop_hab_sampled, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled - 2016")
# dev.off()
# 
# # 2017
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_methods_2017.pdf"))
# ggplot(data = prop_hab_for_plotting %>% filter(year == 2017), aes(x=site, y=prop_hab_sampled, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled - 2017")
# dev.off()
# 
# # 2018
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_methods_2018.pdf"))
# ggplot(data = prop_hab_for_plotting %>% filter(year == 2018), aes(x=site, y=prop_hab_sampled, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled - 2018")
# dev.off()
# 
# # all years, metal tags as total anems
# pdf(file = here("Plots/TotalAnemsandPropHabSampled", "Prop_hab_sampled_allyears_metaltagsTA.pdf"))
# ggplot(data = prop_hab_for_plotting %>% filter(year == 2018), aes(x=site, y=prop_hab_sampled, fill=method)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 1) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   theme(text = element_text(size=12)) +
#   ggtitle("Estimates of proportion habitat sampled - 2018")
# dev.off()
# 
# #################### Save output: ####################
# save(total_anems_by_site, file=here("Data", "total_anems_by_site.RData"))
# save(anems_table, file=here("Data", "anem_sampling_table.RData"))
# 
#################### Old code: ####################
# # and reformat proportion hab sampled estimates for easier plotting
# # this is an annoying method that I'm not going to use but was an indexing victory so I'm leaving it here because I don't want to delete it...
# prop_hab_methods <- c("high TA", "low TA", "mean TA", "metal tags TA") # list of methods (really total anem denomintators) used to estimate prop hab sampled
# prop_hab_for_plotting <- data_frame(site = rep(site_vec, length(prop_hab_methods)*length(years_sampled))) #set up data frame
# prop_hab_for_plotting$year <- rep(NA, length(prop_hab_for_plotting$site))
# prop_hab_for_plotting$prop_hab_sampled <- rep(NA, length(prop_hab_for_plotting$site))
# prop_hab_for_plotting$method <- rep(NA, length(prop_hab_for_plotting$site))
# 
# # fill in the year and method sampled
# for(i in 1:length(years_sampled)) {
#   start_index = (i-1)*length(site_vec)*length(prop_hab_methods)+1
#   end_index = i*length(site_vec)*length(prop_hab_methods)
#   prop_hab_for_plotting$year[start_index:end_index] <- rep(years_sampled[i], length(site_vec)*length(prop_hab_methods)) #fill in the year
#   
#   for(j in 1:length(prop_hab_methods)) {
#     start_index_method <- start_index + (j-1)*length(site_vec)  
#     end_index_method <- start_index_method + length(site_vec)-1
#     prop_hab_for_plotting$method[start_index_method:end_index_method] = rep(prop_hab_methods[j], length(site_vec))
#   }
# }

##### MORE OLD CODE
# # Load relevant libraries
# library(RCurl)
# library(dplyr)
# library(lubridate)
# #library(varhandle)
# library(ggplot2)
# library(cowplot)
# library(here)
# 
# month_list <- c(3,4,5,6,7,8) #months when clownfish were caught (so can filter out winter 2015 trip)


