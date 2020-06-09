# Calculate total number of anemones at each site and proportion of habitat sampled each year based on anemone numbers

#################### Set-up: ####################

# Pull data, functions, constants 
source(here::here('Code', 'Constants_database_common_functions.R'))

# Load packages
library(stringr)
library(ggplot2)

#################### Functions: ####################
# Find the number of tagged anemones visited by each site, each year 
pull_anems_by_year <- function(anemsdf, allanems, divesdf, year_i, site_i, survey_months, dive_types) {  # anemsdf is to find tagged anemones (use anems_Processed); allanems is to find sampled anemones in 2012, before anems were tagged (use all_APCL_anems - anems associated with APCL, regardless if tagge); divesdf is all dives; year_i is year of interest; site_i is site; survey_months is list of months when clownfish were sampled
  
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
    
    # If it is 2012, use anems sampled without specifying that they are tagged (didn't consider dive types)
    anems_visited = as.integer(sampled_anems_from_dives)
  }

  # Put anem counts together in a list (tagged_anems for tagged ones, sampled_anems for sampled, visited_anems for best est of those visited)
  out <- list(tagged_anems = as.integer(tagged_anems_from_dives), sampled_anems = as.integer(sampled_anems_from_dives), visited_anems = anems_visited)
  
  return(out) 
}

#################### Running things: ####################
##### Pull out anems associated with fish from clownfish table (to use to find visited anems in 2012, before anems were tagged)
# Pull out fish
APCL_caught <- fish_db %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, anem_table_id, recap, tag_id, color, size) %>%
  filter(fish_spp == "APCL")

# and their associated anems 
APCL_anems <- anem_db %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_obs_time, anem_id, anem_obs, old_anem_id) %>%
  filter(anem_table_id %in% APCL_caught$anem_table_id)

# and their associated dives
APCL_dives <- dives_db_processed %>%
  select(dive_table_id, dive_type, date, site, gps, year, month) %>%
  filter(dive_table_id %in% APCL_anems$dive_table_id)

# Join together anem and dive info
all_APCL_anems <- left_join(APCL_anems, APCL_dives, by="dive_table_id")  # this is about 250 fewer than anems_Processed... not sure why...

# Remove intermediate data frames for neatness
rm(APCL_caught, APCL_anems, APCL_dives)

##### Find total habitat at each site, as indicated by anemones, via four methods of estimating total anemones at a site (metal tags, all types of tags, anemones seen in at least two years, anemones seen in the 2015 anemone survey)
# Create data frames to store the output for various methods of assessing total anemones
methods = c("metal tags", "all tags", "seen twice", "2015 survey")

total_anems_by_site_metal <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
  mutate(method = "metal tags",
         n_total_anems = NA)
total_anems_by_site_alltags <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
  mutate(method = "all tags",
         n_total_anems = NA)
total_anems_by_site_seen2x <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
  mutate(method = "seen twice",
         n_total_anems = NA)
total_anems_by_site_2015W <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
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
# Cycle through years and sites to find the number of anems visited
anems_visited_dfs <- list()  # set up list to hold output dataframes

for(i in 1:length(years_sampled)) {  # years_sampled is defined in Constants_database_common_constants.R
  # Set up data frame for this year
  av_df <- data.frame(site = site_vec_order$site_name, stringsAsFactors = FALSE) %>%
    mutate(year = years_sampled[i],
           n_tagged_anems = NA,
           n_sampled_anems = NA,
           n_anems = NA)
  
  # Then cycle through the sites
  for(j in 1:length(site_vec_order$site_name)) {
    av_df$n_tagged_anems[j] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db_processed, years_sampled[i], site_vec_order$site_name[j], spring_months, clown_sample_dive_types))$tagged_anems  # use anems_Processed for anemsdf (finding tagged anems visited) and all_APCL_anems for allanems (finding sampled anems even if not tagged) to capture anems sampled for APCL in 2012, when there were no anemone tags
    av_df$n_sampled_anems[j] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db_processed, years_sampled[i], site_vec_order$site_name[j], spring_months, clown_sample_dive_types))$sampled_anems
    av_df$n_anems[j] = (pull_anems_by_year(anems_Processed, all_APCL_anems, dives_db_processed, years_sampled[i], site_vec_order$site_name[j], spring_months, clown_sample_dive_types))$visited_anems
  }
  anems_visited_dfs[[i]] <- av_df
}

# Bind the years together
anems_visited_by_year <- anems_visited_dfs[[1]]
for(i in 2:length(years_sampled)) {
  anems_visited_by_year <- rbind(anems_visited_by_year, anems_visited_dfs[[i]])
}

# Add in total anems per site
anems_visited_by_year <- left_join(anems_visited_by_year, total_anems_by_site, by = "site")

##### Calculate proportion habitat sampled and create a tidied version where NaN, Inf, and values >1 are edited
anems_visited_by_year <- anems_visited_by_year %>%
  mutate(prop_hab_sampled = n_anems/n_total_anems) %>%  # for all methods of determining total number of tags, calculate prop_hab_sampled
  mutate(n_anems_tidied = case_when(n_anems <= n_total_anems ~ n_anems,
                                   n_anems > n_total_anems ~ n_total_anems)) %>%  # to avoid prop_hab_sampled > 1, restrict n_anems in each year and site to the total at that site
  mutate(prop_hab_sampled_tidied = case_when(n_anems_tidied/n_total_anems <= 1 ~ n_anems_tidied/n_total_anems,  
                                              n_anems_tidied/n_total_anems == Inf | is.na(n_anems_tidied/n_total_anems) ~ 0))  # make a tidier version of prop_hab_sampled (no Inf, NaN, > 1 values), using n_anems_tidied


##### Make appendix table with prop hab sampled by site and year and overall, using metal tags as total anems method
# Except, Caridad Proper, Sitio Tugas, and Sitio Lonas were only sampled in 1-2 years and have no metal-tagged anemones. Add their regular-tagged anemone totals in.
sites_with_metal_tags <- c("Palanas","Wangag","N. Magbangon","S. Magbangon","Cabatoan","Caridad Cemetery","Hicgop South","Elementary School",
                           "San Agustin","Poroc San Flower","Poroc Rose","Visca","Gabas","Tamakin Dacot","Haina","Sitio Baybayon")
sites_without_metal_tags <- c("Caridad Proper","Sitio Lonas","Sitio Tugas")

prop_hab_appendix_table <- rbind(anems_visited_by_year %>%  # metal-tagged anems for sites that have them
                                   filter(method == "metal tags", site %in% sites_with_metal_tags) %>%
                                   select(site, year, n_total_anems, n_anems_tidied, prop_hab_sampled_tidied),
                                 anems_visited_by_year %>%  # all tagged anems for sites without metal tags
                                   filter(method == "all tags", site %in% sites_without_metal_tags) %>%
                                   select(site, year, n_total_anems, n_anems_tidied, prop_hab_sampled_tidied))

# Find overall proportion habitat across all sites sampled in each year to add to table
anems_visited_by_year_all_sites <- prop_hab_appendix_table %>%
  group_by(year) %>%
  summarize(region_total_anems = sum(n_total_anems), region_total_sampled = sum(n_anems_tidied)) %>%
  mutate(prop_hab_sampled_tidied_region = region_total_sampled/region_total_anems) %>%
  ungroup()

# And add back into prop_hab_appendix_table
prop_hab_appendix_table <- rbind(prop_hab_appendix_table,
                                 anems_visited_by_year_all_sites %>%
                                   mutate(site = "overall",
                                          n_total_anems = region_total_anems, 
                                          n_anems_tidied = region_total_sampled,
                                          prop_hab_sampled_tidied = prop_hab_sampled_tidied_region) %>%
                                   select(site, year, n_total_anems, n_anems_tidied, prop_hab_sampled_tidied))

##### Do a similar thing by site, finding cumulative proportion habitat sampled by site
time_frames <- c("2012", "2012-2013", "2012-2014", "2012-2015", "2012-2016", "2012-2017", "2012-2018")

# Make the time frames into a list of the right order and length for the site-level data frame (just one method)
time_frame_list_site <- rep(time_frames[1], length(site_vec_order$site_name))
for(i in 2:length(time_frames)) {
  time_frame_list_site <- c(time_frame_list_site, rep(time_frames[i], length(site_vec_order$site_name)))
}

# Set up data frame, just using metal tags tidied method (not other anem methods or area)
cumulative_prop_hab_sampled_by_site <- data.frame(site = rep(site_vec_order$site_name, length(time_frames)),
                                                  method = "metal tags", stringsAsFactors = FALSE) %>%
  mutate(time_frame = time_frame_list_site,
         end_year = str_sub(time_frame, -4, -1),
         total_possible_sample_anems = NA,  # total number of anems that would have been possible to sample in the time frame (sum of total at this site x number of years of sampling)
         total_anems_sampled = NA,  # sum of total anems sampled at this site in that time period
         total_anems_sampled_tidied = NA,  # sum of total anems sampled at that site in that time period not including anems that exceed the total at the site
         total_prop_hab_sampled_anems = NA,
         total_prop_hab_sampled_anems_tidied = NA)
                                                               
# Find total anems that could have been sampled, adding a sample year each time (to be used by Katrina in parentage analysis work)
for(i in 1:length(cumulative_prop_hab_sampled_by_site$method)) {
  end_year <- as.integer(cumulative_prop_hab_sampled_by_site$end_year[i])
  years_to_include <- seq(2012, end_year, by=1)
  site_i = cumulative_prop_hab_sampled_by_site$site[i]
   
  # Total number of anems that could have been sampled (total at site x number of years) 
  cumulative_prop_hab_sampled_by_site$total_possible_sample_anems[i] = ((total_anems_by_site %>% filter(site == site_i, method == "metal tags"))$n_total_anems)*length(years_to_include)

  # Total number sampled (raw value from this method)
  cumulative_prop_hab_sampled_by_site$total_anems_sampled[i] = sum(anems_visited_by_year %>% 
                                                                  filter(site == site_i, year %in% years_to_include, method ==  "metal tags") %>% 
                                                                  select(n_sampled_anems))
  
  # Total number sampled, tidied (so doesn't exceed total at this site)
  cumulative_prop_hab_sampled_by_site$total_anems_sampled_tidied[i] = sum(anems_visited_by_year %>% 
                                                                     filter(site == site_i, year %in% years_to_include, method ==  "metal tags") %>% 
                                                                     select(n_anems_tidied))
}

# Find cumulative proportion habitat sampled for each site
cumulative_prop_hab_sampled_by_site <- cumulative_prop_hab_sampled_by_site %>%
  mutate(total_prop_hab_sampled_anems = total_anems_sampled/total_possible_sample_anems,
         total_prop_hab_sampled_anems_tidied = total_anems_sampled_tidied/total_possible_sample_anems)

##### Do another version of cumulative proportion habitat sampled through time where Caridad Proper, Sitio Tugas and Sitio Lonas's non-metal-tagged anemones are included
# Set up data frame, just using metal tags tidied method for the sites that have metal tags (most) and all tagged anems for those that don't (Caridad Proper, Sitio Tugas, and Sitio Lonas)
cumulative_prop_hab_sampled_by_site_v2 <- data.frame(site = rep(site_vec_order$site_name, length(time_frames)), stringsAsFactors = FALSE) %>%
  mutate(time_frame = time_frame_list_site,
         end_year = str_sub(time_frame, -4, -1),
         total_possible_sample_anems = NA,  # total number of anems that would have been possible to sample in the time frame (sum of total at this site x number of years of sampling)
         total_anems_sampled_tidied = NA,  # sum of total anems sampled at that site in that time period not including anems that exceed the total at the site
         total_prop_hab_sampled_anems_tidied = NA)

# Find total anems that could have been sampled, adding a sample year each time (to be used by Katrina in parentage analysis work)
for(i in 1:length(cumulative_prop_hab_sampled_by_site_v2$site)) {
  end_year <- as.integer(cumulative_prop_hab_sampled_by_site_v2$end_year[i])
  years_to_include <- seq(2012, end_year, by=1)
  site_i = cumulative_prop_hab_sampled_by_site_v2$site[i]
  
  # Total number of anems that could have been sampled (total at site x number of years) 
  cumulative_prop_hab_sampled_by_site_v2$total_possible_sample_anems[i] = ((prop_hab_appendix_table %>% filter(site == site_i, year == 2012))$n_total_anems)*length(years_to_include)
  
  # Total number sampled, tidied (so doesn't exceed total at this site)
  cumulative_prop_hab_sampled_by_site_v2$total_anems_sampled_tidied[i] = sum(prop_hab_appendix_table %>% 
                                                                            filter(site == site_i, year %in% years_to_include) %>% 
                                                                            select(n_anems_tidied))
}

# Find cumulative proportion habitat sampled for each site
cumulative_prop_hab_sampled_by_site_v2 <- cumulative_prop_hab_sampled_by_site_v2 %>%
  mutate(total_prop_hab_sampled_anems_tidied = total_anems_sampled_tidied/total_possible_sample_anems)

#################### Plots: #################### (already have versions of these in this folder from the other script but now can compare)
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

#################### Saving output: ####################
save(anems_visited_by_year, file=here::here("Data/Script_outputs", "anems_visited_by_year.RData"))  # file with total number of anems and prop hab sampled by method
save(cumulative_prop_hab_sampled_by_site, file=here::here("Data/Script_outputs", "cumulative_prop_hab_sampled_by_site.RData"))  # summary of prop hab sampled cumulatively through time by site, method is "metal tags"
save(cumulative_prop_hab_sampled_by_site_v2, file=here::here("Data/Script_outputs","cumulative_prop_hab_sampled_by_site_v2.RData"))  # same as above but using all tagged anems for Caridad Proper, Sitio Tugas, and Sitio Lonas b/c they don't have metal-tagged anems
save(prop_hab_appendix_table, file=here::here("Data/Script_outputs", "prop_hab_appendix_table.RData"))  # numbers for appendix table with total anems and prop habitat sampled by year and site
