# Calculate total number of anemones at each site

#################### Set-up: ####################
# Load relevant libraries
library(RCurl)
library(dplyr)
library(lubridate)
#library(varhandle)
library(ggplot2)
library(cowplot)
library(here)

month_list <- c(3,4,5,6,7,8) #months when clownfish were caught (so can filter out winter 2015 trip)

#################### Functions: ####################
# Load helpful functions from other scripts
script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
# script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
# eval(parse(text = script))

# Functions from Michelle's GitHub helpers script
#helper functions - do various tasks w/database (like assigning dates and site to fish and such)
script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# function to attach 2018 anems to anems previously seen (since right now, anem_obs hasn't been extended to all of them), from AnemLocations.R script
attach2018anems <- function(anemdf) {
  
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
         anems2018$anem_id_unq2[i] = matchanem$anem_id_unq[1] #if there are multiple records from the past that match, take the first one
      } else {
        print(paste("Site does not match for anem_id", testid)) #print message if site doesn't match
      }
      # if not, does an old_anem_id match an anem_id from the past?
    } else if (length(matcholdanem$anem_id) > 0) { #if the old_anem_id matches one from the past
      # if so, does the site match?
      if (matcholdanem$site[1] == anems2018$site[i]) { #check site
        anems2018$anem_id_unq2[i] = matcholdanem$anem_id_unq[1]
      } else {
        print(paste("Site does not match for old_anem_id", testoldid, "and anem_id", testid))
      }
    } else {
      anems2018$anem_id_unq2[i] = anems2018$anem_id_unq[i]
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

# Find the number of tagged anemones visited by each site, each year 
pull_tagged_anems_by_year <- function(anemsdf, year_i) {
  
  # pull all APCL seen
  allfish_fish <- leyte %>% 
    tbl("clownfish") %>%
    select(fish_table_id, anem_table_id, fish_spp, sample_id, cap_id, anem_table_id, recap, tag_id, color, size) %>%
    collect() %>%
    filter(fish_spp == "APCL")
  
  # and anems associated with that
  allfish_anems <- leyte %>%
    tbl("anemones") %>%
    select(anem_table_id, dive_table_id, anem_obs, obs_time, anem_id, anem_obs, old_anem_id) %>%
    collect() %>%
    filter(anem_table_id %in% allfish_fish$anem_table_id)
  
  allfish_dives <- leyte %>%
    tbl("diveinfo") %>%
    select(dive_table_id, dive_type, date, site, gps) %>%
    collect() %>%
    filter(dive_table_id %in% allfish_anems$dive_table_id)
  
  all_anems <- left_join(allfish_anems, allfish_dives, by="dive_table_id") 
  all_anems <- all_anems %>% mutate(year = as.integer(substring(all_anems$date,1,4))) # this line not in the code KC sent me
  
  # for each site in the passed dataframe
  for(i in 1:length(anemsdf$site)) {
    
    # pull out all the dives at that site in the selected year that were clownfish dives - none listed for 2012
    # changed to pull out all dives except those that are R - should check and see how this might affect things...
    dives <- alldives %>%
      mutate(month = as.integer(substring(alldives$date,6,7))) %>%
      filter(site == anemsdf$site[i]) %>%
      filter(year == year_i) %>%
      filter(month %in% month_list) %>% #filter out 2015 winter dives (just anem survey, no clownfish caught)
      #filter(dive_type == "C")
      filter(dive_type != "R") 
    
    # test different ways of calculating anems visited 
    # pull out anems seem on those dives, using tagged anems (this doesn't find any from 2012 - even though there are some tagged at Visca - dive_table_id = 20)
    anems_from_dives <- allanems_2018obsfixed %>%
      filter(dive_table_id %in% dives$dive_table_id) %>% #just the anems on those dives
      distinct(anem_id_unq2, .keep_all = TRUE) #only keep one row of each anem_id_unq2 (which is obs if had it, id if not, updated for 2018 anems) (could probably just use anem_id here too)
    
    anems_total <- sum(count(anems_from_dives, anem_id_unq2)$n)
    
    anemsdf$n_anems[i] = anems_total
    
    # now just try anems sampled without specifying that they are tagged (this is how KC was doing it, built off of her emailed code, in Testing_KC_anems_sampled.R)
    n_anemones_sampled <- all_anems %>%
      filter(site == anemsdf$site[i]) %>%
      filter(year == year_i) %>%
      summarise(n_anem_sampled = n())
    
    anemsdf$n_anems_2[i] = n_anemones_sampled$n_anem_sampled
          #anems_total <- as.integer(anems_from_dives %>% summarise(n())) #count up the number, save it as a number rather than a dataframe
  }
  
  return(anemsdf)
}

#################### Running things: ####################
##### Pull out info from database
leyte <- read_db("Leyte")

# filter out all anems that have an anem_id (which means they were an APCL anem at some point, right?)
allanems <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
  collect() %>%
  filter(!is.na(anem_id))

# filter out all dives of any type
alldives <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, dive_type, date, site, gps) %>%
  collect() %>%
  mutate(year = as.integer(substring(date,1,4)))

# join dives in with anems so can get years associated with each anem (for adding anem_obs to 2018 anems) (some of this code originated in AnemLocations.R)
allanemswithdives <- left_join(allanems, alldives, by="dive_table_id") %>%
  #mutate(year = as.integer(substring(date,1,4))) %>%
  mutate(anem_id = as.numeric(anem_id)) %>% #make anem_id numeric to make future joings/merges easier
  mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), paste("obs", anem_obs, sep=""))) %>% #add unique anem id so can track anems easier (has obs or id in front of whichever number it is reporting: eg. obs192)
  mutate(anem_id_unq2 = anem_id_unq) #add another anem_id_unq to use for matching up the anems seen in 2018 (since anem_obs hasn't been run for the new data yet

# match up anems from 2018 with those seen in previous years - anem_id_unq2 now the column that should identify anems, with same anems having the same string
allanems_2018obsfixed <- attach2018anems(allanemswithdives) #this will give warning messages and print out a couple of errors - just means one anem is off but not a big deal for now

##### Decide the "total" number of anemones at a site (even though, in reality it changes some across years)
total_anems_by_site <- data.frame(site = site_vec) #initialize a data frame with list of sites
total_anems_by_site$total_anems <- rep(NA, length(total_anems_by_site$site)) #add column for total anems (Strategy 1)
total_anems_by_site$total_anems_metal_tags <- rep(NA, length(total_anems_by_site$site)) #column for total anems (Strategy 2: just metal tags)
total_anems_by_site$total_anems_2015 <- rep(NA, length(total_anems_by_site$site)) #column for total anems (Strategy 3: just 2015 anem survey tags)

### Strategy 1: Pull all tagged anems from all dives for each site, then weed out duplicates

for(i in 1:length(site_vec)) { #go through the sites
  
  # pull out all dives at that site
  dives <- alldives %>%
    filter(site == site_vec[i])
  
  # pull out all anems seen on those dives
  anems_from_dives <- allanems_2018obsfixed %>%
    filter(dive_table_id %in% dives$dive_table_id) %>% #just the anems on those dives
    distinct(anem_id_unq2, .keep_all = TRUE) #only keep one row of each anem_id_unq2 (which is obs if had it, id if not, updated for 2018 anems)

  anems_total <- as.integer(anems_from_dives %>% summarise(n())) #count up the number, save it as a number rather than a dataframe
  
  total_anems_by_site$site[i] = site_vec[i]
  total_anems_by_site$total_anems[i] = anems_total 
}

### Strategy 2: Pull all metal tags from all dives for each site
# Find first metal anem tag (placed in May 2015, according to Malin notes in data doc)
metal_tags_2015 <- allanems_2018obsfixed %>%
  mutate(month = as.integer(substring(date,6,7))) %>%
  filter(year == 2015) %>% 
  filter(month == 5) %>%
  arrange(anem_id) # based on looking through these, think the first metal tag is 2001 (but there are a lot of other tags noted before and after that one that aren't in the 2000s and don't have old_anem_ids - emailed Michelle/Malin to clarify what tagging that season looked like)

#first_metal_tag <- min(metal_tags_2015$anem_id) #this didn't work b/c there are still pre-2000s tags scene (or placed?) in May 2015
  


### Strategy 3: 


##### Find the number of tagged anems visited each year in each site
anems_visited_by_site_and_year <- data.frame(site = rep(site_vec, 7)) #initialize a data frame with list of sites for each year
anems_visited_by_site_and_year$year <- c(rep(2012, length(site_vec)), rep(2013, length(site_vec)), rep(2014, length(site_vec)), rep(2015, length(site_vec)), rep(2016, length(site_vec)), rep(2017, length(site_vec)), rep(2018, length(site_vec)))

# Trying by year instead...
anems_visited_2012 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
anems_visited_2012$year <- rep(2012, length(site_vec))
anems_visited_2012$n_anems <- rep(NA, length(site_vec))
anems_visited_2012$n_anems_2 <- rep(NA, length(site_vec))
anems_visited_2012 <- pull_tagged_anems_by_year(anems_visited_2012, 2012)

anems_visited_2013 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
anems_visited_2013$year <- rep(2013, length(site_vec))
anems_visited_2013$n_anems <- rep(NA, length(site_vec))
anems_visited_2013$n_anems_2 <- rep(NA, length(site_vec))
anems_visited_2013 <- pull_tagged_anems_by_year(anems_visited_2013, 2013)

anems_visited_2014 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
anems_visited_2014$year <- rep(2014, length(site_vec))
anems_visited_2014$n_anems <- rep(NA, length(site_vec))
anems_visited_2014$n_anems_2 <- rep(NA, length(site_vec))
anems_visited_2014 <- pull_tagged_anems_by_year(anems_visited_2014, 2014)

anems_visited_2015 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
anems_visited_2015$year <- rep(2015, length(site_vec))
anems_visited_2015$n_anems <- rep(NA, length(site_vec))
anems_visited_2015$n_anems_2 <- rep(NA, length(site_vec))
anems_visited_2015 <- pull_tagged_anems_by_year(anems_visited_2015, 2015)

anems_visited_2016 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
anems_visited_2016$year <- rep(2016, length(site_vec))
anems_visited_2016$n_anems <- rep(NA, length(site_vec))
anems_visited_2016$n_anems_2 <- rep(NA, length(site_vec))
anems_visited_2016 <- pull_tagged_anems_by_year(anems_visited_2016, 2016)

anems_visited_2017 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
anems_visited_2017$year <- rep(2017, length(site_vec))
anems_visited_2017$n_anems <- rep(NA, length(site_vec))
anems_visited_2017$n_anems_2 <- rep(NA, length(site_vec))
anems_visited_2017 <- pull_tagged_anems_by_year(anems_visited_2017, 2017)

anems_visited_2018 <- data.frame(site = site_vec, stringsAsFactors = FALSE)
anems_visited_2018$year <- rep(2018, length(site_vec))
anems_visited_2018$n_anems <- rep(NA, length(site_vec))
anems_visited_2018$n_anems_2 <- rep(NA, length(site_vec))
anems_visited_2018 <- pull_tagged_anems_by_year(anems_visited_2018, 2018)

# and put all the years together
anems_table <- rbind(anems_visited_2012, anems_visited_2013, anems_visited_2014, anems_visited_2015,
                     anems_visited_2016, anems_visited_2017, anems_visited_2018)

# add back in total anems per site
anems_table <- left_join(anems_table, total_anems_by_site, by = "site")

# and calculate proportion habitat sampled
anems_table$prop_hab_sampled = anems_table$n_anems/anems_table$total_anems
anems_table$prop_hab_sampled_2 = anems_table$n_anems_2/anems_table$total_anems

# #################### Plots: ####################
# ## AAAAAGGGGGHHHHH - Got very frustrated with ggplot and gave up on this but seems like it should be an easy plot to make
# ## Ideally, would like a bar plot for each year, with the two ways of estimating prop hab sampled by site
# ggplot(data = anems_table, aes(site, prop_hab_sampled, year)) +
#   geom_bar(stat="identity") + 
#   facet_wrap(.~year)

#################### Save output: ####################
save(total_anems_by_site, file=here("Data", "total_anems_by_site.RData"))
save(anems_table, file=here("Data", "anem_sampling_table.RData"))
