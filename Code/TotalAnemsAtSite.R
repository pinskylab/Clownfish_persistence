# Calculate total number of anemones at each site

#################### Set-up: ####################
# Load relevant libraries
library(RCurl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(here)

tag1 <- 2938 #first new anem tag in 2018 

#################### Functions: ####################
# Load helpful functions from other scripts
script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Functions from Michelle's GitHub helpers script
#helper functions - do various tasks w/database (like assigning dates and site to fish and such)
script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# function to attach 2018 anems to anems previously seen (since right now, anem_obs hasn't been extended to all of them), from AnemLocations.R script
attach2018anems <- function(anemdf) {
  
  #filter out anems that might be repeat visits (2018, tag_id less than the new tags for 2018 or new tag and old tag present) and don't already have an anem_obs
  anems2018 <- anemdf %>% #646
    filter(year == 2018) %>%
    filter(anem_id < tag1 | (anem_id >= tag1 & !is.na(old_anem_id))) %>% #filter out the ones that could could have other anem_ids (sighting of existing anem_id tag or new tag with old_anem_id filled in), checked this (below) and think it is filtering out correctly...
    filter(is.na(anem_obs)) #some anem_obs already filled in so take those out...
  
  #other 2018 anems that aren't candidates for repeat obs (so can add back in later)
  anems2018_2 <- anemdf %>% #270
    filter(year == 2018) %>%
    filter(anem_id >= tag1 & is.na(old_anem_id))
  
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
  collect() 

# join dives in with anems so can get years associated with each anem (for adding anem_obs to 2018 anems) (some of this code originated in AnemLocations.R)
allanemswithdives <- left_join(allanems, alldives, by="dive_table_id") %>%
  mutate(year = as.integer(substring(date,1,4))) %>%
  mutate(anem_id = as.numeric(anem_id)) %>% #make anem_id numeric to make future joings/merges easier
  mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), paste("obs", anem_obs, sep=""))) %>% #add unique anem id so can track anems easier (has obs or id in front of whichever number it is reporting: eg. obs192)
  mutate(anem_id_unq2 = anem_id_unq) #add another anem_id_unq to use for matching up the anems seen in 2018 (since anem_obs hasn't been run for the new data yet

# match up anems from 2018 with those seen in previous years - anem_id_unq2 now the column that should identify anems, with same anems having the same string
allanems_2018obsfixed <- attach2018anems(allanemswithdives) #this will give warning messages and print out a couple of errors - just means one anem is off but not a big deal for now

##### Strategy 1: Pull all anems from all dives for each site, then weed out duplicates
total_anems_by_site <- data.frame(site = site_vec) #initialize a data frame with list of sites
total_anems_by_site$total_anems <- rep(NA, length(total_anems_by_site$site)) #add column for total anems

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

# save the output
save(total_anems_by_site, file=here("Data", "total_anems_by_site.RData"))
