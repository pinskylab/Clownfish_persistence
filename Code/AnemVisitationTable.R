#Make a table of anemones by current ID, site, and binary visited/not by year
#File started: 10/27/17

rm(list=ls())

#################### Set-up: ####################
#Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(lubridate)
library(dbplyr)
library(ggplot2)

#################### Functions: ####################
##### Functions from Michelle's GitHub
#helper functions - do various tasks w/database (like assigning dates and site to fish and such)
script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#################### Running things! ####################
leyte <- read_db("Leyte") 

#select all the unique anem_ids, to initialize output data frame
anem.IDs <- leyte %>% tbl("anemones") %>% select(anem_id, anem_obs) %>% distinct(anem_id) %>% collect() #pull out list of unique anem_ids
anem.IDs <- anem.IDs %>% filter(!is.na(anem_id)) %>% collect() #remove the NAs
anem.IDs <- anem.IDs %>% filter(anem_id != "") %>% collect() #remove blank anem_ids
anem.IDs <- anem
anem.IDs <- anem.IDs[-c(1), ] #remove the -9999 value (which was in the first row)

#initiatlize output data frame with list of unique anem_ids
anem.Visits <- data.frame(anem.IDs) #initialize output data frame with list of unique anem_ids
anem.Visits$anem_id <- as.numeric(anem.Visits$anem_id)

anemIDSum1 <- sum(anem.Visits$anem_id)

#find first and last year visited
dive.Dates <- leyte %>% tbl("diveinfo") %>% select(date) %>% collect() #pull out list of dive dates
first.Year <- min(as.integer(substring(dive.Dates$date, 1, 4))) #find the earliest dive year (by pulling the year out of the date and finding the minimum)
last.Year <- max(as.integer(substring(dive.Dates$date, 1, 4))) #find the latest dive year 

sample.Years <- seq(first.Year, last.Year, 1) #still not sure how to make this work... Michelle got it to work w/separate, put x and y as inputs for the into input

#create visited/not columns for each sampling year (figure out how to do this in a function or for loop!)
anem.Visits$visited.2012 <- rep(0,length(anem.Visits$anem_id))
anem.Visits$visited.2013 <- rep(0,length(anem.Visits$anem_id))
anem.Visits$visited.2014 <- rep(0,length(anem.Visits$anem_id))
anem.Visits$visited.2015 <- rep(0,length(anem.Visits$anem_id))
anem.Visits$visited.2016 <- rep(0,length(anem.Visits$anem_id))
anem.Visits$visited.2017 <- rep(0,length(anem.Visits$anem_id))
#FIGURE OUT HOW TO DO THE ABOVE IN A FUNCTION OR FOR LOOP!!!
# for (i in 1:length(sample.Years)) {
#   x <- paste("visited", as.character(sample.Years[i]), sep=".")
#   anem.Visits <- mutate(anem.Visits, sample.Years[i] = anem.IDs*0)
# }
  
#pull out the site of each anemone
#anem.Info <- leyte %>% tbl("anemones") %>% select(dive_table_id, anem_table_id, anem_id, anem_spp) %>% collect() #pull out anem_id, dive_table_id, anem_table_id to match up with dive info
anem.Info <- leyte %>% tbl("anemones") %>% select(dive_table_id, anem_table_id, anem_id) %>% collect() #pull out anem_id, dive_table_id, anem_table_id to match up with dive info
dive.Sites <- leyte %>% tbl("diveinfo") %>% select(site, dive_table_id) %>% collect() %>% filter(dive_table_id %in% anem.Info$dive_table_id) #pull out dive site and dive_table_id
anem.Info <- left_join(anem.Info, dive.Sites, by='dive_table_id') %>% filter(!is.na(anem_id)) #match up site to anemones using dive_table_id, filter out resulting data frame for just anemones with anem_ids (that aren't NA)
anem.Info$anem_id <- as.numeric(anem.Info$anem_id)
#keeps <- c("anem_id", "site", "anem_spp") #just select site and anem_id columns for merging w/anem.Visits below
keeps <- c("anem_id", "site") #just select site and anem_id columns for merging w/anem.Visits below - not sure why for some reason does funky things when have species in there - lists some anemones twice even if site and species are the same...
anem.Sites <- distinct(anem.Info[keeps])
anemIDSum2 <- sum(anem.Sites$anem_id)

sleuthing <- anem.Info %>% filter(anem_id == 2485) %>% collect()

n_occur <- data.frame(table(anem.Sites$anem_id)) #make a table of the anem_ids and the frequency with which they occur (to see why anem.Sites has more than anem.IDs)
n_occurMultiple <- n_occur[n_occur$Freq > 1,] #find the ones that occur more than once
sleuthing2 <- anem.Sites %>% filter(anem_id == 791) %>% collect()


anem.Sites[anem.Sites$anem_id > 1,]

n_occur[n_occur$Freq > 1,]

tells you which ids occurred more than once.

vocabulary[vocabulary$id %in% n_occur$Var1[n_occur$Freq > 1],]

#add the site to the main data frame - for some reason does funky things when have species in there - lists some anemones twice even if site and species are the same...
anem.Visits <- left_join(anem.Visits, anem.Sites, by="anem_id")

#now, need to say whether was visited or not each year
anem.Years <- leyte %>% tbl("anemones") %>% group_by(anem_id) %>% select(anem_id, dive_table_id, anem_spp) #this thought/line is incomplete

######### STOPPED EDITING HERE

  if(nrow(recap) != 0){
    uni <- fish %>% 
      filter(is.na(cap_id)) %>%  # don't want to remove fish that are capid
      distinct(tag_id) %>% 
      filter(!is.na(tag_id))
  leyte <- read_db("Leyte")
anem <- leyte %>% 
  tbl("anemones") %>% 
  filter(anem_table_id == anem_tbl_id) %>% 
  select(dive_table_id, anem_table_id) %>% 
  
  leyte <- read_db("Leyte")
anem <- leyte %>% 
  tbl("anemones") %>% 
  filter(anem_table_id == anem_tbl_id) %>% 
  select(dive_table_id, anem_table_id) %>% 
  collect()

day <- leyte %>% 
  tbl("diveinfo") %>%
  select(date, dive_table_id) %>%
  collect() %>% 
  filter(dive_table_id %in% anem$dive_table_id)

day <- left_join(day, anem, by ="dive_table_id")
return(day)
  }

  collect()

day <- leyte %>% 
  tbl("diveinfo") %>%
  select(date, dive_table_id) %>%
  collect() %>% 
  filter(dive_table_id %in% anem$dive_table_id)

day <- left_join(day, anem, by ="dive_table_id")
