#Survival estimates for clownfish populations
#Last updated: 7/31/17

rm(list=ls())

#################### Set-up: ####################
#Set path for mark executable so that RMark can find it
MarkPath = "~/Box Sync/Rutgers postdoc/Clownfish_persistence/Code"

#Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(RMark)
library(lubridate)
library(dbplyr)

#Set parameters
min_tag_size <- 6.0 #set minimum size for tagging (6.0cm - but this fluctuated in time so some fish were tagged smaller than that)
min_clip_size <- 3.5 #set minimum size for fin-clipping (3.5cm - need to check this with Michelle!)

#################### Functions: ####################
##### Functions from Michelle's GitHub

#helper functions - do various tasks w/database (like assigning dates and site to fish and such)
script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#datefishcap - "which date was a sample_id captured?"
# script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/date_fish_cap.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
datefishcap <- function(x){ #x = sample_id - the id of the fish for which you are trying to get a date
  leyte <- conleyte()
  fish <- leyte %>% 
    tbl("clownfish") %>% 
    select(anem_table_id, sample_id) %>% 
    collect() %>% 
    filter(sample_id %in% x)
  
  anem <- leyte %>% 
    tbl("anemones") %>% 
    select(dive_table_id, anem_table_id) %>% 
    collect() %>% 
    filter(anem_table_id %in% fish$anem_table_id) 
  
  day <- leyte %>% 
    tbl("diveinfo") %>% 
    select(date, dive_table_id) %>% 
    collect() %>% 
    filter(dive_table_id %in% anem$dive_table_id)
  
  day <- left_join(day, anem, by = "dive_table_id")
  day <- left_join(day, fish, by = "anem_table_id") %>% 
    select(sample_id, date)
  
  return(day)
}


##### Functions written by me functions written by me 
#Pull out all fish from clownfish table at a particular site, appends metadata from dive and anemome tables (might want to not do this by site if recaptured fish ever switch sites...)
FishBySite <- function(site.name){ 
  leyte <- read_db("Leyte") 
  
  # select all dives at a given site (from Michelle function site_recap)
  dives <- leyte %>%
    tbl("diveinfo") %>%
    filter(site == site.name) %>%
    collect()
  
  dives.info <- dives %>%
    select(dive_table_id, dive_type, date, site, gps) %>%
    collect()
  
  # select all anemones from those dives (from Michelle function site_recap) 
  anems <- leyte %>%
    tbl("anemones") %>%
    filter (dive_table_id %in% dives$dive_table_id) %>%
    select(anem_table_id, anem_obs, anem_id, old_anem_id, dive_table_id) %>%
    collect()
  
  # add dive_type, date, site, gps to anemones (so can add to fish later)
  anems <- left_join(anems, dives.info, by="dive_table_id")
  
  # select all fish that are on those anemones (from Michelle function site_recap) - make this capture the year too somehow? and choose the option of which kind of recap to get (capid, tagid, or both)
  fish <- leyte %>%
    tbl("clownfish") %>%
    filter(anem_table_id %in% anems$anem_table_id) %>%
    select(sample_id, cap_id, anem_table_id, recap, tag_id, color, size) %>%
    collect()
  
  # add in the dive_type, date, site, gps
  fish <- left_join(fish, anems, by="anem_table_id")
  
  # pull out just the year and put that in a separate column
  year <- as.integer(substring(fish$date,1,4))
  fish$year <- year
  
  return(fish)
}

# For now, encounter histories for fish just based on year (not on multiple sampling events within a year)
CreateYearlyEncounterHist <- function(start.year, end.year, tagged.fish) {
  
  sample.years <- seq(start.year, end.year, 1) #make a vector of years the fish could have been seen
  encounters <- list(); #initialize an empty list to store the various encounter data frames by year
  # encounters <- data.frame()
  
  for (i in 1:length(sample.years)) { #pull out encounter vector by tag for each year, store each as a data frame of tag ids and binary encounters in the encounter list
    
    year.name <- sample.years[i] #get year 
    var.name <- paste("sighted", as.character(sample.years[i]), sep=".") #create dynamic column names for encounters - sighted.[samplingyear]
    
    encounters[[i]] <- tagged.fish %>% group_by(tag_id) %>% summarise(!!var.name := ifelse(sum(year == sample.years[i])>0, 1, 0)) #create data frames for each year that have a vector of tag ids and a vector of encountered (1) or didn't (0) in 
    
  }
  
  return(encounters)
  
}

#Creates the summarized encounter history by tag id: output is a data frame with 2 columns - tag id and summarized encounter history (i.e. 0010)
CreateEncounterSummary <- function(start.year, end.year, tagged.fish) {
  
  sample.years <- seq(start.year, end.year, 1) #make a vector of years the fish could have been seen
  encounters <- list(); #initialize an empty list to store the various encounter data frames by year
  # encounters <- data.frame()
  
  for (i in 1:length(sample.years)) { #pull out encounter vector by tag for each year, store each as a data frame of tag ids and binary encounters in the encounter list
    
    year.name <- sample.years[i] #get year 
    var.name <- paste("sighted", as.character(sample.years[i]), sep=".") #create dynamic column names for encounters - sighted.[samplingyear]
    
    encounters[[i]] <- tagged.fish %>% group_by(tag_id) %>% summarise(!!var.name := ifelse(sum(year == sample.years[i])>0, 1, 0)) #create data frames for each year that have a vector of tag ids and a vector of encountered (1) or didn't (0) in 
    
  }

  encounters.out <- as.data.frame(encounters[[1]]) #seed summary data frame with list of tag ids and encounter in 1st year
  colnames(encounters.out)<- c("tag_id","encounter.hist") #rename the columns so the encounter.hist column can be pasted to iteratively in next step
  
  for (i in 2:length(sample.years)) {
    encounters.out$encounter.hist <- paste(encounters.out$encounter.hist, encounters[[i]][[2]], sep="") #paste on the other encounter 1/0s so get overall encounter histories as strings
  }
  
  return(encounters.out)
}

#################### Running things! ####################
leyte <- read_db("Leyte") #new function call that replaces conleyte

#could be useful for checking that number of tags that come out in encounter histories is the same as the number of unique tags in the input
  # n.tags <- n_distinct(tagged.fish$tag_id)
  
# Generate a cleaned-up encounter history for Palanas
fishPal <- FishBySite('Palanas') #pull out all the fish from Palanas, with meta-data (like date) appended
fishPalTagged <- fishPal %>% filter(!is.na(tag_id)) #just the tagged fish

#create a cleaned up encounter history data frame, where all binary vectors are together and one column has the combined encounter history as a string
encountersPal <- CreateYearlyEncounterHist(2015,2017,fishPalTagged)
encountersPalClean <- encountersPal[[1]]
encountersPalClean$sighted.2016 <- encountersPal[[2]][["sighted.2016"]] #FIGURE OUT HOW TO DO THIS AS A FUNCTION!
encountersPalClean$sighted.2017 <- encountersPal[[3]][["sighted.2017"]]
encountersPalClean$encounter.hist <- paste(encountersPalClean$sighted.2015, encountersPalClean$sighted.2016, encountersPalClean$sighted.2017, sep="") #FIGURE OUT HOW TO DO THIS AS A FUNCTION!!
  
encountersPal2 <- CreateEncounterSummary(2015, 2017, fishPalTagged) #simpler way of getting encounter history like above using function
  
  


#Next up: make encounter histories just for fish tagged in a particular year, say 2014
#Need to think more broadly about how to make encounter histories when have both within and between year resample events (include both in the same?)
#How to do encounter histories when new fish are being tagged at each re-sample event? Do those go in their own or are they just included with preceeding 0s for the sampling events before they were tagged?
  
hist(fishPalTagged$year)
table(fishPalTagged$year)
table(fishPalTagged$tag_id)


#################### Visualizing things! ####################





# Think about some data visualization plots
# number of fish tagged by year by site
# number of recaptured fish caught by year and by site (and proportion too... but of what? total tagged fish in all previous years? just the year before?)
  # select all sample_ids for which there is not a duplicated tag_id
  if(nrow(recap) != 0){
    uni <- fish %>% 
      filter(is.na(cap_id)) %>%  # don't want to remove fish that are capid
      distinct(tag_id) %>% 
      filter(!is.na(tag_id))
    
    # remove all fish that have a distinct tag_id (were not recaptured)
    fish <- anti_join(fish, uni, by = "tag_id")
  }

# select all fish that have a tag_id or capid (because it is a recap, or are genetically recap)
fish <- fish %>% 
  filter(!is.na(tag_id) | !is.na(cap_id))


#Pull out recaptured fish from Palanas
palRecapFish <- site_recap('Palanas')

testPal <- FishBySite('Palanas')



#Put together a table of total captures (of tagging size) and recaptures by site and year
#Pull out relevant parts of clownfish table: just select APCL, only keep relevant columns
fish <- leyte %>% tbl('clownfish') %>% 
  filter(fish_spp == "APCL") %>% 
  select(anem_table_id,size,sample_id,color,cap_id,recap,tag_id) %>%
  collect()


anem.table <- leyte %>% tbl('anemones') %>% collect()
fish.table <- leyte %>% tbl('clownfish') %>% collect()
dive.table <- leyte %>% tbl('diveinfo') %>% collect()

#Attach dates (based on anem, not fish since not all fish have sample ids - taken from Michelle vonbert.R code (https://github.com/stuartmichelle/Growth/blob/master/code/vonbert.R))
date <- dateanem(fish$anem_table_id)
date$date <- as.Date(date$date)
fish <- left_join(fish,date,by='anem_table_id')

#Attach sites (taken from Michelle vonbert.R code)
site <- siteanem(fish$anem_table_id)
site$dive_table_id <- NULL # remove column so it is not duplicated in the join 
fish <- left_join(fish, site, by = "anem_table_id")

#Pull out just the year and put that in a separate column
year <- as.integer(substring(fish$date,1,4))
fish$year <- year

#Pull out just the fish from Palanas
fish_Pal <- fish %>% filter(name == 'Palanas')

#Pull out just the fish larger than 6.0cm - do I need to do this? R
fish_Pal_tagsize <- fish_Pal %>% filter(size >= min_tag_size)

Paltagfish_summary <- fish_Pal %>% group_by(year) %>% count(recap == 'Y')

SLfish <- fish %>% filter(name == 'Sitio Lonas')

fish_summary

hist(fish$size)
hist(fish$year)
hist(fish$recap)

  
  nfish_hatch <- dataf %>% group_by(hatch_ID) %>% count(hatch_ID)

summarise(cip.completed= sum(code == "a"))
##### Stopped editing here
# get a list of recaptured fish ####
fish <- leyte %>% tbl(clownfish) %>% 
  filter(!is.na(capid) | recap == "Y") %>% 
  select(capid, recap, tagid, size, sample_id, anem_table_id) %>% 
  collect()

# get the first capture of tag recaptured fish ####
tags <- leyte %>% tbl("clownfish") %>% 
  filter(tagid %in% fish$tagid) %>% 
  select(capid, recap, tagid, size, sample_id, anem_table_id) %>% 
  collect()

tags$tagid <- as.character(tags$tagid)
fish <- rbind(fish, tags)

# attach dates based on anem not fish because all fish don't have sample_ids ####
date <- dateanem(fish$anem_table_id)
date$date <- as.Date(date$date)
fish <- left_join(fish, date, by = "anem_table_id")




#run some of Michelle's vonbert code in vonbertMichelle (from the Growth folder on GitHub, copied over 7/3/17 from https://github.com/stuartmichelle/Growth/blob/master/code/vonbert.R)
#use the recap data frame that comes out of that

#make a new column with just the year (1901 means no data)
recap$Year <- format(as.Date(recap$date),"%Y")

#start with Palanas - 35 recaps
Palfish <- recap %>% filter(name=="Palanas")
table(recap$name)

Palfish2 <- fish %>% filter(name=="Palanas")
Palfish2$Year <- format(as.Date(Palfish2$date),"%Y")

table(Palfish2$capid)
table(Palfish2$tagid)



stray_hatch <- dataf %>% group_by(hatch_ID) %>% summarize(nstray=sum(stray_relandHbasins))
nfish_hatch <- dataf %>% group_by(hatch_ID) %>% count(hatch_ID)

df

test = dateanem(2069)

leyte %>% tbl("anemones")

#Code from Michelle git hub growth project
source("../code/conleyte.R")

# connect to db
leyte <- conleyte()

# pull in all fish that have been recaptured over multiple years

# genetically 
capid <- leyte %>% tbl("clownfish") %>% select(sample_id, col, size, tagid, capid, fish_table_id, anem_table_id) %>% filter(!is.na(capid)) %>% collect()

#   # by tag id - skipping because so far these are all within-year recaptures
# recap <- leyte %>% tbl("clownfish") %>% filter(!is.na(sample_id) & recap == "Y") %>% collect()

# determine capid range
z <- max(capid$capid)


# for each capture event, assign a year
capid$year <- substr(capid$sample_id, 5,6)
capid$year <- paste("20", capid$year, sep = "")

#################### Old code from some of Michelle's functions that might be useful later ####################
#conleyte function - connect to leyte database (old, use read_db which is in helpers now...)
# script <- getURL("https://raw.githubusercontent.com/stuartmichelle/Phil_code/master/code/conleyte.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
conleyte <- function(){
  suppressMessages(library(dplyr))
  leyte <- src_mysql(dbname = "Leyte", default.file = path.expand("~/myconfig.cnf"), port = 3306, create = F, host = NULL, user = NULL, password = NULL)
  return(leyte)
}

#date_anem function, updated version from helper scripts: "which date was an anemone observed based on clownfish table?"
#script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/anem_date.R", ssl.verifypeer = FALSE)
#eval(parse(text = script))
date_anem <- function(x){
  leyte <- conleyte()
  anem <- leyte %>% 
    tbl("anemones") %>% 
    filter(anem_table_id == x) %>% 
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

#dive_anem: attaches select dive info to anem data
# script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/dive_anem.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
dive_anem <- function(x){ #x = anem for which to get data
  leyte <- conleyte()
  
  # select anemones and fields to pull from db
  if (length(x) == 1){
    anems <- leyte %>% 
      tbl("anemones") %>% 
      filter(anem_id == x) %>%
      select(dive_table_id, anem_id, old_anem_id, anem_obs, anem_spp, anem_dia, notes, obs_time) %>%
      collect()
  }else{
    anems <- leyte %>% 
      tbl("anemones") %>% 
      filter(anem_id %in% x) %>%
      collect()
  }
  
  # retrieve date and site info for those anemones
  dive <- leyte %>% 
    tbl("diveinfo") %>%
    filter(dive_table_id %in% anems$dive_table_id) %>%
    collect()
  
  # join dive info to anem info
  anems <- left_join(anems, dive, by = "dive_table_id")
  rm(dive)
  return(anems)
}

###THESE FUNCTIONS ARE OLD, BEFORE MICHELLE RE-DID THE DATABASE - NOT SURE DUPLICATED IN ONES ABOVE, THOUGH, SO MIGHT NEED TO EDIT AND USE
#siteanem function - "which site was an anemone observed based on clownfish table?", from: https://github.com/stuartmichelle/Phil_code/blob/master/code/siteanem.R
siteanem <- function(x){
  #source("~/Documents/Philippines/Phil_code/conleyte.R")
  #leyte <- conleyte()
  # connect anem ids to dive ids
  anem <- leyte %>% tbl(anemones) %>% filter(anem_table_id %in% x) %>% select(dive_table_id, anem_table_id) %>% collect()
  # get site
  suppressWarnings(dive <- leyte %>% tbl(diveinfo) %>% filter(id %in% anem$dive_table_id) %>% select(id, name) %>% collect())
  site <- left_join(anem, dive, by = c("dive_table_id" = "id"))
  return(site)
}

#sitefish function, from: https://github.com/stuartmichelle/Phil_code/blob/master/code/sitefish.R
sitefish <- function(x){
  #source("~/Documents/Philippines/Phil_code/conleyte.R")
  #leyte <- conleyte()
  # connect anem ids to dive ids
  anem <- leyte %>% tbl(anemones) %>% filter(anem_table_id %in% x) %>% select(dive_table_id, anem_table_id) %>% collect()
  # get site
  suppressWarnings(dive <- leyte %>% tbl(diveinfo) %>% filter(id %in% anem$dive_table_id) %>% select(id, name) %>% collect())
  site <- left_join(anem, dive, by = c("dive_table_id" = "id"))
  return(site)
}

