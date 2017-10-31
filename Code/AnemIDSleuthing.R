#Anemones with same anem_id and multiple species
#Portion of code in AmenVisitationTable

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

#make a list of anemones (by anem_id) and their sites
anem.Info <- leyte %>% tbl("anemones") %>% select(dive_table_id, anem_table_id, anem_id, anem_spp) %>% collect() #pull out anem_id, dive_table_id, anem_table_id, anem_spp to match up with dive info
dive.Sites <- leyte %>% tbl("diveinfo") %>% select(site, dive_table_id) %>% collect() %>% filter(dive_table_id %in% anem.Info$dive_table_id) #pull out dive site and dive_table_id
anem.Info <- left_join(anem.Info, dive.Sites, by='dive_table_id') %>% filter(!is.na(anem_id)) #match up site to anemones using dive_table_id, filter out resulting data frame for just anemones with anem_ids (that aren't NA)
anem.Info$anem_id <- as.numeric(anem.Info$anem_id) #convert to numeric for joining with other data frames later
keeps <- c("anem_id", "site", "anem_spp") #just select site and anem_id columns for merging w/anem.Visits below - not sure why for some reason does funky things when have species in there - lists some anemones twice even if site and species are the same...
anem.Sites <- distinct(anem.Info[keeps]) #just pull out the relevant values and only distinct rows
anem.Sites <- anem.Sites %>% filter(anem_id != "-9999") %>% collect() #remove -9999 anem_id

#find the anem_ids that pull up at multiple sites
n_occur <- data.frame(table(anem.Sites$anem_id)) #make a table of the anem_ids and the frequency with which they occur (to see why anem.Sites has more than anem.IDs)
n_occurMultiple <- n_occur[n_occur$Freq > 1,] #find the ones that occur more than once

#check out a couple
sleuthing <- anem.Sites %>% filter(anem_id == 2610) %>% collect()

# 11 (at Visca): STME and ENQD
# some are just not listed: like 2610 (Palanas): HECR and <NA>


