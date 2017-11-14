#Make a table of anemones by current ID, site, and binary visited/not by year
#File started: 10/27/17

rm(list=ls())

#################### Set-up: ####################
setwd("~/Box Sync/Rutgers postdoc/Clownfish_persistence/Code") #set working directory

#Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(lubridate)
library(dbplyr)
library(ggplot2)
library(cowplot)

#################### Functions: ####################
##### Functions from Michelle's GitHub
#helper functions - do various tasks w/database (like assigning dates and site to fish and such)
script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#function to find the lat and long for an anem_id (based off of Michelle's sample_latlon function)
anemid_latlong <- function(anem.table.id, latlondata) { #anem.table.id is one anem_table_id value, latlondata is table of GPX data from database (rather than making the function call it each time); will need to think a bit more clearly about how to handle different locations read for different visits to the same anem_id (or different with same anem_obs); for now, just letting every row in anem.Info get a lat-long
  
  #find the dive info and time for this anem observation
  dive <- leyte %>%
    tbl("anemones") %>%
    select(anem_table_id, obs_time, dive_table_id, anem_id) %>%
    collect() %>%
    filter(anem_table_id %in% anem.table.id)
  
  # find the date info and gps unit for this anem observation
  date <- leyte %>% 
    tbl("diveinfo") %>% 
    select(dive_table_id, date, gps, site) %>% 
    collect() %>% 
    filter(dive_table_id %in% dive$dive_table_id)
  
  #join with anem info, format obs time
  anem <- left_join(dive, date, by = "dive_table_id") %>% 
    separate(obs_time, into = c("hour", "minute", "second"), sep = ":") %>% #this line and the next directly from Michelle's code
    mutate(gpx_hour = as.numeric(hour) - 8)
  
  # find the lat long for this anem observation
  latloninfo <- latlondata %>%
    filter(date %in% anem$date) %>% 
    separate(time, into = c("hour", "minute", "second"), sep = ":") %>% 
    filter(as.numeric(hour) == anem$gpx_hour & as.numeric(minute) == anem$minute) 
  
  latloninfo$lat <- as.numeric(latloninfo$lat)
  latloninfo$lon <- as.numeric(latloninfo$lon)
  
  #often get multiple records for each anem_table_id (like if sit there for a while) - so multiple GPS points for same visit to an anemone, not differences across visits
  dups_lat <- which(duplicated(latloninfo$lat)) #vector of positions of duplicate values
  dups_lon <- which(duplicated(latloninfo$lon))
  
  #either take the mean of the lat/lon readings or the duplicated values, depending if there are duplicate points
  if(length(dups_lat) == 0) { #if all latitude points are different
    anem$lat <- round(mean(latloninfo$lat), digits = 5) #take the mean of the latitude values (digits = 5 b/c that is what Michelle had)
    anem$lon <- round(mean(latloninfo$lon), digits = 5) #take the mean of the longitude values
    #lat <- round(mean(latloninfo$lat), digits = 5) 
    #lon <- round(mean(latloninfo$lon), digits = 5)
  }else{
    anem$lat <- latloninfo$lat[dups_lat[1]] #if are duplicates, take the value of the first duplicated point
    anem$lon <- latloninfo$lon[dups_lon[1]]
    #lat <- latloninfo$lat[dups_lat[1]] #if are duplicates, take the value of the first duplicated point
    #lon <- latloninfo$lon[dups_lon[1]]
  }
  
  return(anem)
  
}

#################### Running things! ####################
leyte <- read_db("Leyte") 

#pull out all lat/lon info so can feed into function above (rather than having to pull out from database each time)
alllatlons <- leyte %>%
  tbl("GPX") %>%
  select(lat, lon, time, unit) %>%
  collect(n = Inf) %>% 
  separate(time, into = c("date", "time"), sep = " ") 

#select all the unique anem_ids to initialize output data frame, also grab anem_obs if have it, should double-check that this is gettin the same list as before....
anem.IDs <- leyte %>% tbl("anemones") %>% select(anem_id, anem_obs) %>% collect()
anem.IDs <- anem.IDs %>% filter(!is.na(anem_id)) %>% collect() #remove the NAs
anem.IDs <- anem.IDs %>% filter(anem_id != "") %>% collect() #remove blank anem_ids
anem.IDs <- anem.IDs %>% filter(anem_id != "-9999") %>% collect() #remove -9999 anem_id
anem.IDs <- distinct(anem.IDs, anem_id, .keep_all = TRUE) #just keep distinct observations of anem_id (so don't have multiple observations listed for each)

#initiatlize output data frame with list of unique anem_ids
anem.Obs <- data.frame(anem.IDs) #initialize data frame with list of unique anem_ids and their associated anem_obs
anem.Obs$anem_id <- as.numeric(anem.Obs$anem_id) #make anem_ids numeric for joining in data frame later

#pull out relevant info from dive table and anemone table and merge so can match each anem_id to site and dates visited 
dive.Info <- leyte %>% tbl("diveinfo") %>% select(date, dive_table_id, dive_type, site) %>% collect() #pull out list of dive dates, sites, and types
dive.Info$year <- as.integer(substring(dive.Info$date, 1, 4)) #add a year column
anem.Info <- leyte %>% tbl("anemones") %>% select(dive_table_id, anem_table_id, anem_id, anem_obs) %>% collect() #pull out anem_id, dive_table_id, anem_table_id to match up with dive info

#merge into one dataframe by dive_table_id (to assign date, year, site, dive type to each anem), filter out anem_ids that are NA
anem.AllInfo <- left_join(anem.Info, dive.Info, by="dive_table_id") %>% filter(!is.na(anem_id)) #9853 rows when NAs left in, 4056 when filtered out
anem.AllInfo$anem_id <- as.numeric(anem.AllInfo$anem_id) #make anem_id numeric to make future joins/merges easier
anem.AllInfo <- anem.AllInfo %>% filter(anem_id != "-9999") %>% collect() #remove -9999 anem_id
anem.AllInfo <- anem.AllInfo %>% filter(anem_id != "") %>% collect() #remove blank anem_id

#####  Data frame and analysis treating all anem_ids as separate anemones, not filtering visits by season/field year or dive type at all
anem.Visits <- data.frame(anem.IDs) #initialize output data frame with list of unique anem_ids
anem.Visits$anem_id <- as.numeric(anem.Visits$anem_id) #make anem_ids numeric for joining in data frame later

#add site into main dataframe
keeps <- c("anem_id", "site") #just select site and anem_id columns for merging w/anem.Visits below
anem.Sites <- distinct(anem.AllInfo[keeps])
anem.Visits <- left_join(anem.Visits, anem.Sites, by="anem_id") #add sites into data frame with list of anem_ids

#now, need to say whether anemone was visited each year or not - these create dataframes with a list of anem_ids and 1/0 for that year
#surely there is a way to do this in a function.... will think about this - might be able to do with encounter hist code if join with anem.Vists after each year
anems.2012 <- anem.AllInfo %>% group_by(anem_id) %>% summarize(visited.2012 := ifelse(sum(year == 2012)>0, 1, 0))
anems.2013 <- anem.AllInfo %>% group_by(anem_id) %>% summarize(visited.2013 := ifelse(sum(year == 2013)>0, 1, 0))
anems.2014 <- anem.AllInfo %>% group_by(anem_id) %>% summarize(visited.2014 := ifelse(sum(year == 2014)>0, 1, 0))
anems.2015 <- anem.AllInfo %>% group_by(anem_id) %>% summarize(visited.2015 := ifelse(sum(year == 2015)>0, 1, 0))
anems.2016 <- anem.AllInfo %>% group_by(anem_id) %>% summarize(visited.2016 := ifelse(sum(year == 2016)>0, 1, 0))
anems.2017 <- anem.AllInfo %>% group_by(anem_id) %>% summarize(visited.2017 := ifelse(sum(year == 2017)>0, 1, 0))

#join with anem.Visits
anem.Visits <- left_join(anem.Visits, anems.2012, by="anem_id")
anem.Visits <- left_join(anem.Visits, anems.2013, by="anem_id")
anem.Visits <- left_join(anem.Visits, anems.2014, by="anem_id")
anem.Visits <- left_join(anem.Visits, anems.2015, by="anem_id")
anem.Visits <- left_join(anem.Visits, anems.2016, by="anem_id")
anem.Visits <- left_join(anem.Visits, anems.2017, by="anem_id")

#could plot this in some useful way....
#here, number of anemones visited each year by site
anem.Visits.Sum <- anem.Visits %>% group_by(site) %>% summarize(n.2012 = sum(visited.2012), n.2013 = sum(visited.2013), 
                                                                n.2014 = sum(visited.2014), n.2015 = sum(visited.2015), 
                                                                n.2016 = sum(visited.2016), n.2017 = sum(visited.2017))

##### Grouping anem_ids by anem_obs values (so anemones we already know from tag surveys are actually the same one)
anem.AllInfo <- anem.AllInfo %>% mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), anem_obs)) #add unique anem id to anem.AllInfo

anem.Obs <- left_join(anem.Obs, anem.Sites, by="anem_id") #add sites into data frame with list of anem_ids + anem_obs (same anem.Sites as above with anem.Visits work)
anem.Obs <- anem.Obs %>% mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), anem_obs)) #add anem_id_unq to anem.Obs too (could probably streamline this whole process)


#same as above, need to say whether anemone was visited each year or not - these create dataframes with a list of anem_id_unq and 1/0 for that year
anems.2012.unqID <- anem.AllInfo %>% group_by(anem_id_unq) %>% summarize(visited.2012 := ifelse(sum(year == 2012)>0, 1, 0))
anems.2013.unqID <- anem.AllInfo %>% group_by(anem_id_unq) %>% summarize(visited.2013 := ifelse(sum(year == 2013)>0, 1, 0))
anems.2014.unqID <- anem.AllInfo %>% group_by(anem_id_unq) %>% summarize(visited.2014 := ifelse(sum(year == 2014)>0, 1, 0))
anems.2015.unqID <- anem.AllInfo %>% group_by(anem_id_unq) %>% summarize(visited.2015 := ifelse(sum(year == 2015)>0, 1, 0))
anems.2016.unqID <- anem.AllInfo %>% group_by(anem_id_unq) %>% summarize(visited.2016 := ifelse(sum(year == 2016)>0, 1, 0))
anems.2017.unqID <- anem.AllInfo %>% group_by(anem_id_unq) %>% summarize(visited.2017 := ifelse(sum(year == 2017)>0, 1, 0))

#join anem visted/not data frames from above with anem.Obs
anem.Obs <- left_join(anem.Obs, anems.2012.unqID, by="anem_id_unq")
anem.Obs <- left_join(anem.Obs, anems.2013.unqID, by="anem_id_unq")
anem.Obs <- left_join(anem.Obs, anems.2014.unqID, by="anem_id_unq")
anem.Obs <- left_join(anem.Obs, anems.2015.unqID, by="anem_id_unq")
anem.Obs <- left_join(anem.Obs, anems.2016.unqID, by="anem_id_unq")
anem.Obs <- left_join(anem.Obs, anems.2017.unqID, by="anem_id_unq")

#here, number of anemones visited each year by site
anem.Obs.Sum <- anem.Obs %>% group_by(site) %>% summarize(n.2012 = sum(visited.2012), n.2013 = sum(visited.2013), 
                                                                n.2014 = sum(visited.2014), n.2015 = sum(visited.2015), 
                                                                n.2016 = sum(visited.2016), n.2017 = sum(visited.2017))

#re-configure the data frame so it's easier to plot each site as a different subplot, with years on x-axis and anems visited on y
site_list <- anem.Obs.Sum$site #pull out sites

#make a list with no spaces in the site names (plays better with labeling in ggplot facet)
site_list_no_spaces <- site_list
site_list_no_spaces[2] = "Caridad~Cemetery"
site_list_no_spaces[3] = "Caridad~Proper"
site_list_no_spaces[4] = "Elementary~School"
site_list_no_spaces[7] = "Hicgop~South"
site_list_no_spaces[10] = "Poroc~Rose"
site_list_no_spaces[11] = "Poroc~San~Flower"
site_list_no_spaces[12] = "San~Agustin"
site_list_no_spaces[13] = "Sitio~Baybayon"
site_list_no_spaces[14] = "Sitio~Lonas"
site_list_no_spaces[15] = "Sitio~Tugas"
site_list_no_spaces[16] = "Tamakin~Dacot"

nsites <- length(site_list) #count number of sites
start.Year <- 2012 #first year of sampling
end.Year <- 2017 #last year of sampling
nyears <- end.Year - start.Year + 1 #number of years sampled
site_vec <- rep(site_list_no_spaces, nyears) #repeat sites so can have a data frame with sites, year, nsampled as the columns so easier to plot grouped by site or year in ggplot
year_vec <- c(rep(2012, nsites), rep(2013, nsites), rep(2014, nsites), rep(2015, nsites), rep(2016, nsites), rep(2017, nsites))

visits_vec <- c(anem.Obs.Sum$n.2012, anem.Obs.Sum$n.2013, anem.Obs.Sum$n.2014, anem.Obs.Sum$n.2015, anem.Obs.Sum$n.2016, anem.Obs.Sum$n.2017)

anem.Obs.SumSite <- as.data.frame(site_vec)
anem.Obs.SumSite <- rename(anem.Obs.SumSite, site = site_vec)
anem.Obs.SumSite$year <- year_vec
anem.Obs.SumSite$visits <- visits_vec

# #make a unique id for each anemone - anem_obs if it has it, anem_id w/"id" at front if not (b/c set of numbers used for anem_ids and anem_obs overlap)
# anem.Obs <- anem.Obs %>% mutate(anem_id_unq = ifelse(is.na(anem_obs), paste("id", anem_id, sep=""), anem_obs))
# #group by unique anem id (anem_id_unq), combine visitation histories 
# test3 <- anem.Obs %>% group_by(anem_id_unq)
# test <- anem.AllInfo %>% mutate(anem_id_unq = paste("id", anem_id, sep=""), anem_obs_unq = paste("ob", anem_obs, sep=""))

##### Grouping anem_ids by proximity 
# find the lat/lon for each anem observation (could probably do this without a for loop but then the averaging across multiple gps records per observation got hairy when trying to do it all w/dplyr as vectors...)
anem.AllInfo$lat <- rep(NA, length(anem.AllInfo$anem_table_id))
anem.AllInfo$lon <- rep(NA, length(anem.AllInfo$anem_table_id))

#go through anem_table_ids and find the lat/lon for each anem observation
for (i in 1:length(anem.AllInfo$anem_table_id)) {
  outlatlon <- anemid_latlong(anem.AllInfo$anem_table_id[i], alllatlons) 
  anem.AllInfo$lat[i] <- outlatlon$lat
  anem.AllInfo$lon[i] <- outlatlon$lon
  print(i)
}

#save output, since takes several minutes to run the above for loop
save(anem.AllInfo, file='AnemAllInfowLatLon_14Nov2017.RData')

#find how many times each anem marked w/anem_obs is visited (ultimate goal is here to find the distances between repeat observations of each anemone)
anem.repeatobs <- anem.AllInfo %>% filter(!is.na(anem_obs)) %>% group_by(anem_obs) %>% summarize(nvisits = n()) #just using anem_obs (b/c if not have anem_obs, not visited multiple times by this definition)

#just a visualization of how many times anems are visited (using anem_obs as indicator of revisitation)
hist(anem.repeatobs$nvisits)

#grouping by anem_id_unq (so anem_obs if has one, otherwise anem_id w/"id" in front), finding mean and sd of lat and lon measurements and nvisits
#anem.LatLon <- anem.AllInfo %>% group_by(anem_id_unq) %>% summarize(meanlat = mean(lat, na.rm = TRUE), sdlat = sd(lat, na.rm = TRUE), meanlon = mean(lon, na.rm = TRUE), sdlon = sd(lon, na.rm = TRUE), nvisits = n())
#not all visits have a GPS measurement so find the number of visits that do...
anem.LatLon <- anem.AllInfo %>% filter(!is.na(lat)) %>% group_by(anem_id_unq) %>% summarize(meanlat = mean(lat), sdlat = sd(lat), meanlon = mean(lon, na.rm = TRUE), sdlon = sd(lon, na.rm = TRUE), ngps = n())
anem.LatLon <- left_join(anem.LatLon, select(anem.AllInfo, anem_id_unq, site), by = 'anem_id_unq') #add in site

#testAnem <- anem.AllInfo %>% filter(anem_id_unq == 103) %>% collect() #for looking at particular case studies of examples...

#################### Plots ####################

##### Summary of anemone visits by site and year (counting each anem_id separately)
#2012
#pdf(file='AnemVisits_2012_31Oct2017.pdf')
hist2012 <- ggplot(data=anem.Visits.Sum, aes(site, n.2012)) +   
  geom_bar(position ="dodge", stat="identity") +
  xlab("site") + ylab("tagged anems visited") + ggtitle("2012") +
  theme_bw() +
  theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle =  90, vjust = 1, hjust = 1, size=8))
#dev.off()

#2013
#pdf(file='AnemVisits_2013_31Oct2017.pdf')
hist2013 <- ggplot(data=anem.Visits.Sum, aes(site, n.2013)) +   
  geom_bar(position ="dodge", stat="identity") +
  xlab("site") + ylab("tagged anems visited") + ggtitle("2013") +
  theme_bw() +
  theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle =  90, vjust = 1, hjust = 1, size=8))
#dev.off()

#2014
#pdf(file='AnemVisits_2014_31Oct2017.pdf')
hist2014 <- ggplot(data=anem.Visits.Sum, aes(site, n.2014)) +   
  geom_bar(position ="dodge", stat="identity") +
  xlab("site") + ylab("tagged anems visited") + ggtitle("2014") +
  theme_bw() +
  theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle =  90, vjust = 1, hjust = 1, size=8))
#dev.off()

#2015
#pdf(file='AnemVisits_2015_31Oct2017.pdf')
hist2015 <- ggplot(data=anem.Visits.Sum, aes(site, n.2015)) +   
  geom_bar(position ="dodge", stat="identity") +
  xlab("site") + ylab("tagged anems visited") + ggtitle("2015") +
  theme_bw() +
  theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle =  90, vjust = 1, hjust = 1, size=8))
#dev.off()

#2016
#pdf(file='AnemVisits_2016_31Oct2017.pdf')
hist2016 <- ggplot(data=anem.Visits.Sum, aes(site, n.2016)) +   
  geom_bar(position ="dodge", stat="identity") +
  xlab("site") + ylab("tagged anems visited") + ggtitle("2016") +
  theme_bw() +
  theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size=8))
#dev.off()

#2017
#pdf(file='AnemVisits_2017_31Oct2017.pdf')
hist2017 <- ggplot(data=anem.Visits.Sum, aes(site, n.2017)) +   
  geom_bar(position ="dodge", stat="identity") +
  xlab("site") + ylab("tagged anems visited") + ggtitle("2017") +
  theme_bw() +
  theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle =  90, vjust = 1, hjust = 1, size=8))
#dev.off()

#arrange in a grid, using cowplot (not so much b/c these plots are useful, more to test cowplot)
theme_set(theme_cowplot(font_size=12)) # reduce default font size
#pdf(file='AnemVisits_byYear_usinganem_id_03Nov2017.pdf')
#pdf(file='AnemVisitsbyYearandSite_Nov2017.pdf')
allsubplots <- plot_grid(hist2012, hist2013, hist2014, hist2015, hist2016, hist2017)
#dev.off()

save_plot("AnemVisitsbyYearandSite.pdf", allsubplots,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 3, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
)

##### Summary of anemone visits by site and year (using anem_obs to match anem_ids known to be the same)
pdf(file="AnemObs_VisitsbyYearSite.pdf", width=10, height=3)
ggplot(data=anem.Obs.SumSite, aes(year, visits)) +   
  geom_bar(position ="dodge", stat="identity") +
  facet_grid(.~site_wrap, labeller=label_parsed) +
  xlab("year") + ylab("tagged anems visited") +
  theme_bw() +
  theme(text = element_text(size=8)) +
  theme(axis.text.x = element_text(angle =  90, vjust = 1, hjust = 1, size=8))
dev.off()

##### Visualizing the standard deviation among estimates of anemone positions (just using anem_obs and anem_id as an estimate of repeat anem visits)
pdf(file="AnemGPS_Estimates.pdf")
ggplot(data=anem.LatLon, aes(sdlat, sdlon)) +
  geom_point(aes(color=site, size=ngps)) +
  #facet_grid(.~site_wrap, labeller=label_parsed) +
  xlab("latitude standard deviation") + ylab("longitude standard deviation") + ggtitle("Standard deviation of anem positions across visits") +
  theme_bw()
dev.off()

#for the way the gps measures things, 5th decimal point is meters (4th is 10m, 3rd is 100m, 2nd is 1000m, 1st is 10000m) (talked to Michelle 11/14/17)

# #attempts to get the site names to wrap from here: https://stackoverflow.com/questions/37174316/how-to-fit-long-text-into-ggplot2-facet-titles, didn't work...
# swr = function(string, nwrap=15) {
#   paste(strwrap(string, width=nwrap), collapse="\n")
# }
# swr = Vectorize(swr)
# 
# anem.Obs.SumSite$site_wrap <- swr(anem.Obs.SumSite$site)
# 
# #another idea from here: https://stackoverflow.com/questions/16654691/how-to-dynamically-wrap-facet-label-using-ggplot2
# facet_wrap(~groupwrap, labeller = labeller(groupwrap = label_wrap_gen(10)))





#####

#what about by dive type? could filter first by dive type, then do the same thing as above

#from talking to Michelle 10/31/17 at end of day
#in 2013, 2014, summer 2015, 2016 -> pretty much the same methodology, C dives are what A+C were in 2017
#just anemone surveys in winter 2016 (Jan-March), no clownfish caught - should filter those dates out for clownfish recapture purposes
#2017, would probably want to include both A and C (b/c if Katrina saw an anemone but Gerry didn't catch any fish, that would show up as A and not C)

#many of the anems with different anem_ids are probably actually the same anemone (like if had an old tag that we didn't see and got replaced) - need to connect through anem_obs but kind of messed up right now since a bunch of anem_ids changed yesterday to fix the multiple site issues - talk to Michelle more about how to deal with/use

# ##### LESS PRETTY, UNSORTED, OLD CODE BELOW ############
# #has some useful stuff for determining if anem_id maps to multiple sites or species
# # #Check to see if any anem_ids show up at multiple sites
# # n_occur <- data.frame(table(anem.Sites$anem_id)) #make a table of the anem_ids and the frequency with which they occur (to see why anem.Sites has more than anem.IDs)
# # n_occurMultiple <- n_occur[n_occur$Freq > 1,] #find the ones that occur more than once
# # sleuthing2 <- anem.Sites %>% filter(anem_id == 791) %>% collect()
# 
# #check to see if any anem_ids show up as multiple species - yes, so will add species back in later at some point (lines to do so commented out above)
# n_occur <- data.frame(table(anem.Sites$anem_id)) #make a table of the anem_ids and the frequency with which they occur (to see why anem.Sites has more than anem.IDs)
# n_occurMultiple <- n_occur[n_occur$Freq > 1,] #find the ones that occur more than once
# sleuthing2 <- anem.Sites %>% filter(anem_id == "1061") %>% collect()
# 
# 
# #has code that could be useful for making a function to do the mapping by year
# 
# #now, need to say whether was visited or not each year
# for (i in 1:length(sample.Years)) {
#   year.name <- sample.Years[i]
#   var.name <- paste("visited", as.character(sample.Years[i]), sep=".")
#   encounters[[i]] <- anem.tags %>% group_by(tag_id) %>% summarise(!!var.name := ifelse(sum(year == sample.years[i])>0, 1, 0)) #create data frames for each year that have a vector of tag ids and a vector of encountered (1) or didn't (0) in 
#   
#   
# 
# 
# dive.Dates <- leyte %>% tbl("diveinfo") %>% select(date, dive_table_id, dive_type) %>% collect() #pull out list of dive dates and types
# dive.Dates$year <- as.integer(substring(dive.Dates$date, 1, 4)) #add a year column
# #anem.Info <- leyte %>% tbl("anemones") %>% select(dive_table_id, anem_table_id, anem_id, anem_spp) %>% collect() #pull out anem_id, dive_table_id, anem_table_id to match up with dive info
# anem.Info <- leyte %>% tbl("anemones") %>% select(dive_table_id, anem_table_id, anem_id) %>% collect() #pull out anem_id, dive_table_id, anem_table_id to match up with dive info
# dive.Sites <- leyte %>% tbl("diveinfo") %>% select(site, dive_table_id) %>% collect() %>% filter(dive_table_id %in% anem.Info$dive_table_id) #pull out dive site and dive_table_id
# anem.Info <- left_join(anem.Info, dive.Sites, by='dive_table_id') %>% filter(!is.na(anem_id)) #match up site to anemones using dive_table_id, filter out resulting data frame for just anemones with anem_ids (that aren't NA)
# anem.Info$anem_id <- as.numeric(anem.Info$anem_id)
# #keeps <- c("anem_id", "site", "anem_spp") #just select site and anem_id columns for merging w/anem.Visits below
# keeps <- c("anem_id", "site") #just select site and anem_id columns for merging w/anem.Visits below - not sure why for some reason does funky things when have species in there - lists some anemones twice even if site and species are the same...
# anem.Sites <- distinct(anem.Info[keeps])
# anem.Sites <- anem.Sites %>% filter(anem_id != "-9999") %>% collect() #remove -9999 anem_id
# 
# #create visited/not columns for each sampling year (figure out how to do this in a function or for loop!)
# anem.Visits$visited.2012 <- rep(0,length(anem.Visits$anem_id))
# anem.Visits$visited.2013 <- rep(0,length(anem.Visits$anem_id))
# anem.Visits$visited.2014 <- rep(0,length(anem.Visits$anem_id))
# anem.Visits$visited.2015 <- rep(0,length(anem.Visits$anem_id))
# anem.Visits$visited.2016 <- rep(0,length(anem.Visits$anem_id))
# anem.Visits$visited.2017 <- rep(0,length(anem.Visits$anem_id))
# #FIGURE OUT HOW TO DO THE ABOVE IN A FUNCTION OR FOR LOOP!!!
# # for (i in 1:length(sample.Years)) {
# #   x <- paste("visited", as.character(sample.Years[i]), sep=".")
# #   anem.Visits <- mutate(anem.Visits, sample.Years[i] = anem.IDs*0)
# # }
#   
# 
# #find first and last year visited
# first.Year <- min(as.integer(substring(dive.Dates$date, 1, 4))) #find the earliest dive year (by pulling the year out of the date and finding the minimum)
# last.Year <- max(as.integer(substring(dive.Dates$date, 1, 4))) #find the latest dive year 
# 
# sample.Years <- seq(first.Year, last.Year, 1) #still not sure how to make this work... Michelle got it to work w/separate, put x and y as inputs for the into input
# 
# 
