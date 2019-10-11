# Pulls data from database, saves files for use in analysis

#################### Install packages: ####################
library(RCurl)  # allows running R scripts from GitHub
library(RMySQL)  # connect to SQL database
library(dplyr)
library(here)

#################### Functions: ####################
# # Functions from Michelle's GitHub helpers script so can connect to database
# script <- getURL("https://raw.githubusercontent.com/pinskylab/field/master/scripts/field_helpers.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))

# Function with myconfig file info in it, for some reason new version of R/RStudio can't find the database...
read_db <- function(x) {
  db <- src_mysql(dbname = x,
                  port = 3306,
                  create = F,
                  host = "amphiprion.deenr.rutgers.edu",
                  user = "allisond",
                  password = "fish!NM?717")
  return(db)
}

#################### Pull out database info: ####################
leyte <- read_db("Leyte") 

anem_db <- leyte %>% tbl("anemones") %>% collect()  # anemone table
fish_db <- leyte %>% tbl("clownfish") %>% collect()  # fish caught and seen on clownfish dives table
fish_seen_db <- leyte %>% tbl('clown_sightings') %>% collect()  # fish seen by 'clownfish diver' table
dives_db <- leyte %>% tbl("diveinfo") %>% collect()  # dive info table
gps_db <- leyte %>% tbl("GPX") %>% collect()  # gps info table

#################### Pull data from GitHub: ####################
# Fish obs table linking together individual fish marked and recaught by both genetic and tag methods
download.file(url = "https://github.com/pinskylab/genomics/blob/master/data/fish-obs.RData?raw=true", destfile = "fish-obs.RData")
fish_obs <- readRDS("fish-obs.RData")

#################### Save files: ####################
save(anem_db, file = here::here("Data/Data_from_database", "anem_db.RData"))
save(fish_db, file = here::here("Data/Data_from_database", "fish_db.RData"))
save(fish_seen_db, file = here::here("Data/Data_from_database", "fish_seen_db.RData"))
save(dives_db, file = here::here("Data/Data_from_database", "dives_db.RData"))
save(gps_db, file = here::here("Data/Data_from_database", "gps_db.RData"))
save(fish_obs, file = here::here("Data/From_other_analyses", "fish-obs.RData"))

# saveRDS() - use this so that people can load things and choose their own names rather than being stuck with the names I choose

