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

#################### Pull data from GitHub from other analyses: ####################

# # Fish obs table linking together individual fish marked and recaught by both genetic and tag methods
# download.file(url = "https://github.com/pinskylab/genomics/blob/master/data/fish-obs.RData?raw=true", destfile = "fish-obs.RData")
# fish_obs <- readRDS("fish-obs.RData")

# Fish obs table linking together individual fish marked and recaught by both genetic and tag methods
fish_obs <- read.csv(text = getURL("https://raw.githubusercontent.com/pinskylab/genomics/master/data/fish-obs.csv"), header = TRUE, stringsAsFactors = FALSE)

# Site centers (eyeballed from QGIS by KC)
site_centers <- read.csv(file = here::here("Data/From_other_analyses", "site_centroids.csv"), header = TRUE, stringsAsFactors = FALSE)

# Kernel fits
kernel_summary <- read.csv(text = getURL("https://raw.githubusercontent.com/katcatalano/parentage/master/kernel_fitting/1340_loci/final_results/tables/kernel_fitting_summary.csv"), header = T, stringsAsFactors = F)
k_theta_allyear_95CI_values <- read.csv(text = getURL("https://raw.githubusercontent.com/katcatalano/parentage/master/kernel_fitting/1340_loci/final_results/likelihood_profiles_grid_search/profile95CI_AllYears.csv"), header = T, stringsAsFactors = F)
#k_theta_allyear_95CI_values <- read.csv(text = getURL("https://raw.githubusercontent.com/katcatalano/parentage/master/kernel_fitting/1340_loci/final_results/likelihood_profiles_grid_search/AllYearParams95CI.csv"), header = T, stringsAsFactors = F)

# List of all parents put into parentage analysis (so can match to site)
all_parents <- read.table(text = getURL("https://raw.githubusercontent.com/katcatalano/parentage/master/colony2/20190523_1340loci/input/all_parents_corrected.txt"), header = T, stringsAsFactors = F)

# List of all offspring put into parentage analysis
all_offspring <- read.table(text = getURL("https://raw.githubusercontent.com/katcatalano/parentage/master/colony2/20190523_1340loci/input/all_offspring_corrected.txt"), header = T, stringsAsFactors = F)

# Size-fecundity model from Adam
load(file=here::here('Data', 'loglogFecunditySizeModel.RData'))  # size-fecundity output for best-fit model from Adam, called length_count8llEA

#################### Save files: ####################
# Save database files
save(anem_db, file = here::here("Data/Data_from_database", "anem_db.RData"))
save(fish_db, file = here::here("Data/Data_from_database", "fish_db.RData"))
save(fish_seen_db, file = here::here("Data/Data_from_database", "fish_seen_db.RData"))
save(dives_db, file = here::here("Data/Data_from_database", "dives_db.RData"))
save(gps_db, file = here::here("Data/Data_from_database", "gps_db.RData"))

# Save files from other analyses
save(fish_obs, file = here::here("Data/From_other_analyses", "fish-obs.RData"))  # matches up multiple observations of fish individuals (from Michelle, sourced via GitHub)
save(site_centers, file = here::here("Data/From_other_analyses", "site_centers.RData"))
save(kernel_summary, file = here::here("Data/From_other_analyses", "kernel_summary.RData"))  # kernel fits (from Katrina, sourced via GitHub)
save(all_parents, file = here::here("Data/From_other_analyses", "all_parents.RData"))  # all parents put into parentage analysis (from Katrina, sourced via GitHub)
save(all_offspring, file = here::here("Data/From_other_analyses", "all_offspring.RData"))  # all offspring put into parentage analysis (from Katrina, sourced via GitHub)
save(k_theta_allyear_95CI_values, file = here::here("Data/From_other_analyses", "k_theta_allyear_95CI_values.RData"))  # 95% CI values for k and theta fits
save(length_count8llEA, file= here::here("Data/From_other_analyses", "length_count8llEA.RData"))  # size-fecundity model from Adam

# saveRDS() - use this so that people can load things and choose their own names rather than being stuck with the names I choose

