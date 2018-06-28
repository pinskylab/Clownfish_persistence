# Get an estimate for recapture probability from recapture dives

## To-do: give Cabatoan and Magbangon a 1-m vis for 1 dive during 2017
## Try making a scatter plot, showing all dives and their vis
## Helpful to put the Baybay sites in there too, just for reference?

#################### Set-up: ####################
# Load relevant libraries
library(RCurl)
library(dplyr)
#library(tidyr)
#library(tidyverse)
library(lubridate)
#library(dbplyr)
library(ggplot2)
library(here)

#################### Functions: ####################
# Load helpful functions from other scripts
script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Functions from Michelle's GitHub helpers script
#helper functions - do various tasks w/database (like assigning dates and site to fish and such)
script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Lincoln-Peterson estimate (biased for small sample sizes)
mr_LincolnPeterson <- function(x11, x10, x01) { # x11 = individuals captured in both sampling sessions, x10 = individuals captured only in the first sampling session, x01 = individuals captured only in the second sampling session
  n1 <- x11 + x10     #animals caught + marked in first sampling period (and total number of marked individuals) 
  n2 <- x11 + x01     #number of animals caught in the second sampling event
  m2 <- x11           #number of animals caught in both periods
  r <- n1 + n2 - m2   #number of distinct animals caught during study
  
  Nest <- (n1*n2)/m2
  varN <- ((n1 + 1)*(n2 + 1)*(n1 - m2)*(n2 - m2))/(((m2 + 1)^2)*(m2 + 2))
  prob_r <- m2/n2
  prob_rv2 <- m2/n1
  
  out <- data.frame(Nest = Nest, varN = varN, prob_r = prob_r, prob_rv2 = prob_rv2)
  # # n = total individuals caputured (in second capture event), m = marked individuals captured (in second capture event), M = total marked individuals
  # N <- (n/m)*M  # estimated total population size
  # pr <- m/M    # estimated probablility of recapture
  # out <- data.frame(n = n, prob_recap = pr)
  return(out)
}

# Lincoln-Peterson estimate (biased for small sample sizes), alternate version (THIS IS THE ONE I'VE BEEN USING)
mr_LincolnPeterson_alt <- function(n1, n2, m2) { # x11 = individuals captured in both sampling sessions, x10 = individuals captured only in the first sampling session, x01 = individuals captured only in the second sampling session
  #n1 <- x11 + x10     #animals caught + marked in first sampling period (and total number of marked individuals) 
  #n2 <- x11 + x01     #number of animals caught in the second sampling event
  #m2 <- x11           #number of animals caught in both periods
  r <- n1 + n2 - m2   #number of distinct animals caught during study
  
  Nest <- (n1*n2)/m2
  varN <- ((n1 + 1)*(n2 + 1)*(n1 - m2)*(n2 - m2))/(((m2 + 1)^2)*(m2 + 2))
  prob_r <- m2/n2
  prob_rv2 <- m2/n1
  prop_caught <- r/Nest
  prop_caught_n1 <- n1/Nest
  
  out <- data.frame(Nest = Nest, varN = varN, prob_r = prob_r, prob_rv2 = prob_rv2, prop_caught = prop_caught, prop_caught_n1 = prop_caught_n1)
  # # n = total individuals caputured (in second capture event), m = marked individuals captured (in second capture event), M = total marked individuals
  # N <- (n/m)*M  # estimated total population size
  # pr <- m/M    # estimated probablility of recapture
  # out <- data.frame(n = n, prob_recap = pr)
  return(out)
}

# Chapman estimated (unbiased for large N)
mr_Chapman <- function(n, m, M) {
  N <- (((M+1)*(n+1))/(m+1)) - 1 # estimated total population size
  pr <- m/M  # estimated probability of recapture
  varN <- (N^2*(n-m))/((n+1)*(m+2)) # variance of N
  out <- data.frame(n = n, prob_recap = pr, varN = varN)
  return(out)
}

#################### Running things: ####################
##### Pull out info from database
leyte <- read_db("Leyte")

# Pull out all the fish, filter for just APCL
allfish_fish <- leyte %>% 
  tbl("clownfish") %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, cap_id, anem_table_id, recap, tag_id, color, size) %>%
  collect() %>%
  filter(fish_spp == "APCL")

# Pull out the anems associated with those fish
allfish_anems <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
  collect() %>%
  filter(anem_table_id %in% allfish_fish$anem_table_id)

# Pull out the dives associated with those fish
allfish_dives <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, dive_type, date, site, gps) %>%
  collect() %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id)

# Pull out just the year and put that in a separate column
allfish_dives$year <- as.integer(substring(allfish_dives$date,1,4))

# Join together
allfish <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
allfish <- left_join(allfish, allfish_dives, by="dive_table_id")

##### Try looking at recap dives for Poroc Rose
# recap_dives <- allfish %>% 
#   filter(dive_type == "R")
# 
# # Try Poroc Rose b/c small site - probably don't need to worry about coverage
# PorocRose_2015dives <- allfish %>% #no recap dives in Poroc Rose in 2015...
#   filter(site == "Poroc Rose", year == 2015)
# PorocRose_2016dives <- allfish %>%
#   filter(site == "Poroc Rose", year == 2016)
# PorocRose_2017dives <- allfish %>%
#   filter(site == "Poroc Rose", year == 2017)

# Find the number of fish marked during C dives (sampling session 1) -- this assumes no tags are encountered twice... will need to figure out how to deal with that
#### THIS IS THE PART THAT WOULD NEED TO GET ALTERED FOR AREA - WOULD ONLY WANT THE FISH TAGGED (EITHER NEWLY OR ALREADY) FROM THE AREA COVERED BY THE RECAP DIVE
n_tagged_fish <- allfish %>%
  filter(dive_type == "C", !is.na(tag_id)) %>%
  group_by(site, year) %>%
  summarise(ntags = n())

# Find the number of fish seen during R dives (sampling session 2)
n_fish_seen_R <- allfish %>%
  filter(dive_type == "R") %>%
  group_by(site, year) %>%
  summarise(nfish = n())

# Find the number of previously-tagged fish caught during R dives (sampling session 2)
n_recaptured_fish_R <- allfish %>%
  filter(dive_type == "R", !is.na(tag_id), recap == "Y") %>%
  distinct(site, year, tag_id) %>%
  group_by(site, year) %>%
  summarise(ntags = n())

# just looking...
allfish %>% filter(site == "Poroc San Flower", year == 2016, dive_type == "R")

##### Try running for a few sites just to see how it looks...
n_PSF <- n_fish_seen_R %>% filter(site == "Poroc San Flower", year == 2016)
estimates_PorocSF <- mr_LincolnPeterson((n_fish_seen_R %>% filter(site == "Poroc San Flower", year == 2016))$nfish, (n_recaptured_fish_R %>% filter(site == "Poroc San Flower", year == 2016))$ntags, (n_tagged_fish %>% filter(site == "Poroc San Flower",  year == 2016))$ntags)
estimates_PorocSF2017 <- mr_LincolnPeterson((n_fish_seen_R %>% filter(site == "Poroc San Flower", year == 2017))$nfish, (n_recaptured_fish_R %>% filter(site == "Poroc San Flower", year == 2017))$ntags, (n_tagged_fish %>% filter(site == "Poroc San Flower",  year == 2017))$ntags)
estimates_PorocSF2016_C <- mr_Chapman((n_fish_seen_R %>% filter(site == "Poroc San Flower", year == 2016))$nfish, (n_recaptured_fish_R %>% filter(site == "Poroc San Flower", year == 2016))$ntags, (n_tagged_fish %>% filter(site == "Poroc San Flower",  year == 2016))$ntags)
estimates_PorocSF2017_C <- mr_Chapman((n_fish_seen_R %>% filter(site == "Poroc San Flower", year == 2017))$nfish, (n_recaptured_fish_R %>% filter(site == "Poroc San Flower", year == 2017))$ntags, (n_tagged_fish %>% filter(site == "Poroc San Flower",  year == 2017))$ntags)

##### Try running for all sites with recap dives, without filtering by area at all (not a great plan, will way underestimate proportion sampled)
# initialize data frame
#prop_sampled <- data.frame(site=as.character(site_vec), stringsAsFactors = FALSE)
recap_dive_sites_2016 <- c("Cabatoan","Palanas","Poroc Rose","Poroc San Flower","Wangag")
recap_dive_sites_2017 <- c("Gabas", "Palanas", "Poroc Rose", "Poroc San Flower", "Visca")
recap_dive_sites_2018 <- c("Hicgop South", "Palanas", "San Agustin", "Sitio Baybayon")
prop_sampled <- data.frame(site=c(recap_dive_sites_2016,recap_dive_sites_2017,recap_dive_sites_2018), stringsAsFactors = FALSE)
prop_sampled$year <- c(rep(2016, length(recap_dive_sites_2016)), rep(2017, length(recap_dive_sites_2017)), rep(2018, length(recap_dive_sites_2018)))
prop_sampled$Nest <- rep(NA, length(prop_sampled$year))
prop_sampled$Nvar <- rep(NA, length(prop_sampled$year))
prop_sampled$prop_sampled_noQGIS_total <- rep(NA, length(prop_sampled$year))
prop_sampled$prop_sampled_noQGIS_n1 <- rep(NA, length(prop_sampled$year))

# run through sites and years
for(i in 1:length(prop_sampled$year)){
  siteval <- prop_sampled$site[i]
  y <- prop_sampled$year[i]
  
  n1 = (n_tagged_fish %>% filter(site == siteval, year == y))$ntags
  n2 = (n_fish_seen_R %>% filter(site == siteval, year == y))$nfish
  m2 = (n_recaptured_fish_R %>% filter(site == siteval, year == y))$ntags
  
  # if( (n1 == 0 | n2 == 0 | m2 == 0) == TRUE) {
  #   prop_sampled$Nest[index] = NA
  #   prop_sampled$prop_sampled_noQGIS_total[index] = NA
  #   prop_sampled$prop_sampled_noQGIS_n1[index] = NA
  # } else {
    out = mr_LincolnPeterson_alt(n1, n2, m2)
    prop_sampled$Nest[i] = out$Nest
    prop_sampled$Nvar[i] = out$varN
    prop_sampled$prop_sampled_noQGIS_total[i] = out$prop_caught
    prop_sampled$prop_sampled_noQGIS_n1[i] = out$prop_caught_n1
 # }
}

### Save output file so can access without having to re-run script
save(prop_sampled, file=here("Data", "prop_sampled_RecapDives.RData"))




