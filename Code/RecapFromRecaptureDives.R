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

# Peterson estimate (biased for small sample sizes)
mr_Peterson <- function(n, m, M) { # n = total individuals caputured (in second capture event), m = marked individuals captured (in second capture event), M = total marked individuals
  N <- (n/m)*M  # estimated total population size
  pr <- m/M    # estimated probablility of recapture
  out <- data.frame(n = n, prob_recap = pr)
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

# Pull out the dives associated with thsoe fish
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

##### Try running for a few sites...
n_PSF <- n_fish_seen_R %>% filter(site == "Poroc San Flower", year == 2016)
estimates_PorocSF <- mr_Peterson((n_fish_seen_R %>% filter(site == "Poroc San Flower", year == 2016))$nfish, (n_recaptured_fish_R %>% filter(site == "Poroc San Flower", year == 2016))$ntags, (n_tagged_fish %>% filter(site == "Poroc San Flower",  year == 2016))$ntags)
estimates_PorocSF2017 <- mr_Peterson((n_fish_seen_R %>% filter(site == "Poroc San Flower", year == 2017))$nfish, (n_recaptured_fish_R %>% filter(site == "Poroc San Flower", year == 2017))$ntags, (n_tagged_fish %>% filter(site == "Poroc San Flower",  year == 2017))$ntags)
estimates_PorocSF2016_C <- mr_Chapman((n_fish_seen_R %>% filter(site == "Poroc San Flower", year == 2016))$nfish, (n_recaptured_fish_R %>% filter(site == "Poroc San Flower", year == 2016))$ntags, (n_tagged_fish %>% filter(site == "Poroc San Flower",  year == 2016))$ntags)
estimates_PorocSF2017_C <- mr_Chapman((n_fish_seen_R %>% filter(site == "Poroc San Flower", year == 2017))$nfish, (n_recaptured_fish_R %>% filter(site == "Poroc San Flower", year == 2017))$ntags, (n_tagged_fish %>% filter(site == "Poroc San Flower",  year == 2017))$ntags)
