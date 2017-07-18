#Survival estimates for clownfish populations
#Last updated: 7/18/17

rm(list=ls())

##### Set-up:
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
#library(RMark)

##### Functions:
## Pull in functions from Michelle's GitHub
#conleyte function - connect to leyte database
script <- getURL("https://raw.githubusercontent.com/stuartmichelle/Phil_code/master/code/conleyte.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#dateanemobs - "which date was an anemone observed based on clownfish table? ;This function allows you to find the date based on the anem_table_id"
script <- getURL("https://raw.githubusercontent.com/stuartmichelle/Phil_code/master/code/dateanemobs.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#datefishcap - "which date was a sample_id captured?"
script <- getURL("https://raw.githubusercontent.com/stuartmichelle/Phil_code/master/code/datefishcap.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#siteanem function - "which site was an anemone observed based on clownfish table?", from: https://github.com/stuartmichelle/Phil_code/blob/master/code/siteanem.R
siteanem <- function(x){
  #source("~/Documents/Philippines/Phil_code/conleyte.R")
  leyte <- conleyte()
  # connect anem ids to dive ids
  anem <- leyte %>% tbl("anemones") %>% filter(anem_table_id %in% x) %>% select(dive_table_id, anem_table_id) %>% collect()
  # get site
  suppressWarnings(dive <- leyte %>% tbl("diveinfo") %>% filter(id %in% anem$dive_table_id) %>% select(id, name) %>% collect())
  site <- left_join(anem, dive, by = c("dive_table_id" = "id"))
  return(site)
}

#sitefish function, from: https://github.com/stuartmichelle/Phil_code/blob/master/code/sitefish.R
sitefish <- function(x){
  #source("~/Documents/Philippines/Phil_code/conleyte.R")
  leyte <- conleyte()
  # connect anem ids to dive ids
  anem <- leyte %>% tbl("anemones") %>% filter(anem_table_id %in% x) %>% select(dive_table_id, anem_table_id) %>% collect()
  # get site
  suppressWarnings(dive <- leyte %>% tbl("diveinfo") %>% filter(id %in% anem$dive_table_id) %>% select(id, name) %>% collect())
  site <- left_join(anem, dive, by = c("dive_table_id" = "id"))
  return(site)
}

##### Running things!
leyte = conleyte()
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

