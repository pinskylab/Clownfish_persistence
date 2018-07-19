# Mark-recap analysis to estimate survival, recapture probability, and population size based on tags and genetic IDs
# Somewhat cleaned-up version of code from ClownfishSurvivalEstimates.R

#################### Set-up: ####################
#Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(RMark)
library(lubridate)
#library(dbplyr)
library(ggplot2)
library(here)
#library(rethinking)

#Load data input
load(file=here("Data", "fish_Tagged.RData")) #file with distances appended

#################### Functions: ####################
# Functions and constants from my GitHub function/constant collection
# script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQAXExRItr2tnepmkRr7NRt4hylZrks5aciBtwA%3D%3D", ssl.verifypeer = FALSE)
# eval(parse(text = script))
script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Functions from Michelle's GitHub helpers script
#helper functions - do various tasks w/database (like assigning dates and site to fish and such)
script <- getURL("https://raw.githubusercontent.com/mstuart1/helpers/master/scripts/helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Finds the real parameter estimate from the logit estimate
logit_recip <- function(logitval) {
  recip = (exp(logitval))/(1 + exp(logitval))
  return(recip)
}

# Creates the summarized encounter history by tag id: output is a data frame with 2 columns - tag id and summarized encounter history (i.e. 0010)
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

# # Function to calculate likelihood for 3 years of sampling (check this!!) with fish tagged in year 2 included too
# Like3SampleYears4vars <- function(n111, n110, n101, n100, n011, n010, s1, s2, p2, p3) {
#   #like <- ((s1*p2*s2*p3)^n111)*((s1*p2*(1-s2*p3))^n110)*((s1*(1-p2)*s2*p3)^n101)*((1-s1*p2-s1*(1-p2)*s2*p3)^n100)
#   x2 <- 1 - s2 + s2*(1 - p3) #probability of not being seen after the second sampling session
#   x1 <- 1 - s1 + s1*(1-p2)*x2 #probability of not being seen after the first sampling session
#   
#   p111 <- s1*p2*s2*p3 #probability of encounter history 111
#   p110 <- s1*p2*x2 #probability of encounter history 110
#   p101 <- s1*(1 - p2)*s2*p3 #probability of encounter history 101
#   p100 <- x1 #probability of encounter history 100
#   p011 <- s2*p3 #probability of encounter history 011
#   p010 <- x2 #probability of encounter history 010
#   
#   like <- n111*log(p111) + n110*log(p110) + n101*log(p101) + n100*log(p100) + n011*log(p011) + n010*log(p010)
#   return(like)
# }
# 
# # Same as above but for only two variables (assuming that survival and recapture probs are same through time)
# Like3SampleYears2vars <- function(n111, n110, n101, n100, n011, n010, s, p) {
#   #like <- ((s1*p2*s2*p3)^n111)*((s1*p2*(1-s2*p3))^n110)*((s1*(1-p2)*s2*p3)^n101)*((1-s1*p2-s1*(1-p2)*s2*p3)^n100)
#   x2 <- 1 - s + s*(1 - p) #probability of not being seen after the second sampling session
#   x1 <- 1 - s + s*(1-p)*x2 #probability of not being seen after the first sampling session
#   
#   p111 <- s*p*s*p #probability of encounter history 111
#   p110 <- s*p*x2 #probability of encounter history 110
#   p101 <- s*(1 - p)*s*p #probability of encounter history 101
#   p100 <- x1 #probability of encounter history 100
#   p011 <- s*p #probability of encounter history 011
#   p010 <- x2 #probability of encounter history 010
#   
#   like <- n111*log(p111) + n110*log(p110) + n101*log(p101) + n100*log(p100) + n011*log(p011) + n010*log(p010)
#   return(like)
# }
# 
# # Just cycles through vectors of possible survival and recapture probabilities and calculates the likelihood
# ParamRangeCheck3SY2P_dfoutput <- function(s_vec, p_vec, nvals) { 
#   #set index values for various encouter histories (this is the way R orders them using table)
#   N010 <- 2
#   N011 <- 3
#   N100 <- 4
#   N101 <- 5
#   N110 <- 6
#   N111 <- 7
#   
#   #get numbers of each type of encounter history from the inputs
#   n111 <- nvals$Freq[N111]
#   n110 <- nvals$Freq[N110]
#   n101 <- nvals$Freq[N101]
#   n100 <- nvals$Freq[N100]
#   n011 <- nvals$Freq[N011]
#   n010 <- nvals$Freq[N010]
#   
#   survival <- data.frame(survival=rep(s_vec,length(p_vec)))
#   out <- survival
#   capture_list <- rep(p_vec[1], length(s_vec))
#   for (i in 2:length(p_vec)) {
#     cap_toadd <- rep(p_vec[i], length(s_vec))
#     capture_list <- c(capture_list, cap_toadd)
#   }
#   out$capture <- capture_list
#   like = NULL
#   
#   for (i in 1:dim(out)[1]) {
#     like1 <- Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, out$survival[i], out$capture[i])
#     like <- c(like, like1)
#   }
#   out$like=like
#   return(out)
# }
# 
# # Gradient descent for 3 sample years, estimating 4 parameters (doesn't really seem to work, oscillates between all params essentially at 0 or all at 1)
# GradientDescent3SY4Params <- function(nvals, params0, nsteps, eps, del, highVal, lowVal) {
#   param1 <- rep(NA,1,nsteps)
#   out <- as.data.frame(param1)
#   out$param2 <- rep(NA,1,nsteps)
#   out$param3 <- rep(NA,1,nsteps)
#   out$param4 <- rep(NA,1,nsteps)
#   out$like <- rep(NA,1,nsteps)
#   
#   params <- params0
#   
#   #set index values for various encouter histories (this is the way R orders them using table) (001 is 1)
#   N010 <- 2
#   N011 <- 3
#   N100 <- 4
#   N101 <- 5
#   N110 <- 6
#   N111 <- 7
#   
#   #get numbers of each type of encounter history from the inputs
#   n111 <- nvals$Freq[N111]
#   n110 <- nvals$Freq[N110]
#   n101 <- nvals$Freq[N101]
#   n100 <- nvals$Freq[N100]
#   n011 <- nvals$Freq[N011]
#   n010 <- nvals$Freq[N010]
#   
#   for (i in 1:nsteps) { #doesn't work with Palanas data b /c get NANs
#     #store parameter values and likelihood
#     out$param1[i] <- params[1]
#     out$param2[i] <- params[2]
#     out$param3[i] <- params[3]
#     out$param4[i] <- params[4]
#     out$like[i] <- Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4])
#     
#     #calculate the gradient and update the parameters
#     df <- rep(NA, length(params))
#     df[1] <- (Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1]+eps, params[2], params[3], params[4]) - Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4]))/eps
#     df[2] <- (Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2]+eps, params[3], params[4]) - Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4]))/eps
#     df[3] <- (Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3]+eps, params[4]) - Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4]))/eps
#     df[4] <- (Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4]+eps) - Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4]))/eps
#     
#     params <- params+del*df
#     
#     #if parameters have become <= 0 or >= 1, set to low val or high val
#     for (i in 1:length(params)) {
#       if (params[i] >= 1) {
#         params[i] <- highVal
#       }
#       if (params[i] <= 0) {
#         params[i] <- lowVal
#       }
#     }
#     #print(Like3SampleYears(n111, n110, n101, n100, params[1], params[2], params[3], params[4]))
#   }
#   return(out)
# }
# 
# # Gradient Descent 3 years, 2 params
# GradientDescent3SY2Params <- function(nvals, params0, nsteps, eps, del, highVal, lowVal) {
#   survival_p <- rep(NA,1,nsteps)
#   out <- as.data.frame(survival_p)
#   out$recapture_p <- rep(NA,1,nsteps)
#   #out$param3 <- rep(NA,1,nsteps)
#   #out$param4 <- rep(NA,1,nsteps)
#   out$like <- rep(NA,1,nsteps)
#   
#   params <- params0
#   
#   #set index values for various encouter histories (this is the way R orders them using table) (001 is 1)
#   N010 <- 2
#   N011 <- 3
#   N100 <- 4
#   N101 <- 5
#   N110 <- 6
#   N111 <- 7
#   
#   #get numbers of each type of encounter history from the inputs
#   n111 <- nvals$Freq[N111]
#   n110 <- nvals$Freq[N110]
#   n101 <- nvals$Freq[N101]
#   n100 <- nvals$Freq[N100]
#   n011 <- nvals$Freq[N011]
#   n010 <- nvals$Freq[N010]
#   
#   for (i in 1:nsteps) { #doesn't work with Palanas data b /c get NANs
#     #store parameter values and likelihood
#     out$survival_p[i] <- params[1]
#     out$recapture_p[i] <- params[2]
#     out$like[i] <- Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, params[1], params[2])
#     
#     #calculate the gradient and update the parameters 
#     df <- rep(NA, length(params))
#     df[1] <- (Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, params[1]+eps, params[2]) - Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, params[1], params[2]))/eps
#     df[2] <- (Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, params[1], params[2]+eps) - Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, params[1], params[2]))/eps
#     
#     params <- params+del*df
#     
#     #if parameters have become <= 0 or >= 1, set to low val or high val
#     for (i in 1:length(params)) {
#       if (params[i] >= 1) {
#         params[i] <- highVal
#       }
#       if (params[i] <= 0) {
#         params[i] <- lowVal
#       }
#     }
#   }
#   return(out)
# }


#################### Running things: data manipulation and setup ####################
leyte <- read_db("Leyte")

allfish_fish <- leyte %>% 
  tbl("clownfish") %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, cap_id, anem_table_id, recap, tag_id, color, size) %>%
  collect() #%>%
  #filter(!is.na(tag_id)) 

allfish_anems <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
  collect() %>%
  filter(anem_table_id %in% allfish_fish$anem_table_id)

allfish_dives <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, dive_type, date, site, gps) %>%
  collect() %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id)

# pull out just the year and put that in a separate column
allfish_dives$year <- as.integer(substring(allfish_dives$date,1,4))

#join together
allfish <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
allfish <- left_join(allfish, allfish_dives, by="dive_table_id")

allfish$size <- as.numeric(allfish$size) #make size numeric (rather than a chr) so can do means and such

##### Trying out MARK
#Create input file for MARK
# Pull out all the tags and their encounter history from 2015-2018
encounters_all <- CreateEncounterSummary(2015, 2018, allfish) #simpler way of getting encounter history like above using function

# Remove the NA tag id at the end
encounters_all <- encounters_all[1:(length(encounters_all$tag_id)-1),] 

# Find site and add it in
site_info <- allfish %>% filter(!is.na(tag_id)) %>%
  filter(tag_id %in% encounters_all$tag_id) %>%
  group_by(tag_id) %>%
  summarize(site = site[1]) # this assumes that the fish was only seen/caught at one site - should probably check that at some point...

encounters_all <- left_join(encounters_all, site_info, by="tag_id")

# Rename to fit mark input expectations
encounters_all <- encounters_all %>% rename(ch = encounter.hist)

# Add in tagging size
tag_size <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  group_by(tag_id) %>%
  summarize(tag_size = min(size)) # this assumes it was tagged at the smallest size recorded for it... will not be true when include genetic recaps (and might not be true now anyway...)
  
encounters_all <- left_join(encounters_all, tag_size, by="tag_id")

# Add in tail color when tagged
tail_color <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  group_by(tag_id) %>%
  summarize(tail_color = color[1])

encounters_all <- left_join(encounters_all, tail_color, by="tag_id")

## Code that adds in tail color by year, problem with that is that we only have it for fish that are caught, which isn't allowed for time-varying individual covariates (have to have a value for each fish at each time)
##might be able to get around that by just doing the color at the last time it was caught but not sure that's worth it so for now, commenting out...
# Add in tail color (could probably do these all together in one command... look into that later..)
# tail_color_2015 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#   filter(year == 2015) %>%
#   group_by(tag_id) %>%
#   summarize(tail_color_2015 = color[1])
# 
# tail_color_2016 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#   filter(year == 2016) %>%
#   group_by(tag_id) %>%
#   summarize(tail_color_2016 = color[1])
# 
# tail_color_2017 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#   filter(year == 2017) %>%
#   group_by(tag_id) %>%
#   summarize(tail_color_2017 = color[1])
# 
# tail_color_2018 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#   filter(year == 2018) %>%
#   group_by(tag_id) %>%
#   summarize(tail_color_2018 = color[1])
# 
# #encounters[[i]] <- tagged.fish %>% group_by(tag_id) %>% summarise(!!var.name := ifelse(sum(year == sample.years[i])>0, 1, 0)) #create data frames for each year that have a vector of tag ids and a vector of encountered (1) or didn't (0) in 
# 
# encounters_all <- left_join(encounters_all, tail_color_2015, by="tag_id")
# encounters_all <- left_join(encounters_all, tail_color_2016, by="tag_id")
# encounters_all <- left_join(encounters_all, tail_color_2017, by="tag_id")
# encounters_all <- left_join(encounters_all, tail_color_2018, by="tag_id")
#

## Same issue with size in each year as with color above - don't have a value for the fish that weren't recaptured
#Might be able to add some in based on growth curves using the relationships Michelle has, could get into that at some point but for now ignoring
# Add in size in each of the years
# size_2015 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#   filter(year == 2015) %>%
#   group_by(tag_id) %>% 
#   summarize(size_2015 = mean(size, rm.na = TRUE))
# 
# size_2016 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#   filter(year == 2016) %>%
#   group_by(tag_id) %>% 
#   summarize(size_2016 = mean(size, rm.na = TRUE))
# 
# size_2017 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#   filter(year == 2017) %>%
#   group_by(tag_id) %>% 
#   summarize(size_2017 = mean(size, rm.na = TRUE))
# 
# size_2018 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#   filter(year == 2018) %>%
#   group_by(tag_id) %>% 
#   summarize(size_2018 = mean(size, rm.na = TRUE))
# 
# encounters_all <- left_join(encounters_all, size_2015, by="tag_id")
# encounters_all <- left_join(encounters_all, size_2016, by="tag_id")
# encounters_all <- left_join(encounters_all, size_2017, by="tag_id")
# encounters_all <- left_join(encounters_all, size_2018, by="tag_id")

# Add in distances to anem (from fish.Tagged, loaded above)
encounters_all <- left_join(encounters_all, (fish.Tagged %>% select(tag_id, year_tagged, dist_2016, dist_2017, dist_2018)), by="tag_id") 

# can't have NAs in the covariate columns so changing to mean of distances in each year (should also try replacing with 0, since it shouldn't matter)
#mean2015 <- mean(encounters_all$dist_2015, na.rm=TRUE)
mean2016 <- mean(encounters_all$dist_2016, na.rm=TRUE)
mean2017 <- mean(encounters_all$dist_2017, na.rm=TRUE)
mean2018 <- mean(encounters_all$dist_2018, na.rm=TRUE)

# replace the NAs (for years before the fish was tagged) with the mean distance that year
encounters_means <- encounters_all %>%
  rename(dist2016 = dist_2016, dist2017 = dist_2017, dist2018 = dist_2018) %>% #first, rename distance columns so easier for MARK to find
  mutate(dist2016 = replace(dist2016, is.na(dist2016), mean2016)) %>%
  mutate(dist2017 = replace(dist2017, is.na(dist2017), mean2017)) %>%
  mutate(dist2018 = replace(dist2018, is.na(dist2018), mean2018))

# replace the NAs (for years before the fish was tagged) with 0 - could also just calculate distance to that anem, based on the first tag number in the year the fish was tagged in the original script that finds the distance
encounters_0 <- encounters_all %>%
  rename(dist2016 = dist_2016, dist2017 = dist_2017, dist2018 = dist_2018) %>% #first, rename distance columns so easier for MARK to find
  mutate(dist2016 = replace(dist2016, is.na(dist2016), 0)) %>%
  mutate(dist2017 = replace(dist2017, is.na(dist2017), 0)) %>%
  mutate(dist2018 = replace(dist2018, is.na(dist2018), 0))

#################### Running things: MARK models ####################

#### Models with distance (with 0s put in for distance pre-fish first caught), size-at-tagging, site, tail color-at-tagging
# data
eall <- encounters_0 %>%
  select(ch, site, tag_size, tail_color, dist2016, dist2017, dist2018) %>% #select relevant columns
  mutate(site = as.factor(site))

eall <- eall[complete.cases(eall),] #for using tag size, need no NAs...

# process data and make ddl
eall.processed = process.data(eall, model="CJS", begin.time=2015, groups="site")
eall.ddl = make.design.data(eall.processed)

eall.processed2 = process.data(eall, model="CJS", begin.time=2015, groups="tail_color")
eall.ddl2 = make.design.data(eall.processed2)


# set models for Phi
Phi.dot = list(formula=~1, link="logit")
Phi.time = list(formula=~time, link="logit")
Phi.site = list(formula=~site, link="logit")
Phi.size = list(formula=~tag_size, link="logit")
Phi.color = list(formula=~tail_color, link="logit")

# set models for p
p.dot = list(formula=~1, link="logit")
p.time = list(formula=~time, link="logit")
p.dist = list(formula=~dist, link="logit")
p.time.plus.dist = list(formula=~time+dist, link="logit")
p.site = list(formula=~site, link="logit")

# run some models
eall.Phi.dot.p.dot = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
eall.Phi.dot.p.time = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.dot, p=p.time)) #has issues running, can't estimate p in 2018
eall.Phi.dot.p.dist = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist)) #think this ran!!!
eall.Phi.size.p.dot = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.size, p=p.dot))
eall.Phi.size.p.site = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.size, p=p.site))
eall.Phi.size.p.dist = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.size, p=p.dist))
eall.Phi.site.p.site = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.site, p=p.site)) #this one has issues converging for some of the sites
eall.Phi.color.p.dot = mark(eall.processed2, eall.ddl2, model.parameters=list(Phi=Phi.color, p=p.dot))
eall.Phi.site.p.dot = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.site, p=p.dot)) #hmmm, something seems funny with this one...
eall.Phi.site.p.dist = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.site, p=p.dist))

#edist.Phi.dot.p.dist.plus.time = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.dot, p=p.time.plus.dist))
# edist.Phi.time.p.dot = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.dot))
# edist.Phi.time.p.time = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.time))
# edist.Phi.time.p.dist = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.dist))
# edist.Phi.time.p.time.plus.dist = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.time.plus.dist))

# organize output to get ready to make plots
###### Constant survival and recapture prob, no covariates (eall.Phi.dot.p.dot)
eall.constant = as.data.frame(eall.Phi.dot.p.dot$results$beta) %>% 
  mutate(param = c("Phi","p")) %>% #add a parameters column to make it easier to plot in ggplot
  mutate(upper = logit_recip(ucl), lower = logit_recip(lcl), est = logit_recip(estimate)) #do the reciprocal transform on upper and lower confidence limits (need to check what those are - 95? SE? and that transforming them just straight up is the right way to go)

###### Survival constant, p varies by distance (eall.Phi.dot.p.dist)
mindist = min(eall$dist2016,eall$dist2017,eall$dist2018) #this is 0, obv...
maxdist1 = max(eall$dist2016,eall$dist2017,eall$dist2018) #seems kind of high and probably due to an error
maxdist2 = 200 #encompasses highest values in 2016, 2018 
#tail(sort(eall$dist2016),10)
#tail(sort(eall$dist2017),10) #really 2017 that has the high distances...
#tail(sort(eall$dist2018),10)
dist.values1 = mindist+(0:30)*(maxdist1-mindist)/30
dist.values2 = mindist+(0:30)*(maxdist2-mindist)/30

# Not getting this to work - it doesn't change with distance... maybe/probably I'm specifying the index incorrectly?
#pbydist1 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values),indices=c(7))
#pbydist2 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values2),indices=c(1))
#pbydist3 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values3),indices=c(1))

# Set up dataframe 
eall.Phi.dot.p.dist.results = as.data.frame(eall.Phi.dot.p.dist$results$beta) 

# Large range of distances (up to 5000m)
pbydist_longrange= data.frame(dist = dist.values1) %>%
  mutate(p_logit = (eall.Phi.dot.p.dist.results$estimate[2]) + (eall.Phi.dot.p.dist.results$estimate[3]*dist)) %>%
  mutate(p_lcl_logit = eall.Phi.dot.p.dist.results$lcl[2] + (eall.Phi.dot.p.dist.results$lcl[3]*dist)) %>%
  mutate(p_ucl_logit = eall.Phi.dot.p.dist.results$ucl[2] + (eall.Phi.dot.p.dist.results$ucl[3]*dist)) %>%
  mutate(p = logit_recip(p_logit)) %>%
  mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
  mutate(p_ucl = logit_recip(p_ucl_logit)) 

# Short range of distances (up to 200m)
pbydist_shortrange = data.frame(dist = dist.values2) %>%
  mutate(p_logit = (eall.Phi.dot.p.dist.results$estimate[2]) + (eall.Phi.dot.p.dist.results$estimate[3]*dist)) %>%
  mutate(p_lcl_logit = eall.Phi.dot.p.dist.results$lcl[2] + (eall.Phi.dot.p.dist.results$lcl[3]*dist)) %>%
  mutate(p_ucl_logit = eall.Phi.dot.p.dist.results$ucl[2] + (eall.Phi.dot.p.dist.results$ucl[3]*dist)) %>%
  mutate(p = logit_recip(p_logit)) %>%
  mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
  mutate(p_ucl = logit_recip(p_ucl_logit)) 

###### Survival varies by tagging size, p constant (eall.Phi.size.p.dot)
minsize = min(eall$tag_size) #right now, this is 1.6 (b/c of a typo.... Michelle fixed - need to re-call data from db and run distance calcs again)
maxsize = max(eall$tag_size)
size.values = minsize+(0:30)*(maxsize-minsize)/30
Phibysize = covariate.predictions(eall.Phi.size.p.dot,data=data.frame(tag_size=size.values),indices=c(1)) #should do this the way I did distance too, make sure get the same thing...

###### Survival varies by tail color, p constant (eall.Phi.color.p.dot)
eall.color = as.data.frame(eall.Phi.color.p.dot$results$real) %>%
  mutate(param = c(rep("Phi",6),"p")) %>%
  mutate(color = c("BW","O","W","Y","YP","YR","all")) %>%
  mutate(row = seq(1,7,1)) #makes plotting easier

###### Survival varies by tagging size, p varies by distance (eall.Phi.size.p.dist)
minsize = min(eall$tag_size) #right now, this is 2.6 (prob still a typo? Michelle fixed the 1.6 one but are there more that are too small?)
maxsize = max(eall$tag_size)
size.values = minsize+(0:30)*(maxsize-minsize)/30
Phibysize = covariate.predictions(eall.Phi.size.p.dist,data=data.frame(tag_size=size.values),indices=c(1)) #should do this the way I did distance too, make sure get the same thing...

mindist = min(eall$dist2016,eall$dist2017,eall$dist2018) #this is 0, obv...
maxdist1 = max(eall$dist2016,eall$dist2017,eall$dist2018) #seems kind of high and probably due to an error
maxdist2 = 200 #encompasses highest values in 2016, 2018 
dist.values1 = mindist+(0:30)*(maxdist1-mindist)/30
dist.values2 = mindist+(0:30)*(maxdist2-mindist)/30

eall.Phi.size.p.dist.results = as.data.frame(eall.Phi.size.p.dist$results$beta) 

# Size
pbydistPhisize_size= data.frame(tag_size = size.values) %>%
  mutate(Phi_logit = (eall.Phi.size.p.dist.results$estimate[1]) + (eall.Phi.size.p.dist.results$estimate[2]*tag_size)) %>%
  mutate(Phi_lcl_logit = eall.Phi.size.p.dist.results$lcl[1] + (eall.Phi.size.p.dist.results$lcl[2]*tag_size)) %>%
  mutate(Phi_ucl_logit = eall.Phi.size.p.dist.results$ucl[1] + (eall.Phi.size.p.dist.results$ucl[2]*tag_size)) %>%
  mutate(Phi = logit_recip(Phi_logit)) %>%
  mutate(Phi_lcl = logit_recip(Phi_lcl_logit)) %>%
  mutate(Phi_ucl = logit_recip(Phi_ucl_logit)) 

# Large range of distances (up to 5000m)
pbydistPhisize_longrange= data.frame(dist = dist.values1) %>%
  mutate(p_logit = (eall.Phi.size.p.dist.results$estimate[3]) + (eall.Phi.size.p.dist.results$estimate[4]*dist)) %>%
  mutate(p_lcl_logit = eall.Phi.size.p.dist.results$lcl[3] + (eall.Phi.size.p.dist.results$lcl[4]*dist)) %>%
  mutate(p_ucl_logit = eall.Phi.size.p.dist.results$ucl[3] + (eall.Phi.size.p.dist.results$ucl[4]*dist)) %>%
  mutate(p = logit_recip(p_logit)) %>%
  mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
  mutate(p_ucl = logit_recip(p_ucl_logit)) 

# Short range of distances (up to 200m)
pbydistPhisize_shortrange = data.frame(dist = dist.values2) %>%
  mutate(p_logit = (eall.Phi.size.p.dist.results$estimate[3]) + (eall.Phi.size.p.dist.results$estimate[4]*dist)) %>%
  mutate(p_lcl_logit = eall.Phi.size.p.dist.results$lcl[3] + (eall.Phi.size.p.dist.results$lcl[4]*dist)) %>%
  mutate(p_ucl_logit = eall.Phi.size.p.dist.results$ucl[3] + (eall.Phi.size.p.dist.results$ucl[4]*dist)) %>%
  mutate(p = logit_recip(p_logit)) %>%
  mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
  mutate(p_ucl = logit_recip(p_ucl_logit)) 

###### Survival varies by site, p constant (eall.Phi.site.p.dot)
eall.Phisite = as.data.frame(eall.Phi.site.p.dot$results$real) %>%
  mutate(param = c(rep("Phi",15),"p")) %>%
  mutate(site = c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
                 "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag","all")) %>%
  mutate(row = seq(1,16,1))

###### Survival varies by site, p by distance (eall.Phi.site.p.dist) - doesn't look like this one parsed distance correctly?
eall.Phisitepdist_Phi = as.data.frame(eall.Phi.site.p.dist$results$real[1:15,]) %>%
  mutate(param = c(rep("Phi",15))) %>%
  mutate(site = c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
                  "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag")) %>%
  mutate(row = seq(1,15,1))

###### Survival varies by site, p varies by site (eall.Phi.site.p.site)
eall.site = as.data.frame(eall.Phi.site.p.site$results$beta) %>%
  mutate(param = c(rep("Phi",15),rep("p",15))) %>% #add a parameter column to make it easier to plot in ggplot
  mutate(site = c(rep(c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
                        "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag"),2))) #add site column (Cabatoan is intercept)

eall.site.real = as.data.frame(eall.Phi.site.p.site$results$real) %>%
  mutate(param = c(rep("Phi",15),rep("p",15))) %>% #add a parameter column to make it easier to plot in ggplot
  mutate(site = c(rep(c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
                        "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag"),2))) %>% #add site column (Cabatoan is intercept)
  mutate(row = seq(1,30,1)) #add row numbers to make plotting easier

# make some plots
###### Constant survival and recapture prob, no covariates (eall.Phi.dot.p.dot)
pdf(file = here("Plots/PhiandpEstimates", "eall_Phiandp_constant.pdf"))
ggplot(data = eall.constant, aes(param, est, color=param)) +
  geom_pointrange(aes(ymin=lower, ymax=upper), size=1) +
  xlab("parameter") + ylab("estimate") + ggtitle("Constant p and Phi (eall.Phi.dot.p.dot)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

###### Survival constant, p varies by distance (eall.Phi.dot.p.dist)
# all distances included, up to 5000m from anem (even though large ones are probably an error...)
pdf(file = here("Plots/PhiandpEstimates", "eall_Phidot_pdist_alldists.pdf"))
ggplot(data = pbydist_longrange, aes(dist, p)) +
  geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
  geom_line(color="black") +
  xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("Constant Phi, p by distance (eall.Phi.dot.p.dist))") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

# most distances included (all from 2016, 2018, not some of the higher ones from 2017)
pdf(file = here("Plots/PhiandpEstimates", "eall_Phidot_pdist_shortdists.pdf"))
ggplot(data = pbydist_shortrange, aes(dist, p)) +
  geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
  geom_line(color="black") +
  xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("Constant Phi, p by distance (eall.Phi.dot.p.dist))") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

###### Survival varies by tagging size, p constant (eall.Phi.size.p.dot)
pdf(file = here("Plots/PhiandpEstimates", "eall_Phisize.pdf"))
ggplot(data = Phibysize$estimates, aes(covdata, estimate)) +
  geom_ribbon(aes(ymin=lcl,ymax=ucl),color="light blue",fill="light blue") +
  geom_line(color="black") +
  xlab("size at tagging (cm)") + ylab("Phi estimate") + ggtitle("Phi varies by tagging size, p constant (eall.Phi.size.p.dot)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

###### Survival varies by tail color, p constant (eall.Phi.color.p.dot)
pdf(file = here("Plots/PhiandpEstimates", "eall_Phicolor.pdf"))
ggplot(data = eall.color, aes(row, estimate, color=color, shape=param)) +
  geom_pointrange(aes(ymin=lcl, ymax=ucl), size=1) +
  xlab("parameter") + ylab("estimate") + ggtitle("Phi by tail color, p constant (eall.Phi.color.p.dot)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

###### Survival varies by tagging size, p varies by distance (eall.Phi.size.p.dist)
# all distances included, up to 5000m from anem (even though large ones are probably an error...)
pdf(file = here("Plots/PhiandpEstimates", "eall_Phisize_pdist_alldists.pdf"))
ggplot(data = pbydistPhisize_longrange, aes(dist, p)) +
  geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
  geom_line(color="black") +
  xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("p: Phi by size, p by distance (eall.Phi.size.p.dist))") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

# most distances included (all from 2016, 2018, not some of the higher ones from 2017)
pdf(file = here("Plots/PhiandpEstimates", "eall_Phisize_pdist_shortdists.pdf"))
ggplot(data = pbydistPhisize_shortrange, aes(dist, p)) +
  geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
  geom_line(color="black") +
  xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("p: Phi by size, p by distance (eall.Phi.size.p.dist))") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

# Phi by tag_size
pdf(file = here("Plots/PhiandpEstimates", "eall_Phisize_pdist_Phi.pdf"))
ggplot(data = pbydistPhisize_size, aes(tag_size, Phi)) +
  geom_ribbon(aes(ymin=Phi_lcl,ymax=Phi_ucl),color="light blue",fill="light blue") +
  geom_line(color="black") +
  xlab("size at tagging (cm)") + ylab("Phi estimate") + ggtitle("Phi: Phi by size, p by distance (eall.Phi.size.p.dist))") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

###### Survival varies by site, p constant (eall.Phi.site.p.dot)
pdf(file = here("Plots/PhiandpEstimates", "eall_Phisite.pdf"))
ggplot(data = eall.Phisite, aes(row, estimate, color=site, shape=param)) +
  geom_pointrange(aes(ymin=lcl, ymax=ucl), size=1) +
  xlab("parameter") + ylab("Phi estimate") + ggtitle("Phi by site, p constant (eall.Phi.site.p.dot)") +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

###### Survival varies by site, p by distance (eall.Phi.site.p.dist) (doesn't look like this one did distance right - has a couple of params for it...)
pdf(file = here("Plots/PhiandpEstimates", "eall_Phisitepdist_Phi.pdf"))
ggplot(data = eall.Phisitepdist_Phi, aes(row, estimate, color=site)) +
  geom_pointrange(aes(ymin=lcl, ymax=ucl), size=1) +
  xlab("parameter") + ylab("Phi estimate") + ggtitle("Phi: Phi by site, p by dist (eall.Phi.site.p.dist)") +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

###### Survival varies by site, p varies by site (eall.Phi.site.p.site)
pdf(file = here("Plots/PhiandpEstimates", "eall_Phiandp_site.pdf"))
ggplot(data = eall.site.real, aes(row, estimate, color=site, shape=param)) +
  geom_pointrange(aes(ymin=lcl, ymax=ucl), size=1) +
  xlab("parameter") + ylab("estimate") + ggtitle("Phi and p by site (eall.Phi.site.p.site)") +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
dev.off()

###### Comparing Phi and p across models

#################### Saving files and output ####################
save(encounters_all, file=here("Data", "encounters_all.RData"))
save(encounters_means, file=here("Data", "encounters_means.RData"))
save(encounters_0, file=here("Data", "encounters_0.RData"))
save(allfish, file=here("Data", "allfish.RData"))

####### Recap on what KC needs:


#### TO DOS
# Do runs with size for Phi, site for Phi?, dist for p
# Make plots to look at those
# Compare models (AIC? goodness of fit?) and summarize in some way
# Check the outlier distances, re-run with any errors fixed there and a fresh pull from the db
# Move on to metrics - write-up and placeholder script for those!


#################### OLD CODE - BOTH MARK AND NOT ####################
######## ORIGINAL WORK WITH SIZE AS A COVARIATE
# ##### Try running a MARK model with size-at-tagging as an individual constraint, site as grouping
# site_size_at_tagging <- function() {
#   # select relevant data
#   e2 <- encounters_all %>%
#     select(ch, tag_size, site) %>%
#     mutate(site = as.factor(site))
#   
#   e2_2 <- e2[complete.cases(e2),] #for using tag size, need no NAs...
# 
#   
#   # process data, make ddl
#   e2.processed = process.data(e2, begin.time=2015, groups="site")
#   e2.ddl = make.design.data(e2.processed)
#   
#   # process data, make ddl for complete cases
#   e2_2.processed = process.data(e2_2, begin.time=2015, groups="site")
#   e2_2.ddl = make.design.data(e2_2.processed)
#   
#   
#   # define models for Phi and p
#   Phi.dot = list(formula=~1, link="logit") #one survival for all sizes and sites and times
#   Phi.site = list(formula=~site, link="logit") #survival depends on site
#   Phi.size = list(formula=~tag_size, link="logit") #survival depends on size
#   Phi.size.x.site = list(formula=~tag_size*site, link="logit") #survival depends on site and site-dependent slope for size
#   Phi.time = list(formula=~time, link="logit") #survival varies with time, same for all sites
#   Phi.time.plus.site = list(formula=~time+site, link="logit")
#   
#   p.dot = list(formula=~1, link="logit") #one recap prob for all sizes, sites, and time
#   p.site = list(formula=~site, link="logit")
#   p.time = list(formula=~time, link="logit")
#   
#   # create model list
#   #e2.cml = create.model.list("CJS")
#   
#   #run and return models
#   #return(mark.wrapper(e2.cml, data=e2.processed, ddl=e2.ddl))
#   
#   #create models
#   e2.Phi.dot.p.dot = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
#   e2.Phi.site.p.dot = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.site, p=p.dot))
#   e2.Phi.size.p.dot = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.size, p=p.dot))
#   e2.Phi.size.p.site = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.size, p=p.site))
#   e2.Phi.size.x.site.p.site = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.size.x.site, p=p.site))
#   e2.Phi.site.p.site = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.site, p=p.site))
#   
#   e2_2.Phi.size.p.dot = mark(e2_2.processed, e2_2.ddl, model.parameters=list(Phi=Phi.size, p=p.dot))
#   e2_2.Phi.size.p.site = mark(e2_2.processed, e2_2.ddl, model.parameters=list(Phi=Phi.size, p=p.site))
#   e2_2.Phi.size.x.site.p.site = mark(e2_2.processed, e2_2.ddl, model.parameters=list(Phi=Phi.size.x.site, p=p.site))
#   
# }
# 
# e2.results = site_size_at_tagging()
# 
# 
# # Make some plots
# # e2.Phi.dot.p.dot (constant survival and recapture prob the whole time)
# e2.constant = as.data.frame(e2.Phi.dot.p.dot$results$beta) %>% 
#   mutate(param = c("Phi","p")) %>% #add a parameters column to make it easier to plot in ggplot
#   mutate(upper = logit_recip(ucl), lower = logit_recip(lcl), est = logit_recip(estimate)) #do the reciprocal transform on upper and lower confidence limits (need to check what those are - 95? SE? and that transforming them just straight up is the right way to go)
#   
# e2.site = as.data.frame(e2.Phi.site.p.site$results$beta) %>%
#   mutate(param = c(rep("Phi",15),rep("p",15))) %>% #add a parameter column to make it easier to plot in ggplot
#   mutate(site = c(rep(c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
#                         "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag"),2))) #add site column (Cabatoan is intercept)
# 
# e2.site.real = as.data.frame(e2.Phi.site.p.site$results$real) %>%
#   mutate(param = c(rep("Phi",15),rep("p",15))) %>% #add a parameter column to make it easier to plot in ggplot
#   mutate(site = c(rep(c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
#                         "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag"),2))) %>% #add site column (Cabatoan is intercept)
#   mutate(row = seq(1,30,1)) #add row numbers to make plotting easier
# 
# #Think this is handled for now, should check again to be sure...
# #NEED TO ADD THE EFFECT OF EACH SITE TO THE INTERCEPT
# #### NEED TO SWITCH IT BACK FROM THE LOGIT SCALE FIRST!! AND THAT AFFECTS THE CONFIDENCE INTERVALS, RIGHT? CHECK THAT AGAIN....
# #### ALSO SHOULD PLOT +- SE TOO, SEE HOW THAT COMPARES/MATCHES UP TO LCL and UCL
# 
# # Constant Phi and p across time and site, no covariates (e2.Phi.dot.p.dot)
# pdf(file = here("Plots/PhiandpEstimates", "e2_Phiandp_constant.pdf"))
# ggplot(data = e2.constant, aes(param, est, color=param)) +
#   geom_pointrange(aes(ymin=lower, ymax=upper), size=1) +
#   xlab("parameter") + ylab("estimate") + ggtitle("e2.Phi.dot.p.dot (constant p and Phi)") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme_bw()
# dev.off()
#   
# # Constant Phi and p across time, variable across site, no covariates (e2.Phi.site.p.site)
# pdf(file = here("Plots/PhiandpEstimates", "e2_Phiandp_site.pdf"))
# ggplot(data = e2.site.real, aes(row, estimate, color=site, shape=param)) +
#   geom_pointrange(aes(ymin=lcl, ymax=ucl), size=1) +
#   xlab("parameter") + ylab("estimate") + ggtitle("e2.Phi.site.p.site (p and Phi by site)") +
#   #scale_y_continuous(limits = c(0, 1)) +
#   theme_bw()
# dev.off()
# 
# # Phi has tagging size as a covariate, not dependent on site or time, p constant (e2_2.Phi.size.p.dot)
# minsize = min(e2_2$tag_size) #right now, this is 1.6 (prob b/c of a typo.... Michelle is looking at it)
# maxsize = max(e2_2$tag_size)
# size.values = minsize+(0:30)*(maxsize-minsize)/30
# Phibysize = covariate.predictions(e2_2.Phi.size.p.dot,data=data.frame(tag_size=size.values),indices=c(1))
# 
# pdf(file = here("Plots/PhiandpEstimates", "e2_2Phisize.pdf"))
# ggplot(data = Phibysize$estimates, aes(covdata, estimate)) +
#   geom_ribbon(aes(ymin=lcl,ymax=ucl),color="light blue",fill="light blue") +
#   geom_line(color="black") +
#   xlab("size at tagging (cm)") + ylab("Phi estimate") + ggtitle("e2_2.Phi.size.p.dot") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme_bw()
# dev.off()

###### ORIGINAL WORK WITH DISTANCE AS A TIME-VARYING INDIVIDUAL CONSTRAINT
# ##### Try running time-varying individual constraints - distance
# # data
# edist <- fish.Tagged %>%
#   select(ch, site, dist_2016, dist_2017, dist_2018) %>% #select just site and distance (not including 2015, first potential tagging year)
#   rename(dist2016 = dist_2016, dist2017 = dist_2017, dist2018 = dist_2018) %>% 
#   mutate(dist2016 = replace(dist2016, is.na(dist2016), 0)) %>% #replace NAs (before a fish was tagged) with 0s
#   mutate(dist2017 = replace(dist2017, is.na(dist2017), 0)) %>%
#   mutate(dist2018 = replace(dist2018, is.na(dist2018), 0))
# 
# # process data and make ddl
# edist.processed = process.data(edist, model="CJS", begin.time=2015)
# edist.ddl = make.design.data(edist.processed)
# 
# # set models for Phi
# Phi.dot = list(formula=~1, link="logit")
# Phi.time = list(formula=~time, link="logit")
# 
# # set models for p
# p.dot = list(formula=~1, link="logit")
# p.time = list(formula=~time, link="logit")
# p.dist = list(formula=~dist, link="logit")
# p.time.plus.dist = list(formula=~time+dist, link="logit")
# 
# # run some models
# edist.Phi.dot.p.dot = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
# edist.Phi.dot.p.time = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.dot, p=p.time)) #has issues running, can't estimate p in 2018
# edist.Phi.dot.p.dist = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist)) #think this ran!!!
# # edist.Phi.dot.p.dist.plus.time = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.dot, p=p.time.plus.dist))
# # edist.Phi.time.p.dot = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.dot))
# # edist.Phi.time.p.time = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.time))
# # edist.Phi.time.p.dist = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.dist))
# # edist.Phi.time.p.time.plus.dist = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.time.plus.dist))
# 
# # make some plots
# 
# # Phi constant, not dependent on site or time, p depends on distance to anem (edist.Phi.dot.p.dist)
# mindist = min(edist$dist2016,edist$dist2017,edist$dist2018) 
# maxdist = max(edist$dist2016,edist$dist2017,edist$dist2018) #seems kind of high and probably due to an error
# maxdist2 = 400 #encompasses all but the highest value (which is an order of magnitude higher - 4874 compared to next highest of 352, both in 2017)
# maxdist3 = 200 #encompasses highest values in 2016, 2018 
# #tail(sort(edist$dist2016),10)
# #tail(sort(edist$dist2017),10) #really 2017 that has the high distances...
# #tail(sort(edist$dist2018),10)
# dist.values = mindist+(0:30)*(maxdist-mindist)/30
# dist.values2 = mindist+(0:30)*(maxdist2-mindist)/30
# dist.values3 = mindist+(0:30)*(maxdist3-mindist)/30
# 
# # Not getting this to work - it doesn't change with distance... maybe/probably I'm specifying the index incorrectly?
# #pbydist1 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values),indices=c(7))
# #pbydist2 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values2),indices=c(1))
# #pbydist3 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values3),indices=c(1))
# 
# edist.Phi.dot.p.dist.results = as.data.frame(edist.Phi.dot.p.dist$results$beta) 
# 
# # Large range of distances (up to 5000m)
# pbydist_range1 = data.frame(dist = dist.values) %>%
#   mutate(p_logit = (edist.Phi.dot.p.dist.results$estimate[2]) + (edist.Phi.dot.p.dist.results$estimate[3]*dist)) %>%
#   mutate(p_lcl_logit = edist.Phi.dot.p.dist.results$lcl[2] + (edist.Phi.dot.p.dist.results$lcl[3]*dist)) %>%
#   mutate(p_ucl_logit = edist.Phi.dot.p.dist.results$ucl[2] + (edist.Phi.dot.p.dist.results$ucl[3]*dist)) %>%
#   mutate(p = logit_recip(p_logit)) %>%
#   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#   mutate(p_ucl = logit_recip(p_ucl_logit)) 
# 
# # Mid range of distances (up to 400m)
# pbydist_range2 = data.frame(dist = dist.values2) %>%
#   mutate(p_logit = (edist.Phi.dot.p.dist.results$estimate[2]) + (edist.Phi.dot.p.dist.results$estimate[3]*dist)) %>%
#   mutate(p_lcl_logit = edist.Phi.dot.p.dist.results$lcl[2] + (edist.Phi.dot.p.dist.results$lcl[3]*dist)) %>%
#   mutate(p_ucl_logit = edist.Phi.dot.p.dist.results$ucl[2] + (edist.Phi.dot.p.dist.results$ucl[3]*dist)) %>%
#   mutate(p = logit_recip(p_logit)) %>%
#   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#   mutate(p_ucl = logit_recip(p_ucl_logit)) 
#   
# # Short range of distances (up to 200m)
# pbydist_range3 = data.frame(dist = dist.values3) %>%
#   mutate(p_logit = (edist.Phi.dot.p.dist.results$estimate[2]) + (edist.Phi.dot.p.dist.results$estimate[3]*dist)) %>%
#   mutate(p_lcl_logit = edist.Phi.dot.p.dist.results$lcl[2] + (edist.Phi.dot.p.dist.results$lcl[3]*dist)) %>%
#   mutate(p_ucl_logit = edist.Phi.dot.p.dist.results$ucl[2] + (edist.Phi.dot.p.dist.results$ucl[3]*dist)) %>%
#   mutate(p = logit_recip(p_logit)) %>%
#   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#   mutate(p_ucl = logit_recip(p_ucl_logit)) 
# 
# # Long dist values (up to 5000m from anem)
# pdf(file = here("Plots/PhiandpEstimates", "edist_Phidot_pdist_distvalues1.pdf"))
# ggplot(data = pbydist_range1, aes(dist, p)) +
#   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#   geom_line(color="black") +
#   xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("edist.Phi.dot.p.dist (constant Phi, p by dist)") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme_bw()
# dev.off()
# 
# # Mid range dist values (up to 400m from anem)
# pdf(file = here("Plots/PhiandpEstimates", "edist_Phidot_pdist_distvalues2.pdf"))
# ggplot(data = pbydist_range2, aes(dist, p)) +
#   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#   geom_line(color="black") +
#   xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("edist.Phi.dot.p.dist (constant Phi, p by dist)") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme_bw()
# dev.off()
# 
# # Short range dist values (up to 200m from anem)
# pdf(file = here("Plots/PhiandpEstimates", "edist_Phidot_pdist_distvalues3.pdf"))
# ggplot(data = pbydist_range3, aes(dist, p)) +
#   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#   geom_line(color="black") +
#   xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("edist.Phi.dot.p.dist (constant Phi, p by dist)") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme_bw()
# dev.off()


##### Try running MARK (again), this time with time-varying individual constraints!
# Trying one without site as a group, just with distance to anem as the time-varying individual constraint
# data
e1 <- encounters_means %>% 
  select(ch, dist2016, dist2017, dist2018)

# process data and make ddl
fish.processed = process.data(e1, model="CJS", begin.time=2015)
fish.ddl = make.design.data(fish.processed)

# set models for Phi 
Phi.dot = list(formula=~1, link="logit")
Phi.time = list(formula=~time, link="logit")

# set models for p
p.dot = list(formula=~1, link="logit")
p.time = list(formula=~time, link="logit")
p.dist = list(formula=~dist, link="logit")
p.dist.plus.time = list(formula=~dist+time, link="logit")

# run models individually
e1.Phi.dot.p.dot = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
e1.Phi.dot.p.time = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.dot, p=p.time))
e1.Phi.dot.p.dist = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist))
e1.Phi.dot.p.dist.plus.time = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist.plus.time))
e1.Phi.time.p.dot = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.time, p=p.dot))
e1.Phi.time.p.time = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.time, p=p.time))
e1.Phi.time.p.dist = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.time, p=p.dist))
e1.Phi.time.p.dist.plus.time = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.time, p=p.dist.plus.time))

# run models using clm and wrapper
e1.cml = create.model.list("CJS")
e1.model.wrap = mark.wrapper(e1.cml, data=fish.processed, ddl=fish.ddl)

e1.model.wrap[[1]]$design.matrix # this is one of the ones with distance... not convinced it worked right..
e1.model.wrap[[2]]$design.matrix #
e1.model.wrap[[5]]$design.matrix

e1.collect = collect.models(table=TRUE)
# make plots to compare the models

# compare using AIC?

# run models using site as a group

# run models with tagging size as covariate

# run models with actual size as a covariate (can fill in missing data for un-recaught fish using Michelle's growth? or just assume it stays the same if not ever caught again, mean between two sizes on either size if one caputure missed, not sure what would do for fish pre-capture)

# run models with distance where NAs were replaced with 0 and where actual distance to anem was calculated whether or not fish had been tagged yet, compare all three methods to make sure it doesn't matter

# make some sort of data frame with estimate, 95% confidence bounds, with different assumptions (time, distance, site, etc. matters) shown, so can use as want in persistence metrics and analyses

##### Now, try running MARK! (original attempt)
# specifying site as group
# make process and design data
#with site
tagged.site.process = process.data(encounters_all, model="CJS", begin.time=2015, groups="site")
tagged.site.ddl = make.design.data(tagged.site.process)

#without site
tagged.process = process.data(encounters_all, model="CJS", begin.time=2015)
tagged.ddl = make.design.data(tagged.process)

# set different relationships
Phi.dot = list(formula=~1)
Phi.time = list(formula=~time)
Phi.site = list(formula=~site)
Phi.siteplustime = list(formua=~site+time)

p.dot = list(formula=~1)
p.time = list(formula=~time)
p.site = list(formula=~site)
p.siteplustime = list(formua=~site+time) # does this still pool sites at some level, assuming they come from the same dist and are related (like would in rethinking) or just break them out completely separately?

# run some models
tagged.phi.dot.p.dot = mark(tagged.process, tagged.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
tagged.site.phi.dot.p.dot = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
tagged.site.phi.time.p.dot = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.time, p=p.dot))
tagged.site.phi.dot.p.time = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.dot, p=p.time))
tagged.site.phi.time.p.time = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.time, p=p.time))
tagged.site.phi.site.p.site = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.site, p=p.site))
tagged.site.phi.time.p.siteplustime = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.site, p=p.siteplustime))

#these don't run... must have something wrong? 
tagged.site.phi.siteplustime.p.time = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.siteplustime, p=p.time))
tagged.site.phi.siteplustime.p.siteplustime = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.siteplustime, p=p.siteplustime))

#tagged.site <- mark(data=encounters_all, model = "CJS")

# Using the output
## Note from MARK book - the Beta estimates at the top relate to the link function - thethey're the values actually estimated before the actual paramter values are re-constituted; they're estimates on the logit scale
## So for survival and recap, the Beta estimates are the logit values (ln(p/(1-p))), for linear model constraints, they are the coefficients (slopes) for each term in the linear model
## See pp.226-229 in MARK book for details
#################################### Testing out doing mark-recap in rethinking ############# 




#################################### Old stuff below ############# 

##### Previous work using own likelihood functions...
# Now, create encounter histories for fish at all sites combined
encountersAll <- CreateEncounterSummary(2015, 2017, allfish) #simpler way of getting encounter history like above using function
nvals_allfish <- as.data.frame(table(encountersAll$encounter.hist))

# Try gradient descent on this
params0 <- c(0.5, 0.5, 0.5, 0.5) #s1, s2, p2, p3
params0_2 <- c(0.5, 0.5) #s, p
highVal <- 0.999
lowVal <- 0.001
eps <- 0.001
del <- 0.001
p_vec <- seq(0.01, 0.99, 0.01) #recapture prob
s_vec <- seq(0.01, 0.99, 0.01) #survival prob
Val1 <- c("001","010","011","100","101","110","111")


all_GradDesc <- GradientDescent3SY4Params(nvals_allfish, params0, 50, 0.001, 0.001, 0.999, 0.001) #works! except param 2 and param 4 are exactly the same...
#param1 (s1) = 0.15, param2 (s2) = 0.36, param3 (p2) = 0.35, param4 (p3) = 0.36, like = -345.4786
all_GradDesc_1000tsteps <- GradientDescent3SY4Params(nvals_allfish, params0, 1000, 0.001, 0.001, 0.999, 0.001) #param 2 and 4 stll exactly the same... seems weird
all_GradDesc_1000tsteps_2p <- GradientDescent3SY2Params(nvals_allfish, params0_2, 1000, 0.001, 0.001, 0.999, 0.001) #param 2 and 4 stll exactly the same... seems weird

# Now try for different sites
# Cabatoan
encountersCabatoan <- CreateEncounterSummary(2015, 2017, (allfish %>% filter(site == site_vec[Cabatoan]))) #simpler way of getting encounter history like above using function
nvals_Cabatoan <- as.data.frame(table(encountersCabatoan$encounter.hist))
nvals_Cabatoan_clean <- c(0, 15, 0, 31, 0, 0, 0)
nvals_Cabat <- as.data.frame(Val1)
nvals_Cabat$Freq <- nvals_Cabatoan_clean

Cabat_ParamRangeScan <- ParamRangeCheck3SY2P_dfoutput(s_vec, p_vec, nvals_Cabat)
smaxCabat <- Cabat_ParamRangeScan$survival[which.max(Cabat_ParamRangeScan$like)]
pmaxCabat <- Cabat_ParamRangeScan$capture[which.max(Cabat_ParamRangeScan$like)]
Cabat_GradDesc_1000tsteps <- GradientDescent3SY4Params(nvals_Cabat, params0, 1000, 0.001, 0.001, 0.999, 0.001) #param 2 and 4 stll exactly the same... seems weird


N010 <- 2
N011 <- 3
N100 <- 4
N101 <- 5
N110 <- 6
N111 <- 7
# Try gradient descent on this
params0 <- c(0.5, 0.5, 0.5, 0.5) #s1, s2, p2, p3
highVal <- 0.999
lowVal <- 0.001
eps <- 0.001
del <- 0.001
all_GradDesc <- GradientDescent3SY4Params(nvals_allfish, params0, 50, 0.001, 0.001, 0.999, 0.001) #works! except param 2 and param 4 are exactly the same...
#param1 (s1) = 0.15, param2 (s2) = 0.36, param3 (p2) = 0.35, param4 (p3) = 0.36, like = -345.4786
all_GradDesc_1000tsteps <- GradientDescent3SY4Params(nvals_allfish, params0, 1000, 0.001, 0.001, 0.999, 0.001) #param 2 and 4 stll exactly the same... seems weird


# Caridad Cemetery

# Caridad Proper

# Elementary School

# Gabas

# Haina

# Hicgop South

# Magbangon

# Palanas

# Poroc Rose

# Poroc San Flower


# Pop size
Pal <- allfish %>% filter(site == "Palanas")
Pal_pop <- as.data.frame(table(Pal$year))

Wan <- allfish %>% filter(site == "Wangag")
Wan_pop <- as.data.frame(table(Wan$year))

SB <- allfish %>% filter(site == "Sitio Baybayon")
SB_pop <- as.data.frame(table(SB$year))

Vis <- allfish %>% filter(site == "Visca")
Vis_pop <- as.data.frame(table(Vis$year))



#################################### R MARK work ############# 
# need to get data in format for R Mark
# column ch - has cenounter history (like 1011001)
# then columns with characteristics - should have tag_id, size (at first capture), site, could have distance to any anem it was caught it to that year's track for each year...
#encounters[[i]] <- tagged.fish %>% group_by(tag_id) %>% summarise(!!var.name := ifelse(sum(year == sample.years[i])>0, 1, 0)) #create data frames for each year that have a vector of tag ids and a vector of encountered (1) or didn't (0) in 


#Load practice data sets
data(dipper)

#Example from the RMark pdf, just to test mark
dipper.Phidot.pdot=mark(dipper,threads=1)

data(example.data)


