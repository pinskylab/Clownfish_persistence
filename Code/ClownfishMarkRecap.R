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


#################### Running things: ####################
leyte <- read_db("Leyte")

# fish <- leyte %>%
#   tbl("clownfish")
# 
# anem <- leyte %>%
#   tbl("anemones")
# 
# dives <- leyte %>%
#   tbl("diveinfo")
# save(fish, file='fish_db_05-07-2018')
# save(anem, file='anem_db_05-07-2018')
# save(dives, file='dives_db_05-07-2018')

# # Pull out all the tagged fish (here, using saved database tables, not leyte)
# allfish_fish <- leyte %>% 
#   tbl("clownfish") %>%
#   select(fish_table_id, anem_table_id, fish_spp, sample_id, cap_id, anem_table_id, recap, tag_id, color, size) %>%
#   collect() %>%
#   filter(!is.na(tag_id)) 
# 
# allfish_anems <- leyte %>%
#   tbl("anemones") %>%
#   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
#   collect() %>%
#   filter(anem_table_id %in% allfish_fish$anem_table_id)
# 
# allfish_dives <- leyte %>%
#   tbl("diveinfo") %>%
#   select(dive_table_id, dive_type, date, site, gps) %>%
#   collect() %>%
#   filter(dive_table_id %in% allfish_anems$dive_table_id)
# 
# # pull out just the year and put that in a separate column
# allfish_dives$year <- as.integer(substring(allfish_dives$date,1,4))
# 
# #join together
# allfish <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
# allfish <- left_join(allfish, allfish_dives, by="dive_table_id")
# 
# # 
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

# Add in tail color (could probably do these all together in one command... look into that later..)
tail_color_2015 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  filter(year == 2015) %>%
  group_by(tag_id) %>%
  summarize(tail_color_2015 = color[1])

tail_color_2016 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  filter(year == 2016) %>%
  group_by(tag_id) %>%
  summarize(tail_color_2016 = color[1])

tail_color_2017 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  filter(year == 2017) %>%
  group_by(tag_id) %>%
  summarize(tail_color_2017 = color[1])

tail_color_2018 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  filter(year == 2018) %>%
  group_by(tag_id) %>%
  summarize(tail_color_2018 = color[1])

#encounters[[i]] <- tagged.fish %>% group_by(tag_id) %>% summarise(!!var.name := ifelse(sum(year == sample.years[i])>0, 1, 0)) #create data frames for each year that have a vector of tag ids and a vector of encountered (1) or didn't (0) in 

encounters_all <- left_join(encounters_all, tail_color_2015, by="tag_id")
encounters_all <- left_join(encounters_all, tail_color_2016, by="tag_id")
encounters_all <- left_join(encounters_all, tail_color_2017, by="tag_id")
encounters_all <- left_join(encounters_all, tail_color_2018, by="tag_id")

# Add in size in each of the years
size_2015 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  filter(year == 2015) %>%
  group_by(tag_id) %>% 
  summarize(size_2015 = mean(size, rm.na = TRUE))

size_2016 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  filter(year == 2016) %>%
  group_by(tag_id) %>% 
  summarize(size_2016 = mean(size, rm.na = TRUE))

size_2017 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  filter(year == 2017) %>%
  group_by(tag_id) %>% 
  summarize(size_2017 = mean(size, rm.na = TRUE))

size_2018 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
  filter(year == 2018) %>%
  group_by(tag_id) %>% 
  summarize(size_2018 = mean(size, rm.na = TRUE))

encounters_all <- left_join(encounters_all, size_2015, by="tag_id")
encounters_all <- left_join(encounters_all, size_2016, by="tag_id")
encounters_all <- left_join(encounters_all, size_2017, by="tag_id")
encounters_all <- left_join(encounters_all, size_2018, by="tag_id")

save(encounters_all, file=here("Data", "encounters_all.RData"))

##### Now, try running MARK!
# specifying site as group
# make process and design data
tagged.site.process = process.data(encounters_all, model="CJS", begin.time=2015, groups="site")
tagged.site.ddl = make.design.data(tagged.site.process)

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
