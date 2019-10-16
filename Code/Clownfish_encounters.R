# Prep for mark-recap analysis - creating set of encountered fish, projecting their growth
# Code taken from ClownfishMarkRecap.R and modified

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))

# Load files from new, simple growth analysis (Growth_analysis.R)
load(file = here::here("Data/Script_outputs", "recap_pairs_year.RData"))
load(file = here::here("Data/Script_outputs", "growth_info_estimate.RData"))

# Find mean Linf and K from those runs
Linf_mean <- mean(growth_info_estimate$Linf_est)  # 10.71 (old: 10.58)
k_mean <- mean(growth_info_estimate$k_est)  # 0.864 (old: 0.928)

#################### Functions: ####################

# Calculates second length, t_i is in terms of years - from Fabers section of Hampton et al
Fabers_model <- function(Linf, L_release, K, t_i, e_i) {
  growth <- (Linf - L_release)*(1 - exp(-(K*t_i))) + e_i
  length_out <- L_release + growth
  return(length_out)
}

# Growth-increment VBL
growthIncrementVBL <- function(Linf, L_release, K) {
  length_out = L_release + (Linf - L_release)*(1 - exp(-K))
  return(length_out)
}

# Creates the summarized encounter history by fish_id: output is a data frame with 2 columns - fish_id (either capXX or tagXXXXX or sampleXXX) and summarized encounter history (i.e. 0010)
CreateEncounterSummary <- function(start.year, end.year, tagged.fish) { #start.year is first year a fish could have been tagged (either PIT or genetically), end.year is last year of sampling, tagged.fish is dataframe with all fish with cap_ids or tag_ids
  
  sample.years <- seq(start.year, end.year, 1) #make a vector of years the fish could have been seen
  encounters <- list(); #initialize an empty list to store the various encounter data frames by year
  
  for (i in 1:length(sample.years)) { #pull out encounter vector by tag for each year, store each as a data frame of tag ids and binary encounters in the encounter list
    year.name <- sample.years[i] #get year 
    var.name <- paste("sighted", as.character(sample.years[i]), sep=".") #create dynamic column names for encounters - sighted.[samplingyear]
    
    encounters[[i]] <- tagged.fish %>% #create data frames for each year that have a vector of fish_ids and a vector of encountered (1) or didn't (0) in 
      group_by(fish_indiv) %>% 
      summarise(!!var.name := ifelse(sum(year == sample.years[i])>0, 1, 0)) #encounter history
  }
  
  encounters.out <- as.data.frame(encounters[[1]]) #seed summary data frame with list of tag ids and encounter in 1st year
  colnames(encounters.out)<- c("fish_indiv","encounter.hist") #rename the columns so the encounter.hist column can be pasted to iteratively in next step
  
  for (i in 2:length(sample.years)) {
    encounters.out$encounter.hist <- paste(encounters.out$encounter.hist, encounters[[i]][[2]], sep="") #paste on the other encounter 1/0s so get overall encounter histories as strings
  }
  
  return(encounters.out)
}

#################### Running things: ####################
# Pull out all fish marked in some way, either by PIT tag or by gen_id 
marked_fish <- allfish_caught %>%
  filter(!is.na(fish_indiv)) 
saveRDS(marked_fish, file = here::here("Data/Script_outputs", "marked_fish.RData"))

mf_orig <- readRDS(file = here::here("Data/Script_outputs", "marked_fish.RData"))

fI_1 <- data.frame(mf_orig %>% select(fish_indiv))
fI_2 <- fish_obs %>% select(fish_indiv) %>% mutate(fish_indiv = as.character(fish_indiv))

n_FI_1 <- fI_1 %>% distinct(fish_indiv)  # 3052

test_diff_in1 <- dplyr::setdiff(fI_1, fI_2)
test_diff_in2 <- dplyr::setdiff(fI_2, fI_1)

fI_1_table <- as.data.frame(table(fI_1$fish_indiv))
fI_2_table <- as.data.frame(table(fI_2$fish_indiv))

#################### Plots: ####################

#################### Saving output: ####################