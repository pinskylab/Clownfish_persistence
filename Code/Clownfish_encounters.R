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

# Fill in missing sizes (either fill in mean for that year or 0 for pre-capture years or project size using growth curve)
projectSize <- function(year, size_year_vec, first_cap_vec, prev_size_vec, mean_year_size, Linf_mean, k_mean) {
  out_size_vec <- size_year_vec
  for(i in 1:length(size_year_vec)) {  # for each fish
    size_test <- size_year_vec[i]  # pull out the size for fish i in this year
    first_cap <- first_cap_vec[i]  # pull out the year of first capture
    prev_size <- prev_size_vec[i]  # pull out the size in the previous year
    
    if(is.na(size_test)) {  # if there isn't a size recorded
      if(first_cap > year) {  # and this is a year before the first capture of this fish
        out_size_vec[i] <- mean_year_size  # fill in the mean size for that year (or 0)
      }
      else {
        out_size_vec[i] <- growthIncrementVBL(Linf_mean, prev_size, k_mean)  # project growth from the previous size
      }
    }
  }
  return(out_size_vec)
}

#################### Running things: ####################
# Pull out all fish marked in some way, either by PIT tag or by gen_id 
marked_fish <- allfish_caught %>%
  filter(!is.na(fish_indiv)) 
#saveRDS(marked_fish, file = here::here("Data/Script_outputs", "marked_fish.RData"))

# Prep data for MARK by making encounter histories
encounters_list <- CreateEncounterSummary(2012, 2018, marked_fish)  # 3053 fish on 10/16/19 - same as distinct fish_indivs in current fish-obs table so that's good

# Find fish traits: site, tail color, size, life stage (all at first time captured), join with encounter histories 
trait_info <- marked_fish %>% 
  group_by(fish_indiv) %>%
  arrange(year) %>%
  summarize(site = site[1],
            first_capture_year = min(year), #year this fish was first captured
            cap_size = size[1], #earlier did min(size) and didn't arrange by year, could go either way
            cap_color = color[1],
            cap_stage = sex[1])

# Join trait info with encounter histories
encounters_list <- left_join(encounters_list, trait_info, by="fish_indiv")

# Find sizes for captured fish, mean size of captured fish in each year, and join size in each year with encounters_all data frame
encounters_size <- encounters_list
fish_capture_sizes <- list()
mean_size <- rep(NA, length(years_sampled))

for(i in 1:length(years_sampled)) {
  
  # Find mean size in year i for each fish
  size_df <- marked_fish %>% 
    #filter(fish_indiv %in% encounters_all$fish_indiv) %>%  # don't think this should pull out any fish, right? hold over from when I was doing tags and genetics separately
    filter(year == years_sampled[i]) %>%
    filter(!is.na(size)) %>% 
    group_by(fish_indiv) %>%
    summarize(fish_size = mean(size, rm.na = TRUE))   
  
  mean_size[i] <- mean(size_df$fish_size, rm.na = TRUE)  # find the mean size of fish for that year
  
  names(size_df)[2] = paste("size", years_sampled[i], sep="")  # rename the size column to be sizeYEAR for MARK (e.g. size2012)
  
  fish_capture_sizes[[i]] <- size_df
  
  encounters_size <- left_join(encounters_size, size_df, by = "fish_indiv")
}

# Mean size across years
#mean_size_overall <- mean(c(encounters_size$size2012, encounters_size$size2013, encounters_size$size2014, encounters_

# Fill in missing sizes with either projected sizes or mean size from each year if pre- first capture year
encounters_size_means_by_year <- encounters_size %>% mutate(size2012 = if_else(is.na(size2012), mean_size[1], size2012))  # replace 2012 size NAs
encounters_size_means_by_year$size2013 <- projectSize(2013, encounters_size$size2013, encounters_size$first_capture_year, encounters_size_means_by_year$size2012, mean_size[2], Linf_mean, k_mean)
encounters_size_means_by_year$size2014 <- projectSize(2014, encounters_size$size2014, encounters_size$first_capture_year, encounters_size_means_by_year$size2013, mean_size[3], Linf_mean, k_mean)
encounters_size_means_by_year$size2015 <- projectSize(2015, encounters_size$size2015, encounters_size$first_capture_year, encounters_size_means_by_year$size2014, mean_size[4], Linf_mean, k_mean)
encounters_size_means_by_year$size2016 <- projectSize(2016, encounters_size$size2016, encounters_size$first_capture_year, encounters_size_means_by_year$size2015, mean_size[5], Linf_mean, k_mean)
encounters_size_means_by_year$size2017 <- projectSize(2017, encounters_size$size2017, encounters_size$first_capture_year, encounters_size_means_by_year$size2016, mean_size[6], Linf_mean, k_mean)
encounters_size_means_by_year$size2018 <- projectSize(2018, encounters_size$size2018, encounters_size$first_capture_year, encounters_size_means_by_year$size2017, mean_size[7], Linf_mean, k_mean)

# Fill in missing sizes with either projected sizes or 0 if pre- first capture year
encounters_size_0 <- encounters_size %>% mutate(size2012 = if_else(is.na(size2012), 0, size2012))  # replace 2012 size NAs
encounters_size_0$size2013 <- projectSize(2013, encounters_size$size2013, encounters_size$first_capture_year, encounters_size_0$size2012, 0, Linf_mean, k_mean)
encounters_size_0$size2014 <- projectSize(2014, encounters_size$size2014, encounters_size$first_capture_year, encounters_size_0$size2013, 0, Linf_mean, k_mean)
encounters_size_0$size2015 <- projectSize(2015, encounters_size$size2015, encounters_size$first_capture_year, encounters_size_0$size2014, 0, Linf_mean, k_mean)
encounters_size_0$size2016 <- projectSize(2016, encounters_size$size2016, encounters_size$first_capture_year, encounters_size_0$size2015, 0, Linf_mean, k_mean)
encounters_size_0$size2017 <- projectSize(2017, encounters_size$size2017, encounters_size$first_capture_year, encounters_size_0$size2016, 0, Linf_mean, k_mean)
encounters_size_0$size2018 <- projectSize(2018, encounters_size$size2018, encounters_size$first_capture_year, encounters_size_0$size2017, 0, Linf_mean, k_mean)


#################### Saving output: ####################
save(marked_fish, file = here::here("Data/Script_outputs", "marked_fish.RData"))  # all obs of fish with a fish_indiv, plus dive, anem, etc. info
save(encounters_list, file = here::here("Data/Script_outputs", "encounters_list.RData"))  # list of encounter histories
save(encounters_size_means, file = here::here("Data/Script_outputs", "encounters_size_means.RData"))  # encounter histories with sizes and missing ones filled in with projections or means
save(encounters_size_0, file = here::here("Data/Script_outputs", "encounters_size_0.RData"))  # encounter histories with sizes and missing ones filled in with projections or 0s