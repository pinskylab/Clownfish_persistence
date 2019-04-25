# Looking at observed population size (and size-structure?) across sampling years - do pops seem to be persistent?

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))

# Load libraries
library(ggplot2)

# Load in site_areas
load(file=here::here('Data','site_areas.RData'))

# Load file with proportion habitat sampled estimates (from TotalAnemsAtSite.R)
#load(file=here("Data",'anem_sampling_table.RData')) #file with anems, after 2018 anems matched

# Load in new proportion habitat sampled file
load(file=here::here("Data/Script_outputs", "anems_visited_by_year.RData"))

# Set parameters (here from EggRecruitRelationship.R)
sample_months <- c(3,4,5,6,7,8) #set months to avoid winter 2015 anem surveys
#dive_list <- c("0","C","D","E")
min_female_size <- 5  # not sure what this should be, just setting it as 5 for now (not based on any particular data, except that when finding breeding size for uncertainty metrics, looked like more fish marked YP and <6cm than I would have expected)

# For time series, only consider sites with more than one sampling year (so not Caridad Proper, Sitio Lonas, or Sitio Tugas)
site_vec_timeseries <- c("Cabatoan", "Caridad Cemetery", "Elementary School", "Gabas", "Haina", "Hicgop South",
              "N. Magbangon", "S. Magbangon", "Palanas", "Poroc Rose", "Poroc San Flower", "San Agustin", "Sitio Baybayon", 
              "Tamakin Dacot", "Visca", "Wangag")

#################### Functions: ####################


#################### Running things: ####################
# Goal - what does the pop look like over time?:
# Metric 1: number of females - at site level and overall
# Metric 2: number of males - at site level and overall
# Metric 3: size structure through time - at site level and overall

##### Number of females (just from clownfish table for now - could change this later)
## Find raw number of females in database, now using sex column instead of color/size - weirdly this doesn't have as many site/year combos - should look into this more... sticking with the other method for now
females_df_F <- allfish_caught %>%
  #filter(dive_type %in% dive_list) %>%  # not filtering by dives b/c most of the fish in 2017 were getting removed? but there are duplicates of tagged fish and such?? not sure how to deal with this...
  #filter(color %in% c('YP','Y')) %>%
  #filter(size >= min_female_size) %>%
  filter(sex == "F") %>%
  group_by(year,site) %>%
  summarize(nFemalesRaw = n())

## Find raw number of females in database, using color/size
females_df <- allfish_caught %>%
  #filter(dive_type %in% dive_list) %>%  # not filtering by dives b/c most of the fish in 2017 were getting removed? but there are duplicates of tagged fish and such?? not sure how to deal with this...
  filter(color %in% c('YP','Y')) %>%
  filter(size >= min_female_size) %>%
  #filter(sex == "F") %>%
  group_by(year,site) %>%
  summarize(nFemalesRaw = n())

## Scale by proportion of habitat sampled and probability of capturing a fish
# Join in proportion of habitat sampled in each year at each site, then make proportions that are Inf or >1 (b/c of metal tag issues) equal to 1 
females_df <- left_join(females_df, anems_visited_by_year %>% 
                          filter(method == "metal tags") %>% 
                          select(site, year, prop_hab_sampled_tidied), by = c('year', 'site'))  
# females_df_F <- left_join(females_df_F, anems_visited_by_year %>% 
#                           filter(method == "metal tags") %>% 
#                           select(site, year, prop_hab_sampled_tidied), by = c('year', 'site'))  

# Join in whether or not a site was sampled in a particular year
females_df <- left_join(site_visits, females_df, by = c('year', 'site'))

# Find average probability of capturing a fish (prob_r) from KC code
prob_r_avg <- mean(prob_r)

# Add in prob_r, scale up females by proportion habitat sampled and probability of capturing a fish
females_df <- females_df %>%
  mutate(prob_r = prob_r_avg) %>%
  mutate(nFemalesScaled = nFemalesRaw/(prob_r*prop_hab_sampled_tidied)) #%>%
  #mutate(sampled = if_else(is.na(sampled), 0, 1))

# Fix sites that are sampled but with 0 proportion habitat sampled (S. Magbangon in 2017)
females_df <- females_df %>%
  mutate(nFemalesEstimated = ifelse(is.infinite(nFemalesScaled), NA, nFemalesScaled))

# # Add an 'all sites' factor
# females_df_allsites <- females_df %>% 
#   group_by(year) %>%
#   summarize(prop_hab_sampled_as = sum(sampled))
#             prop_hab_sampled = )


# Fit linear models for each site and for overall
females_df_models <- data.frame(site = site_vec_timeseries) %>%
  mutate(intercept = NA,
         coeff =  NA,
         intercept_se = NA,
         coeff_se = NA)

for(i in 1:length(site_vec_timeseries)) {
  site_val = site_vec_timeseries[i]
  df <- females_df %>% 
    filter(site == site_val) %>%
    filter(!is.na(nFemalesEstimated))
  
  testm <- lm(nFemalesEstimated ~ year, data=df)
  testm_S <- summary(testm)
    
  females_df_models$intercept[i] = testm$coefficients[1]
  females_df_models$coeff[i] = testm$coefficients[2]
  females_df_models$intercept_se[i] = testm_S$coefficients[1,2]
  females_df_models$coeff_se[i] = testm_S$coefficients[2,2]
}


#################### Plots: ####################
# Scaled-up number of females at each site through time
pdf(file = here('Plots/DataCharacteristics', 'Scaled_F_through_time_by_site.pdf'))
ggplot(data = females_df, aes(x=year, y=nFemalesEstimated)) +
  #geom_bar(stat='identity') +
  geom_point() +
  geom_abline(data = females_df_models, aes(intercept = intercept, slope = coeff)) +
  #geom_ribbon(data = females_df_models, aes(ymin = ))
  xlab('year') + ylab('# scaled females') + ggtitle('Females by site through time') +
  facet_wrap(~site) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

# Scaled-up number of females at each site through time - attempts at figures
pdf(file = here('Plots/FigureDrafts', 'Time_series_scaled_F_by_site_with_lines.pdf'))
ggplot(data = females_df %>% filter(site %in% site_vec_timeseries), aes(x=year, y=nFemalesEstimated)) +
  #geom_bar(stat='identity') +
  geom_point() +
  geom_abline(data = females_df_models, aes(intercept = intercept, slope = coeff)) +
  #geom_ribbon(data = females_df_models, aes(ymin = ))
  xlab('year') + ylab('# scaled females') + #ggtitle('Females by site through time') +
  facet_wrap(~site) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

# Scaled-up number of females at each site through time - attempts at figures
pdf(file = here('Plots/FigureDrafts', 'Time_series_scaled_F_by_site_no_lines.pdf'))
ggplot(data = females_df, aes(x=year, y=nFemalesEstimated)) +
  #geom_bar(stat='identity') +
  geom_point() +
  #geom_abline(data = females_df_models, aes(intercept = intercept, slope = coeff)) +
  #geom_ribbon(data = females_df_models, aes(ymin = ))
  xlab('year') + ylab('# scaled females') + ggtitle('Females by site through time') +
  facet_wrap(~site) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()


#################### Saving output: ####################

#################### Old code: ####################

#################### Code from elsewhere that might be useful: ####################
# ##### Find number of breeding females
# # Find number of breeding females seen or captured
# breedingF <- allAPCL %>% 
#   group_by(site,year) %>%
#   filter(color %in% c("Y","YP")) %>% #just filter out YP fish (and Y for 2012, when there was no YP)
#   filter(size >= min_size_breedingF) %>% #but also make sure they're big enough to be a breeder (sometimes small fish get marked YP too)
#   summarise(NbreedingF = n())
# 
# # Find number of breeding females tagged
# breedingF_tagged <- allAPCL %>%
#   group_by(year,site) %>%
#   filter(color %in% c("Y","YP")) %>%
#   filter(!is.na(tag_id)) %>%
#   summarise(Ntagged_breedingF = n()) 
# 
# # Combine - since didn't tag all years, might not have any females for some site/years
# breedingF_info <- left_join(breedingF, breedingF_tagged, by = c("site" = "site", "year" = "year")) # no 2012 obs here - why?
# 
# # Find average prob_r from KC code
# prob_r_avg <- mean(prob_r)
# 
# 
# demog_info_recruits <- recruits_info %>%
#   group_by(site) %>%
#   summarize(mean_est_R = mean(totalR_est_metalTA, na.rm=TRUE),
#             high_est_R = max(totalR_est_metalTA, na.rm=TRUE),
#             low_est_R = min(totalR_est_metalTA, na.rm=TRUE),
#             mean_raw_R = mean(Nrecruits, na.rm=TRUE),
#             high_raw_R = max(Nrecruits, na.rm=TRUE),
#             low_raw_R = min(Nrecruits, na.rm=TRUE))
# 
# demog_info_eggs <- breedingF_info %>%
#   group_by(site) %>%
#   summarize(mean_est_eggs = mean(est_eggs_metalTA, na.rm=TRUE),
#             high_est_eggs = max(est_eggs_metalTA, na.rm=TRUE),
#             low_est_eggs = min(est_eggs_metalTA, na.rm=TRUE),
#             mean_raw_F = mean(NbreedingF_combo, na.rm=TRUE),
#             high_raw_F = max(NbreedingF_combo, na.rm=TRUE),
#             low_raw_F = min(NbreedingF_combo, na.rm=TRUE))
# 
