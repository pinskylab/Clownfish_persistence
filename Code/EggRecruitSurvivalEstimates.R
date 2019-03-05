# Estimate egg-recruit survival (here, using Johnson et al method) - was in PersistenceMetrics script, pulling it out here

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))

# Load in site_areas
load(file=here::here('Data','site_areas.RData'))

# Load file with proportion habitat sampled estimates (from TotalAnemsAtSite.R)
load(file=here("Data",'anem_sampling_table.RData')) #file with anems, after 2018 anems matched

# Load metrics on size distribution by color
load(file=here('Data', 'female_sizes.RData'))  # distribution of sizes of females
load(file=here('Data', 'size_by_color_metrics.RData'))  # mean and sd of sizes by YR, O, YP

# Load best estimates of metrics
load(file=here('Data', 'best_est_metrics.RData'))  # from Metrics_with_uncertainty.R, run with best estimates of each param
load(file=here('Data', 'LEP_ests.RData'))

# Load all parentage matches (as of Nov 2018, N and S Mag separated but before 2016, 2017, 2018 genotypes are in)
parentage_moms <- read.csv(file=here('Data','20181017colony_migest_mums_allyears.csv'), stringsAsFactors = FALSE)
parentage_dads <- read.csv(file=here('Data','20181017colony_migest_dads_allyears.csv'), stringsAsFactors = FALSE)
parentage_trios <- read.csv(file=here('Data','20181017colony_migest_trios_allyears.csv'), stringsAsFactors = FALSE)

# Raw mean egg counts
egg_counts_AY_data <- c(479, 590, 586, 305, 679, 683, 387, 720, 427, 688, 169, 655, 414, 352, 1102, 265, 1886, 904,
                        851, 160, 648, 766, 1060, 670, 351, 557)  # from egg_data2018f.csv in Adam's repository
mean_eggs_per_clutch_from_counts <- mean(egg_counts_AY_data)

#################### Functions: ####################

#################### Running things: ####################
# Combine parentage files (mums, dads, trios) - first rename columns so they match across the files, add a column for match type, then rbind
parentage_dads <- parentage_dads %>% 
  dplyr::rename(parent_site = par2_site, nmatches = n_dad, offspring_site = offs_site) %>% 
  mutate(match_type = rep('dad', dim(parentage_dads)[1]))
parentage_moms <- parentage_moms %>% 
  dplyr::rename(parent_site = par1_site, nmatches = n_mum, offspring_site = offs_site) %>%
  mutate(match_type = rep('mom', dim(parentage_moms)[1]))
parentage_trios <- parentage_trios %>% 
  dplyr::rename(parent_site = par1_site, nmatches = n_trios, offspring_site = offs_site) %>%
  mutate(match_type = rep('trio', dim(parentage_trios)[1]))

parentage_matches_raw <- rbind(parentage_dads, parentage_moms, parentage_trios)

########## Assess survival from egg-recruit in a new way - use total parentage matches/N parents genotypes
# Total number of offspring identified
n_offspring_parentage <- sum(parentage_matches_raw$nmatches)  # all offspring identified via parentage

# Total number of parents identified - should have a gen_id for all captures, not just the times(s) it was clipped, right?
n_parents_parentage_df <- allfish_caught %>%
  filter(!is.na(gen_id)) %>%
  filter(size >= min_breeding_M_size | color == 'YP' | color == 'O') %>%  # should think more about what this min size to be
  #filter(color %in% c('YP','O')) %>%
  distinct(gen_id, .keep_all = TRUE) 
n_parents_parentage <- length(n_parents_parentage_df$gen_id)

# could try a mix of denominators?
# 1311 if filter just color in YP, O - no size requirement
# 1664 if filter bigger than 6 or YP or O
# 1515 if filter bigger than 7 or YP or O
# investigating fish that aren't YP or O - pretty evenly divided among years, mix of YR, W, Y colors
# test <- n_parents_parentage_df %>% filter(color != 'YP' & color != 'O')

# For our genotyped 'parents', what fraction are males vs. females?
#n_parents_parentage_F <- 

# Estimate survival male-female

# Not going to worry about that for now, just going to multiply all parents by LEP for a fish that reaches breeding female 
# Or could do LEP for their size? Or could just do LEP for min clip size?
# For now, just doing LEP of 6cm size

# Sum up total site area (all site areas times all years sampled) - total possible sampling area
sites_for_total_areas <- c('Cabatoan', 'Caridad Cemetery', 'Elementary School', 'Gabas', 
                           'Haina', 'Hicgop South', 'N. Magbangon', 'Palanas', 'Poroc Rose',
                           'Poroc San Flower', 'San Agustin', 'Sitio Baybayon', 'Tamakin Dacot',
                           'Visca', 'Wangag', 'S. Magbangon')  # excluding Caridad Proper, Sitio Lonas, Sitio Tugas from this (should I exclude the Sitio Lonas match too?)
site_areas_modified <- site_areas %>% filter(site %in% sites_for_total_areas)

total_area_all_years <- sum(site_areas_modified$kmsq_area)*length(years_sampled)

# Find proportion of habitat sampled overall - sum of area sampled in each
sampled_area_each_year <- left_join(site_areas_modified, anems_table %>% select(site, year, prop_hab_sampled_metal_TA), by = 'site')
sampled_area_each_year <- sampled_area_each_year %>%
  mutate(area_sampled = prop_hab_sampled_metal_TA*kmsq_area)
total_area_sampled <- sum(sampled_area_each_year$area_sampled)

total_prop_hab_sampled <- total_area_sampled/total_area_all_years

# Estimate average number of eggs produced by an individual already having reached 

# Estimate survival from adult-recruit
prop_F_M <- 0.5  # saying 50% of the "adults" we clip are males that won't make it to females -- reasonable? could check this. But LEP takes that into account, right?
tagged_offspring_3.5cm <- n_parents_parentage*LEP_ests$LEP_3.5cm
tagged_offspring_6cm <- n_parents_parentage*LEP_ests$LEP_6cm
recruited_offspring <- n_offspring_parentage/(total_prop_hab_sampled*mean(prob_r))  # scale up by proportion of habitat sampled and probability of catching a fish
surv_egg_recruit <- recruited_offspring/tagged_offspring_3.5cm



############ STOPPED EDITING/RUNNING HERE! ####################


#How does that compare to what we had from the estimated egg-recruit relationship?
surv_RfromE <- findRfromE(slope_mod1, 1, intercept_mod1)


#################### Plots: ####################

#################### Saving output: ####################
save(surv_egg_recruit, file=here('Data','surv_egg_recruit_est.RData'))

#################### Old code: ####################




