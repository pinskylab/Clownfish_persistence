# Estimating persistence metrics 

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))  # Pull in common constants, functions, saved data files, and processed data files 

##### Load libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

##### Load files from other scripts within this repository 

# Load file with proportion habitat sampled estimates (from Code/Total_anems_proportion_hab_sampled.R)
load(file = here::here("Data/Script_outputs", "anems_visited_by_year.RData"))  # has total anems at each site and proportion habitat sampled at each site in each year
load(file = here::here("Data/Script_outputs", "total_area_sampled_through_time.RData"))  # has total area sampled across time (for egg-recruit survival estimate)

# Load file with site widths and distances between sites and coordinates of N/S edges of sampling area (from Code/Site_widths_and_distances.R)
load(file = here::here("Data/Script_outputs", "site_width_info.RData"))
load(file = here::here("Data/Script_outputs", "site_dist_info.RData"))
load(file = here::here("Data/Script_outputs", "sampling_area_edges.RData"))
load(file = here::here("Data/Script_outputs", "site_buffer_info.RData"))

# Load density dependence estimates (from Code/Density_dependence_scaling.R)
load(file=here::here("Data/Script_outputs", "anems_APCL_and_not.RData"))
perc_APCL_val = (anems_APCL_and_not %>% filter(perc_hab == "APCL"))$value
perc_UNOC_val = (anems_APCL_and_not %>% filter(perc_hab == "UNOC"))$value

# Load simple VBL growth analysis (from Code/Growth_analysis.R)
load(file = here::here("Data/Script_outputs", "growth_info_estimate.RData"))
load(file = here::here("Data/Script_outputs", "recap_pairs_year.RData"))  # all recap pairs a year apart, for plotting purposes

# Load survival and recap from MARK models (from Code/ClownfishMarkRecap.R)
load(file = here::here("Data/Script_outputs", "best_fit_model_dfs.RData"))  # load best MARK model
Phi_size_pos = 17  # placement of size effect for survival
p_int_pos = 18  # placement of intercept for recap prob
p_size_pos = 19  # placement of size effect for recap prob
p_dist_pos = 20  # placement of distance effect for recap prob
# survival_output = readRDS(file=here::here("Data/Script_outputs", "eall_mean_Phi_size_p_size_plus_dist.RData"))  # MARK output (lowest AICc model)
# Phi_int_pos = 1  # placement of intercept for survival
# Phi_size_pos = 2  # placement of size effect for survival
# p_int_pos = 3  # placement of intercept for recap prob
# p_size_pos = 4  # placement of size effect for recap prob
# p_dist_pos = 5  # placement of distance effect for recap prob

##### NEED TO CHECK THE TimeSeriesPersistence SCRIPT #####
# Load output from abundance trend (for plotting) (from Code/TimeSeriesPersistence.R)
load(file = here::here("Data/Script_outputs", "site_trends_all.RData"))
load(file = here::here("Data/Script_outputs", "site_trends_time.RData"))
########################################################## 

# # Figure out where these outputs came from so can source those scripts too!
# # Should have two options: source the files that create these outputs or load them from Data folder
# #load(file=here('Data', 'female_sizes.RData'))  # sizes of females from data
# load(file=here::here('Data', 'eall_mean_Phi_size_p_size_plus_dist.RData'))  # MARK output (lowest AICc model) - not sure this has been udpated since I changed the growth in it? Should check!

##### Load input from other analyses outside this repository

# Size transition info (Michelle analysis in genomics repo) - switch this to just females from males (is that reasonable?)
recap_first_male = readRDS(file=here::here("Data/From_other_analyses", "recap_first_male.RData"))
recap_first_female = readRDS(file=here::here("Data/From_other_analyses", "recap_first_female.RData"))

##### Set parameters (for running IPM, for calculating connectivity, for uncertainty runs, etc.)

# Fecundity (for LEP) 
size_fecundity_model = length_count8llEA  # best fit model from Adam fecundity work, both length and eggs on log scale
eggs_intercept_log = size_fecundity_model$coefficients[1]  # on log-scale
eggs_slope_log = size_fecundity_model$coefficients[2]
eyed_effect = size_fecundity_model$coefficients[3]  # effect of eggs having visible eyes

# Set params for IPM structure for calculating LEP
n_bins = 100
n_tsteps = 100
# start_recruit_size = 3.5  # size of recruit that starts out the IPM for LEP
start_recruit_sd = 0.1  # should estimate this somehow! Or do sensitivity to it?

# Parameters for incorporating growth
size_around_mean_L1_size = 0.1  # how much buffer to give for finding fish of that size (so searching for the sd of L2 of fish within 0.2 cm range)

# Parent size (for "tagged" eggs produced by parents)
parent_size = 6.0

# Parameters for what-ifs
n_sites = length(site_vec_order$site_name)  # number of sites for habitat, region width, and larval navigation simulations
perc_hab_vals <- seq(from=0.05,to=1.0, by=0.05)  # set range of percent habitat to test for sensitivity to amount of habitat in region

##### Parameter info (candidates for uncertainty)

# # Growth (for LEP)
# # s = exp(-0.0148)  # what is this?? goes into dnorm for growth part... sd around the mean size? Not sure where this estimate came from... should update it with my estimates
# k_growth_mean = mean(growth_info_estimate$k_est)  # from Growth_analysis growth work (very simple)
# Linf_growth_mean = mean(growth_info_estimate$Linf_est)  # from Growth_analysis growth work (very simple)
# 
# # estimating s - sd of size in a year from a starting size in year 1 (using mean size of L1)
# mean_L1_size <- mean(recap_pairs_year$L1)
# size_around_mean_L1_size = 0.1  # how much buffer to give for finding fish of that size (so searching for the sd of L2 of fish within 0.2 cm range)
# s <- sd((recap_pairs_year %>% filter(mean_L1_size - size_around_mean_L1_size <= L1 & L1 <= mean_L1_size + size_around_mean_L1_size))$L2)

#k_growth_mean = 0.9447194  # lowest AIC model
#Linf_growth_mean = 10.50670  # lowest AIC model
#Linf_growth_sd = sqrt(1.168163)  # from variance for Linf in lowest AIC model

# ##### Find total number of metal-tagged anemones
# n_total_metal_anems <- anems_visited_by_year %>% filter(method == "metal tags") %>% filter(year == "2018") %>% summarize(n_total_anems_all_sites = sum(n_total_anems))
# 
# ##### Parameters for what-ifs
# n_sites = length(site_vec_order$site_name)  # number of sites for habitat, region width, and larval navigation simulations
# perc_hab_vals <- seq(from=0.05,to=1.0, by=0.05)  # set range of percent habitat to test for sensitivity to amount of habitat in region


#################### Functions: ####################
# # Find probability of dispersing distance d with all-years fit (old version, where theta=0.5)
# disp_kernel_all_years <- function(d, k, theta) {  # theta = 0.5, equation for p(d) in eqn. 6c in Bode et al. 2018
#   z = exp(k)
#   disp = (z/2)*exp(-(z*d)^(theta))
#   return(disp)
# }
#
# # Dispersal kernel with best-fit params as a function of d - think about where to put this now that egg-recruit survival will depend on dispersal kernel
# disp_allyears_d <- function(d) {  # theta = 0.5, equation for p(d) in eqn. 6c in Bode et al. 2018
#   z = exp(k_allyears)
#   disp = (z/2)*exp(-(z*d)^(theta_allyears))
#   return(disp)
# }


# For flexible theta and k, function of d, k, and theta
disp_kernel_all_years <- function(d, k, theta) {  # generalization of equation for p(d) in eqn. 6c in Bode et al. 2018
  z = exp(k)
  z_front = (z*theta)/gamma(1/theta)
  disp = (z_front/2)*exp(-(z*d)^(theta))
  return(disp)
}

# For flexible theta and k, function of just d (takes k_d and theta from the environment)
disp_allyears_d <- function(d) {  # generalization of equation for p(d) in eqn. 6c in Bode et al. 2018
  z = exp(k_allyears)
  z_front = (z*theta_allyears)/gamma(1/theta_allyears)
  disp = (z_front/2)*exp(-(z*d)^(theta_allyears))
  return(disp)
}

# # SHOULD FIND A WAY OF CHOOSING BASED ON THETA!
# # Find probability of dispersing distance d with all-years fit (new version, where theta=1)
# disp_kernel_all_years <- function(d, k, theta) {  # theta = 1, equation for p(d) in eqn. 6c in Bode et al. 2018
#   z = exp(k)
#   disp = z*exp(-(z*d)^(theta))
#   return(disp)
# }
#
# # Dispersal kernel with best-fit params as a function of d - think about where to put this now that egg-recruit survival will depend on dispersal kernel
# disp_allyears_d <- function(d) {  # theta = 1, equation for p(d) in eqn. 6c in Bode et al. 2018
#   z = exp(k_allyears)
#   disp = z*exp(-(z*d)^(theta_allyears))
#   return(disp)
# }

# Find beta distribution parameters from mean and variance (mostly for prob_r uncertainty)
findBetaDistParams <- function(mu, var) {
  alpha = (((1-mu)/var) - (1/mu))*mu^2
  beta = alpha*(1/mu -1)

  out = list(alpha = alpha, beta = beta)
  return(out)
}

# Growth
VBL_growth <- function(Linf, k_growth, length) {
  Ls = Linf - (Linf - length)*exp(-k_growth)
  return(Ls)
}

# Find eggs by fish size (eyed eggs) - should double check used log (ln) and not log10 (base 10 log)
findEggs = function(fish_size, egg_size_intercept, egg_size_slope, eyed_effect) {
  count_logged = egg_size_intercept + egg_size_slope*log(fish_size) + eyed_effect
  raw_eggs = exp(count_logged)
  return(raw_eggs)
}

# Find scaled number of tagged recruits we would expect to have found if we sampled the whole area and caught all the fish
scaleTaggedRecruits = function(offspring_assigned, total_prop_hab_sampled, prob_capture, prop_total_disp_area_sampled, prop_hab) {
  rto <- offspring_assigned/(total_prop_hab_sampled*prob_capture*prop_total_disp_area_sampled*prop_hab)
  #recruited_tagged_offspring_total_oursites = offspring_assigned_to_parents/(total_prop_hab_sampled*prob_capture)  # scale by proportion of habitat in our sites we sampled and by prob of catching a fish
  #recruited_tagged_offspring_total = recruited_tagged_offspring_oursites/(prop_total_disp_area_sampled*prop_hab)
  return(rto)
}

# Just consider recruits at our sites
scaleTaggedRecruits_local = function(offspring_assigned, total_prop_hab_sampled, prob_capture) {
  rto <- offspring_assigned/(total_prop_hab_sampled*prob_capture)
  #recruited_tagged_offspring_total_oursites = offspring_assigned_to_parents/(total_prop_hab_sampled*prob_capture)  # scale by proportion of habitat in our sites we sampled and by prob of catching a fish
  #recruited_tagged_offspring_total = recruited_tagged_offspring_oursites/(prop_total_disp_area_sampled*prop_hab)
  return(rto)
}

# Find scaled number of tagged recruits by open habitat, density dependence
scaleTaggedRecruitsDD = function(recruited_tagged_offspring, perc_UNOC, perc_APCL) {
  recruited_tagged_offspring_total_DD = ((perc_APCL+perc_UNOC)/perc_UNOC)*recruited_tagged_offspring  # scale as if currently-occupied APCL habitat is open
  return(recruited_tagged_offspring_total_DD)
}

# Find recruits per egg (using Johnson et al. like method)
findRecruitsPerTaggedEgg = function(tagged_recruits, tagged_eggs) {
  tagged_recruits/tagged_eggs
}

# Find LEP (SHOULD CHECK, UPDATE THIS!)
findLEP = function(min_size, max_size, n_bins, t_steps, Sint, Sl, s, Linf, k_growth, eggs_per_clutch, clutches_per_year,
                   breeding_size, start_recruit_size, start_recruit_sd, egg_size_slope, egg_size_intercept, eyed_effect) {

  # Create vector of lengths
  lengths_vec = seq(min_size, max_size, length.out = n_bins)
  dx = diff(lengths_vec[1:2])

  # Make matrix
  xmat = matrix(rep(lengths_vec, n_bins), nrow=n_bins, ncol=n_bins, byrow=TRUE)  # 100 rows of lengths_vec
  ymat = t(xmat)

  # Survival probability (based on size)
  S = logit_recip(Sint + lengths_vec*Sl)  # survival probability (based on size) - make sure doing this right, with the whole logit_recip thing...
  S[is.na(S)] = 1  # replace the NaNs (coming from Inf/Inf in logit_recip) with 1
  Smat = matrix(rep(S,n_bins), nrow=n_bins, ncol=n_bins, byrow = TRUE)  # survival part of kernel

  # Growth probability
  Ls = Linf - (Linf - xmat)*exp(-k_growth)  # expected length (z') in next time step based on current length (z)
  Lmat = dnorm(ymat,Ls,s)  # turn that into a matrix of probability densities

  # Make the kernel with the survival + growth inputs
  K = Lmat*Smat
  K = dx*K

  # Iterate with 1 subadult to start  # 3.5, 0.1
  N0 = dnorm(lengths_vec, start_recruit_size, start_recruit_sd)  # create a size distribution for one recruit

  # Integrate and scale to have it integrate to 1
  N0 = N0/sum(N0*dx)

  # Iterate through time
  N = matrix(rep(NA,n_bins*n_tsteps), nrow=n_bins, ncol=n_tsteps)
  N[,1] = N0  # initial conditions
  # iterate through time
  for (t in 2:n_tsteps) {
    N[,t] = K %*% N[, t-1]
  }

  # Make fecundity matrix
  # Start with vector of whether or not a fish is reproducing (has it reached female or not yet)
  Fvec = rep(0,n_bins)  # start with 0 (so no reproduction) for all sizes
  Fvec[which(lengths_vec > breeding_size)] = 1

  # # Vector of eggs produced per clutch, dependent on female size - old, non-log-transformed version
  # eggs_by_size_vec = egg_size_intercept + lengths_vec*egg_size_slope
  # eggs_by_size_vec[eggs_by_size_vec < 0] = 0  # make negative values of eggs 0

  # Vector of eyed eggs produced per clutch, dependent on female size - log-transformed version
  eggs_by_size_vec = findEggs(lengths_vec, egg_size_intercept, egg_size_slope, eyed_effect)
  eggs_by_size_vec[eggs_by_size_vec < 0] = 0  # make negative values of eggs 0

  # Multiply breeding-or-not vec by eggs fec to get eggs-per-clutch at each size, then multiply by clutches per year
  Fec_by_size = Fvec*eggs_by_size_vec*clutches_per_year

  # And put this into a matrix, like ymat, where lengths are down the rows and columns are replicated
  Fmat = matrix(rep(Fec_by_size, n_tsteps), nrow=n_bins, ncol=n_tsteps, byrow=FALSE)  # 100 columns of Fec_by_size

  # Now multiply the size-structure produced in each year by the fecundity-by-length
  egg_out = N*Fmat

  # And integrate over all the eggs that get produced
  LEP = sum(egg_out*dx)

  # # Compare to the LEP estimate without size-dependent fecundity
  # # Find the number of breeding adults one recruit produces
  # breeding_Adults <- colSums(N[lengths_vec>breeding_size,]*dx)
  #
  # # Eggs produced by one breeding adult in one year
  # eggs <- eggs_per_clutch*clutches_per_year
  #
  # # Combine to get LEP
  # LEP_nossF <- sum(breeding_Adults*eggs)
  #
  # out = list(LEP=LEP, LEP_nossF = LEP_nossF)

  return(LEP)
}

# Run through a metric calculation with one set of parameters
calcMetrics <- function(param_set, site_based_surv_sets, sites_and_dists, site_vec, DD) {

  sites <- site_vec$site_name

  # Define function with the right parameters
  disp_allyears <- function(d) {  # theta = 1, equation for p(d) in eqn. 6c in Bode et al. 2018
    z = exp(param_set$k_connectivity)
    z_front = (z*param_set$theta_connectivity)/gamma(1/param_set$theta_connectivity)
    disp = (z_front/2)*exp(-(z*d)^(param_set$theta_connectivity))
    return(disp)
  }
  
  # Create connectivity matrix
  Cmat <- sites_and_dists %>% select(org_site, dest_site, d1_km, d2_km, org_alpha_order, org_geo_order, dest_alpha_order, dest_geo_order)
  for(i in 1:length(Cmat$org_site)) {
    Cmat$prob_disp[i] <- integrate(disp_allyears, Cmat$d1_km[i], Cmat$d2_km[i])$value
  }

  # Find LEP (in terms of eggs) by site
  LEP_by_site <- data.frame(site = (site_vec %>% arrange(geo_order))$site_name, LEP = NA, LEP_parents = NA)
  for(i in 1:length(LEP_by_site$site)) {
    Sint_set <- (site_based_surv_sets %>% filter(site == LEP_by_site$site[i]))$Sint
    Sl_set <- (site_based_surv_sets %>% filter(site == LEP_by_site$site[i]))$Sl

    LEP_by_site$LEP[i] = findLEP(param_set$min_size, param_set$max_size, param_set$n_bins, param_set$t_steps, Sint_set, Sl_set,
                                 param_set$s, param_set$Linf, param_set$k_growth, param_set$eggs_per_clutch, param_set$clutches_per_year,
                                 param_set$breeding_size, param_set$start_recruit_size, param_set$start_recruit_sd,
                                 param_set$egg_size_slope, param_set$egg_size_intercept, param_set$eyed_effect)
    LEP_by_site$LEP_parents[i] = findLEP(param_set$min_size, param_set$max_size, param_set$n_bins, param_set$t_steps, Sint_set, Sl_set,
                                                 param_set$s, param_set$Linf, param_set$k_growth, param_set$eggs_per_clutch, param_set$clutches_per_year,
                                                 param_set$breeding_size, parent_size, param_set$start_recruit_sd,
                                                 param_set$egg_size_slope, param_set$egg_size_intercept, param_set$eyed_effect)
  }

  # # Find LEP (in terms of eggs) - COULD MAKE MORE OF THESE PARAMETERS PULLED FROM A DISTRIBUTION!
  # LEP = findLEP(param_set$min_size, param_set$max_size, param_set$n_bins, param_set$t_steps, param_set$Sint, param_set$Sl,
  #               param_set$s, param_set$Linf, param_set$k_growth, param_set$eggs_per_clutch, param_set$clutches_per_year,
  #               param_set$breeding_size, param_set$start_recruit_size, param_set$start_recruit_sd,
  #               param_set$egg_size_slope, param_set$egg_size_intercept, param_set$eyed_effect)

  # Find egg-recruit survival (recruits/egg)
  tagged_recruits_val = scaleTaggedRecruits(param_set$offspring_assigned_to_parents, param_set$total_prop_hab_sampled,
                                            param_set$prob_r, param_set$prop_total_disp_area_sampled, param_set$prop_hab)

  tagged_recruits_local = scaleTaggedRecruits_local(param_set$offspring_assigned_to_parents, param_set$total_prop_hab_sampled,
                                                    param_set$prob_r)

  # # Find LEP_parents by site
  # LEP_parents = findLEP(param_set$min_size, param_set$max_size, param_set$n_bins, param_set$t_steps, param_set$Sint, param_set$Sl,
  #                       param_set$s, param_set$Linf, param_set$k_growth, param_set$eggs_per_clutch, param_set$clutches_per_year,
  #                       param_set$breeding_size, 6.0, param_set$start_recruit_sd,
  #                       param_set$egg_size_slope, param_set$egg_size_intercept, param_set$eyed_effect)

 LEP_parents_mean <- mean(LEP_by_site$LEP_parents)

 tagged_eggs_val = param_set$n_parents*LEP_parents_mean  # not doing this by site for now... just averaging LEP parents across sites

  if(DD == "TRUE"){
    tagged_recruits_val_scaled = scaleTaggedRecruitsDD(tagged_recruits_val, param_set$perc_UNOC, param_set$perc_APCL)
    tagged_recruits_local_scaled = scaleTaggedRecruitsDD(tagged_recruits_local, param_set$perc_UNOC, param_set$perc_APCL)
  } else {
    tagged_recruits_val_scaled = tagged_recruits_val
    tagged_recruits_local_scaled = tagged_recruits_local
  }
  #recruits_per_egg = findRecruitsPerTaggedEgg(tagged_recruits_val, tagged_eggs_val)
  recruits_per_egg = tagged_recruits_val_scaled/tagged_eggs_val

  # # Find egg-recruit survival (recruits/egg) - RIGHT NOW, USING JOHNSON-LIKE ESTIMATE BUT COULD MAKE THIS A DISTRIBUTION TOO
  # recruits_per_egg = param_set$recruits_per_egg

  # Find LEP in terms of recruits
  LEP_by_site <- LEP_by_site %>%
    mutate(LEP_R = LEP*recruits_per_egg,
           LEP_R_local = LEP*(tagged_recruits_local_scaled/tagged_eggs_val))

  LEP_R_mean <- mean(LEP_by_site$LEP_R)

  LEP_R_local_mean = mean(LEP_by_site$LEP)*(tagged_recruits_local_scaled/tagged_eggs_val)

  # # Find LEP in terms of recruits
  # LEP_R = LEP*recruits_per_egg
  # LEP_R_local = LEP*(tagged_recruits_local_scaled/tagged_eggs_val)

  # Find connectivity matrix - EVENTUALLY, WILL USE CONFIDENCE INTERVALS AROUND DISPERSAL KERNELS TO DO THIS - FOR ALL-YEARS ONE? NOT SURE...
  #conn_matrix = Cmatrix
  conn_matrix <- matrix(NA,ncol=max(Cmat$org_geo_order), nrow=max(Cmat$org_geo_order))
  for(i in 1:length(Cmat$org_site)) {
    column = Cmat$org_geo_order[i]  # column is origin
    row = Cmat$dest_geo_order[i]  # row is destination
    conn_matrix[row, column] = Cmat$prob_disp[i]
  }

  # Make realized connectivity matrix, both in the matrix form and dataframe form
  conn_matrixR = conn_matrix
  for(i in 1:length(LEP_by_site$site)) {
    conn_matrixR[,i] <- conn_matrix[,i]*LEP_by_site$LEP_R[i]
  }
  Cmat <- left_join(Cmat, LEP_by_site %>% select(site, LEP_R), by=c("org_site" = "site"))
  Cmat <- Cmat %>%
    mutate(prob_disp_R = prob_disp*LEP_R)

  # # Make realized connectivity matrix, both in the matrix form and dataframe form
  # conn_matrixR = conn_matrix*LEP_R  # matrix form (for eigenvalues)
  # Cmat <- Cmat %>%
  #   mutate(prob_disp_R = prob_disp*LEP_R)  # dataframe form (for plotting)

  # Assess network persistence
  eig_cR = eigen(conn_matrixR)

  # Pull out self-persistence values (diagonals of realized connectivity matrix)
  SP_values = data.frame(site = sites, stringsAsFactors = FALSE) %>%
    mutate(SP_value = NA, org_geo_order = NA)
  for(i in 1:length(sites)) {
    # SP_values$SP_value[i] = conn_matrixR[i,i]  # pull out diagonal entries of realized connectivity matrix - I think something is going wrong here, since Poroc Rose is always 0 in SP but doesn't seem to be elsewhere
    site_val = sites[i]
    SP_values$SP_value[i] = (Cmat %>% filter(org_site == site_val & dest_site == site_val))$prob_disp_R
    SP_values$org_geo_order[i] = (Cmat %>% filter(org_site == site_val & dest_site == site_val))$org_geo_order
  }

  # Put outputs together into one list
  out = list(NP = Re(eig_cR$values[1]), SP = SP_values, LEP_by_site = LEP_by_site, LEP_R_mean = LEP_R_mean, LEP_R_local_mean = LEP_R_local_mean, recruits_per_egg = recruits_per_egg,
             conn_matrix = conn_matrix, conn_matrixR = conn_matrixR, Cmat = Cmat)
}

# Calculate metrics across many runs (this is slow, if have time, should really try to write this without a for loop....)
calcMetricsAcrossRuns <- function(n_runs, param_sets, site_based_surv_sets, site_dist_info, site_vec_order, set_name, DD) {

  LEP_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP', inout = 'input', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_min_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP', inout = 'input', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_max_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP', inout = 'input', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_R_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP_R', inout = 'input', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_R_min_out_df <-  data.frame(value = rep(NA, n_runs), metric = 'LEP_R_max', inout = 'input', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_R_max_out_df <-  data.frame(value = rep(NA, n_runs), metric = 'LEP_R_max', inout = 'input', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_R_local_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP_R_local', inout = 'input', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  RperE_out_df <- data.frame(value = rep(NA, n_runs), metric = 'recruits per egg', inout = 'input', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  NP_out_df <- data.frame(value = rep(NA, n_runs), metric = 'NP', inout = 'output', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)

  metric_vals <- data.frame(run = seq(1:n_runs), LEP = NA, LEP_R = NA, LEP_R_local = NA, recruits_per_egg = NA, NP = NA, uncertainty_type = set_name, stringsAsFactors = FALSE)

  # Create vector of sites for SP dataframe
  runsrepped = rep(1, length(site_vec_order$site_name))
  for(i in 2:n_runs) {
    runsrepped = c(runsrepped, rep(i, length(site_vec_order$site_name)))
  }

  SP_out_df <- data.frame(value = rep(NA, length(site_vec_order$site_name)*n_runs), metric = 'SP', run = NA,
                          site = NA, org_geo_order = NA, uncertainty_type = set_name)
  LEP_by_site_out_df <- data.frame(value = rep(NA, length(site_vec_order$site_name)*n_runs), metric = "LEP", run = NA,
                                               site = NA, uncertainty_type = set_name)
  LEP_R_by_site_out_df <- data.frame(value = rep(NA, length(site_vec_order$site_name)*n_runs), metric = "LEP_R", run = NA,
                                                 site = NA, uncertainty_type = set_name)

  # Calculate the metrics for each parameter set, fill into the data frames - would be faster if this wasn't a for loop (and doesn't need to be, since the runs don't depend on each other at all...)
  for(i in 1:n_runs) {

    # Select parameter set
    params <- param_sets[i,]

    surv_params <- data.frame(site = c(no_space_sites_revisited, sites_not_revisited), Sint = NA, Sl = NA, stringsAsFactors = FALSE)
    for(j in 1:length(c(no_space_sites_revisited, sites_not_revisited))) {
      surv_params$site[j] = c(no_space_sites_revisited, sites_not_revisited)[j]
      surv_params$Sint[j] = site_based_surv_sets[[j]]$Sint[i]
      surv_params$Sl[j] = site_based_surv_sets[[j]]$Sl[i]
    }

    # Do the run
    metrics_output = calcMetrics(params, surv_params, site_dist_info, site_vec_order, DD)

    print(i)

    # Fill in the metrics
    LEP_out_df$value[i] = mean(metrics_output$LEP_by_site$LEP)
    LEP_min_out_df$value[i] = min(metrics_output$LEP_by_site$LEP)
    LEP_max_out_df$value[i] = max(metrics_output$LEP_by_site$LEP)
    LEP_R_out_df$value[i] = metrics_output$LEP_R_mean
    LEP_R_min_out_df$value[i] = min(metrics_output$LEP_by_site$LEP_R)
    LEP_R_max_out_df$value[i] = max(metrics_output$LEP_by_site$LEP_R)
    LEP_R_local_out_df$value[i] = metrics_output$LEP_R_local_mean
    RperE_out_df$value[i] = metrics_output$recruits_per_egg
    NP_out_df$value[i] = metrics_output$NP

    metric_vals$LEP[i] = mean(metrics_output$LEP_by_site$LEP)
    metric_vals$LEP_R[i] = metrics_output$LEP_R_mean
    metric_vals$LEP_R_local[i] = metrics_output$LEP_R_local_mean
    metric_vals$recruits_per_egg[i] = metrics_output$recruits_per_egg
    metric_vals$NP[i] = metrics_output$NP

    # Pull out the SP metrics
    start_index = (i-1)*length(site_vec_order$site_name)+1
    end_index = i*length(site_vec_order$site_name)

    SP_out_df$site[start_index:end_index] = metrics_output$SP$site
    SP_out_df$value[start_index:end_index] = metrics_output$SP$SP_value
    SP_out_df$run[start_index:end_index] = rep(i, length(site_vec_order$site_name))
    SP_out_df$org_geo_order[start_index:end_index] = metrics_output$SP$org_geo_order

    LEP_by_site_out_df$site[start_index:end_index] = metrics_output$LEP_by_site$site
    LEP_by_site_out_df$value[start_index:end_index] = metrics_output$LEP_by_site$LEP
    LEP_by_site_out_df$run[start_index:end_index] = rep(i, length(site_vec_order$site_name))

    LEP_R_by_site_out_df$site[start_index:end_index] = metrics_output$LEP_by_site$site
    LEP_R_by_site_out_df$value[start_index:end_index] = metrics_output$LEP_by_site$LEP*metrics_output$recruits_per_egg
    LEP_R_by_site_out_df$run[start_index:end_index] = rep(i, length(site_vec_order$site_name))

  }

  # Add some of the changing parameters in, so can look at in plots
  metric_vals_with_params <- metric_vals %>%
    mutate(breeding_size = param_sets$breeding_size,
           Linf = param_sets$Linf,
           k_growth = param_sets$k_growth,
           #Sint = param_sets$Sint,
           k_connectivity = param_sets$k_connectivity,
           theta_connectivity = param_sets$theta_connectivity,
           recruits_per_egg = param_sets$recruits_per_egg,
           prob_r = param_sets$prob_r,
           assigned_offspring = param_sets$offspring_assigned_to_parents,
           start_recruit_size = param_sets$start_recruit_size,
           total_prop_hab_sampled = param_sets$total_prop_hab_sampled)

  SP_vals_with_params <- left_join(SP_out_df, metric_vals_with_params, by='run') %>%
    dplyr::rename(SP = value)

  out = list(LEP_out_df = LEP_out_df, LEP_min_out_df = LEP_min_out_df, LEP_max_out_df = LEP_max_out_df, LEP_R_out_df = LEP_R_out_df,
             LEP_R_min_out_df = LEP_R_min_out_df, LEP_R_max_out_df = LEP_R_max_out_df, LEP_R_local_out_df = LEP_R_local_out_df,
             RperE_out_df = RperE_out_df, NP_out_df = NP_out_df, LEP_by_site_out_df = LEP_by_site_out_df, LEP_R_by_site_out_df = LEP_R_by_site_out_df,
             metric_vals_with_params = metric_vals_with_params, SP_vals_with_params = SP_vals_with_params, params_in = param_sets,
             uncertainty_type = set_name)

  return(out)
}

##### Functions for what-ifs
# Find width of sites and distance between them for different percent habitats (and sampling region lengths)
find_site_locations <- function(nsites, perc_hab_val, total_region) {
  
  site_info <- data.frame(org_site = seq(from=1, to=nsites, by=1))
  site_width <- (total_region*perc_hab_val)/nsites
  hab_break <- (total_region*(1-perc_hab_val))/(nsites-1)
  
  for(i in 1:nsites){
    site_info$S_edge_org[i] = (i-1)*(site_width + hab_break)
    site_info$N_edge_org[i] = (i-1)*(site_width + hab_break) + site_width
    site_info$center_org[i] = (i-1)*(site_width + hab_break) + (site_width/2)
    site_info$width_org[i] = site_info$N_edge_org[i] - site_info$S_edge_org[i]
  }
  return(site_info)
}

# Turn into a bigger data frame that can be used in metric functions
make_output_with_dist <- function(nsites, perc_hab_val, total_region) {
  
  # Find locations of each site
  org_df <- find_site_locations(nsites, perc_hab_val, total_region)
  
  # Repeat each row nsites times
  org_df_out <- org_df %>% slice(rep(1:n(), each = nsites))
  
  # Rename to be for dest
  dest_df <- org_df %>%
    dplyr::rename(dest_site = org_site, S_edge_dest = S_edge_org, 
                  N_edge_dest = N_edge_org, center_dest = center_org,
                  width_dest = width_org) 
  
  # Put into a bigger data frame with a row for each org and dest combo
  dest_df_out <- do.call("rbind", replicate(nsites, dest_df, simplify = FALSE))
  
  # Put them together, add placeholder distance columns
  out_df <- cbind(org_df_out, dest_df_out) %>%
    mutate(dist_mid_to_S = NA, dist_mid_to_N = NA, d1_km = NA, d2_km = NA)
  
  # Find distances between each site of sites
  for(i in 1:length(out_df$org_site)) {
    out_df$dist_mid_to_N[i] = abs(out_df$center_org[i] - out_df$N_edge_dest[i])
    out_df$dist_mid_to_S[i] = abs(out_df$center_org[i] - out_df$S_edge_dest[i])
  }
  
  # Find d1 and d2 (min and max distances to each site)
  for(i in 1:length(out_df$org_site)) {
    out_df$d1_km[i] = min(c(out_df$dist_mid_to_S[i], out_df$dist_mid_to_N[i]))
    out_df$d2_km[i] = max(c(out_df$dist_mid_to_S[i], out_df$dist_mid_to_N[i]))
  }
  
  # Fix the self-selfs (will just need to make d1_km = 0 when fix the normalization, I think!)
  out_df <- out_df %>%
    mutate(d1_km = case_when(org_site != dest_site ~ d1_km,
                             org_site == dest_site ~ 0),  # min distance is 0 for self-self distances
           d2_km = case_when(org_site != dest_site ~ d2_km,
                             org_site == dest_site ~ width_dest/2))  # max distance is site width for self-self distances (now 1/2 width)
  
  # Add the random other columns the calcMetrics function needs
  out_df <- out_df %>% 
    mutate(org_alpha_order = org_site,
           org_geo_order = org_site,
           dest_alpha_order = dest_site,
           dest_geo_order = dest_site)
  
  return(out_df)
}


#################### Running things: ####################

########## Find number of tagged fish, genotyped fish, fish that are marked somehow (to report in paper)
n_fish_genotyped <- allfish_caught %>%
  filter(!is.na(gen_id)) %>%
  distinct(gen_id) %>%
  summarize(nfish = n())

n_fish_PIT_tagged <- allfish_caught %>%
  filter(!is.na(tag_id)) %>%
  distinct(tag_id) %>%
  summarize(nfish = n())

n_fish_marked <- allfish_caught %>%
  filter(!is.na(fish_indiv)) %>%
  distinct(fish_indiv) %>%
  summarize(nfish = n())

#################### Find "best estimates" of input parameters if additional calculations needed before use in finding metrics ####################

##### Growth (for LEP)
k_growth_mean = mean(growth_info_estimate$k_est)  # mean of various estimates from different pairs of recaptures for fish caught multiple times
Linf_growth_mean = mean(growth_info_estimate$Linf_est)  # mean of various estimates from different pairs of recaptures for fish caught multiple times

mean_L1_size <- mean(recap_pairs_year$L1)  # mean size of first recapture of a fish for those caught a year apart, use as the starting size to estimate the sd of the size the fish will be in a year
s <- sd((recap_pairs_year %>% filter(mean_L1_size - size_around_mean_L1_size <= L1 & L1 <= mean_L1_size + size_around_mean_L1_size))$L2)  # find the sd of the L2 size of fish that start at size mean_L1_size, within a buffer of +- size_around_mean_L1_size (defined above in Set-up parameters section)

##### Breeding size (used in LEP)
breeding_size_mean <- mean(recap_first_female$size)

##### Recruit size - use mean size of potential (genotyped) offspring
mean_sampled_offspring_size <- mean(all_offspring$size, na.rm = TRUE)  
#start_recruit_size = mean_sampled_offspring_size

# ##### Find total number of metal-tagged anemones
# n_total_metal_anems <- anems_visited_by_year %>% filter(method == "metal tags") %>% filter(year == "2018") %>% summarize(n_total_anems_all_sites = sum(n_total_anems))

##### Scaling factors recruits to find egg-recruit survival

### Probability of capturing a fish (Pc)
prob_r_mean <- mean(prob_r)  # average value of prob r from each recap dive

### Proportion of sampling area that is habitat (Ps) (to scale egg-recruit survival to avoid counting mortality from dispersal to non-habitat twice)
# Find total length of sampling area from N to S (assumes sites are all in a line)
total_range_of_sampling_area <- geosphere::distHaversine(c((sampling_area_edges %>% filter(edge == "north"))$lon, (sampling_area_edges %>% filter(edge == "north"))$lat),
                                                         c((sampling_area_edges %>% filter(edge == "south"))$lon, (sampling_area_edges %>% filter(edge == "south"))$lat))
# Total area covered by sites (sum of site widths)
total_sum_of_site_widths <- sum(site_width_info$width_m)

# (Rough) proportion of sampling area that is habitat
prop_sampling_area_habitat <- total_sum_of_site_widths/total_range_of_sampling_area
Ps <- prop_sampling_area_habitat

### (Pd) How much of the dispersal kernel area from each site did we actually sample? (So can scale up "tagged" recruits found to account for areas they might have gone that weren't in our sites)
# Find the number and proportion of total parents at each site       
all_parents_site <- all_parents_by_site %>%
  group_by(site) %>%
  summarize(nparents = n()) %>%
  mutate(prop_parents = nparents/sum(nparents))

# Add in site info
all_parents_site <- left_join(all_parents_site, site_width_info %>% select(site, site_geo_order, dist_to_N_edge_km, dist_to_S_edge_km), by = "site")

# Find area of dispersal kernel from middle of each site within sampling area to the north and south of each site
all_parents_site <- all_parents_site %>%
  mutate(disp_area_N_within_sites = NA,
         disp_area_S_within_sites = NA)

for(i in 1:length(all_parents_site$site)) {
  all_parents_site$disp_area_N_within_sites[i] = integrate(disp_allyears_d, 0, all_parents_site$dist_to_N_edge_km[i])$value
  all_parents_site$disp_area_S_within_sites[i] = integrate(disp_allyears_d, 0, all_parents_site$dist_to_S_edge_km[i])$value
}

# Find proportion of total area under dispersal kernel
all_parents_site <- all_parents_site %>%
  mutate(total_disp_area_within_sites = disp_area_N_within_sites + disp_area_S_within_sites,
         prop_disp_area_within_sites = total_disp_area_within_sites,  # now, dispersal kernel normalized to 0.5 so total to N and S combined is 1
         total_parent_area_sampled = total_disp_area_within_sites*nparents)

all_parents_site_summarized <- all_parents_site %>%
  summarize(total_parent_kernel_area = sum(nparents),  
            sampled_parent_kernel_area = sum(total_parent_area_sampled),
            prop_parent_kernel_area_sampled = sampled_parent_kernel_area/total_parent_kernel_area)

prop_total_disp_area_sampled_best_est <- all_parents_site_summarized$prop_parent_kernel_area_sampled
Pd <- prop_total_disp_area_sampled_best_est

### (Ph) What proportion of habitat at our sites did we sample over time? (to scale up to account for habitat areas we didn't sample)
total_prop_hab_sampled_through_time <- (total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area  # use proportion of metal-tagged anemones sampled at each site
Ph <- total_prop_hab_sampled_through_time

##### Survival by site (here, using median site survival for missing sites)

# Set up data frame for best-estimate survival parameters by site
site_surv_best_est <- data.frame(site = no_space_sites_revisited, Sint = NA, 
                                 lcl = NA, ucl = NA,
                                 Sl = best_fit_model_dfs$results$estimate[Phi_size_pos], 
                                 stringsAsFactors = FALSE)

# Cabatoan is the first intercept
site_surv_best_est$Sint[1] = best_fit_model_dfs$results$estimate[1]
site_surv_best_est$lcl[1] = best_fit_model_dfs$results$lcl[1]
site_surv_best_est$ucl[1] = best_fit_model_dfs$results$ucl[1]

# fill in rest of sites, which are additions to the intercept for Cabatoan
for(i in 2:length(no_space_sites_revisited)) {
  site_surv_best_est$Sint[i] = best_fit_model_dfs$results$estimate[1] + best_fit_model_dfs$results$estimate[i]
  site_surv_best_est$lcl[i] = best_fit_model_dfs$results$lcl[i]
  site_surv_best_est$ucl[i] = best_fit_model_dfs$results$ucl[i]
}

# Fill in median of estimated sites for sites not revisited so not estimated (Caridad Proper, Sitio Lonas, Sitio Tugas)
median_site_pos = 3
site_surv_notrevisited = data.frame(site = sites_not_revisited, Sint = (best_fit_model_dfs$results$estimate[1] + best_fit_model_dfs$results$estimate[median_site_pos]),
                                    lcl = best_fit_model_dfs$results$lcl[median_site_pos], ucl = best_fit_model_dfs$results$ucl[median_site_pos],
                                    Sl = best_fit_model_dfs$results$estimate[Phi_size_pos], stringsAsFactors = FALSE)

# Bind them together
site_surv_best_est <- rbind(site_surv_best_est, site_surv_notrevisited)

##### Parentage assignment rate
assignment_rate = n_offspring_matched/n_offspring_genotyped  # proportion of genotyped offspring that were assigned to parents in parentage analysis


########## Estimate some of the "best estimates" of metrics/parameters

# ##### How big are the offspring? What size should we use for recruits?
# mean_sampled_offspring_size <- mean(all_offspring$size, na.rm = TRUE)  # might be duplicate observations of fish in here, but I don't think so - should check
# #mean_sampled_offspring_size <- mean(n_offspring_genotypes_df$size, rm.na = TRUE)  # this is just one of the obs of each of these fish... not sure how many duplicates there are, should really check...
# 
# start_recruit_size = mean_sampled_offspring_size

# ### Best-estimate survival parameters by site
# site_surv_best_est <- data.frame(site = no_space_sites_revisited, Sint = NA, 
#                                  lcl = NA, ucl = NA,
#                                  Sl = best_fit_model_dfs$results$estimate[Phi_size_pos], 
#                                  stringsAsFactors = FALSE)
# 
# # Cabatoan
# site_surv_best_est$Sint[1] = best_fit_model_dfs$results$estimate[1]
# site_surv_best_est$lcl[1] = best_fit_model_dfs$results$lcl[1]
# site_surv_best_est$ucl[1] = best_fit_model_dfs$results$ucl[1]
# 
# # fill in rest of sites
# for(i in 2:length(no_space_sites_revisited)) {
#   site_surv_best_est$Sint[i] = best_fit_model_dfs$results$estimate[1] + best_fit_model_dfs$results$estimate[i]
#   site_surv_best_est$lcl[i] = best_fit_model_dfs$results$lcl[i]
#   site_surv_best_est$ucl[i] = best_fit_model_dfs$results$ucl[i]
# }
# 
# # Fill in median of estimated sites for sites not revisited
# median_site_pos = 3
# site_surv_notrevisited = data.frame(site = sites_not_revisited, Sint = (best_fit_model_dfs$results$estimate[1] + best_fit_model_dfs$results$estimate[median_site_pos]),
#                                     lcl = best_fit_model_dfs$results$lcl[median_site_pos], ucl = best_fit_model_dfs$results$ucl[median_site_pos],
#                                     Sl = best_fit_model_dfs$results$estimate[Phi_size_pos], stringsAsFactors = FALSE)
# 
# # # Fill in average of other sites for sites not revisited and estimated (Sitio Lonas, Sitio Tugas, Caridad Proper)
# # avg_Sint = sum(best_fit_model_dfs$results$estimate[1]*16, best_fit_model_dfs$results$estimate[2:16])/16  
# # avg_Sint_withoutCC = sum(best_fit_model_dfs$results$estimate[1]*15, best_fit_model_dfs$results$estimate[3:16])/15  # avg Sint of all sites except Caridad Cemetery  
# # site_surv_notrevisited = data.frame(site = sites_not_revisited, Sint = avg_Sint, Sl = best_fit_model_dfs$results$estimate[Phi_size_pos], stringsAsFactors = FALSE)
# 
# # Bind them together
# site_surv_best_est <- rbind(site_surv_best_est, site_surv_notrevisited)

# ##### LEP, from a few different starting sizes
# LEP_different_sizes <- site_surv_best_est %>%
#   mutate(LEP_mean_breeding_size = NA,
#          LEP_6cm = NA,
#          LEP_3.5cm = NA)
# 
# for(i in 1:length(LEP_different_sizes$site)) {
#   
#   # Starting from mean transition to female size
#   LEP_different_sizes$LEP_mean_breeding_size[i] <- findLEP(min_size, max_size, n_bins, t_steps, LEP_different_sizes$Sint[i], LEP_different_sizes$Sl[i],
#                                                            s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts, 
#                                                            clutches_per_year_mean, breeding_size_mean, breeding_size_mean, start_recruit_sd, 
#                                                            eggs_slope_log, eggs_intercept_log, eyed_effect)
#   # Starting at tagging size (6cm)
#   LEP_different_sizes$LEP_6cm[i] <- findLEP(min_size, max_size, n_bins, t_steps, LEP_different_sizes$Sint[i], LEP_different_sizes$Sl[i],
#                                             s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts, 
#                                             clutches_per_year_mean, breeding_size_mean, 6, start_recruit_sd, 
#                                             eggs_slope_log, eggs_intercept_log, eyed_effect)
#   
#   # Starting at fin-clip size (3.5cm)
#   LEP_different_sizes$LEP_3.5cm[i] <- findLEP(min_size, max_size, n_bins, t_steps, LEP_different_sizes$Sint[i], LEP_different_sizes$Sl[i],
#                                               s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts, 
#                                               clutches_per_year_mean, breeding_size_mean, 3.5, start_recruit_sd, 
#                                               eggs_slope_log, eggs_intercept_log, eyed_effect)
# }
# 
# ##### Find estimated number of tagged eggs, estimated number of tagged recruits (scaled up), the egg-recruit survival
# # How many potential offspring (eggs) were produced by those potential (genotyped) parents ("tagged" adults b/c genetically marked)?
# tagged_eggs_6cm <- n_parents_genotyped*mean(LEP_different_sizes$LEP_6cm)
# tagged_eggs_3.5cm <- n_parents_genotyped*mean(LEP_different_sizes$LEP_3.5cm)
# # # tagged_eggs_6cm <- n_parents_genotyped*LEP_6cm  # Use LEP from what size here? How to avoid double-counting if those parents mated together?
# # # tagged_eggs_3.5cm <- n_parents_genotyped*LEP_3.5cm
# 
# # Scale up the number of "tagged offspring" by the probability of catching a fish (prob_r), proportion of site area we sampled (Ph), proportion of dispersal kernel area we sampled (Pd), and proportion of sample area that is habitat (Ps)
# recruited_tagged_offspring_scaled <- n_offspring_matched/(prob_r_mean*Ps*Pd*Ph)
# recruited_tagged_offspring_oursites <- n_offspring_matched/(prob_r_mean*Ph)  # just scaling up to recruits arriving back to our sites
# 
# # Scale up by DD
# recruited_tagged_offspring_scaled_DD <- scaleTaggedRecruitsDD(recruited_tagged_offspring_scaled, perc_UNOC_val, perc_APCL_val)
# recruited_tagged_offspring_oursites_DD <- scaleTaggedRecruitsDD(recruited_tagged_offspring_oursites, perc_UNOC_val, perc_APCL_val)
# 
# # Assignment rate
# assignment_rate = n_offspring_matched/n_offspring_genotyped  # proportion of genotyped offspring that were assigned to parents in parentage analysis
# 
# ##### Estimate egg-recruit survival
# # Estimate survival from eggs-recruits by seeing how many "tagged" offspring we found out of eggs "tagged" parents produced
# recruits_per_egg_best_est <- recruited_tagged_offspring_scaled/tagged_eggs_6cm
# 
# # New best est with DD
# recruits_per_egg_best_est_DD <- recruited_tagged_offspring_scaled_DD/tagged_eggs_6cm
# 
# # Recruits per egg just arriving to our sites (for an LRP for our sites that includes dispersal mortality in the egg-recruit surv)
# recruits_per_egg_oursites <- recruited_tagged_offspring_oursites/tagged_eggs_6cm
# recruits_per_egg_oursites_DD <- scaleTaggedRecruitsDD(recruited_tagged_offspring_oursites, perc_UNOC_val, perc_APCL_val)/tagged_eggs_6cm
# 
# 

# ##### Sensitivity test to different max sizes for LEP
# max_size_sens <- c(9,9.5,10,10.5,11,11.5,12,12.5,13)  # max sizes to try
# LEP_max_size_sens <- list()
# for(i in 1:length(max_size_sens)) {
#   
# }

########## Estimate survival from egg to recruit (method similar to Johnson et al. 2018):

# ##### Scaling factors for recruits to find egg-recruit survival
# 
# ### (Ps) Find proportion of sampling area that is habitat (to scale egg-recruit survival to avoid counting mortality from dispersal to non-habitat twice)
# # Right now this assumes the whole sampling region and all sites are totally in a line, which they are not...
# # Find total length of sampling area from N to S
# total_range_of_sampling_area <- geosphere::distHaversine(c((sampling_area_edges %>% filter(edge == "north"))$lon, (sampling_area_edges %>% filter(edge == "north"))$lat),
#                                                          c((sampling_area_edges %>% filter(edge == "south"))$lon, (sampling_area_edges %>% filter(edge == "south"))$lat))
# # Total area covered by sites (sum of site widths)
# total_sum_of_site_widths <- sum(site_width_info$width_m)
# 
# # (Rough) proportion of sampling area that is habitat
# prop_sampling_area_habitat <- total_sum_of_site_widths/total_range_of_sampling_area
# Ps <- prop_sampling_area_habitat
# 
# ### (Pd) How much of the dispersal kernel area from each site did we actually sample? (So can scale up "tagged" recruits found to account for areas they might have gone that weren't in our sites)
# # Find the number of parents at each site       # OLD COMMENT, RIGHT? (eventually, this parent file pull will go in Constants_database_common_functions). Just putting it here for now b/c going to use original 913, with same distribution as current parents
# all_parents_site <- all_parents_by_site %>%
#   group_by(site) %>%
#   summarize(nparents = n()) %>%
#   mutate(prop_parents = nparents/sum(nparents))
# 
# # #### CONFIRM THAT PARENTS DON'T NEED TO BE MULTIPLIED BY 2 ANYMORE!! Yes, looks like that based on written methods in SI.
# # # Total dispersal kernel area (total parents*2 - total area dispersing north of site is 1 and south is 1 for each parent) - think this isn't true with the new normalization...
# # # total_parent_kernel_area = n_parents_genotyped*2
# # total_parent_kernel_area = n_parents_genotyped
# 
# # Add in site info
# all_parents_site <- left_join(all_parents_site, site_width_info %>% select(site, site_geo_order, dist_to_N_edge_km, dist_to_S_edge_km), by = "site")
# 
# # Find area within sampling area to the north and south of each site
# all_parents_site <- all_parents_site %>%
#   mutate(disp_area_N_within_sites = NA,
#          disp_area_S_within_sites = NA)
# 
# for(i in 1:length(all_parents_site$site)) {
#   all_parents_site$disp_area_N_within_sites[i] = integrate(disp_allyears_d, 0, all_parents_site$dist_to_N_edge_km[i])$value
#   all_parents_site$disp_area_S_within_sites[i] = integrate(disp_allyears_d, 0, all_parents_site$dist_to_S_edge_km[i])$value
# }
# 
# # Find proportion of total area under dispersal kernel
# ### OLD COMMENT - now total area to Inf is 1 (where total area to INF is 2 (1 for each side)) covered within sample sites
# all_parents_site <- all_parents_site %>%
#   mutate(total_disp_area_within_sites = disp_area_N_within_sites + disp_area_S_within_sites,
#          prop_disp_area_within_sites = total_disp_area_within_sites,  # now, dispersal kernel normalized to 0.5 so total to N and S combined is 1
#          total_parent_area_sampled = total_disp_area_within_sites*nparents)
# 
# all_parents_site_summarized <- all_parents_site %>%
#   summarize(total_parent_kernel_area = sum(nparents),  # no longer need to multiply by 2 because sum of kernel to N and kernel to S is 1
#             sampled_parent_kernel_area = sum(total_parent_area_sampled),
#             prop_parent_kernel_area_sampled = sampled_parent_kernel_area/total_parent_kernel_area)
# 
# prop_total_disp_area_sampled_best_est <- all_parents_site_summarized$prop_parent_kernel_area_sampled
# Pd <- prop_total_disp_area_sampled_best_est
# 
# ### (Ph) What proportion of habitat at our sites did we sample over time?
# # Find the total prop habitat sampled over time - need to update this to have the other sites...
# #total_prop_hab_sampled_through_time <- (total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area*mean(prob_r)  # scale up by proportion of habitat sampled and probability of catching a fish
# total_prop_hab_sampled_through_time <- (total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area  # scale up by proportion of habitat sampledP
# Ph <- total_prop_hab_sampled_through_time

# ##### Find estimated number of tagged eggs, estimated number of tagged recruits (scaled up), the egg-recruit survival
# # How many potential offspring (eggs) were produced by those potential (genotyped) parents ("tagged" adults b/c genetically marked)?
# tagged_eggs_6cm <- n_parents_genotyped*mean(LEP_different_sizes$LEP_6cm)
# tagged_eggs_3.5cm <- n_parents_genotyped*mean(LEP_different_sizes$LEP_3.5cm)
# # # tagged_eggs_6cm <- n_parents_genotyped*LEP_6cm  # Use LEP from what size here? How to avoid double-counting if those parents mated together?
# # # tagged_eggs_3.5cm <- n_parents_genotyped*LEP_3.5cm
# 
# # Scale up the number of "tagged offspring" by the probability of catching a fish (prob_r), proportion of site area we sampled (Ph), proportion of dispersal kernel area we sampled (Pd), and proportion of sample area that is habitat (Ps)
# recruited_tagged_offspring_scaled <- n_offspring_matched/(prob_r_mean*Ps*Pd*Ph)
# recruited_tagged_offspring_oursites <- n_offspring_matched/(prob_r_mean*Ph)  # just scaling up to recruits arriving back to our sites
# 
# # Scale up by DD
# recruited_tagged_offspring_scaled_DD <- scaleTaggedRecruitsDD(recruited_tagged_offspring_scaled, perc_UNOC_val, perc_APCL_val)
# recruited_tagged_offspring_oursites_DD <- scaleTaggedRecruitsDD(recruited_tagged_offspring_oursites, perc_UNOC_val, perc_APCL_val)
# 
# # Assignment rate
# assignment_rate = n_offspring_matched/n_offspring_genotyped  # proportion of genotyped offspring that were assigned to parents in parentage analysis
# 
# ##### Estimate egg-recruit survival
# # Estimate survival from eggs-recruits by seeing how many "tagged" offspring we found out of eggs "tagged" parents produced
# recruits_per_egg_best_est <- recruited_tagged_offspring_scaled/tagged_eggs_6cm
# 
# # New best est with DD
# recruits_per_egg_best_est_DD <- recruited_tagged_offspring_scaled_DD/tagged_eggs_6cm
# 
# # Recruits per egg just arriving to our sites (for an LRP for our sites that includes dispersal mortality in the egg-recruit surv)
# recruits_per_egg_oursites <- recruited_tagged_offspring_oursites/tagged_eggs_6cm
# recruits_per_egg_oursites_DD <- scaleTaggedRecruitsDD(recruited_tagged_offspring_oursites, perc_UNOC_val, perc_APCL_val)/tagged_eggs_6cm


#################### Estimate metrics ####################

# Put best-estimate intput parameters into one data frame (where starting recruit size is mean size of actual offspring collected ~ 4.38)
param_best_est_mean_collected_offspring <- data.frame(t_steps = n_tsteps) %>%
  mutate(min_size = min_size, max_size = max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year_mean = clutches_per_year_mean,  # average fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect =  eyed_effect,  # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean, 
         #Sl = Sl_mean, Sint = Sint_mean,
         breeding_size = breeding_size_mean, #recruits_per_egg = recruits_per_egg_best_est,
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Calculate the metrics for the best estimates, without and with (_DD) density-dependence compensated for
best_est_metrics_mean_offspring <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info, site_vec_order, "FALSE")
best_est_metrics_mean_offspring_DD <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info, site_vec_order, "TRUE")

# # Do I need any of these? maybe not...
# LEP_best_est <- mean(best_est_metrics_mean_offspring$LEP_by_site$LEP)  # how to deal with site variation?
# LEP_best_est_min <- min(best_est_metrics_mean_offspring$LEP_by_site$LEP)
# LEP_best_est_max <- max(best_est_metrics_mean_offspring$LEP_by_site$LEP)
# LEP_R_best_est_min <- min(best_est_metrics_mean_offspring$LEP_by_site$LEP_R)
# LEP_R_best_est_max <- max(best_est_metrics_mean_offspring$LEP_by_site$LEP_R)
# LEP_R_best_est <- best_est_metrics_mean_offspring$LEP_R_mean
# NP_best_est <- best_est_metrics_mean_offspring$NP
# SP_best_est <- best_est_metrics_mean_offspring$SP
# LEP_R_local_best_est <- best_est_metrics_mean_offspring$LEP_R_local_mean
# 
# LEP_best_est_DD <- mean(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP)  # how to deal with, show site variation?
# LEP_best_est_min_DD <- min(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP)
# LEP_best_est_max_DD <- max(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP)
# LEP_R_best_est_min_DD <- min(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP_R)
# LEP_R_best_est_max_DD <- max(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP_R)
# LEP_R_best_est_DD <- best_est_metrics_mean_offspring_DD$LEP_R_mean
# NP_best_est_DD <- best_est_metrics_mean_offspring_DD$NP
# SP_best_est_DD <- best_est_metrics_mean_offspring_DD$SP
# LEP_R_local_best_est_DD <- best_est_metrics_mean_offspring_DD$LEP_R_local_mean
# 
# LEP_R_oursites <- recruits_per_egg_oursites*LEP_best_est
# LEP_R_oursites_DD <- recruits_per_egg_oursites_DD*LEP_best_est

#################### Generate sets of parameters for uncertainty: ####################

##### Growth - NEED TO DECIDE HOW DOING GROWTH UNCERTAINTY!
Linf_set = growth_info_estimate$Linf_est  # growth (Linf)
k_growth_set = growth_info_estimate$k_est  # growth (k)

### Create better growth estimate parameter sets
growth_set_intercept_est_mean <- mean(growth_info_estimate$intercept_est)  # mean of the intercept estimates from the fits with different growth pairs for fish with multiple
growth_set_intercept_se_mean <- mean(growth_info_estimate$intercept_se)  # mean of intercept se estimates
growth_set_slope_est_mean <- mean(growth_info_estimate$slope_est)  # mean of slope estimates from the fits with different growth pairs randomly chosen for fish with multiple
growth_set_slope_se_mean <- mean(growth_info_estimate$slope_se)

# Select 1000 slopes and intercepts using the mean and se from the fits, then find k and Linf from that
growth_set_params <- data.frame(slope_est = runif(n_runs, growth_set_slope_est_mean-growth_set_slope_se_mean, growth_set_slope_est_mean+growth_set_slope_se_mean),  # choose slope values from within mean+-SE range
                                intercept_est = runif(n_runs, growth_set_intercept_est_mean-growth_set_intercept_se_mean, growth_set_intercept_est_mean+growth_set_intercept_se_mean))  # choose intercept values from within mean+-SE range
growth_set_params <- growth_set_params %>%
  mutate(k_est = -log(slope_est),
         Linf_est = intercept_est/(1-slope_est))

# mutate(k_est = -log(slope_est),
#        Linf_est = intercept_est/(1 - slope_est))


#Sint_set = rnorm(n_runs, mean = Sint_mean, sd = Sint_se)  # adult survival 
#Sl_set = rnorm(n_runs, mean = Sl_mean, sd = Sl_se)  # adult size-survival relationship
#k_connectivity_set = k_connectivity_values$V1

##### Dispersal
# Select values from the 95% CI, weighted by log-likelihood
min_LL <- min(k_theta_allyear_95CI_values$log_like)
dispersal_sample_set <- k_theta_allyear_95CI_values %>%
  mutate(deviation = min_LL/log_like) %>%
  sample_n(size = n_runs, weight = deviation, replacement = T) 
k_connectivity_set <- dispersal_sample_set$k_eval
theta_connectivity_set <- dispersal_sample_set$theta_eval

# # KC weighting code
# minall <- min(mddall_df$log_like)
# mddall_weighted <- mddall_df %>%
#   mutate(deviation=minall/log_like) %>%
#   sample_n(size=1000, weight=deviation, replacement=T) %>%
#   mutate(Year="2012-18")
# 
# # Old method
# k_theta_indices = sample(1:length(k_theta_allyear_95CI_values$k_eval), size=n_runs, replace=FALSE)
# k_connectivity_set = k_theta_allyear_95CI_values$k_eval[k_theta_indices]
# theta_connectivity_set = k_theta_allyear_95CI_values$theta_eval[k_theta_indices]

##### Create survival parameter sets, sampling from within 95% confidence bounds for Cabatoan intercept, addition for each site, and size effect 
site_surv_param_sets <- list()

# Cabatoan first
site_surv_df <- data.frame(Sint_C = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[1], max=best_fit_model_dfs$results$ucl[1]), 
                           Sint_site = 0, 
                           Sl = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[Phi_size_pos], max=best_fit_model_dfs$results$ucl[Phi_size_pos])) %>%
  mutate(Sint = Sint_C + Sint_site)
site_surv_param_sets[[1]] <- site_surv_df

# then rest of sites where survival was estimated
for(i in 2:length(no_space_sites_revisited)) {
  site_surv_df <- data.frame(Sint_C = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[1], max=best_fit_model_dfs$results$ucl[1]), 
                             Sint_site = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[i], max=best_fit_model_dfs$results$ucl[i]), 
                             Sl = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[Phi_size_pos], max=best_fit_model_dfs$results$ucl[Phi_size_pos])) %>%
    mutate(Sint = Sint_C + Sint_site)
  site_surv_param_sets[[i]] <- site_surv_df
}

# then for the three that use average (Caridad Proper, Sitio Lonas, Sitio Tugas)
median_site_Sint_lcl = best_fit_model_dfs$results$lcl[median_site_pos]
median_site_Sint_ucl = best_fit_model_dfs$results$ucl[median_site_pos]

for(i in length(no_space_sites_revisited)+1:length(site_vec_order$site_name)) {
  site_surv_df <- data.frame(Sint_C = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[1], max=best_fit_model_dfs$results$ucl[1]), 
                             Sint_site = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[median_site_pos], max=best_fit_model_dfs$results$ucl[median_site_pos]), 
                             Sl = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[Phi_size_pos], max=best_fit_model_dfs$results$ucl[Phi_size_pos])) %>%
    mutate(Sint = Sint_C + Sint_site)
  site_surv_param_sets[[i]] <- site_surv_df
}

# avg_Sint_lcl = sum(best_fit_model_dfs$results$lcl[1]*16, best_fit_model_dfs$results$lcl[2:16])/16  # avg Sint lcl, all sites included
# avg_Sint_ucl = sum(best_fit_model_dfs$results$ucl[1]*16, best_fit_model_dfs$results$ucl[2:16])/16  # avg Sint ucl, all sites included
# avg_Sint_lcl_withoutCC = sum(best_fit_model_dfs$results$lcl[1]*15, best_fit_model_dfs$results$lcl[3:16])/15  # avg Sint lcl of all sites except Caridad Cemetery  
# avg_Sint_ucl_withoutCC = sum(best_fit_model_dfs$results$ucl[1]*15, best_fit_model_dfs$results$ucl[3:16])/15  # avg Sint lcl of all sites except Caridad Cemetery  
# 
# for(i in length(no_space_sites_revisited)+1:length(site_vec_order$site_name)) {
#   site_surv_df <- data.frame(Sint_C = NA, 
#                              Sint_site = NA, 
#                              Sint = runif(n=n_runs, min=avg_Sint_lcl_withoutCC, max=avg_Sint_ucl_withoutCC),
#                              Sl = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[Phi_size_pos], max=best_fit_model_dfs$results$ucl[Phi_size_pos])) 
#   site_surv_param_sets[[i]] <- site_surv_df
# }

# Make one with median Sint for all sites to use for what-if simulations
site_surv_median_Sint_param_sets <- list()
for(i in 1:length(site_vec_order$site_name)) {
  site_surv_df <- data.frame(Sint_C = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[1], max=best_fit_model_dfs$results$ucl[1]), 
                             Sint_site = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[median_site_pos], max=best_fit_model_dfs$results$ucl[median_site_pos]), 
                             Sl = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[Phi_size_pos], max=best_fit_model_dfs$results$ucl[Phi_size_pos])) %>%
    mutate(Sint = Sint_C + Sint_site)
  site_surv_median_Sint_param_sets[[i]] <- site_surv_df
}

# # Make one with average Sint (without CC) for all sites
# site_surv_avg_Sint_param_sets <- list()
# for(i in 1:length(site_vec_order$site_name)) {
#   site_surv_df <- data.frame(Sint_C = NA, Sint_site = NA, Sl = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[Phi_size_pos], max=best_fit_model_dfs$results$ucl[Phi_size_pos]),
#                              Sint = runif(n=n_runs, min=avg_Sint_lcl_withoutCC, max=avg_Sint_ucl_withoutCC))
#   site_surv_avg_Sint_param_sets[[i]] <- site_surv_df
# }

# Make one for best ests too so can use same code
site_surv_best_est_sets <- list()
# Cabatoan first
site_surv_df <- data.frame(Sint_C = rep(best_fit_model_dfs$results$estimate[1], n_runs),
                           Sint_site = rep(0, n_runs), 
                           Sl = rep(best_fit_model_dfs$results$estimate[Phi_size_pos], n_runs)) %>%
  mutate(Sint = Sint_C + Sint_site)
site_surv_best_est_sets[[1]] <- site_surv_df
# then rest of sites that are revisited
for(i in 2:length(no_space_sites_revisited)) {
  site_surv_df <- data.frame(Sint_C = rep(best_fit_model_dfs$results$estimate[1], n_runs),
                             Sint_site = rep(best_fit_model_dfs$results$estimate[i], n_runs), 
                             Sl = rep(best_fit_model_dfs$results$estimate[Phi_size_pos], n_runs)) %>%
    mutate(Sint = Sint_C + Sint_site)
  site_surv_best_est_sets[[i]] <- site_surv_df
}
# then three that use average
for(i in length(no_space_sites_revisited)+1:length(site_vec_order$site_name)) {
  site_surv_df <- data.frame(Sint_C = rep(best_fit_model_dfs$results$estimate[1], n_runs),
                             Sint_site = rep(best_fit_model_dfs$results$estimate[median_site_pos], n_runs), 
                             Sl = rep(best_fit_model_dfs$results$estimate[Phi_size_pos], n_runs)) %>%
    mutate(Sint = Sint_C + Sint_site)
  site_surv_best_est_sets[[i]] <- site_surv_df
}

# # Make one for best ests too so can use same code
# site_surv_best_est_sets <- list()
# # Cabatoan first
# site_surv_df <- data.frame(Sint_C = rep(best_fit_model_dfs$results$estimate[1], n_runs),
#                            Sint_site = rep(0, n_runs), 
#                            Sl = rep(best_fit_model_dfs$results$estimate[Phi_size_pos], n_runs)) %>%
#   mutate(Sint = Sint_C + Sint_site)
# site_surv_best_est_sets[[1]] <- site_surv_df
# # then rest of sites that are revisited
# for(i in 2:length(no_space_sites_revisited)) {
#   site_surv_df <- data.frame(Sint_C = rep(best_fit_model_dfs$results$estimate[1], n_runs),
#                              Sint_site = rep(best_fit_model_dfs$results$estimate[i], n_runs), 
#                              Sl = rep(best_fit_model_dfs$results$estimate[Phi_size_pos], n_runs)) %>%
#     mutate(Sint = Sint_C + Sint_site)
#   site_surv_best_est_sets[[i]] <- site_surv_df
# }
# # then three that use average
# for(i in length(no_space_sites_revisited)+1:length(site_vec_order$site_name)) {
#   site_surv_df <- data.frame(Sint_C = NA, 
#                              Sint_site = NA, 
#                              Sint = rep(avg_Sint, n_runs),
#                              Sl = rep(best_fit_model_dfs$results$estimate[Phi_size_pos], n_runs))
#   site_surv_best_est_sets[[i]] <- site_surv_df
# }

##### Breeding size (transition to female)
#breeding_size_set = sample(female_sizes$size, n_runs, replace=TRUE)  # transition to female size (replace should be true, right?)
breeding_size_set = sample(recap_first_female$size, n_runs, replace = TRUE)  # transition to female size, pulled from first-observed sizes at F for recaught fish, shoud I make this more of a distribution?

##### Probability of capturing a fish (Pc)
prob_r_beta_params = findBetaDistParams(mean(prob_r), var(prob_r))  # find beta distribution parameters for prob r distrubtion from normal mean and variance
prob_r_set_fodder = rbeta(n_runs, prob_r_beta_params$alpha, prob_r_beta_params$beta, 0)  # should the non-centrality parameter be 0?
prob_r_set_truncated = prob_r_set_fodder[prob_r_set_fodder >= min(prob_r)]  # get rid of any values lower than the lowest observed value, then re-sample from that truncated vector... removes about 100 obs (down to 898 in one case)
prob_r_set = sample(prob_r_set_truncated, n_runs, replace = TRUE)  # sample from that truncated set

##### Recruit size (SHOULD I PULL FROM THE ACTUAL OFFSPRING SIZES, LIKE FOR FEMALE SIZE?)
start_recruit_size_set <- runif(n_runs, min = 3.5, max = 6.0)  # pull from range of offspring sizes
# start_recruit_size_options <- data.frame(recruit_size = c('3.5cm', '4.75cm', '6.0cm', 'mean offspring'),
#                                          size = c(3.5, 4.75, 6.0, mean_sampled_offspring_size), stringsAsFactors = FALSE)

###### Uncertainty in the number of offspring that get matched through parentage analysis (feeds into uncertainty in recruits-per-egg)
n_offspring_parentage_set <- rbinom(n_runs, n_offspring_genotyped, assignment_rate)  # number of assigned offspring using just uncertainty in binomial (assigned/not), for recruits-per-egg est

#################### Generate sets of parameters for uncertainty and do uncertainty runs: ####################
### Make parameter set with different kinds of uncertainty included 

# Uncertainty in start recruit size only
param_set_start_recruit <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = start_recruit_size_set, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,
         # Sl = Sl_mean, Sint = Sint_mean,
         breeding_size = breeding_size_mean, 
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in growth only (both Linf and k)
param_set_growth <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         #k_growth = k_growth_set, s = s, Linf = Linf_set, 
         k_growth = growth_set_params$k_est, s=s, Linf = growth_set_params$Linf_est,
         #Sl = Sl_mean, Sint = Sint_mean,
         breeding_size = breeding_size_mean, 
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in survival only
param_set_survival <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean, 
         #Sl = Sl_set, , Sint = Sint_set,
         breeding_size = breeding_size_mean, 
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in breeding size only
param_set_breeding_size <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean, 
         #Sl = Sl_mean, Sint = Sint_mean,
         breeding_size = breeding_size_set, 
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in offspring assigned to parents (affects recruits-per-egg)
param_set_offspring_assigned <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,
         #Sl = Sl_mean, Sint = Sint_mean,
         breeding_size = breeding_size_mean,
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_parentage_set, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in probability catching a fish (affects recruits-per-egg)
param_set_prob_r <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean, 
         #Sl = Sl_mean, Sint = Sint_mean,
         breeding_size = breeding_size_mean,
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_set, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in both offspring assigned to parents and probability catching a fish (affects recruits-per-egg)
param_set_prob_r_offspring_assigned <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,
         #Sl = Sl_mean, Sint = Sint_mean,
         breeding_size = breeding_size_mean,
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_set, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_parentage_set, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in dispersal only
param_set_dispersal <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,
         #Sl = Sl_mean, Sint = Sint_mean,
         breeding_size = breeding_size_mean, 
         k_connectivity = k_connectivity_set, theta_connectivity = theta_connectivity_set,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in all parameters included for now 
param_set_full <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = start_recruit_size_set, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = growth_set_params$k_est, s = s, Linf = growth_set_params$Linf_est,
         #k_growth = k_growth_set, s = s, Linf = Linf_set, 
         #Sl = Sl_set,  Sint = Sint_set,
         breeding_size = breeding_size_set,
         k_connectivity = k_connectivity_set, theta_connectivity = theta_connectivity_set,  # dispersal kernel parameters
         prob_r = prob_r_set, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_parentage_set, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat)  

### Run metrics for a bunch of different types of uncertainty - No DD
output_uncert_start_recruit <- calcMetricsAcrossRuns(n_runs, param_set_start_recruit, site_surv_best_est_sets, site_dist_info, site_vec_order, "start recruit size", FALSE)
output_uncert_growth <- calcMetricsAcrossRuns(n_runs, param_set_growth, site_surv_best_est_sets, site_dist_info, site_vec_order, "growth", FALSE)
output_uncert_survival <- calcMetricsAcrossRuns(n_runs, param_set_survival, site_surv_param_sets, site_dist_info, site_vec_order, "survival", FALSE)
output_uncert_breeding_size <- calcMetricsAcrossRuns(n_runs, param_set_breeding_size, site_surv_best_est_sets, site_dist_info, site_vec_order, "breeding size", FALSE)
output_uncert_offspring_assigned <- calcMetricsAcrossRuns(n_runs, param_set_offspring_assigned, site_surv_best_est_sets, site_dist_info, site_vec_order, "assigned offspring", FALSE)
output_uncert_prob_r <- calcMetricsAcrossRuns(n_runs, param_set_prob_r, site_surv_best_est_sets, site_dist_info, site_vec_order, "prob r", FALSE)
#output_uncert_prob_r_and_offspring_assigned <- calcMetricsAcrossRuns(n_runs, param_set_prob_r_offspring_assigned, site_surv_best_est_sets, site_dist_info, site_vec_order, "assigned offspring and prob r", FALSE)
output_uncert_dispersal <- calcMetricsAcrossRuns(n_runs, param_set_dispersal, site_surv_best_est_sets, site_dist_info, site_vec_order, "dispersal k", FALSE)
output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info, site_vec_order, "all", FALSE)

# output_uncert_growth <- calcMetricsAcrossRuns(n_runs, param_set_growth, site_surv_best_est_sets, site_dist_info, site_vec_order, "growth", FALSE)
# output_uncert_growth_DD <- calcMetricsAcrossRuns(n_runs, param_set_growth, site_surv_best_est_sets, site_dist_info, site_vec_order, "growth", TRUE)
# output_uncert_all_DD <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info, site_vec_order, "all", TRUE)


# # Join them together for plotting purposes - DO I NEED THIS?
# LEP_uncert <- rbind(output_uncert_start_recruit$LEP_out_df, output_uncert_growth$LEP_out_df,
#                  output_uncert_survival$LEP_out_df, output_uncert_breeding_size$LEP_out_df,
#                  output_uncert_offspring_assigned$LEP_out_df,output_uncert_prob_r$LEP_out_df,
#                  #output_uncert_prob_r_and_offspring_assigned$LEP_out_df, 
#                  output_uncert_dispersal$LEP_out_df, output_uncert_all$LEP_out_df)
# 
# LEP_by_site_uncert <- rbind(output_uncert_start_recruit$LEP_by_site_out_df, output_uncert_growth$LEP_by_site_out_df,
#                                output_uncert_survival$LEP_by_site_out_df, output_uncert_breeding_size$LEP_by_site_out_df,
#                                output_uncert_offspring_assigned$LEP_by_site_out_df, output_uncert_prob_r$LEP_by_site_out_df,
#                                output_uncert_dispersal$LEP_by_site_out_df, output_uncert_all$LEP_by_site_out_df)
# 
# LEP_R_uncert <- rbind(output_uncert_start_recruit$LEP_R_out_df, output_uncert_growth$LEP_R_out_df,
#                     output_uncert_survival$LEP_R_out_df, output_uncert_breeding_size$LEP_R_out_df,
#                     output_uncert_offspring_assigned$LEP_R_out_df, output_uncert_prob_r$LEP_R_out_df,
#                     #output_uncert_prob_r_and_offspring_assigned$LEP_R,
#                     output_uncert_dispersal$LEP_R_out_df, output_uncert_all$LEP_R_out_df)
# 
# LEP_R_by_site_uncert <- rbind(output_uncert_start_recruit$LEP_R_by_site_out_df, output_uncert_growth$LEP_R_by_site_out_df,
#                                  output_uncert_survival$LEP_R_by_site_out_df, output_uncert_breeding_size$LEP_R_by_site_out_df,
#                                  output_uncert_offspring_assigned$LEP_R_by_site_out_df, output_uncert_prob_r$LEP_R_by_site_out_df,
#                                  output_uncert_dispersal$LEP_R_by_site_out_df, output_uncert_all$LEP_R_by_site_out_df)
# 
# RperE_uncert <- rbind(output_uncert_start_recruit$RperE_out_df, output_uncert_growth$RperE_out_df,
#                     output_uncert_survival$RperE_out_df, output_uncert_breeding_size$RperE_out_df,
#                     output_uncert_offspring_assigned$RperE_out_df, output_uncert_prob_r$RperE_out_df,
#                     #output_uncert_prob_r_and_offspring_assigned$RperE_out_df,
#                     output_uncert_dispersal$RperE_out_df, output_uncert_all$RperE_out_df)
# 
# NP_uncert <- rbind(output_uncert_start_recruit$NP_out_df, output_uncert_growth$NP_out_df,
#                    output_uncert_survival$NP_out_df, output_uncert_breeding_size$NP_out_df,
#                    output_uncert_offspring_assigned$NP_out_df, output_uncert_prob_r$NP_out_df,
#                    #output_uncert_prob_r_and_offspring_assigned$NP_out_df,
#                    output_uncert_dispersal$NP_out_df, output_uncert_all$NP_out_df)

###### Now run with DD
### Run metrics for a bunch of different types of uncertainty
output_uncert_start_recruit_DD <- calcMetricsAcrossRuns(n_runs, param_set_start_recruit, site_surv_best_est_sets, site_dist_info, site_vec_order, "start recruit size", TRUE)
output_uncert_growth_DD <- calcMetricsAcrossRuns(n_runs, param_set_growth, site_surv_best_est_sets, site_dist_info, site_vec_order, "growth", TRUE)
output_uncert_survival_DD <- calcMetricsAcrossRuns(n_runs, param_set_survival, site_surv_param_sets, site_dist_info, site_vec_order, "survival", TRUE)
output_uncert_breeding_size_DD <- calcMetricsAcrossRuns(n_runs, param_set_breeding_size, site_surv_best_est_sets, site_dist_info, site_vec_order, "breeding size", TRUE)
output_uncert_offspring_assigned_DD <- calcMetricsAcrossRuns(n_runs, param_set_offspring_assigned, site_surv_best_est_sets, site_dist_info, site_vec_order, "assigned offspring", TRUE)
output_uncert_prob_r_DD <- calcMetricsAcrossRuns(n_runs, param_set_prob_r, site_surv_best_est_sets, site_dist_info, site_vec_order, "prob r", TRUE)
#output_uncert_prob_r_and_offspring_assigned_DD <- calcMetricsAcrossRuns(n_runs, param_set_prob_r_offspring_assigned, site_surv_best_est_sets, site_dist_info, site_vec_order, "assigned offspring and prob r", TRUE)
output_uncert_dispersal_DD <- calcMetricsAcrossRuns(n_runs, param_set_dispersal, site_surv_best_est_sets, site_dist_info, site_vec_order, "dispersal k", TRUE)
output_uncert_all_DD <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info, site_vec_order, "all", TRUE)

#################### Analyze output, prepare for plotting: ####################
# Join together DD runs with different types of uncertainty for plotting purposes (FIGURE OUT WHICH ONES TO KEEP!)
LEP_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_out_df, output_uncert_growth_DD$LEP_out_df,
                    output_uncert_survival_DD$LEP_out_df, output_uncert_breeding_size_DD$LEP_out_df,
                    output_uncert_offspring_assigned_DD$LEP_out_df, output_uncert_prob_r_DD$LEP_out_df,
                    #output_uncert_prob_r_and_offspring_assigned_DD$LEP_out_df, 
                    output_uncert_dispersal_DD$LEP_out_df, output_uncert_all_DD$LEP_out_df)

LEP_by_site_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_by_site_out_df, output_uncert_growth_DD$LEP_by_site_out_df,
                       output_uncert_survival_DD$LEP_by_site_out_df, output_uncert_breeding_size_DD$LEP_by_site_out_df,
                       output_uncert_offspring_assigned_DD$LEP_by_site_out_df, output_uncert_prob_r_DD$LEP_by_site_out_df,
                       output_uncert_dispersal_DD$LEP_by_site_out_df, output_uncert_all_DD$LEP_by_site_out_df)

LEP_R_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_R_out_df, output_uncert_growth_DD$LEP_R_out_df,
                      output_uncert_survival_DD$LEP_R_out_df, output_uncert_breeding_size_DD$LEP_R_out_df,
                      output_uncert_offspring_assigned_DD$LEP_R_out_df, output_uncert_prob_r_DD$LEP_R_out_df,
                      #output_uncert_prob_r_and_offspring_assigned_DD$LEP_R,
                      output_uncert_dispersal_DD$LEP_R_out_df, output_uncert_all_DD$LEP_R_out_df)

LEP_R_by_site_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_R_by_site_out_df, output_uncert_growth_DD$LEP_R_by_site_out_df,
                               output_uncert_survival_DD$LEP_R_by_site_out_df, output_uncert_breeding_size_DD$LEP_R_by_site_out_df,
                               output_uncert_offspring_assigned_DD$LEP_R_by_site_out_df, output_uncert_prob_r_DD$LEP_R_by_site_out_df,
                               output_uncert_dispersal_DD$LEP_R_by_site_out_df, output_uncert_all_DD$LEP_R_by_site_out_df)

RperE_uncert_DD <- rbind(output_uncert_start_recruit_DD$RperE_out_df, output_uncert_growth_DD$RperE_out_df,
                      output_uncert_survival_DD$RperE_out_df, output_uncert_breeding_size_DD$RperE_out_df,
                      output_uncert_offspring_assigned_DD$RperE_out_df, output_uncert_prob_r_DD$RperE_out_df,
                      #output_uncert_prob_r_and_offspring_assigned_DD$RperE_out_df,
                      output_uncert_dispersal_DD$RperE_out_df, output_uncert_all_DD$RperE_out_df)

NP_uncert_DD <- rbind(output_uncert_start_recruit_DD$NP_out_df, output_uncert_growth_DD$NP_out_df,
                   output_uncert_survival_DD$NP_out_df, output_uncert_breeding_size_DD$NP_out_df,
                   output_uncert_offspring_assigned_DD$NP_out_df, output_uncert_prob_r_DD$NP_out_df,
                   #output_uncert_prob_r_and_offspring_assigned_DD$NP_out_df,
                   output_uncert_dispersal_DD$NP_out_df, output_uncert_all_DD$NP_out_df)

##### Find 95% range
## NP with DD
NP_decimal_points <- 2  # bin size of 0.01
lower_edge <- 0.025
upper_edge <- 0.975
lower_index = 25
upper_index = 975

NP_values_rounded <- data.frame(NP_raw = output_uncert_all_DD$NP_out_df$value) %>% 
  mutate(NP_bin = round(NP_raw, 2))
NP_range <- data.frame(table(NP_values_rounded$NP_bin)) %>%
  mutate(rel_freq = Freq/sum(Freq),
         cum_freq = cumsum(rel_freq))
NP_lower <- NP_range$Var1[NP_range$cum_freq >= 0.025][1]

NP_vec <- sort(output_uncert_all_DD$NP_out_df$value)
NP_lower <- round(NP_vec[25], 2)
NP_upper <- round(NP_vec[975], 2)

## LEP averaged across sites
LEP_avg_vec <- sort(output_uncert_all_DD$LEP_out_df$value)
LEP_avg_lower <- round(LEP_avg_vec[lower_index])
LEP_avg_upper <- round(LEP_avg_vec[upper_index])

## LEP at each site
LEP_site_limits <- data.frame(site = site_vec_order$site_name, lower = NA, upper = NA)
for(i in 1:n_sites) {
  LEP_site_vec <- sort((output_uncert_all_DD$LEP_by_site_out_df %>% filter(site == i))$value) 
  LEP_site_limits$lower[i] <- LEP_site_vec[lower_index]
  LEP_site_limits$upper[i] <- LEP_site_vec[upper_index]
}

## Egg-recruit survival
egg_recruit_survival_vec <- sort(output_uncert_all_DD$RperE_out_df$value)
egg_recruit_survival_lower <- egg_recruit_survival_vec[lower_index]
egg_recruit_survival_upper <- egg_recruit_survival_vec[upper_index]

## LEP R averaged across sites
LRP_average_vec <- sort(output_uncert_all_DD$LEP_R_out_df$value)
LRP_average_lower <- LRP_average_vec[lower_index]
LRP_average_upper <- LRP_average_vec[upper_index]

## LRP local averaged across sites
LRP_local_average_vec <- sort(output_uncert_all_DD$LEP_R_local_out_df$value)
LRP_local_average_lower <- LRP_local_average_vec[lower_index]
LRP_local_average_upper <- LRP_local_average_vec[upper_index]

# NP_vec <- output_uncert_all_DD$NP_out_df$value
# NP_breaks <- seq(0, 1, by = 0.)
# > breaks = seq(1.5, 5.5, by=0.5) 
# > duration.cut = cut(duration, breaks, right=FALSE) 
# > duration.freq = table(duration.cut) 

##### Find best estimates
LEP_best_est <- mean(best_est_metrics_mean_offspring$LEP_by_site$LEP)  # mean LEP across sites
LEP_R_best_est <- best_est_metrics_mean_offspring$LEP_R_mean
LEP_R_local_best_est <- best_est_metrics_mean_offspring$LEP_R_local_mean
NP_best_est <- best_est_metrics_mean_offspring$NP
SP_best_est <- best_est_metrics_mean_offspring$SP

LEP_R_best_est_DD <- best_est_metrics_mean_offspring_DD$LEP_R_mean  # mean LRP across sites
LEP_R_local_best_est_DD <- best_est_metrics_mean_offspring_DD$LEP_R_local_mean
NP_best_est_DD <- best_est_metrics_mean_offspring_DD$NP
SP_best_est_DD <- best_est_metrics_mean_offspring_DD$SP

# # Do I need any of these? maybe not...
# LEP_best_est <- mean(best_est_metrics_mean_offspring$LEP_by_site$LEP)  # how to deal with site variation?
# LEP_best_est_min <- min(best_est_metrics_mean_offspring$LEP_by_site$LEP)
# LEP_best_est_max <- max(best_est_metrics_mean_offspring$LEP_by_site$LEP)
# LEP_R_best_est_min <- min(best_est_metrics_mean_offspring$LEP_by_site$LEP_R)
# LEP_R_best_est_max <- max(best_est_metrics_mean_offspring$LEP_by_site$LEP_R)
# LEP_R_best_est <- best_est_metrics_mean_offspring$LEP_R_mean
# NP_best_est <- best_est_metrics_mean_offspring$NP
# SP_best_est <- best_est_metrics_mean_offspring$SP
# LEP_R_local_best_est <- best_est_metrics_mean_offspring$LEP_R_local_mean
# 
# LEP_best_est_DD <- mean(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP)  # how to deal with, show site variation?
# LEP_best_est_min_DD <- min(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP)
# LEP_best_est_max_DD <- max(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP)
# LEP_R_best_est_min_DD <- min(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP_R)
# LEP_R_best_est_max_DD <- max(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP_R)
# LEP_R_best_est_DD <- best_est_metrics_mean_offspring_DD$LEP_R_mean
# NP_best_est_DD <- best_est_metrics_mean_offspring_DD$NP
# SP_best_est_DD <- best_est_metrics_mean_offspring_DD$SP
# LEP_R_local_best_est_DD <- best_est_metrics_mean_offspring_DD$LEP_R_local_mean
# 
# LEP_R_oursites <- recruits_per_egg_oursites*LEP_best_est
# LEP_R_oursites_DD <- recruits_per_egg_oursites_DD*LEP_best_est

# ##### Find range of uncertainty for LEP average and LRP average (with DD)
# LRP_best_est_avg_DD_min <- min(output_uncert_all_DD$LEP_R_out_df$value)
# LRP_best_est_avg_DD_max <- max(output_uncert_all_DD$LEP_R_out_df$value)
# 
# LEP_best_est_avg_min <- min(output_uncert_all_DD$LEP_out_df$value)
# LEP_best_est_avg_max <- max(output_uncert_all_DD$LEP_out_df$value)
# 
# ##### Find range of LR estimates with uncertainty (with DD)
# LR_DD_min <- min(output_uncert_all_DD$LEP_R_local_out_df$value)
# LR_DD_max <- max(output_uncert_all_DD$LEP_R_local_out_df$value)
# 
# #### Find range of LR with all recruits estimates with uncertainty (with DD)
# LR_all_recruits_DD_min <- min(output_uncert_all_offspring_all_DD$LEP_R_local_out_df$value)
# LR_all_recruits_DD_max <- max(output_uncert_all_offspring_all_DD$LEP_R_local_out_df$value)
# 
# # ##### SP estimate range and % > 1 for Haina and Wangag
# SP_Haina_min <- min(output_uncert_all_DD$SP_vals_with_params %>% filter(site == "Haina") %>% select(SP))
# SP_Haina_max <- max(output_uncert_all_DD$SP_vals_with_params %>% filter(site == "Haina") %>% select(SP))
# SP_Wangag_min <- min(output_uncert_all_DD$SP_vals_with_params %>% filter(site == "Wangag") %>% select(SP))
# SP_Wangag_max <- max(output_uncert_all_DD$SP_vals_with_params %>% filter(site == "Wangag") %>% select(SP))
# 
# SP_Haina_perc_above_1 <- output_uncert_all_DD$SP_vals_with_params %>% filter(site == "Haina") %>% filter(SP >= 0.5) %>% summarize(n_estimates = n()) 
# SP_Wangag_perc_above_1 <- output_uncert_all_DD$SP_vals_with_params %>% filter(site == "Wangag") %>% filter(SP >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)
# 
# ##### Find % of metrics above persistence thresholds
# LRP_est_avg_DD_above_1 <- output_uncert_all_DD$LEP_R_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)
# LR_est_DD_above_1 <- output_uncert_all_DD$LEP_R_local_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)
# NP_est_DD_above_1 <- output_uncert_all_DD$NP_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)
# 
# LR_all_recruits_est_DD_above_1 <- output_uncert_all_offspring_all_DD$LEP_R_local_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)

# ########## Find range and % above threshold of values without DD
# # Find range of uncertainty for LEP average and LRP average (without DD)
# LRP_best_est_avg_min <- min(output_uncert_all$LEP_R_out_df$value)
# LRP_best_est_avg_max <- max(output_uncert_all$LEP_R_out_df$value)
# 
# # Find range of Se (without DD)
# RperE_min <- min(output_uncert_all$RperE_out_df$value)
# RperE_max <- max(output_uncert_all$RperE_out_df$value)
# 
# # Find range of LR estimates with uncertainty (without DD)
# LR_min <- min(output_uncert_all$LEP_R_local_out_df$value)
# LR_max <- max(output_uncert_all$LEP_R_local_out_df$value)
# 
# # Find range of NP
# NP_min <- min(output_uncert_all$NP_out_df$value)
# NP_max <- max(output_uncert_all$NP_out_df$value)
# 
# # Find range of LR with all recruits estimates with uncertainty (with DD)
# LR_all_recruits_min <- min(output_uncert_all_offspring_all$LEP_R_local_out_df$value)
# LR_all_recruits_max <- max(output_uncert_all_offspring_all$LEP_R_local_out_df$value)
# 
# # SP estimate range and % > 1 for Haina 
# SP_Haina_min_noDD <- min(output_uncert_all$SP_vals_with_params %>% filter(site == "Haina") %>% select(SP))
# SP_Haina_max_noDD <- max(output_uncert_all$SP_vals_with_params %>% filter(site == "Haina") %>% select(SP))
# 
# SP_Haina_perc_above_1_noDD <- output_uncert_all$SP_vals_with_params %>% filter(site == "Haina") %>% filter(SP >= 0.5) %>% summarize(n_estimates = n()) 
# 
# # Find % of metrics above persistence thresholds
# LRP_est_avg_above_1 <- output_uncert_all$LEP_R_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)
# LR_est_above_1 <- output_uncert_all$LEP_R_local_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)
# NP_est_above_1 <- output_uncert_all$NP_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)
# 
# LR_all_recruits_est_above_1 <- output_uncert_all_offspring_all$LEP_R_local_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)

####### NEED TO EDIT GROWTH UNCERTAINTY INPUTS IN PARAMS AND RE-RUN THESE WHAT-IFS WITH THOSE!!
#################### What-if calculations: ####################
##### What-if calculation 1) what if all genotyped offspring came from the population?

# # Find egg-recruit survival if all offspring we genotyped are included (say they all came from this pop)
# recruits_per_egg_all_offspring <- n_offspring_genotyped/tagged_eggs_6cm

# start at mean size of actual offspring collected (about 4.45)
param_best_est_mean_collected_offspring_all_offspring <- data.frame(t_steps = n_tsteps) %>%
  mutate(min_size = min_size, max_size = max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year_mean = clutches_per_year_mean,  # average fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect =  eyed_effect,  # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,
         #Sl = Sl_mean, , Sint = Sint_mean,
         breeding_size = breeding_size_mean, #recruits_per_egg = recruits_per_egg_all_offspring,
         k_connectivity = k_allyears, theta_connectivity = theta_allyears, # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = 1,  # here, assuming all the offspring arriving came from these sites so don't need to scale it up for the dispersal kernel area (right?)
         offspring_assigned_to_parents = n_offspring_genotyped, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# # NOT SURE WHY THE ONE ABOVE IS WRONG BUT IT GIVES A LOWER NP THAN THE USUAL BEST ESTIMATE!
# param_best_est_mean_collected_offspring_all_offspring <- data.frame(t_steps = n_tsteps) %>%
#   mutate(min_size = min_size, max_size = max_size, n_bins = n_bins,  # LEP-IPM matrix
#          eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year_mean = clutches_per_year_mean,  # average fecundity info
#          egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect =  eyed_effect,  # size-dependent fecundity info
#          start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
#          k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean, 
#          #Sl = Sl_mean, Sint = Sint_mean,
#          breeding_size = breeding_size_mean, #recruits_per_egg = recruits_per_egg_best_est,
#          k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
#          prob_r = prob_r_mean, offspring_assigned_to_parents = n_offspring_genotyped, n_parents = n_parents_genotyped,
#          total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
#          perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Calculate the metrics for the best estimates
best_est_metrics_mean_offspring_all_offspring <- calcMetrics(param_best_est_mean_collected_offspring_all_offspring, site_surv_best_est, site_dist_info, site_vec_order, FALSE)
best_est_metrics_mean_offspring_all_offspring_DD <- calcMetrics(param_best_est_mean_collected_offspring_all_offspring, site_surv_best_est, site_dist_info, site_vec_order, TRUE)

# Save as separate items, for plotting ease
# LEP_best_est_all_offspring <- data.frame(recruit_size = c("3.5cm", "4.75cm", "6.0cm", "mean offspring"),
#                            LEP = c(best_est_metrics_3.5cm_all_offspring$LEP, best_est_metrics_4.75cm_all_offspring$LEP, best_est_metrics_6.0cm_all_offspring$LEP, best_est_metrics_mean_offspring_all_offspring$LEP), stringsAsFactors = FALSE)
# LEP_R_best_est_all_offspring <- data.frame(recruit_size = c("3.5cm", "4.75cm", "6.0cm", "mean offspring"),
#                              LEP_R = c(best_est_metrics_3.5cm_all_offspring$LEP_R, best_est_metrics_4.75cm_all_offspring$LEP_R, best_est_metrics_6.0cm_all_offspring$LEP_R, best_est_metrics_mean_offspring_all_offspring$LEP_R), stringsAsFactors = FALSE)
# NP_best_est_all_offspring <- data.frame(recruit_size = c("3.5cm", "4.75cm", "6.0cm", "mean offspring"),
#                           NP = c(best_est_metrics_3.5cm_all_offspring$NP, best_est_metrics_4.75cm_all_offspring$NP, best_est_metrics_6.0cm_all_offspring$NP, best_est_metrics_mean_offspring_all_offspring$NP), stringsAsFactors = FALSE)
# SP_best_est_all_offspring <- rbind(best_est_metrics_3.5cm_all_offspring$SP %>% mutate(recruit_size = "3.5cm"),
#                                    best_est_metrics_4.75cm_all_offspring$SP %>% mutate(recruit_size = "4.75cm"),
#                                    best_est_metrics_6.0cm_all_offspring$SP %>% mutate(recruit_size = "6.0cm"),
#                                    best_est_metrics_mean_offspring_all_offspring$SP %>% mutate(recruit_size = "mean offspring"))

# Do uncertainty run with those parameters
# Uncertainty in all parameters included for now 
param_set_full_all_offspring <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = start_recruit_size_set, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         #k_growth = k_growth_set, s = s, Linf = Linf_set,
         k_growth = growth_set_params$k_est, s = s, Linf = growth_set_params$Linf_est,
         #Sl = Sl_set, , Sint = Sint_set,
         breeding_size = breeding_size_set,
         k_connectivity = k_connectivity_set, theta_connectivity = theta_connectivity_set,  # dispersal kernel parameters
         prob_r = prob_r_set, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled_best = 1,
         offspring_assigned_to_parents = n_offspring_genotyped, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat)  

output_uncert_all_offspring_all <- calcMetricsAcrossRuns(n_runs, param_set_full_all_offspring, site_surv_param_sets, site_dist_info, site_vec_order, "all: alloff", FALSE)
output_uncert_all_offspring_all_DD <- calcMetricsAcrossRuns(n_runs, param_set_full_all_offspring, site_surv_param_sets, site_dist_info, site_vec_order, "all: alloff", TRUE)

##### What-if calculation 2) What would LRP need to be for NP to be 1? 
#egg_recruit_survival_vec_WI2 <- seq(from = recruits_per_egg_best_est, to = 0.001, by = 0.00001)
LRP_vec_WI2 = seq(from=0.1, to=10.0, by=0.01)  # used to only go up to 5....
NP_output_vec_WI2  <- rep(NA, length(LRP_vec_WI2))

for(i in 1:length(NP_output_vec_WI2)) {
  conn_matR_WI2 = best_est_metrics_mean_offspring$conn_matrix*LRP_vec_WI2[i]
  eig_CR_WI2 = eigen(conn_matR_WI2)
  #NP_WI2 = Re(eig_cR_WI2$values[1])
  #NP_output_vec_WI2[i] = eig_CR_WI2$values[1]
  NP_output_vec_WI2[i] = Re(eig_CR_WI2$values[1])
}

LRP_for_NP = LRP_vec_WI2[which(NP_output_vec_WI2 >= 1)[1]]  # 8.84 on 12/10/19 (3.99 (earlier), 8.64 for WSN)

##### What-if calculation 3) To get the required LRP for NP, what do egg-recruit survival (for our LEP) and LEP (for our egg-recruit-survival) need to be?
#egg_recruit_survival_for_NP <- LRP_for_NP/best_est_metrics_mean_offspring$LEP
egg_recruit_survival_for_NP <- LRP_for_NP/mean(best_est_metrics_mean_offspring$LEP_by_site$LEP)

LEP_for_NP <- LRP_for_NP/best_est_metrics_mean_offspring$recruits_per_egg

LEP_for_NP_DD <- LRP_for_NP/best_est_metrics_mean_offspring_DD$recruits_per_egg

# % of estimates of LRP (with DD) that are >= LRP_for_NP
LRP_ests_above_LRP_for_NP <- output_uncert_all_DD$LEP_R_out_df %>% filter(value >= LRP_for_NP) %>% summarize(n_estimates = n()) %>% select(n_estimates)

##### What-if calculation 4) - how much habitat would we need for the sampling region to be persistent? (code from Sensitivity_to_site_size.R)
# Set width of region
region_width_km <- total_range_of_sampling_area/1000  # length of region (if it was a straight line...)

# Go through different habitat percentages and find site locations, distances, etc. for equally sized and spaced sites  
site_dist_info_habperc <- list()
for(i in 1:length(perc_hab_vals)) {
  site_dist_out_df <- make_output_with_dist(n_sites, perc_hab_vals[i], region_width_km)
  site_dist_out_df$org_site <- (as.data.frame(site_vec_NS, stringsAsFactors = FALSE) %>% slice(rep(1:n(), each=n_sites)))$site_vec_NS  # replace numeric org_site with names
  site_dist_out_df$dest_site <- rep(site_vec_NS, n_sites)  # replace numeric dest_site with names
  site_dist_info_habperc[[i]] <- site_dist_out_df
}

# Make site-survs use numbers as names to match
median_site_Sint <- (site_surv_best_est %>% filter(site == "Elementary School"))$Sint
site_surv_best_est_med_Sint <- site_surv_best_est %>%
  mutate(Sint = median_site_Sint)
# site_surv_best_est_avg_Sint <- site_surv_best_est %>%
#   mutate(Sint = avg_Sint_withoutCC)
#site_surv_best_est_hab_sens$site <- c(5,6,10,15,18,8,3,1,14,13,12,19,4,17,15,2,7,11,9)  # geographical order no. of sites, in the order they are in site_surv_best_est listing

# Find best estimate and uncertainty metrics for those habitat configurations with real site-specific survs (and including compensation for density-dependence)
# perc_hab_best_ests_realSurvs <- list()
# perc_hab_uncertainty_realSurvs <- list()
# for(i in 1:length(perc_hab_vals)) {
#   perc_hab_best_ests_realSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info_habperc[[i]], site_vec_order, TRUE)
#   perc_hab_uncertainty_realSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info_habperc[[i]], site_vec_order, "all: hab sens", TRUE)
# }
# 
# # Make a data frame to plot - actual site-specific survs
# NP_vec_perc_hab_realSurvs_best_est <- perc_hab_best_ests_realSurvs[[1]]$NP  # best estimate NP
# NP_vec_perc_hab_realSurvs_sd <- sd(perc_hab_uncertainty_realSurvs[[1]]$NP_out_df$value)  # sd of NP values with uncertainty
# NP_vec_perc_hab_realSurvs_min <- min(perc_hab_uncertainty_realSurvs[[1]]$NP_out_df$value)  # min of NP values with uncertainty
# NP_vec_perc_hab_realSurvs_max <- max(perc_hab_uncertainty_realSurvs[[1]]$NP_out_df$value)  # max of NP values with uncertainty
# 
# for(i in 2:length(perc_hab_vals)) {
#   NP_vec_perc_hab_realSurvs_best_est <- c(NP_vec_perc_hab_realSurvs_best_est, perc_hab_best_ests_realSurvs[[i]]$NP)
#   NP_vec_perc_hab_realSurvs_sd <- c(NP_vec_perc_hab_realSurvs_sd, sd(perc_hab_uncertainty_realSurvs[[i]]$NP_out_df$value))
#   NP_vec_perc_hab_realSurvs_min <- c(NP_vec_perc_hab_realSurvs_min, min(perc_hab_uncertainty_realSurvs[[i]]$NP_out_df$value))
#   NP_vec_perc_hab_realSurvs_max <- c(NP_vec_perc_hab_realSurvs_max, max(perc_hab_uncertainty_realSurvs[[i]]$NP_out_df$value))
# }
# 
# NP_by_perc_hab_realSurvs <- data.frame(perc_hab = perc_hab_vals,
#                                        NP = NP_vec_perc_hab_realSurvs_best_est,
#                                        NP_sd = NP_vec_perc_hab_realSurvs_sd,
#                                        NP_min = NP_vec_perc_hab_realSurvs_min,
#                                        NP_max = NP_vec_perc_hab_realSurvs_max, stringsAsFactors = FALSE)
# 
# Find best estimate and uncertainty metrics for those habitat configurations with average-site survs (and including compensation for density-dependence)
perc_hab_best_ests_avgSurvs <- list()
perc_hab_uncertainty_avgSurvs <- list()
for(i in 1:length(perc_hab_vals)) {
  perc_hab_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_med_Sint, site_dist_info_habperc[[i]], site_vec_order, TRUE)
  perc_hab_uncertainty_avgSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_median_Sint_param_sets, site_dist_info_habperc[[i]], site_vec_order, "all: hab sens, avg surv", TRUE)
}

# ###### 
# # Make a data frame to plot - avg survs (all sites have the same surv, pulled from range each run in uncertainty runs)
# NP_vec_perc_hab_avgSurvs_best_est <- perc_hab_best_ests_avgSurvs[[1]]$NP  # best estimate NP
# NP_vec_perc_hab_avgSurvs_sd <- sd(perc_hab_uncertainty_avgSurvs[[1]]$NP_out_df$value)  # sd of NP values with uncertainty
# NP_vec_perc_hab_avgSurvs_min <- min(perc_hab_uncertainty_avgSurvs[[1]]$NP_out_df$value)  # min of NP values with uncertainty
# NP_vec_perc_hab_avgSurvs_max <- max(perc_hab_uncertainty_avgSurvs[[1]]$NP_out_df$value)  # max of NP values with uncertainty
# 
# for(i in 2:(length(perc_hab_vals))) {
#   NP_vec_perc_hab_avgSurvs_best_est <- c(NP_vec_perc_hab_avgSurvs_best_est, perc_hab_best_ests_avgSurvs[[i]]$NP)
#   NP_vec_perc_hab_avgSurvs_sd <- c(NP_vec_perc_hab_avgSurvs_sd, sd(perc_hab_uncertainty_avgSurvs[[i]]$NP_out_df$value))
#   NP_vec_perc_hab_avgSurvs_min <- c(NP_vec_perc_hab_avgSurvs_min, min(perc_hab_uncertainty_avgSurvs[[i]]$NP_out_df$value))
#   NP_vec_perc_hab_avgSurvs_max <- c(NP_vec_perc_hab_avgSurvs_max, max(perc_hab_uncertainty_avgSurvs[[i]]$NP_out_df$value))
# }
# 
# NP_by_perc_hab_avgSurvs <- data.frame(perc_hab = perc_hab_vals,
#                                        NP = NP_vec_perc_hab_avgSurvs_best_est,
#                                        NP_sd = NP_vec_perc_hab_avgSurvs_sd,
#                                        NP_min = NP_vec_perc_hab_avgSurvs_min,
#                                        NP_max = NP_vec_perc_hab_avgSurvs_max, stringsAsFactors = FALSE)
# 
# # Find percent above 1 for what-if
# NP_above1_perc_hab_realSurvs <- data.frame(perc_hab = perc_hab_vals, perc_persistent = NA)
# NP_above1_perc_hab_avgSurvs <- data.frame(perc_hab = perc_hab_vals, perc_persistent = NA)
# 
# for(i in 1:(length(perc_hab_vals))) {
#   #NP_above1_perc_hab_realSurvs$perc_persistent[i] = sum(perc_hab_uncertainty_realSurvs[[i]]$NP_out_df$value >= 1)/n_runs
#   NP_above1_perc_hab_avgSurvs$perc_persistent[i] = sum(perc_hab_uncertainty_avgSurvs[[i]]$NP_out_df$value >= 1)/n_runs
# }

# Make dfs for plot where each run is its own line - perc hab on x axis, NP on y axis
perc_hab_plot_df <- perc_hab_uncertainty_avgSurvs[[1]]$NP_out_df %>% select(value,run) %>% mutate(perc_hab = perc_hab_vals[1])  # one for uncertainty runs
for(i in 2:length(perc_hab_vals)) {
  perc_hab_plot_df <- rbind(perc_hab_plot_df, perc_hab_uncertainty_avgSurvs[[i]]$NP_out_df %>% select(value,run) %>% mutate(perc_hab = perc_hab_vals[i]))
}
perc_hab_best_est_NP_df <- data.frame(perc_hab = rep(NA,length(perc_hab_vals)), value = rep(NA,length(perc_hab_vals)))  # and one for best ests
for(i in 1:length(perc_hab_vals)) {
  perc_hab_best_est_NP_df$perc_hab[i] = perc_hab_vals[i]
  perc_hab_best_est_NP_df$value[i] = perc_hab_best_ests_avgSurvs[[i]]$NP
}


##### What if - same habitat density, wider region - should re-do this so the patches stay in the center, just extend the region!!
region_width_list <- c(region_width_km, 35, 40, 45, 50, 55)

# Go through different region widths for the real habitat percentage and find site locations, distances, etc. for equally sized and spaced sites
site_dist_info_wider_region <- list()
for(i in 1:length(region_width_list)) {
  site_dist_out_df <- make_output_with_dist(n_sites, Ps, region_width_list[i])
  site_dist_out_df$org_site <- (as.data.frame(site_vec_NS, stringsAsFactors = FALSE) %>% slice(rep(1:n(), each=n_sites)))$site_vec_NS  # replace numeric org_site with names
  site_dist_out_df$dest_site <- rep(site_vec_NS, n_sites)  # replace numeric dest_site with names
  site_dist_info_wider_region[[i]] <- site_dist_out_df
}

# Find best estimate and uncertainty metrics for those habitat configurations with average-site survs (and including compensation for density-dependence)
wider_region_best_ests_avgSurvs <- list()
wider_region_uncertainty_avgSurvs <- list()
for(i in 1:length(region_width_list)) {
  wider_region_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_med_Sint, site_dist_info_wider_region[[i]], site_vec_order, TRUE)
  wider_region_uncertainty_avgSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_median_Sint_param_sets, site_dist_info_wider_region[[i]], site_vec_order, "all: region sens, avg surv", TRUE)
}

# # Make a data frame to plot - avg survs (all sites have the same surv, pulled from range each run in uncertainty runs)
# NP_vec_wider_region_avgSurvs_best_est <- wider_region_best_ests_avgSurvs[[1]]$NP  # best estimate NP
# NP_vec_wider_region_avgSurvs_sd <- sd(wider_region_uncertainty_avgSurvs[[1]]$NP_out_df$value)  # sd of NP values with uncertainty
# NP_vec_wider_region_avgSurvs_min <- min(wider_region_uncertainty_avgSurvs[[1]]$NP_out_df$value)  # min of NP values with uncertainty
# NP_vec_wider_region_avgSurvs_max <- max(wider_region_uncertainty_avgSurvs[[1]]$NP_out_df$value)  # max of NP values with uncertainty
# 
# for(i in 2:(length(region_width_list))) {
#   NP_vec_wider_region_avgSurvs_best_est <- c(NP_vec_wider_region_avgSurvs_best_est, wider_region_best_ests_avgSurvs[[i]]$NP)
#   NP_vec_wider_region_avgSurvs_sd <- c(NP_vec_wider_region_avgSurvs_sd, sd(wider_region_uncertainty_avgSurvs[[i]]$NP_out_df$value))
#   NP_vec_wider_region_avgSurvs_min <- c(NP_vec_wider_region_avgSurvs_min, min(wider_region_uncertainty_avgSurvs[[i]]$NP_out_df$value))
#   NP_vec_wider_region_avgSurvs_max <- c(NP_vec_wider_region_avgSurvs_max, max(wider_region_uncertainty_avgSurvs[[i]]$NP_out_df$value))
# }
# 
# NP_by_wider_region_avgSurvs <- data.frame(region_width = region_width_list,
#                                       NP = NP_vec_wider_region_avgSurvs_best_est,
#                                       NP_sd = NP_vec_wider_region_avgSurvs_sd,
#                                       NP_min = NP_vec_wider_region_avgSurvs_min,
#                                       NP_max = NP_vec_wider_region_avgSurvs_max, stringsAsFactors = FALSE)
# 
# # Find percent above 1 for what-if
# NP_above1_wider_region_realSurvs <- data.frame(region_width = region_width_list, perc_persistent = NA)
# NP_above1_wider_region_avgSurvs <- data.frame(region_width = region_width_list, perc_persistent = NA)
# 
# for(i in 1:(length(region_width_list))) {
#   NP_above1_wider_region_avgSurvs$perc_persistent[i] = sum(wider_region_uncertainty_avgSurvs[[i]]$NP_out_df$value >= 1)/n_runs
# }

# Make dfs for plot where each run is it's own line - region width on x axis, NP on y axis
wider_region_plot_df <- wider_region_uncertainty_avgSurvs[[1]]$NP_out_df %>% select(value,run) %>% mutate(region_width = region_width_list[1])  # one for uncertainty runs
for(i in 2:length(region_width_list)) {
  wider_region_plot_df <- rbind(wider_region_plot_df, wider_region_uncertainty_avgSurvs[[i]]$NP_out_df %>% select(value,run) %>% mutate(region_width = region_width_list[i]))
}
wider_region_best_est_NP_df <- data.frame(region_width = rep(NA,length(region_width_list)), value = rep(NA,length(region_width_list)))  # and one for best ests
for(i in 1:length(region_width_list)) {
  wider_region_best_est_NP_df$region_width[i] = region_width_list[i]
  wider_region_best_est_NP_df$value[i] = wider_region_best_ests_avgSurvs[[i]]$NP
}

##### What if - 100% habitat, wider region
# Go through different region widths for the real habitat percentage and find site locations, distances, etc. for equally sized and spaced sites
site_dist_info_wider_region_all_hab <- list()
for(i in 1:length(region_width_list)) {
  site_dist_out_df <- make_output_with_dist(n_sites, 1.0, region_width_list[i])
  site_dist_out_df$org_site <- (as.data.frame(site_vec_NS, stringsAsFactors = FALSE) %>% slice(rep(1:n(), each=n_sites)))$site_vec_NS  # replace numeric org_site with names
  site_dist_out_df$dest_site <- rep(site_vec_NS, n_sites)  # replace numeric dest_site with names
  site_dist_info_wider_region_all_hab[[i]] <- site_dist_out_df
}

# Find best estimate and uncertainty metrics for those habitat configurations with average-site survs (and including compensation for density-dependence)
wider_region_all_hab_best_ests_avgSurvs <- list()
wider_region_all_hab_uncertainty_avgSurvs <- list()
for(i in 1:length(region_width_list)) {
  wider_region_all_hab_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_med_Sint, site_dist_info_wider_region_all_hab[[i]], site_vec_order, TRUE)
  wider_region_all_hab_uncertainty_avgSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_median_Sint_param_sets, site_dist_info_wider_region_all_hab[[i]], site_vec_order, "all: region sens, all hab", TRUE)
}

# # Make a data frame to plot - avg survs (all sites have the same surv, pulled from range each run in uncertainty runs)
# NP_vec_wider_region_all_hab_avgSurvs_best_est <- wider_region_all_hab_best_ests_avgSurvs[[1]]$NP  # best estimate NP
# NP_vec_wider_region_all_hab_avgSurvs_sd <- sd(wider_region_all_hab_uncertainty_avgSurvs[[1]]$NP_out_df$value)  # sd of NP values with uncertainty
# NP_vec_wider_region_all_hab_avgSurvs_min <- min(wider_region_all_hab_uncertainty_avgSurvs[[1]]$NP_out_df$value)  # min of NP values with uncertainty
# NP_vec_wider_region_all_hab_avgSurvs_max <- max(wider_region_all_hab_uncertainty_avgSurvs[[1]]$NP_out_df$value)  # max of NP values with uncertainty
# 
# for(i in 2:(length(region_width_list))) {
#   NP_vec_wider_region_all_hab_avgSurvs_best_est <- c(NP_vec_wider_region_all_hab_avgSurvs_best_est, wider_region_all_hab_best_ests_avgSurvs[[i]]$NP)
#   NP_vec_wider_region_all_hab_avgSurvs_sd <- c(NP_vec_wider_region_all_hab_avgSurvs_sd, sd(wider_region_all_hab_uncertainty_avgSurvs[[i]]$NP_out_df$value))
#   NP_vec_wider_region_all_hab_avgSurvs_min <- c(NP_vec_wider_region_all_hab_avgSurvs_min, min(wider_region_all_hab_uncertainty_avgSurvs[[i]]$NP_out_df$value))
#   NP_vec_wider_region_all_hab_avgSurvs_max <- c(NP_vec_wider_region_all_hab_avgSurvs_max, max(wider_region_all_hab_uncertainty_avgSurvs[[i]]$NP_out_df$value))
# }
# 
# NP_by_wider_region_all_hab_avgSurvs <- data.frame(region_width = region_width_list,
#                                           NP = NP_vec_wider_region_all_hab_avgSurvs_best_est,
#                                           NP_sd = NP_vec_wider_region_all_hab_avgSurvs_sd,
#                                           NP_min = NP_vec_wider_region_all_hab_avgSurvs_min,
#                                           NP_max = NP_vec_wider_region_all_hab_avgSurvs_max, stringsAsFactors = FALSE)
# 
# # Find percent above 1 for what-if
# NP_above1_wider_region_all_hab_realSurvs <- data.frame(region_width = region_width_list, perc_persistent = NA)
# NP_above1_wider_region_all_hab_avgSurvs <- data.frame(region_width = region_width_list, perc_persistent = NA)
# 
# for(i in 1:(length(region_width_list))) {
#   NP_above1_wider_region_all_hab_avgSurvs$perc_persistent[i] = sum(wider_region_all_hab_uncertainty_avgSurvs[[i]]$NP_out_df$value >= 1)/n_runs
# }

# Make dfs for plot where each run is it's own line - region width on x axis, NP on y axis
wider_region_all_hab_plot_df <- wider_region_all_hab_uncertainty_avgSurvs[[1]]$NP_out_df %>% select(value,run) %>% mutate(region_width = region_width_list[1])  # one for uncertainty runs
for(i in 2:length(region_width_list)) {
  wider_region_all_hab_plot_df <- rbind(wider_region_all_hab_plot_df, wider_region_all_hab_uncertainty_avgSurvs[[i]]$NP_out_df %>% select(value,run) %>% mutate(region_width = region_width_list[i]))
}
wider_region_all_hab_best_est_NP_df <- data.frame(region_width = rep(NA,length(region_width_list)), value = rep(NA,length(region_width_list)))  # and one for best ests
for(i in 1:length(region_width_list)) {
  wider_region_all_hab_best_est_NP_df$region_width[i] = region_width_list[i]
  wider_region_all_hab_best_est_NP_df$value[i] = wider_region_all_hab_best_ests_avgSurvs[[i]]$NP
}

##### What if for % habitat and region width combined
# Find the site dist info
perc_hab_vals_long_form <- seq(from=0.01, to=1.0, by = 0.01)
region_width_list_long_form <- seq(from = 20, to = 55, by = 1)

site_dist_info_wider_region_perc_hab <- list()
list_val = 1
for(i in 1:length(region_width_list_long_form)) {
  for(j in 1:length(perc_hab_vals_long_form)) {
    site_dist_out_df <- make_output_with_dist(n_sites, perc_hab_vals_long_form[j], region_width_list_long_form[i])
    site_dist_out_df$org_site <- (as.data.frame(site_vec_NS, stringsAsFactors = FALSE) %>% slice(rep(1:n(), each=n_sites)))$site_vec_NS  # replace numeric org_site with names
    site_dist_out_df$dest_site <- rep(site_vec_NS, n_sites)  # replace numeric dest_site with names
    site_dist_info_wider_region_perc_hab[[list_val]] <- site_dist_out_df
    list_val = list_val + 1
  }
}

# Find best estimate and uncertainty metrics for those habitat configurations with average-site survs (and including compensation for density-dependence)
wider_region_perc_hab_best_ests_avgSurvs <- list()
#wider_region_perc_hab_uncertainty_avgSurvs <- list()
for(i in 1:(length(region_width_list_long_form)*length(perc_hab_vals_long_form))) {
  wider_region_perc_hab_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_med_Sint, site_dist_info_wider_region_perc_hab[[i]], site_vec_order, TRUE)
  #  wider_region_perc_hab_uncertainty_avgSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_avg_Sint_param_sets, site_dist_info_wider_region_perc_hab[[i]], site_vec_order, "all: region sens, perc hab sens, avg surv", TRUE)
  print(i)
}

# Make dfs for plot where region width, perc hab, and NP all listed
wider_region_perc_hab_plot_df <- data.frame(NP = rep(NA, length(region_width_list_long_form)*length(perc_hab_vals_long_form)))
for(i in 1:(length(region_width_list_long_form)*length(perc_hab_vals_long_form))) {
  wider_region_perc_hab_plot_df$NP[i] = wider_region_perc_hab_best_ests_avgSurvs[[i]]$NP
}
wider_region_perc_hab_plot_df$region_width = (as.data.frame(region_width_list_long_form, stringsAsFactors = FALSE) %>% slice(rep(1:n(), each = length(perc_hab_vals_long_form))))$region_width_list_long_form
wider_region_perc_hab_plot_df$perc_hab = rep(perc_hab_vals_long_form, length(region_width_list_long_form))

# Put NP into categories: 0-.2, 0.2-0.4, 0.4-0.6, 0.6-0.8, 0.8-1.0, >1.0
wider_region_perc_hab_plot_df_categories <- wider_region_perc_hab_plot_df %>%
  mutate(NP_category = case_when(NP < 0.2 ~ "< 0.2",
                                 (0.2 <= NP & NP < 0.4) ~ "[0.2 - 0.4)",
                                 (0.4 <= NP & NP < 0.6) ~ "[0.4 - 0.6)",
                                 (0.6 <= NP & NP < 0.8) ~ "[0.6 - 0.8)",
                                 (0.8 <= NP & NP < 1) ~ "[0.8 - 1)",
                                 (NP >= 1) ~ ">= 1")) %>%
  mutate(order_6 = case_when(NP < 0.2 ~ 1,
                           (0.2 <= NP & NP < 0.4) ~ 2,
                           (0.4 <= NP & NP < 0.6) ~ 3,
                           (0.6 <= NP & NP < 0.8) ~ 4,
                           (0.8 <= NP & NP < 1) ~ 5,
                           (NP >= 1) ~ 6)) %>%
  mutate(NP_two_category = case_when(NP < 1 ~ "< 1",
                                     NP >= 1 ~ ">= 1")) %>%
  mutate(order_2 = case_when(NP < 1 ~ 1,
                             NP >= 1 ~ 2)) %>%
  mutate(NP_categories = case_when(NP < 0.25 ~ "< 0.25",
                                 (0.25 <= NP & NP < 0.5) ~ "[0.25 - 0.5)",
                                 (0.5 <= NP & NP < 0.75) ~ "[0.5 - 0.75)",
                                 (0.75 <= NP & NP < 1) ~ "[0.75 - 1)",
                                 (NP >= 1) ~ ">= 1")) %>%
  mutate(order_5 = case_when(NP < 0.25 ~ 1,
                                   (0.25 <= NP & NP < 0.5) ~ 2,
                                   (0.5 <= NP & NP < 0.75) ~ 3,
                                   (0.75 <= NP & NP < 1) ~ 4,
                                   (NP >= 1) ~ 5))

save(wider_region_perc_hab_best_ests_avgSurvs, file=here::here("Data/Script_outputs","wider_region_perc_hab_best_ests_avgSurvs.RData"))
save(wider_region_perc_hab_plot_df, file=here::here("Data/Script_outputs","wider_region_perc_hab_plot_df.RData"))
#save(wider_region_perc_hab_uncertainty_avgSurvs, file=here::here("Data/Script_outputs","wider_region_perc_hab_uncertainty_avgSurvs.RData"))
#save(site_surv_avg_Sint_param_sets, file=here::here("Data/Script_outputs","site_surv_avg_Sint_param_sets.RData"))

##### What if larvae could navigate? Integrate kernel from 500 or 1000m on either side of each patch
# Load generated parameter sets
load(file = here::here("Data/Script_outputs", "param_set_full.RData"))
load(file = here::here("Data/Script_outputs","site_surv_param_sets.RData"))

# # Re-estimate number of offspring, where scaling for loss to dispersal to non-habitat is altered for the amount of navigation
findPropHabWithLarvalNav <- function(nav_buffer, site_buffer_info_df, site_width_info_df, sampling_range) {  # nav_buffer = distance larvae can navigate in m, site_buffer_info_df = site_buffer_info, site_width_info_df = site_width_info, sampling_range = total_range_of_sampling_area
  # total sum of site widths without buffer
  sum_of_site_widths <- sum(site_width_info_df$width_m)
  buffer_to_add_N <- rep(NA, length(site_width_info_df$site))
  buffer_to_add_S <- rep(NA, length(site_width_info_df$site))
  
  # go through sites, check if buffer is less than max amount, then add either buffer or max
  for(i in 1:length(site_width_info$site)){
    site_val <- site_width_info$site[i]
    max_buffer_N <- (site_buffer_info_df %>% filter(site == site_val))$max_buffer_N_m
    max_buffer_S <- (site_buffer_info_df %>% filter(site == site_val))$max_buffer_S_m
    
    # Check north buffer
    if(nav_buffer > max_buffer_N) {
      buffer_to_add_N[i] <- max_buffer_N  # if bigger than max buffer, just add the max
    } else {
      buffer_to_add_N[i] <- nav_buffer  # otherwise, add the navigation buffer
    }
    
    # Check south buffer
    if(nav_buffer > max_buffer_S) {
      buffer_to_add_S[i] <- max_buffer_S  # if bigger than max buffer, just add the max
    } else {
      buffer_to_add_S[i] <- nav_buffer  # otherwise, add the navigation buffer
    }
  }
  
  total_hab <- sum_of_site_widths + sum(buffer_to_add_N) + sum(buffer_to_add_S)
  
  prop_hab <- total_hab/total_range_of_sampling_area
  
  return(prop_hab)
}

# Make new param set
makeParamSetWithSampledHab <- function(param_set, Ps) {  # param_set is param_best_est_mean_collected_offspring, Ps is prop_hab from findPropHabWithLarvalNav
 param_set_out <- param_set
 param_set_out$prop_hab <- Ps
 return(param_set_out)
}

# # (Rough) proportion of sampling area that is habitat
# total_sum_of_site_widths <- sum(site_width_info$width_m)
# prop_sampling_area_habitat <- total_sum_of_site_widths/total_range_of_sampling_area
# Ps <- prop_sampling_area_habitat

# Re-estimate the distances data frames
findSiteDistsWithNavBuffer <- function(nav_buffer, site_dist_info_df, site_buffer_info_df) {  # nav_buffer = distance larvae can navigate in km, site_dist_info_df = site_dist_info, site_buffer_info_df = site_buffer_info
  site_dist_info_out <- site_dist_info_df  # modify existing site_dist_info with buffers built into d1 and d2
  # Go through site combos (one row in site_dist_info)
  # For each set of sites, decide if destination is 1) the same as the origin, 2) north of the origin, or 3) south of the origin
  for(i in 1:length(site_dist_info_df$org_site)) {
    site_org_val <- site_dist_info_df$org_site[i]
    site_dest_val <- site_dist_info_df$dest_site[i]
    site_org_geo <- site_dist_info_df$org_geo_order[i]
    site_dest_geo <- site_dist_info_df$dest_geo_order[i]
    
    # if the origin and destination sites are the same: d1 is unchanged, d2 is d2 + the buffer up to the max buffer in either direction for that site (which leaves a bit on the table on the other side but...)
    if(site_org_val == site_dest_val) {  
      max_buffer_N_or_S <- min((site_buffer_info_df %>% filter(site == site_org_val))$max_buffer_N_km, (site_buffer_info_df %>% filter(site == site_org_val))$max_buffer_S_km)  # find min max buffer to N or S
      if(nav_buffer > max_buffer_N_or_S) {
        site_dist_info_out$d2_km[i] <- site_dist_info_df$d2_km[i] + max_buffer_N_or_S  # if the nav buffer is greater than the max buffer, use max buffer
      } else {
        site_dist_info_out$d2_km[i] <- site_dist_info_df$d2_km[i] + nav_buffer
      }
    } else if(site_org_geo > site_dest_geo) {  # if the destination is north of the origin...
      
      # buffer d1 (distance to the close side of the destination site)
      max_buffer_S <- (site_buffer_info_df %>% filter(site == site_dest_val))$max_buffer_S_km  # find the max buffer to the south site of the destination
      if(nav_buffer > max_buffer_S) {
        site_dist_info_out$d1_km[i] <- site_dist_info_df$d1_km[i] - max_buffer_S  # d1 is d1 - the buffer up to the max buffer S for the destination site
      } else {
        site_dist_info_out$d1_km[i] <- site_dist_info_df$d1_km[i] - nav_buffer
      }
      
      # buffer d2 (distance to the far side of the destination site)
      max_buffer_N <- (site_buffer_info_df %>% filter(site == site_dest_val))$max_buffer_N_km  # find the max buffer to the north site of the destination
      if(nav_buffer > max_buffer_N) {
        site_dist_info_out$d2_km[i] <- site_dist_info_df$d2_km[i] + max_buffer_N  # d2 is d2 + the buffer up to the max buffer N for the destination site
      } else {
        site_dist_info_out$d2_km[i] <- site_dist_info_df$d2_km[i] + nav_buffer
      }
    } else if(site_org_geo < site_dest_geo) {  # if the destination site is to the south of the origin...
      
      # buffer d1 (distance to the close site of the destination site), d1 is d1 - the buffer up to the max buffer N for the destination site
      max_buffer_N <- (site_buffer_info_df %>% filter(site == site_dest_val))$max_buffer_N_km  # find the max buffer to the north site of the destination
      if(nav_buffer > max_buffer_N) {  
        site_dist_info_out$d1_km[i] <- site_dist_info_df$d1_km[i] - max_buffer_N  # d1 is d1 - the buffer up to the max buffer N for the destination site
      } else {
        site_dist_info_out$d1_km[i] <- site_dist_info_df$d1_km[i] - nav_buffer
      }
      
      # buffer d2 (distance to far side of destination site), d2 is d2 + the buffer up to the max buffer S for the destination site
      max_buffer_S <- (site_buffer_info_df %>% filter(site == site_dest_val))$max_buffer_S_km  # find the max buffer to the south site of the destination
      if(nav_buffer > max_buffer_S) {  
        site_dist_info_out$d2_km[i] <- site_dist_info_df$d2_km[i] + max_buffer_S  # d2 is d2 + the buffer up to the max buffer S for the destination site
      } else {
        site_dist_info_out$d2_km[i] <- site_dist_info_df$d2_km[i] + nav_buffer
      }
    }
  }
  # If it's the same as the origin, d1 is unchanged, d2 is d2 + the buffer up to the max buffer in either direction for that site (which leaves a bit on the table on the other side but...)
  # If it's north of the origin, d1 is d1 - the buffer up to the max buffer S for the destination site
  # If it's north of the origin, d2 is d2 + the buffer up to the max buffer N for the destination site
  # If it's south of the origin, d1 is d1 - the buffer up to the max buffer N for the destination site
  # If it's south of the origin, d2 is d2 + the buffer up to the max buffer S for the destination site
  
  return(site_dist_info_out %>% select(-dist_mid_to_S_m, -dist_mid_to_S_km, -dist_mid_to_N_m, -dist_mid_to_N_km))
}

# Go through nav distances and find parameter sets and distance data frames
nav_buffer_vec <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

larv_nav_prop_hab <- rep(NA,length(nav_buffer_vec))
larv_nav_params_list <- list()
larv_nav_site_dists_list <- list()
larv_nav_best_ests <- list()
larv_nav_uncert_all <- list()

for(i in 1:length(nav_buffer_vec)) {
  # Find prop hab for scaling
  larv_nav_prop_hab[i] <- findPropHabWithLarvalNav(nav_buffer_vec[i]*1000, site_buffer_info, site_width_info, total_range_of_sampling_area)
  # Make new param set
  larv_nav_params_list[[i]] <- makeParamSetWithSampledHab(param_best_est_mean_collected_offspring, larv_nav_prop_hab[i])
  # Find distances
  larv_nav_site_dists_list[[i]] <- findSiteDistsWithNavBuffer(nav_buffer_vec[i], site_dist_info, site_buffer_info)
  # Make new param set with uncertainty with larv_nav prop_hab
  param_set_full_larv_nav <- param_set_full
  param_set_full_larv_nav$prop_hab <- rep(larv_nav_prop_hab[i],n_runs)
  # Do best ests
  larv_nav_best_ests[[i]] <- calcMetrics(larv_nav_params_list[[i]], site_surv_best_est, larv_nav_site_dists_list[[i]], site_vec_order, "TRUE")
  larv_nav_uncert_all[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full_larv_nav, site_surv_param_sets, larv_nav_site_dists_list[[i]], site_vec_order, "all", TRUE)
}

# Make dfs for plot - km larv nav on x axis, NP on y axis
# Make dfs
larv_nav_plot_df <- larv_nav_uncert_all[[1]]$NP_out_df %>% select(value,run) %>% mutate(larv_nav = nav_buffer_vec[1])  # one for uncertainty runs
for(i in 2:length(nav_buffer_vec)) {
  larv_nav_plot_df <- rbind(larv_nav_plot_df, larv_nav_uncert_all[[i]]$NP_out_df %>% select(value,run) %>% mutate(larv_nav = nav_buffer_vec[i]))
}
larv_nav_best_est_NP_df <- data.frame(larv_nav = rep(NA,length(nav_buffer_vec)), value = rep(NA,length(nav_buffer_vec)))  # and one for best ests
for(i in 1:length(nav_buffer_vec)) {
  larv_nav_best_est_NP_df$larv_nav[i] = nav_buffer_vec[i]
  larv_nav_best_est_NP_df$value[i] = larv_nav_best_ests[[i]]$NP
}

# # Test new function
# test_site_dist <- findSiteDistsWithNavBuffer(1, site_dist_info, site_buffer_info)
# 
# # Create new site distance data frames
# site_dist_info_200m <- site_dist_info %>%
#   mutate(d1_km = if_else(d1_km-0.2 > 0, d1_km - 0.2, 0),  # near edge of site is 200m closer (but can't be negative...)
#          d2_km = d2_km + 0.2)  # far edge of site is 200m farther away
# site_dist_info_500m <- site_dist_info %>%
#   mutate(d1_km = if_else(d1_km-0.5 > 0, d1_km - 0.5, 0),  # near edge of site is 500m closer (but can't be negative...)
#          d2_km = d2_km + 0.5)  # far edge of site is 500m farther away
# site_dist_info_1000m <- site_dist_info %>%
#   mutate(d1_km = if_else(d1_km - 1 > 0, d1_km - 1, 0),  # near edge of site is 1000m closer
#          d2_km = d2_km + 1)
# 
# # Do runs (with DD)
# best_est_metrics_larv_nav_200m_DD <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info_200m, site_vec_order, "TRUE")
# best_est_metrics_larv_nav_500m_DD <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info_500m, site_vec_order, "TRUE")
# best_est_metrics_larv_nav_1000m_DD <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info_1000m, site_vec_order, "TRUE")
# output_uncert_all_DD_larv_nav_200m <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info_200m, site_vec_order, "all", TRUE)
# output_uncert_all_DD_larv_nav_500m <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info_500m, site_vec_order, "all", TRUE)
# output_uncert_all_DD_larv_nav_1000m <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info_1000m, site_vec_order, "all", TRUE)
# 
# # And without DD
# best_est_metrics_larv_nav_200m <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info_200m, site_vec_order, "FALSE")
# best_est_metrics_larv_nav_500m <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info_500m, site_vec_order, "FALSE")
# best_est_metrics_larv_nav_1000m <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info_1000m, site_vec_order, "FALSE")
# output_uncert_all_larv_nav_200m <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info_200m, site_vec_order, "all", FALSE)
# output_uncert_all_larv_nav_500m <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info_500m, site_vec_order, "all", FALSE)
# output_uncert_all_larv_nav_1000m <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info_1000m, site_vec_order, "all", FALSE)

##### What if our sites retained all of the offspring they produced? Persistent then? Yes, because LRP is > 1, right?
# And makes sense that best est isn't NP>1 even at 100% habitat because still not retaining all of the offspring produced there...

##### Are fish getting evicted from the LEP matrix?
max_size_test_vec <- seq(from=breeding_size_mean, to=20, by=0.1 )
LEP_test_max_size <- data.frame(max_size = max_size_test_vec, LEP = NA)

for(i in 1:length(max_size_test_vec)) {
  LEP_test_max_size$LEP[i] = findLEP(param_best_est_mean_collected_offspring$min_size, max_size_test_vec[i],
                                     param_best_est_mean_collected_offspring$n_bins, param_best_est_mean_collected_offspring$t_steps,
                                     site_surv_best_est_avg_Sint$Sint[1], site_surv_best_est_avg_Sint$Sl[1],
                                     param_best_est_mean_collected_offspring$s, param_best_est_mean_collected_offspring$Linf,
                                     param_best_est_mean_collected_offspring$k_growth, param_best_est_mean_collected_offspring$eggs_per_clutch,
                                     param_best_est_mean_collected_offspring$clutches_per_year_mean,
                                     param_best_est_mean_collected_offspring$breeding_size, param_best_est_mean_collected_offspring$start_recruit_size,
                                     param_best_est_mean_collected_offspring$start_recruit_sd, param_best_est_mean_collected_offspring$egg_size_slope,
                                     param_best_est_mean_collected_offspring$egg_size_intercept, param_best_est_mean_collected_offspring$eyed_effect)
}
  
##### What is LRP if we include all the arriving offspring (not just those locally produced)?
# LRP_all_offspring = LEP*Re (right?) # 10.43057 (best_est_metrics_mean_offspring_all_offspring$LEP_R_mean)


# ##### What-if calculation 4) For the largest patch to be SP (highest p(i,i)), what would LRP need to be? And then what would LEP and egg-recruit survival need to be to acheive that LRP?
# highest_self_disp <- max(best_est)

##### What-if calculation 3) What would egg-recruit survival need to be for one of the patches to be SP? 

##### What-if calculation 4) What would local retention need to be for one of the patches to be SP?

##### What-if calculation 5) If we include the ghost population recruits too, is the population NP persistent?

#################### Metrics and parameters summarized for easy access: ####################
params_summary <- data.frame(param = c("k_disp","k_disp_lcl","k_disp_ucl",
                                        "theta_disp","theta_disp_lcl","theta_disp_ucl",
                                        "Linf","Linf_lcl","Linf_ucl",
                                        "k_growth","k_growth_lcl","k_growth_ucl",
                                        "n_parents_genotyped","n_offspring_genotyped",
                                        "n_offspring_matched","assignment_rate",
                                        "Ph","Pc","Ps","Pd",
                                        "p_APCL","p_UNOC"),
                              value = c(k_allyears, min(k_theta_allyear_95CI_values$k_eval), max(k_theta_allyear_95CI_values$k_eval),
                                        theta_allyears, min(k_theta_allyear_95CI_values$theta_eval), max(k_theta_allyear_95CI_values$theta_eval),
                                        Linf_growth_mean, min(growth_set_params$Linf_est), max(growth_set_params$Linf_est),
                                        k_growth_mean, min(growth_set_params$k_est), max(growth_set_params$k_est),
                                        n_parents_genotyped, n_offspring_genotyped, 
                                        n_offspring_matched, assignment_rate,
                                        Ph, prob_r_mean, Ps, Pd,
                                        perc_APCL_val, perc_UNOC_val))

metrics_summary <- data.frame(metric = c("LEP avg best est", "LEP avg lower", "LEP avg upper",
                                         "LRP_DD avg best est", "LRP_DD avg lower", "LRP_DD avg upper",
                                         "NP_DD best est", "NP_DD lower", "NP_DD upper",
                                         "egg_recruit surv best est", "egg recruit surv lower", "egg recruit surv upper",
                                         "LPR_local avg best est", "LRP_local avg lower", "LRP_local avg upper"),
                              value = c(LEP_best_est, LEP_avg_lower, LEP_avg_upper,
                                        LEP_R_best_est_DD, LRP_average_lower, LRP_average_upper,
                                        NP_best_est_DD, NP_lower, NP_upper,
                                        best_est_metrics_mean_offspring_DD$recruits_per_egg, egg_recruit_survival_lower, egg_recruit_survival_upper,
                                        LEP_R_local_best_est_DD, LRP_local_average_lower, LRP_local_average_upper))

                                               


# metrics_params_summary <- data.frame(metric_param = c("k_disp","theta_disp","k_disp_lcl","k_disp_ucl","theta_disp_lcl","theta_disp_ucl",
#                                                       "L_inf","L_inf_lcl","Linf_ucl","k_growth","k_growth_lcl","k_growth_ucl",
#                                                       "n_parents_genotyped","n_offspring_genotyped","n_offspring_matched", "assignment_rate",
#                                                       "Ph","Pc","Ps","Pd",
#                                                       "best_est_NP_DD","NP_DD_lcl","NP_DD_ucl",
#                                                       "best_est_LEP_avg","best_est_LEP_site_min","best_est_LEP_site_max",
#                                                       "best_est_LRP_avg","best_est_LRP_site_min","best_est_LRP_site_max",
#                                                       "best_est_local_replacement_avg","best_est_local_replacement_site_min","best_est_local_replacement_site_max",
#                                                       "best_est_recruits_per_egg_avg","rperE_avg_lcl","rperE_avg_ucl",
#                                                       "WI_all_offs_NP_DD","WI_all_offs_NP",
#                                                       "WI_LRP_for_NP","WI_LEP_for_NP","WI_LEP_for_NP_DD","WI_RperE_for_NP",
#                                                       "breeding_size_mean", "min_breeding_transition_size", "max_breeding_transition_size",
#                                                       "P_DD"),
#                                      type = c("param","param","param","param","param","param",
#                                               "param","param","param","param","param","param",
#                                               "param","param","param","param",
#                                               "param","param","param","param",
#                                               "metric","metric","metric",
#                                               "metric","metric","metric",
#                                               "metric","metric","metric",
#                                               "metric","metric","metric",
#                                               "metric","metric","metric",
#                                               "metric","metric",
#                                               "metric","metric","metric","metric",
#                                               "param","param","param",
#                                               "param"),
#                                      value = c(k_allyears, theta_allyears, min(k_connectivity_set), max(k_connectivity_set), min(theta_connectivity_set), max(theta_connectivity_set),
#                                                Linf_growth_mean, min(Linf_set), max(Linf_set), k_growth_mean, min(k_growth_set), max(k_growth_set),
#                                                n_parents_genotyped, n_offspring_genotyped, n_offspring_matched, assignment_rate,
#                                                Ph, prob_r_mean, Ps, Pd,
#                                                best_est_metrics_mean_offspring_DD$NP, min(output_uncert_all_DD$NP_out_df$value), max(output_uncert_all_DD$NP_out_df$value),
#                                                LEP_best_est, min(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP), max(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP),
#                                                best_est_metrics_mean_offspring_DD$LEP_R_mean, min(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP_R), max(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP_R),
#                                                best_est_metrics_mean_offspring_DD$LEP_R_local_mean, min(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP_R_local), max(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP_R_local),
#                                                best_est_metrics_mean_offspring_DD$recruits_per_egg, min(output_uncert_all_DD$RperE_out_df$value), max(output_uncert_all_DD$RperE_out_df$value),
#                                                best_est_metrics_mean_offspring_all_offspring_DD$NP, best_est_metrics_mean_offspring_all_offspring$NP,
#                                                LRP_for_NP,LEP_for_NP,LEP_for_NP_DD,egg_recruit_survival_for_NP,
#                                                breeding_size_mean, min(recap_first_female$size), max(recap_first_female$size),
#                                                (perc_APCL_val+perc_UNOC_val)/perc_UNOC_val))

#################### Save output: ####################
# Runs without DD compensation
save(param_best_est_mean_collected_offspring, file=here::here("Data/Script_outputs", "param_best_est_mean_collected_offspring.RData"))
save(best_est_metrics_mean_offspring, file=here::here("Data/Script_outputs", "best_est_metrics_mean_offspring.RData"))
save(param_set_full, file=here::here("Data/Script_outputs", "param_set_full.RData"))
save(output_uncert_start_recruit, file=here::here("Data/Script_outputs", "output_uncert_start_recruit.RData"))
save(output_uncert_growth, file=here::here("Data/Script_outputs", "output_uncert_growth.RData"))
save(output_uncert_survival, file=here::here("Data/Script_outputs", "output_uncert_survival.RData"))
save(output_uncert_breeding_size, file=here::here("Data/Script_outputs", "output_uncert_breeding_size.RData"))
save(output_uncert_offspring_assigned, file=here::here("Data/Script_outputs", "output_uncert_offspring_assigned.RData"))
save(output_uncert_prob_r, file=here::here("Data/Script_outputs", "output_uncert_prob_r.RData"))
#save(output_uncert_prob_r_and_offspring_assigned, file=here::here("Data/Script_outputs", "output_uncert_prob_r_and_assigned_offspring.RData"))
save(output_uncert_dispersal, file=here::here("Data/Script_outputs", "output_uncert_dispersal.RData"))
save(output_uncert_all, file=here::here("Data/Script_outputs", "output_uncert_all.RData"))

# Runs with DD compensation
save(best_est_metrics_mean_offspring_DD, file=here::here("Data/Script_outputs", "best_est_metrics_mean_offspring_DD.RData"))
save(output_uncert_start_recruit_DD, file=here::here("Data/Script_outputs", "output_uncert_start_recruit_DD.RData"))
save(output_uncert_growth_DD, file=here::here("Data/Script_outputs", "output_uncert_growth_DD.RData"))
save(output_uncert_survival_DD, file=here::here("Data/Script_outputs", "output_uncert_survival_DD.RData"))
save(output_uncert_breeding_size_DD, file=here::here("Data/Script_outputs", "output_uncert_breeding_size_DD.RData"))
save(output_uncert_offspring_assigned_DD, file=here::here("Data/Script_outputs", "output_uncert_offspring_assigned_DD.RData"))
save(output_uncert_prob_r_DD, file=here::here("Data/Script_outputs", "output_uncert_prob_r_DD.RData"))
#save(output_uncert_prob_r_and_offspring_assigned_DD, file=here::here("Data/Script_outputs", "output_uncert_prob_r_and_assigned_offspring_DD.RData"))
save(output_uncert_dispersal_DD, file=here::here("Data/Script_outputs", "output_uncert_dispersal_DD.RData"))
save(output_uncert_all_DD, file=here::here("Data/Script_outputs", "output_uncert_all_DD.RData"))

# Params
save(site_surv_param_sets, file=here::here("Data/Script_outputs","site_surv_param_sets.RData"))

# What-ifs
save(best_est_metrics_mean_offspring_all_offspring, file=here::here("Data/Script_outputs","best_est_metrics_mean_offspring_all_offspring.RData"))  # without DD compensation
save(best_est_metrics_mean_offspring_all_offspring_DD, file=here::here("Data/Script_outputs","best_est_metrics_mean_offspring_all_offspring_DD.RData"))  # with DD compensation
#save(perc_hab_best_ests_realSurvs, file=here::here("Data/Script_outputs","perc_hab_best_ests_realSurvs.RData"))
#save(perc_hab_uncertainty_realSurvs, file=here::here("Data/Script_outputs","perc_hab_uncertainty_realSurvs.RData"))
save(perc_hab_best_ests_avgSurvs, file=here::here("Data/Script_outputs","perc_hab_best_ests_avgSurvs.RData"))
save(perc_hab_uncertainty_avgSurvs, file=here::here("Data/Script_outputs","perc_hab_uncertainty_avgSurvs.RData"))
#save(NP_by_perc_hab_realSurvs, file=here::here("Data/Script_outputs","NP_by_perc_hab_realSurvs.RData"))
#save(NP_by_perc_hab_avgSurvs, file=here::here("Data/Script_outputs","NP_by_perc_hab_avgSurvs.RData"))
#save(NP_above1_perc_hab_realSurvs, file=here::here("Data/Script_outputs","NP_above1_perc_hab_realSurvs.RData"))
#save(NP_above1_perc_hab_avgSurvs, file=here::here("Data/Script_outputs","NP_above1_perc_hab_avgSurvs.RData"))

#save(NP_by_wider_region_avgSurvs, file=here::here("Data/Script_outputs","NP_by_wider_region_avgSurvs.RData"))
#save(NP_above1_wider_region_avgSurvs, file=here::here("Data/Script_outputs","NP_above1_wider_region_avgSurvs.RData"))
#save(NP_by_wider_region_all_hab_avgSurvs, file=here::here("Data/Script_outputs","NP_by_wider_region_all_hab_avgSurvs.RData"))
#save(NP_above1_wider_region_all_hab_avgSurvs, file=here::here("Data/Script_outputs","NP_above1_wider_region_all_hab_avgSurvs.RData"))

save(wider_region_best_ests_avgSurvs, file=here::here("Data/Script_outputs","wider_region_best_ests_avgSurvs.RData"))
save(wider_region_uncertainty_avgSurvs, file=here::here("Data/Script_outputs","wider_region_uncertainty_avgSurvs.RData"))
save(wider_region_all_hab_best_ests_avgSurvs, file=here::here("Data/Script_outputs","wider_region_all_hab_best_ests_avgSurvs.RData"))
save(wider_region_all_hab_uncertainty_avgSurvs, file=here::here("Data/Script_outputs","wider_region_all_hab_uncertainty_avgSurvs.RData"))

# Larval navigation sensitivity
save(larv_nav_prop_hab, file=here::here("Data/Script_outputs","larv_nav_prop_hab.RData"))
save(larv_nav_params_list, file=here::here("Data/Script_outputs","larv_nav_params_list.RData"))
save(larv_nav_best_ests, file=here::here("Data/Script_outputs","larv_nav_best_ests.RData"))
save(larv_nav_uncert_all, file=here::here("Data/Script_outputs","larv_nav_uncert_all.RData"))

# # From older version...
# save(best_est_metrics_larv_nav_200m, file=here::here("Data/Script_outputs","best_est_metrics_larv_nav_200m.RData"))  # 200m each side of patch, no DD comp
# save(best_est_metrics_larv_nav_200m_DD, file=here::here("Data/Script_outputs","best_est_metrics_larv_nav_200m_DD.RData"))  # 200m each side of patch, DD comp
# save(best_est_metrics_larv_nav_500m, file=here::here("Data/Script_outputs","best_est_metrics_larv_nav_500m.RData"))  # 500m each side of patch, no DD comp
# save(best_est_metrics_larv_nav_500m_DD, file=here::here("Data/Script_outputs","best_est_metrics_larv_nav_500m_DD.RData"))  # 500m each side of patch, DD comp
# save(best_est_metrics_larv_nav_1000m, file=here::here("Data/Script_outputs","best_est_metrics_larv_nav_1000m.RData"))  # 1000m each side of patch, no DD comp
# save(best_est_metrics_larv_nav_1000m_DD, file=here::here("Data/Script_outputs","best_est_metrics_larv_nav_1000m_DD.RData"))  # 1000m each side of patch, DD comp
# save(output_uncert_all_larv_nav_200m, file=here::here("Data/Script_outputs","output_uncert_all_larv_nav_200m.RData"))  # 200m each side of patch, no DD comp
# save(output_uncert_all_DD_larv_nav_200m, file=here::here("Data/Script_outputs","output_uncert_all_DD_larv_nav_200m.RData"))  # 200m each side of patch, DD comp
# save(output_uncert_all_larv_nav_500m, file=here::here("Data/Script_outputs","output_uncert_all_larv_nav_500m.RData"))  # 500m each side of patch, no DD comp
# save(output_uncert_all_DD_larv_nav_500m, file=here::here("Data/Script_outputs","output_uncert_all_DD_larv_nav_500m.RData"))  # 500m each side of patch, DD comp
# save(output_uncert_all_larv_nav_1000m, file=here::here("Data/Script_outputs","output_uncert_all_larv_nav_1000m.RData"))  # 1000m each side of patch, no DD comp
# save(output_uncert_all_DD_larv_nav_1000m, file=here::here("Data/Script_outputs","output_uncert_all_DD_larv_nav_1000m.RData"))  # 1000m each side of patch, DD comp

# Summary of parameters and metrics, for easy reference while writing
#save(metrics_params_summary, file=here::here("Data/Script_outputs","metrics_params_summary.RData"))
save(params_summary, file=here::here("Data/Script_outputs","params_summary.RData"))
save(metrics_summary, file=here::here("Data/Script_outputs","metrics_summary.RData"))

# ###################################### STUFF TO RUN OVERNIGHT!
# ### Run metrics for a bunch of different types of uncertainty - No DD
# output_uncert_start_recruit <- calcMetricsAcrossRuns(n_runs, param_set_start_recruit, site_surv_best_est_sets, site_dist_info, site_vec_order, "start recruit size", FALSE)
# output_uncert_growth <- calcMetricsAcrossRuns(n_runs, param_set_growth, site_surv_best_est_sets, site_dist_info, site_vec_order, "growth", FALSE)
# output_uncert_survival <- calcMetricsAcrossRuns(n_runs, param_set_survival, site_surv_param_sets, site_dist_info, site_vec_order, "survival", FALSE)
# output_uncert_breeding_size <- calcMetricsAcrossRuns(n_runs, param_set_breeding_size, site_surv_best_est_sets, site_dist_info, site_vec_order, "breeding size", FALSE)
# output_uncert_offspring_assigned <- calcMetricsAcrossRuns(n_runs, param_set_offspring_assigned, site_surv_best_est_sets, site_dist_info, site_vec_order, "assigned offspring", FALSE)
# output_uncert_prob_r <- calcMetricsAcrossRuns(n_runs, param_set_prob_r, site_surv_best_est_sets, site_dist_info, site_vec_order, "prob r", FALSE)
# #output_uncert_prob_r_and_offspring_assigned <- calcMetricsAcrossRuns(n_runs, param_set_prob_r_offspring_assigned, site_surv_best_est_sets, site_dist_info, site_vec_order, "assigned offspring and prob r", FALSE)
# output_uncert_dispersal <- calcMetricsAcrossRuns(n_runs, param_set_dispersal, site_surv_best_est_sets, site_dist_info, site_vec_order, "dispersal k", FALSE)
# output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info, site_vec_order, "all", FALSE)
# 
# # output_uncert_growth <- calcMetricsAcrossRuns(n_runs, param_set_growth, site_surv_best_est_sets, site_dist_info, site_vec_order, "growth", FALSE)
# # output_uncert_growth_DD <- calcMetricsAcrossRuns(n_runs, param_set_growth, site_surv_best_est_sets, site_dist_info, site_vec_order, "growth", TRUE)
# # output_uncert_all_DD <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info, site_vec_order, "all", TRUE)
# 
# 
# # Join them together for plotting purposes
# LEP_uncert <- rbind(output_uncert_start_recruit$LEP_out_df, output_uncert_growth$LEP_out_df,
#                     output_uncert_survival$LEP_out_df, output_uncert_breeding_size$LEP_out_df,
#                     output_uncert_offspring_assigned$LEP_out_df,output_uncert_prob_r$LEP_out_df,
#                     #output_uncert_prob_r_and_offspring_assigned$LEP_out_df, 
#                     output_uncert_dispersal$LEP_out_df, output_uncert_all$LEP_out_df)
# 
# LEP_by_site_uncert <- rbind(output_uncert_start_recruit$LEP_by_site_out_df, output_uncert_growth$LEP_by_site_out_df,
#                             output_uncert_survival$LEP_by_site_out_df, output_uncert_breeding_size$LEP_by_site_out_df,
#                             output_uncert_offspring_assigned$LEP_by_site_out_df, output_uncert_prob_r$LEP_by_site_out_df,
#                             output_uncert_dispersal$LEP_by_site_out_df, output_uncert_all$LEP_by_site_out_df)
# 
# LEP_R_uncert <- rbind(output_uncert_start_recruit$LEP_R_out_df, output_uncert_growth$LEP_R_out_df,
#                       output_uncert_survival$LEP_R_out_df, output_uncert_breeding_size$LEP_R_out_df,
#                       output_uncert_offspring_assigned$LEP_R_out_df, output_uncert_prob_r$LEP_R_out_df,
#                       #output_uncert_prob_r_and_offspring_assigned$LEP_R,
#                       output_uncert_dispersal$LEP_R_out_df, output_uncert_all$LEP_R_out_df)
# 
# LEP_R_by_site_uncert <- rbind(output_uncert_start_recruit$LEP_R_by_site_out_df, output_uncert_growth$LEP_R_by_site_out_df,
#                               output_uncert_survival$LEP_R_by_site_out_df, output_uncert_breeding_size$LEP_R_by_site_out_df,
#                               output_uncert_offspring_assigned$LEP_R_by_site_out_df, output_uncert_prob_r$LEP_R_by_site_out_df,
#                               output_uncert_dispersal$LEP_R_by_site_out_df, output_uncert_all$LEP_R_by_site_out_df)
# 
# RperE_uncert <- rbind(output_uncert_start_recruit$RperE_out_df, output_uncert_growth$RperE_out_df,
#                       output_uncert_survival$RperE_out_df, output_uncert_breeding_size$RperE_out_df,
#                       output_uncert_offspring_assigned$RperE_out_df, output_uncert_prob_r$RperE_out_df,
#                       #output_uncert_prob_r_and_offspring_assigned$RperE_out_df,
#                       output_uncert_dispersal$RperE_out_df, output_uncert_all$RperE_out_df)
# 
# NP_uncert <- rbind(output_uncert_start_recruit$NP_out_df, output_uncert_growth$NP_out_df,
#                    output_uncert_survival$NP_out_df, output_uncert_breeding_size$NP_out_df,
#                    output_uncert_offspring_assigned$NP_out_df, output_uncert_prob_r$NP_out_df,
#                    #output_uncert_prob_r_and_offspring_assigned$NP_out_df,
#                    output_uncert_dispersal$NP_out_df, output_uncert_all$NP_out_df)
# 
# ###### Now run with DD
# ### Run metrics for a bunch of different types of uncertainty
# output_uncert_start_recruit_DD <- calcMetricsAcrossRuns(n_runs, param_set_start_recruit, site_surv_best_est_sets, site_dist_info, site_vec_order, "start recruit size", TRUE)
# output_uncert_growth_DD <- calcMetricsAcrossRuns(n_runs, param_set_growth, site_surv_best_est_sets, site_dist_info, site_vec_order, "growth", TRUE)
# output_uncert_survival_DD <- calcMetricsAcrossRuns(n_runs, param_set_survival, site_surv_param_sets, site_dist_info, site_vec_order, "survival", TRUE)
# output_uncert_breeding_size_DD <- calcMetricsAcrossRuns(n_runs, param_set_breeding_size, site_surv_best_est_sets, site_dist_info, site_vec_order, "breeding size", TRUE)
# output_uncert_offspring_assigned_DD <- calcMetricsAcrossRuns(n_runs, param_set_offspring_assigned, site_surv_best_est_sets, site_dist_info, site_vec_order, "assigned offspring", TRUE)
# output_uncert_prob_r_DD <- calcMetricsAcrossRuns(n_runs, param_set_prob_r, site_surv_best_est_sets, site_dist_info, site_vec_order, "prob r", TRUE)
# #output_uncert_prob_r_and_offspring_assigned_DD <- calcMetricsAcrossRuns(n_runs, param_set_prob_r_offspring_assigned, site_surv_best_est_sets, site_dist_info, site_vec_order, "assigned offspring and prob r", TRUE)
# output_uncert_dispersal_DD <- calcMetricsAcrossRuns(n_runs, param_set_dispersal, site_surv_best_est_sets, site_dist_info, site_vec_order, "dispersal k", TRUE)
# output_uncert_all_DD <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info, site_vec_order, "all", TRUE)
# 
# # Join them together for plotting purposes
# LEP_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_out_df, output_uncert_growth_DD$LEP_out_df,
#                        output_uncert_survival_DD$LEP_out_df, output_uncert_breeding_size_DD$LEP_out_df,
#                        output_uncert_offspring_assigned_DD$LEP_out_df, output_uncert_prob_r_DD$LEP_out_df,
#                        #output_uncert_prob_r_and_offspring_assigned_DD$LEP_out_df, 
#                        output_uncert_dispersal_DD$LEP_out_df, output_uncert_all_DD$LEP_out_df)
# 
# LEP_by_site_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_by_site_out_df, output_uncert_growth_DD$LEP_by_site_out_df,
#                                output_uncert_survival_DD$LEP_by_site_out_df, output_uncert_breeding_size_DD$LEP_by_site_out_df,
#                                output_uncert_offspring_assigned_DD$LEP_by_site_out_df, output_uncert_prob_r_DD$LEP_by_site_out_df,
#                                output_uncert_dispersal_DD$LEP_by_site_out_df, output_uncert_all_DD$LEP_by_site_out_df)
# 
# LEP_R_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_R_out_df, output_uncert_growth_DD$LEP_R_out_df,
#                          output_uncert_survival_DD$LEP_R_out_df, output_uncert_breeding_size_DD$LEP_R_out_df,
#                          output_uncert_offspring_assigned_DD$LEP_R_out_df, output_uncert_prob_r_DD$LEP_R_out_df,
#                          #output_uncert_prob_r_and_offspring_assigned_DD$LEP_R,
#                          output_uncert_dispersal_DD$LEP_R_out_df, output_uncert_all_DD$LEP_R_out_df)
# 
# LEP_R_by_site_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_R_by_site_out_df, output_uncert_growth_DD$LEP_R_by_site_out_df,
#                                  output_uncert_survival_DD$LEP_R_by_site_out_df, output_uncert_breeding_size_DD$LEP_R_by_site_out_df,
#                                  output_uncert_offspring_assigned_DD$LEP_R_by_site_out_df, output_uncert_prob_r_DD$LEP_R_by_site_out_df,
#                                  output_uncert_dispersal_DD$LEP_R_by_site_out_df, output_uncert_all_DD$LEP_R_by_site_out_df)
# 
# RperE_uncert_DD <- rbind(output_uncert_start_recruit_DD$RperE_out_df, output_uncert_growth_DD$RperE_out_df,
#                          output_uncert_survival_DD$RperE_out_df, output_uncert_breeding_size_DD$RperE_out_df,
#                          output_uncert_offspring_assigned_DD$RperE_out_df, output_uncert_prob_r_DD$RperE_out_df,
#                          #output_uncert_prob_r_and_offspring_assigned_DD$RperE_out_df,
#                          output_uncert_dispersal_DD$RperE_out_df, output_uncert_all_DD$RperE_out_df)
# 
# NP_uncert_DD <- rbind(output_uncert_start_recruit_DD$NP_out_df, output_uncert_growth_DD$NP_out_df,
#                       output_uncert_survival_DD$NP_out_df, output_uncert_breeding_size_DD$NP_out_df,
#                       output_uncert_offspring_assigned_DD$NP_out_df, output_uncert_prob_r_DD$NP_out_df,
#                       #output_uncert_prob_r_and_offspring_assigned_DD$NP_out_df,
#                       output_uncert_dispersal_DD$NP_out_df, output_uncert_all_DD$NP_out_df)
# 
# output_uncert_all_offspring_all <- calcMetricsAcrossRuns(n_runs, param_set_full_all_offspring, site_surv_param_sets, site_dist_info, site_vec_order, "all: alloff", FALSE)
# output_uncert_all_offspring_all_DD <- calcMetricsAcrossRuns(n_runs, param_set_full_all_offspring, site_surv_param_sets, site_dist_info, site_vec_order, "all: alloff", TRUE)
# 
# perc_hab_best_ests_avgSurvs <- list()
# perc_hab_uncertainty_avgSurvs <- list()
# for(i in 1:length(perc_hab_vals)) {
#   perc_hab_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_avg_Sint, site_dist_info_habperc[[i]], site_vec_order, TRUE)
#   perc_hab_uncertainty_avgSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_avg_Sint_param_sets, site_dist_info_habperc[[i]], site_vec_order, "all: hab sens, avg surv", TRUE)
# }
# 
# wider_region_best_ests_avgSurvs <- list()
# wider_region_uncertainty_avgSurvs <- list()
# for(i in 1:length(region_width_list)) {
#   wider_region_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_avg_Sint, site_dist_info_wider_region[[i]], site_vec_order, TRUE)
#   wider_region_uncertainty_avgSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_avg_Sint_param_sets, site_dist_info_wider_region[[i]], site_vec_order, "all: region sens, avg surv", TRUE)
# }
# 
# wider_region_all_hab_best_ests_avgSurvs <- list()
# wider_region_all_hab_uncertainty_avgSurvs <- list()
# for(i in 1:length(region_width_list)) {
#   wider_region_all_hab_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_avg_Sint, site_dist_info_wider_region_all_hab[[i]], site_vec_order, TRUE)
#   wider_region_all_hab_uncertainty_avgSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_avg_Sint_param_sets, site_dist_info_wider_region_all_hab[[i]], site_vec_order, "all: region sens, all hab", TRUE)
# }
# 
# wider_region_perc_hab_best_ests_avgSurvs <- list()
# #wider_region_perc_hab_uncertainty_avgSurvs <- list()
# for(i in 1:(length(region_width_list_long_form)*length(perc_hab_vals_long_form))) {
#   wider_region_perc_hab_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_avg_Sint, site_dist_info_wider_region_perc_hab[[i]], site_vec_order, TRUE)
#   #  wider_region_perc_hab_uncertainty_avgSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_avg_Sint_param_sets, site_dist_info_wider_region_perc_hab[[i]], site_vec_order, "all: region sens, perc hab sens, avg surv", TRUE)
#   print(i)
# }
# # Runs without DD compensation
# save(param_best_est_mean_collected_offspring, file=here::here("Data/Script_outputs", "param_best_est_mean_collected_offspring.RData"))
# save(best_est_metrics_mean_offspring, file=here::here("Data/Script_outputs", "best_est_metrics_mean_offspring.RData"))
# save(param_set_full, file=here::here("Data/Script_outputs", "param_set_full.RData"))
# save(output_uncert_start_recruit, file=here::here("Data/Script_outputs", "output_uncert_start_recruit.RData"))
# save(output_uncert_growth, file=here::here("Data/Script_outputs", "output_uncert_growth.RData"))
# save(output_uncert_survival, file=here::here("Data/Script_outputs", "output_uncert_survival.RData"))
# save(output_uncert_breeding_size, file=here::here("Data/Script_outputs", "output_uncert_breeding_size.RData"))
# save(output_uncert_offspring_assigned, file=here::here("Data/Script_outputs", "output_uncert_offspring_assigned.RData"))
# save(output_uncert_prob_r, file=here::here("Data/Script_outputs", "output_uncert_prob_r.RData"))
# #save(output_uncert_prob_r_and_offspring_assigned, file=here::here("Data/Script_outputs", "output_uncert_prob_r_and_assigned_offspring.RData"))
# save(output_uncert_dispersal, file=here::here("Data/Script_outputs", "output_uncert_dispersal.RData"))
# save(output_uncert_all, file=here::here("Data/Script_outputs", "output_uncert_all.RData"))
# 
# # Runs with DD compensation
# save(best_est_metrics_mean_offspring_DD, file=here::here("Data/Script_outputs", "best_est_metrics_mean_offspring_DD.RData"))
# save(output_uncert_start_recruit_DD, file=here::here("Data/Script_outputs", "output_uncert_start_recruit_DD.RData"))
# save(output_uncert_growth_DD, file=here::here("Data/Script_outputs", "output_uncert_growth_DD.RData"))
# save(output_uncert_survival_DD, file=here::here("Data/Script_outputs", "output_uncert_survival_DD.RData"))
# save(output_uncert_breeding_size_DD, file=here::here("Data/Script_outputs", "output_uncert_breeding_size_DD.RData"))
# save(output_uncert_offspring_assigned_DD, file=here::here("Data/Script_outputs", "output_uncert_offspring_assigned_DD.RData"))
# save(output_uncert_prob_r_DD, file=here::here("Data/Script_outputs", "output_uncert_prob_r_DD.RData"))
# #save(output_uncert_prob_r_and_offspring_assigned_DD, file=here::here("Data/Script_outputs", "output_uncert_prob_r_and_assigned_offspring_DD.RData"))
# save(output_uncert_dispersal_DD, file=here::here("Data/Script_outputs", "output_uncert_dispersal_DD.RData"))
# save(output_uncert_all_DD, file=here::here("Data/Script_outputs", "output_uncert_all_DD.RData"))
# 
# # Params
# save(site_surv_param_sets, file=here::here("Data/Script_outputs","site_surv_param_sets.RData"))
# 
# # What-ifs
# save(best_est_metrics_mean_offspring_all_offspring, file=here::here("Data/Script_outputs","best_est_metrics_mean_offspring_all_offspring.RData"))  # without DD compensation
# save(best_est_metrics_mean_offspring_all_offspring_DD, file=here::here("Data/Script_outputs","best_est_metrics_mean_offspring_all_offspring_DD.RData"))  # with DD compensation
# save(perc_hab_best_ests_avgSurvs, file=here::here("Data/Script_outputs","perc_hab_best_ests_avgSurvs.RData"))
# save(perc_hab_uncertainty_avgSurvs, file=here::here("Data/Script_outputs","perc_hab_uncertainty_avgSurvs.RData"))
# 
# save(wider_region_best_ests_avgSurvs, file=here::here("Data/Script_outputs","wider_region_best_ests_avgSurvs.RData"))
# save(wider_region_uncertainty_avgSurvs, file=here::here("Data/Script_outputs","wider_region_uncertainty_avgSurvs.RData"))
# save(wider_region_all_hab_best_ests_avgSurvs, file=here::here("Data/Script_outputs","wider_region_all_hab_best_ests_avgSurvs.RData"))
# save(wider_region_all_hab_uncertainty_avgSurvs, file=here::here("Data/Script_outputs","wider_region_all_hab_uncertainty_avgSurvs.RData"))
# 
# save(wider_region_perc_hab_best_ests_avgSurvs, file=here::here("Data/Script_outputs","wider_region_perc_hab_best_ests_avgsSurvs.RData"))


#################### Plots: ####################

########## Main text figures ##########

##### Figure 1 (schematic - made outside of R)
##### Figure 2 (map + photo - made in SiteMap.R script)

##### Figure 3 (demographic and dispersal inputs: survival curve, growth curve, dispersal kernel, transition size to female)
# Prep work for dispersal kernel figure
distance_vec <- seq(from=0, to=50, by=0.01)
connectivity_est_vec <- disp_kernel_all_years(distance_vec, k_allyears, theta_allyears) 
connectivity_sens <- matrix(ncol = length(distance_vec), nrow = n_runs)
for(i in 1:n_runs) {
  connectivity_sens[i,] = disp_kernel_all_years(distance_vec, k_connectivity_set[i], theta_connectivity_set[i])
}
connectivity_lowest_bound = apply(connectivity_sens, 2, min)
connectivity_highest_bound = apply(connectivity_sens, 2, max)
dispersal_df <- data.frame(distance = distance_vec, kernel_bestfit = connectivity_est_vec, kernel_CI1 = connectivity_lowest_bound, kernel_CI2 = connectivity_highest_bound)

# Dispersal kernel plot
dispersal_kernel_plot <- ggplot(data=dispersal_df, aes(x=distance, y=kernel_bestfit, ymin=kernel_CI1, ymax=kernel_CI2)) +
  geom_line(color='black') +
  geom_ribbon(alpha=0.5, color='gray') +
  xlab('distance (km)') + ylab('dispersal probability') + #ggtitle('Dispersal kernel') +
  theme_bw() 

# Growth curve - plot the data and linear model for fish caught about a year apart, models for only one recapture pair per fish, pairs selected randomly and models fit 1000x
growth_curve_plot <- ggplot(data = recap_pairs_year, aes(x = L1, y = L2)) +
  geom_point(shape = 1, color = "dark gray") +
  #geom_abline(aes(intercept = mean(growth_info_estimate$intercept_est), slope = mean(growth_info_estimate$slope_est)), color = "black", lwd = 1.5) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", lwd = 1.5, linetype = "dashed") +  #  1:1 line
  geom_ribbon(aes(x=seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1)),
                  ymin = (min(growth_info_estimate$intercept_est) + min(growth_info_estimate$slope_est)*seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1))),
                  ymax = (max(growth_info_estimate$intercept_est) + max(growth_info_estimate$slope_est)*seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1)))), fill = "dark gray", alpha = 0.75) +
  geom_abline(aes(intercept = mean(growth_info_estimate$intercept_est), slope = mean(growth_info_estimate$slope_est)), color = "black", lwd = 1.5) +
  xlab("length (cm)") + ylab("length (cm) next year") + #ggtitle('Growth') +
  theme_bw()

# Alternate version if want to do a VBL instead
# growth_df <- data.frame(length1 = seq(min_size, max_size, length.out = n_bins*10)) %>%
#   mutate(length2 = VBL_growth(Linf_growth_mean, k_growth_mean, length1))
# 
# growth_curve_plot <- ggplot(data=growth_df, aes(x=length1, y=length2)) +
#   geom_point(color='black') +
#   xlab('size (cm)') + ylab('size next year') + ggtitle('b) Growth curve') +
#   ylim(0,13) +
#   theme_bw()

# Breeding size distribution
breeding_size_plot <- ggplot(data = recap_first_female, aes(x=size)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept=breeding_size_mean, color='black') +
  xlab('transition size (cm)') + #ggtitle('Female transition size') +
  theme_bw()

# Survival plot (for now doing an example site - Palanas, but could change that...) - now doing Elementary School, which has median survival
survival_output_to_plot <- best_fit_model_dfs$surv_site_size %>% filter(site == "Elementary School")

survival_plot <- ggplot(data = survival_output_to_plot, aes(size, estimate_prob)) +
  geom_ribbon(aes(ymin=lcl_prob,ymax=ucl_prob),color="gray",fill="gray") +
  geom_line(color="black") +
  xlab("size (cm)") + ylab("probability of survival \n (ex. site: Elementary School)") + #ggtitle("Annual survival (ex. Palanas)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

pdf(file = here('Plots/FigureDrafts','Parameter_inputs.pdf'), width=6, height=6)
plot_grid(dispersal_kernel_plot, growth_curve_plot, survival_plot, breeding_size_plot, labels = c("a","b","c","d"), nrow=2)
dev.off()

##### Figure 4 (abundance + replacement metrics)
# LEP - why is this so nuts? (because CP,ST,SL and CC have wide uncertainty around survival, gives wide range of LEP) - here shown without them
LEP_plot <- ggplot(data = output_uncert_all$LEP_by_site_out_df %>% filter(site %in% c(1,4,5,6,7,8,9,10,11,12,13,14,17,18,19)), aes(x=value)) +
  geom_histogram(bins=100, color = 'gray', fill = 'gray') +
  #geom_freqpoly() +
  #geom_vline(data = LEP_best_est, aes(xintercept = LEP, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_best_est %>% filter(recruit_size == "mean offspring"))$LEP, color='black') +
  geom_vline(xintercept = LEP_best_est, color = "black") +
  xlab("lifetime egg production (LEP)") + #ggtitle('LEP') +
  theme_bw()

LEP_plot_freq <- ggplot(data = output_uncert_all$LEP_by_site_out_df %>% filter(site %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)), aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=100, color = 'gray', fill = 'gray') +
  #geom_freqpoly() +
  #geom_vline(data = LEP_best_est, aes(xintercept = LEP, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_best_est %>% filter(recruit_size == "mean offspring"))$LEP, color='black') +
  geom_vline(xintercept = LEP_best_est, color = "black") +
  xlab("lifetime egg production (LEP)") + ylab("relative frequency") + #ggtitle('LEP') +
  theme_bw()

# LEP_plot <- ggplot(data = output_uncert_all$LEP_by_site_out_df %>% filter(site %in% c(1,2,4,5,6,7,8,9,10,11,12,13,14,17,18,19)), aes(x=value)) +
#   geom_histogram(bins=500, color = 'gray', fill = 'gray') +
#   #geom_vline(data = LEP_best_est, aes(xintercept = LEP, color = recruit_size)) +
#   #geom_vline(xintercept = (LEP_best_est %>% filter(recruit_size == "mean offspring"))$LEP, color='black') +
#   geom_vline(xintercept = LEP_best_est, color = "black") +
#   xlab('LEP') + ggtitle('LEP') +
#   theme_bw()
# 
# LEP_plot <- ggplot(data = output_uncert_all$LEP_out_df, aes(x=value)) +
#   geom_histogram(bins=500, color = 'gray', fill = 'gray') +
#   #geom_vline(data = LEP_best_est, aes(xintercept = LEP, color = recruit_size)) +
#   #geom_vline(xintercept = (LEP_best_est %>% filter(recruit_size == "mean offspring"))$LEP, color='black') +
#   geom_vline(xintercept = LEP_best_est, color = "black") +
#   xlab('LEP') + ggtitle('LEP') +
#   theme_bw()

# LRP with DD - this is by site, with CC,CP,SL,ST removed
# LEP_R_plot_DD <- ggplot(data = output_uncert_all_DD$LEP_R_by_site_out_df %>% filter(site %in% c(1,4,5,6,7,8,9,10,11,12,13,14,17,18,19)), aes(x=value)) +
#   geom_histogram(bins=100, color = 'gray', fill = 'gray') +
#   #geom_vline(data = LEP_R_best_est, aes(xintercept = LEP_R, color = recruit_size)) +
#   #geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "mean offspring"))$LEP_R, color = 'black') +
#   geom_vline(xintercept = LEP_R_best_est_DD, color = "black") +
#   xlab('LRP') + #ggtitle('LRP') +
#   theme_bw()

# LRP with DD - this is average across sites
LEP_R_plot_DD <- ggplot(data = output_uncert_all_DD$LEP_R_out_df, aes(x=value)) +
  geom_histogram(bins=50, color = 'gray', fill = 'gray') +
  #geom_vline(data = LEP_R_best_est, aes(xintercept = LEP_R, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "mean offspring"))$LEP_R, color = 'black') +
  geom_vline(xintercept = LEP_R_best_est_DD, color = "black") +
  xlab(bquote("lifetime recruit production (LRP)")) + #ggtitle('LRP') +
  theme_bw()

LEP_R_plot_DD_freq <- ggplot(data = output_uncert_all_DD$LEP_R_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=50, color = 'gray', fill = 'gray') +
  #geom_vline(data = LEP_R_best_est, aes(xintercept = LEP_R, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "mean offspring"))$LEP_R, color = 'black') +
  geom_vline(xintercept = LEP_R_best_est_DD, color = "black") +
  xlab(bquote("lifetime recruit production (LRP)")) + ylab("relative frequency") + #ggtitle('LRP') +
  theme_bw()

# LRP_local with DD
LEP_R_local_plot_DD <- ggplot(data = output_uncert_all_DD$LEP_R_local_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = "gray", fill = "gray") +
  geom_vline(xintercept = LEP_R_local_best_est_DD, color = "black") +
  xlab(bquote("local replacement (LR)")) + #ggtitle("Local replacement") +
  theme_bw()

LEP_R_local_plot_DD_freq <- ggplot(data = output_uncert_all_DD$LEP_R_local_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color = "gray", fill = "gray") +
  geom_vline(xintercept = LEP_R_local_best_est_DD, color = "black") +
  xlab(bquote("local replacement (LR)")) + ylab("relative frequency") + #ggtitle("Local replacement") +
  theme_bw()

# Abundance trend 
Fig4_abundance_plot <- ggplot(data = site_trends_time, aes(x=year, y=mean_nF, group=site)) +
  geom_line(color="grey") +
  geom_line(data=site_trends_all, aes(x=year, y=mean_nF), color = "black", size=1.5) +
  xlab("year") + ylab("# females") + #ggtitle("Estimated abundance through time") +
  scale_x_continuous(breaks=c(2,4,6), labels=c("2013","2015","2017")) +
  theme_bw()

# Put them together
pdf(file = here::here('Plots/FigureDrafts', 'Abundance_LEP_LRP_LocalReplacement.pdf'), width=6, height=6)
plot_grid(Fig4_abundance_plot,LEP_plot, LEP_R_plot_DD, LEP_R_local_plot_DD,
          labels = c("a","b","c","d"), nrow=2)
dev.off()

# Put them together, using frequency plots instead of histograms
pdf(file = here::here('Plots/FigureDrafts', 'Abundance_LEP_LRP_LocalReplacement_FreqPlots.pdf'), width=6, height=6)
plot_grid(Fig4_abundance_plot, LEP_plot_freq, LEP_R_plot_DD_freq, LEP_R_local_plot_DD_freq,
          labels = c("a","b","c","d"), nrow=2)
dev.off()

##### Figure 5 (SP metrics, NP metrics, connectivity matrices)
best_est_metrics_mean_offspring_DD$Cmat$org_site <- replace(best_est_metrics_mean_offspring_DD$Cmat$org_site, 
                                                            best_est_metrics_mean_offspring_DD$Cmat$org_site=="Tamakin Dacot", 
                                                            "Tomakin Dako")
best_est_metrics_mean_offspring_DD$Cmat$dest_site <- replace(best_est_metrics_mean_offspring_DD$Cmat$dest_site, 
                                                            best_est_metrics_mean_offspring_DD$Cmat$dest_site=="Tamakin Dacot", 
                                                            "Tomakin Dako")
output_uncert_all_DD$SP_vals_with_params$site <- replace(output_uncert_all_DD$SP_vals_with_params$site, 
                                                         output_uncert_all_DD$SP_vals_with_params$site=="Tamakin Dacot", 
                                                         "Tomakin Dako")
SP_best_est_DD$site <- replace(SP_best_est_DD$site, SP_best_est_DD$site=="Tamakin Dacot", "Tomakin Dako")

# SP (accounting for DD)
SP_plot_DD <- ggplot(data = output_uncert_all_DD$SP_vals_with_params, aes(x=reorder(site, org_geo_order), y=SP)) +
  geom_violin(fill="grey") +
  geom_point(data = SP_best_est_DD, aes(x = site, y = SP_value), color = "black") +
  xlab("\npatch") + ylab(bquote("self persistence (SP)")) + #ggtitle("Self-persistence") +
  #ylim(c(0,0.65)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  #theme(axis.title.x = element_text(margin = margin(r=500)))

# NP (accounting for DD)
NP_plot_DD <- ggplot(data = output_uncert_all_DD$NP_out_df, aes(x=value)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = NP_best_est_DD, color = "black") +
  xlab(expression(lambda[c])) + #ggtitle('Network persistence') +
  theme_bw() 

NP_plot_DD_freq <- ggplot(data = output_uncert_all_DD$NP_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = NP_best_est_DD, color = "black") +
  xlab(expression(lambda[c])) + ylab("relative frequency") + #ggtitle('Network persistence') +
  theme_bw() 

# realized connectivity matrix
realized_C_plot_DD <- ggplot(data = best_est_metrics_mean_offspring_DD$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_R)) +
  scale_fill_gradient(high='black', low='white', name='Recruits') +
  xlab('\norigin') + ylab('destination') + #ggtitle('Realized connectivity matrix') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(angle=45,hjust=1)) 

# pdf(file = here::here('Plots/FigureDrafts', 'SP_NP_connMatrixR.pdf'), width=10, height=5)  # hacked the color scales being comparable, deal with for real if people like showing both
# plot_grid(SP_plot_DD, realized_C_plot_DD, NP_plot_DD, rel_widths=c(1,1.5,1), labels = c("a","b","c"), nrow=1)
# dev.off()

# All three together - NP as histogram
pdf(file = here::here('Plots/FigureDrafts', 'SP_NP_connMatrixR.pdf'), width=11, height=5)  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(SP_plot_DD, realized_C_plot_DD, NP_plot_DD, rel_widths=c(1.2,1.5,1.2), labels = c("a","b","c"), nrow=1)
dev.off()

# All three together - NP as frequency
pdf(file = here::here('Plots/FigureDrafts', 'SP_NP_connMatrixR_freq.pdf'), width=11, height=5)  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(SP_plot_DD, realized_C_plot_DD, NP_plot_DD_freq, rel_widths=c(1.2,1.5,1.2), labels = c("a","b","c"), nrow=1)
dev.off()

##### Figure 6 (what ifs) - new version where there are the four what-if subplots  a) more habitat, b) wider region at same habitat density, c) wider region with 100% hab, d) larval navigation together
# A: more habitat
perc_hab_plot <- ggplot(data = perc_hab_plot_df %>% filter(run==1), aes(x=perc_hab, y=value)) +
  geom_line(color="dark gray", alpha=0.2) +
  xlab("percent habitat") + ylab(expression(lambda[c])) +
  theme_bw()
for(i in 2:n_runs) {
  perc_hab_plot <- perc_hab_plot +
    geom_line(data=perc_hab_plot_df %>% filter(run==i), color="dark gray", alpha=0.2)
}
perc_hab_plot <- perc_hab_plot +
  geom_line(data = perc_hab_best_est_NP_df, aes(x=perc_hab,y=value),color="black") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = prop_sampling_area_habitat, color = "orange") 

# B: wider region, same habitat density
wider_region_plot <- ggplot(data = wider_region_plot_df %>% filter(run==1), aes(x=region_width, y=value)) +
  geom_line(color="dark gray", alpha=0.2) +
  xlab("region width (km)") + ylab(expression(lambda[c])) +
  theme_bw()
for(i in 2:n_runs) {
  wider_region_plot <- wider_region_plot +
    geom_line(data=wider_region_plot_df %>% filter(run==i), color="dark gray", alpha=0.2)
}
wider_region_plot <- wider_region_plot +
  geom_line(data = wider_region_best_est_NP_df, aes(x=region_width,y=value),color="black") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = region_width_km, color = "orange") 

# C: wider region, perc habitat density
# color schemes are from ColorBrewer: https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=5
lambda_c_color_scheme_6 <- c("#ffffcc","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494")  # 6 colors
lambda_c_color_scheme_5 <- c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494")  # 5 colors
lambda_c_color_scheme_4 <- c("#ffffcc","#a1dab4","#41b6c4","#225ea8")  # 4 colors
lambda_c_color_scheme_2 <- c("#c7e9b4","#253494")  # 2 colors

# Version 1 - kind of like a heat map
wider_region_perc_hab_plot_v1 <- ggplot(wider_region_perc_hab_plot_df, aes(x=region_width, y=perc_hab, color=NP, fill=NP)) +
  geom_point(size=4, shape=15) +
  scale_color_gradient(labels = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_fill_gradient(labels = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  geom_hline(yintercept = Ps, color = "orange") +
  geom_vline(xintercept = region_width_km, color = "orange") +
  xlab("region width (km)") + ylab("percent habitat") + labs(fill = expression(lambda[c]), color = expression(lambda[c])) +
  theme_bw() 

# Version 2 - kind of like a heat map but categories instead of continuous
wider_region_perc_hab_plot_v2 <- ggplot(wider_region_perc_hab_plot_df_categories, aes(x=region_width, y=perc_hab, color=reorder(NP_category,order_6), fill=reorder(NP_category,order_6))) +
  geom_point(size=4, shape=15) +
  scale_color_manual(values = lambda_c_color_scheme_6) +
  scale_fill_manual(values = lambda_c_color_scheme_6) +
  geom_hline(yintercept = Ps, color = "orange") +
  geom_vline(xintercept = region_width_km, color = "orange") +
  xlab("region width (km)") + ylab("percent habitat") + labs(fill = expression(lambda[c]), color = expression(lambda[c])) +
  theme_bw() 

# Version 3 - only 5 categories
wider_region_perc_hab_plot_v3 <- ggplot(wider_region_perc_hab_plot_df_categories, aes(x=region_width, y=perc_hab, color=reorder(NP_categories,order_5), fill=reorder(NP_categories,order_5))) +
  geom_point(size=4, shape=15) +
  scale_color_manual(values = lambda_c_color_scheme_5) +
  scale_fill_manual(values = lambda_c_color_scheme_5) +
  geom_hline(yintercept = Ps, color = "orange") +
  geom_vline(xintercept = region_width_km, color = "orange") +
  theme(legend.title = element_text( size=2), legend.text=element_text(size=2)) +
  theme(legend.key.size = unit(1, "cm")) +
  xlab("region width (km)") + ylab("percent habitat") + labs(fill = expression(lambda[c]), color = expression(lambda[c])) +
  theme_bw() 

# Version 4 - only 2 categories
wider_region_perc_hab_plot_v4 <- ggplot(wider_region_perc_hab_plot_df_categories, aes(x=region_width, y=perc_hab, color=reorder(NP_two_category,order_2), fill=reorder(NP_two_category,order_2))) +
  geom_point(size=4, shape=15) +
  scale_color_manual(values = lambda_c_color_scheme_2) +
  scale_fill_manual(values = lambda_c_color_scheme_2) +
  geom_hline(yintercept = Ps, color = "orange") +
  geom_vline(xintercept = region_width_km, color = "orange") +
  xlab("region width (km)") + ylab("percent habitat") + labs(fill = expression(lambda[c]), color = expression(lambda[c])) +
  theme_bw() 

# Version 5 - 3D - with plotly
wider_region_perc_hab_plot_3D <- plot_ly(y=wider_region_perc_hab_plot_df$region_width, x=wider_region_perc_hab_plot_df$perc_hab, z=wider_region_perc_hab_plot_df$NP, type="scatter3d", mode="markers",color=wider_region_perc_hab_plot_df$NP)

# With gg3d
# wider_region_perc_hab_plot_3D <- ggplot(wider_region_perc_hab_plot_df, aes(x=region_width, y=NP, z=perc_hab, color=NP, fill=NP)) +
#   theme_void() +
#   axes_3D() +
#   stat_3D() +
#   xlab("region width (km)") + ylab(expression(lambda[c])) 

plot_ly(x=temp, y=pressure, z=dtime, type="scatter3d", mode="markers", color=temp)
                                    
# # C: wider region, 100% habitat density (original plot C)
# wider_region_all_hab_plot <- ggplot(data = wider_region_all_hab_plot_df %>% filter(run==1), aes(x=region_width, y=value)) +
#   geom_line(color="dark gray", alpha=0.2) +
#   xlab("region width (km)") + ylab(expression(lambda[c])) +
#   theme_bw()
# for(i in 2:n_runs) {
#   wider_region_all_hab_plot <- wider_region_all_hab_plot +
#     geom_line(data=wider_region_all_hab_plot_df %>% filter(run==i), color="dark gray", alpha=0.2)
# }
# wider_region_all_hab_plot <- wider_region_all_hab_plot +
#   geom_line(data = wider_region_all_hab_best_est_NP_df, aes(x=region_width,y=value),color="black") +
#   geom_hline(yintercept = 1, color = "blue") 

# D: larval navigation
larv_nav_plot <- ggplot(data = larv_nav_plot_df %>% filter(run==1), aes(x=larv_nav, y=value)) +
  geom_line(color="gray", alpha=0.2) +
  xlab("larval navigation (km)") + ylab(expression(lambda[c])) +
  theme_bw()
for(i in 2:n_runs) {
  larv_nav_plot <- larv_nav_plot +
    geom_line(data=larv_nav_plot_df %>% filter(run==i), color="gray", alpha=0.2)
}
larv_nav_plot <- larv_nav_plot +
  geom_line(data = larv_nav_best_est_NP_df, aes(x=larv_nav,y=value),color="black") +
  geom_hline(yintercept = 1, color = "blue") 

# All together
pdf(file=here::here("Plots/FigureDrafts","What_if_4_panels.pdf"), width=6, height=6)
plot_grid(perc_hab_plot, wider_region_plot, wider_region_all_hab_plot, larv_nav_plot, 
          nrow=2, labels=c("a","b","c","d"))
dev.off()

# All together with a sample 3d plot
pdf(file=here::here("Plots/FigureDrafts","What_if_4_panels_3D.pdf"), width=6, height=6)
plot_grid(perc_hab_plot, wider_region_plot, wider_region_perc_hab_plot_v3, larv_nav_plot, 
          nrow=2, labels=c("a","b","c","d"))
dev.off()

# Pull out legend from 3D plot
legend_3D <- get_legend(wider_region_perc_hab_plot_v3)
wider_region_perc_hab_plot_v3 <- wider_region_perc_hab_plot_v3 + theme(legend.position = "none")
wider_region_perc_hab_plot_v3_with_legend <- plot_grid(wider_region_perc_hab_plot_v3, legend_3D, nrow=1, rel_widths = c(0.68, 0.32))
pdf(file=here::here("Plots/FigureDrafts","What_if_4_panels_3D.pdf"), width=6, height=6)
plot_grid(perc_hab_plot, wider_region_plot, wider_region_perc_hab_plot_v3_with_legend, larv_nav_plot, 
          nrow=2, labels=c("a","b","c","d"))
dev.off()

                                                       
                                                       
top_row_what_ifs_plot <- plot_grid(perc_hab_plot, wider_region_plot, nrow=1, labels=c("a","b"))
bottom_row_what_ifs_plot <- plot_grid(wider_region_perc_hab_plot_v3, legend_3D, larv_nav_plot, nrow=1, 
                        rel_widths = c(0.70,0.3,1), labels = c("c",NULL,"d"))
pdf(file=here::here("Plots/FigureDrafts","What_if_4_panels_3D.pdf"), width=6, height=6)
plot_grid(top_row_what_ifs_plot, bottom_row_what_ifs_plot,
          nrow=2)
dev.off()
                        
legend <- get_legend(plot.mpg)

# and replot suppressing the legend
plot.mpg <- plot.mpg + theme(legend.position='none')

# Now plots are aligned vertically with the legend to the right
ggdraw(plot_grid(plot_grid(p1, plot.mpg, ncol=1, align='v'),
                 plot_grid(NULL, legend, ncol=1),
                 rel_widths=c(1, 0.2)))
# Different 3D plot options
pdf(file=here::here("Plots/FigureDrafts","What_if_4_panels_3D_options.pdf"), width=6, height=6)
plot_grid(wider_region_perc_hab_plot_v1, wider_region_perc_hab_plot_v3, wider_region_perc_hab_plot_v4, NULL,
          nrow=2, labels=c("a","b","c",NULL))
dev.off()


#### Also, save them individually...
# Larval navigation
pdf(file=here::here("Plots/FigureDrafts","larv_nav_NP.pdf"))
larv_nav_plot
dev.off()

# Previous versions of Fig 6
##### Figure 6 (what ifs)
## NP by perc hab - actual site-specific survivals
# showing min and max of NP in ribbon
NP_perc_hab_realSurvs_plot_min_max <- ggplot(data = NP_by_perc_hab_realSurvs, aes(x=perc_hab, y=NP, ymin=NP_min, ymax=NP_max)) +
  geom_line(color="black") +
  geom_ribbon(alpha=0.5, color="gray", fill="gray") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = prop_sampling_area_habitat, color = "orange") +
  xlab("proportion habitat") + ylab("NP (site-specific survivals)") +
  theme_bw()

# showing +- sd in ribbon
NP_perc_hab_realSurvs_plot_sd <- ggplot(data = NP_by_perc_hab_realSurvs, aes(x=perc_hab, y=NP, ymin=NP-NP_sd, ymax=NP+NP_sd)) +
  geom_line(color="black") +
  geom_ribbon(alpha=0.5, color="gray", fill="gray") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = prop_sampling_area_habitat, color = "orange") +
  xlab("proportion habitat") + ylab("NP (site-specific survivals)") +
  theme_bw()
  
## NP by perc hab - average site-specific survivals
# showing min and max of NP in ribbon
NP_perc_hab_avgSurvs_plot_min_max <- ggplot(data = NP_by_perc_hab_avgSurvs, aes(x=perc_hab, y=NP, ymin=NP_min, ymax=NP_max)) +
  geom_line(color="black") +
  geom_ribbon(alpha=0.5, color="gray", fill="gray") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = prop_sampling_area_habitat, color = "orange") +
  xlab("proportion habitat") + ylab("NP (average survivals)") +
  theme_bw()

# showing +- sd in ribbon
NP_perc_hab_avgSurvs_plot_sd <- ggplot(data = NP_by_perc_hab_avgSurvs, aes(x=perc_hab, y=NP, ymin=NP-NP_sd, ymax=NP+NP_sd)) +
  geom_line(color="black") +
  geom_ribbon(alpha=0.5, color="gray", fill="gray") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = prop_sampling_area_habitat, color = "orange") +
  #xlab("proportion habitat") + ylab("NP (average survivals)") +
  xlab("proportion habitat") + ylab(bquote("network persistence (NP"[DD] ~")")) +
  theme_bw()

# percent of runs with NP>1, realSurvs
perc_hab_realSurvs_persistent_plot <- ggplot(data = NP_above1_perc_hab_realSurvs, aes(x=perc_hab, y=perc_persistent*100)) +
  geom_line(color="black") +
  geom_vline(xintercept = prop_sampling_area_habitat, color = "orange") +
  xlab("proportion habitat") + ylab("percent runs persistent") +
  theme_bw()

# percent of runs with NP>1, avgSurvs
perc_hab_avgSurvs_persistent_plot <- ggplot(data = NP_above1_perc_hab_avgSurvs, aes(x=perc_hab, y=perc_persistent*100)) +
  geom_line(color="black") +
  geom_vline(xintercept = prop_sampling_area_habitat, color = "orange") +
  xlab("proportion habitat") + ylab("% estimates persistent") +
  theme_bw()

### Put them together
# with sd as ribbon
pdf(file=here::here("Plots/FigureDrafts","NP_by_per_hab_sd_ribbon.pdf"), width=6, height=3)
plot_grid(NP_perc_hab_realSurvs_plot_sd, NP_perc_hab_avgSurvs_plot_sd, nrow=1, labels=c("a","b"))
dev.off()

# with min and max as ribbon
pdf(file=here::here("Plots/FigureDrafts","NP_by_per_hab_min_max_ribbon.pdf"), width=6, height=3)
plot_grid(NP_perc_hab_realSurvs_plot_min_max, NP_perc_hab_avgSurvs_plot_min_max, nrow=1, labels=c("a","b"))
dev.off()

# avgSurvs, ribbon with sd and % runs persistent -THIS IS THE ONE IN THE DRAFT RIGHT NOW!
pdf(file=here::here("Plots/FigureDrafts","NP_perc_hab_perc_persist_sd_ribbon.pdf"), width=6, height=3)
plot_grid(NP_perc_hab_avgSurvs_plot_sd, perc_hab_avgSurvs_persistent_plot, nrow=1, labels=c("a","b"))
dev.off()

# just average survs, sd
pdf(file=here::here("Plots/FigureDrafts","NP_by_perc_hab_avgSurvs_sd.pdf"))
ggplot(data = NP_by_perc_hab_avgSurvs, aes(x=perc_hab, y=NP, ymin=NP-NP_sd, ymax=NP+NP_sd)) +
  geom_line(color="black") +
  geom_ribbon(alpha=0.5, color="gray", fill="gray") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = prop_sampling_area_habitat, color = "orange") +
  xlab("proportion habitat") + ylab(bquote("network persistence (NP"[DD] ~")")) +
  theme_bw()
dev.off()

# just average survs, min max
pdf(file=here::here("Plots/FigureDrafts","NP_by_perc_hab_avgSurvs_min_max.pdf"))
ggplot(data = NP_by_perc_hab_avgSurvs, aes(x=perc_hab, y=NP, ymin=NP_min, ymax=NP_max)) +
  geom_line(color="black") +
  geom_ribbon(alpha=0.5, color="gray", fill="gray") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = prop_sampling_area_habitat, color = "orange") +
  xlab("proportion habitat") + ylab("NP") +
  theme_bw()
dev.off()

##### Fig 6 alternates (other what-ifs) - wider region, same perc hab
# NP across region widths showing +- sd in ribbon
NP_wider_region_avgSurvs_plot_sd <- ggplot(data = NP_by_wider_region_avgSurvs, aes(x=region_width, y=NP, ymin=NP-NP_sd, ymax=NP+NP_sd)) +
  geom_line(color="black") +
  geom_ribbon(alpha=0.5, color="gray", fill="gray") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = region_width_km, color = "orange") +
  #xlab("proportion habitat") + ylab("NP (average survivals)") +
  xlab("region width (km)") + ylab(bquote("network persistence (NP"[DD] ~")")) +
  theme_bw()

# percent of runs with NP>1, avgSurvs
wider_region_avgSurvs_persistent_plot <- ggplot(data = NP_above1_wider_region_avgSurvs, aes(x=region_width, y=perc_persistent*100)) +
  geom_line(color="black") +
  geom_vline(xintercept = region_width_km, color = "orange") +
  xlab("region width (km)") + ylab("% estimates persistent") +
  theme_bw()

# both together
pdf(file=here::here("Plots/FigureDrafts","NP_wider_region_perc_persist_sd_ribbon.pdf"), width=6, height=3)
plot_grid(NP_wider_region_avgSurvs_plot_sd, wider_region_avgSurvs_persistent_plot, nrow=1, labels=c("a","b"))
dev.off()

##### Fig 6 alternates (other what-ifs) - wider region, 100% perc hab
# NP across region widths showing +- sd in ribbon
NP_wider_region_all_hab_avgSurvs_plot_sd <- ggplot(data = NP_by_wider_region_all_hab_avgSurvs, aes(x=region_width, y=NP, ymin=NP-NP_sd, ymax=NP+NP_sd)) +
  geom_line(color="black") +
  geom_ribbon(alpha=0.5, color="gray", fill="gray") +
  geom_hline(yintercept = 1, color = "blue") +
  #xlab("proportion habitat") + ylab("NP (average survivals)") +
  xlab("region width (km)") + ylab(bquote("network persistence (NP"[DD] ~")")) +
  theme_bw()

# percent of runs with NP>1, avgSurvs
wider_region_all_hab_avgSurvs_persistent_plot <- ggplot(data = NP_above1_wider_region_all_hab_avgSurvs, aes(x=region_width, y=perc_persistent*100)) +
  geom_line(color="black") +
  xlab("region width (km)") + ylab("% estimates persistent") +
  theme_bw()

# both together
pdf(file=here::here("Plots/FigureDrafts","NP_wider_region_all_hab_perc_persist_sd_ribbon.pdf"), width=6, height=3)
plot_grid(NP_wider_region_all_hab_avgSurvs_plot_sd, wider_region_all_hab_avgSurvs_persistent_plot, nrow=1, labels=c("a","b"))
dev.off()

##### Fig 6 alternates (other what-ifs) - larval navigation
# 200 meters either side of patch
SP_plot_DD_200m <- ggplot(data = output_uncert_all_DD_larv_nav_200m$SP_vals_with_params, aes(x=reorder(site, org_geo_order), y=SP)) +
  geom_violin(fill="grey") +
  geom_point(data = best_est_metrics_larv_nav_200m_DD$SP, aes(x = site, y = SP_value), color = "black") +
  xlab("\nsite") + ylab(bquote("self persistence (SP"[DD] ~")")) + #ggtitle("Self-persistence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

NP_plot_DD_freq_200m <- ggplot(data = output_uncert_all_DD_larv_nav_200m$NP_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = best_est_metrics_larv_nav_200m_DD$NP, color = "black") +
  xlab(bquote("network persistence (NP"[DD] ~")")) + ylab("relative frequency") + #ggtitle('Network persistence') +
  theme_bw() 

realized_C_plot_DD_200m <- ggplot(data = best_est_metrics_larv_nav_200m_DD$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_R)) +
  scale_fill_gradient(high='black', low='white', name='Recruits') +
  xlab('\norigin') + ylab('destination') + #ggtitle('Realized connectivity matrix') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(angle=45,hjust=1)) 

pdf(file = here::here('Plots/FigureDrafts', 'SP_NP_connMatrixR_freq_200m.pdf'), width=11, height=5)  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(SP_plot_DD_200m, realized_C_plot_DD_200m, NP_plot_DD_freq_200m, rel_widths=c(1.2,1.5,1.2), labels = c("a","b","c"), nrow=1)
dev.off()

# 500 meters either side of patch
SP_plot_DD_500m <- ggplot(data = output_uncert_all_DD_larv_nav_500m$SP_vals_with_params, aes(x=reorder(site, org_geo_order), y=SP)) +
  geom_violin(fill="grey") +
  geom_point(data = best_est_metrics_larv_nav_500m_DD$SP, aes(x = site, y = SP_value), color = "black") +
  xlab("\nsite") + ylab(bquote("self persistence (SP"[DD] ~")")) + #ggtitle("Self-persistence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

NP_plot_DD_freq_500m <- ggplot(data = output_uncert_all_DD_larv_nav_500m$NP_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = best_est_metrics_larv_nav_500m_DD$NP, color = "black") +
  xlab(bquote("network persistence (NP"[DD] ~")")) + ylab("relative frequency") + #ggtitle('Network persistence') +
  theme_bw() 

realized_C_plot_DD_500m <- ggplot(data = best_est_metrics_larv_nav_500m_DD$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_R)) +
  scale_fill_gradient(high='black', low='white', name='Recruits') +
  xlab('\norigin') + ylab('destination') + #ggtitle('Realized connectivity matrix') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(angle=45,hjust=1)) 

pdf(file = here::here('Plots/FigureDrafts', 'SP_NP_connMatrixR_freq_500m.pdf'), width=11, height=5)  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(SP_plot_DD_500m, realized_C_plot_DD_500m, NP_plot_DD_freq_500m, rel_widths=c(1.2,1.5,1.2), labels = c("a","b","c"), nrow=1)
dev.off()

# 1000 meters either side of patch
SP_plot_DD_1000m <- ggplot(data = output_uncert_all_DD_larv_nav_1000m$SP_vals_with_params, aes(x=reorder(site, org_geo_order), y=SP)) +
  geom_violin(fill="grey") +
  geom_point(data = best_est_metrics_larv_nav_1000m_DD$SP, aes(x = site, y = SP_value), color = "black") +
  xlab("\nsite") + ylab(bquote("self persistence (SP"[DD] ~")")) + #ggtitle("Self-persistence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

NP_plot_DD_freq_1000m <- ggplot(data = output_uncert_all_DD_larv_nav_1000m$NP_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = best_est_metrics_larv_nav_1000m_DD$NP, color = "black") +
  xlab(bquote("network persistence (NP"[DD] ~")")) + ylab("relative frequency") + #ggtitle('Network persistence') +
  theme_bw() 

realized_C_plot_DD_1000m <- ggplot(data = best_est_metrics_larv_nav_1000m_DD$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_R)) +
  scale_fill_gradient(high='black', low='white', name='Recruits') +
  xlab('\norigin') + ylab('destination') + #ggtitle('Realized connectivity matrix') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(angle=45,hjust=1)) 

# All three together - NP as frequency
pdf(file = here::here('Plots/FigureDrafts', 'SP_NP_connMatrixR_freq_1000m.pdf'), width=11, height=5)  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(SP_plot_DD_1000m, realized_C_plot_DD_1000m, NP_plot_DD_freq_1000m, rel_widths=c(1.2,1.5,1.2), labels = c("a","b","c"), nrow=1)
dev.off()

########## Appendix figures #########

##### Proportion of kernel sampled - made elsewhere

##### Survival by size and site
surv_by_site_to_plot <- best_fit_model_dfs$surv_site_size %>%  # change site from factor to character so can update spelling of Tomakin Dako
  mutate_if(is.factor, as.character) %>%
  mutate(site = if_else(site == "Tamakin Dacot", "Tomakin Dako", site))

pdf(file = here::here("Plots/FigureDrafts", "APP_FIG_surv_by_size_and_site_Phisiteplussize_psizeplusdist.pdf"))
ggplot(data = surv_by_site_to_plot, aes(x = size, y = estimate_prob, fill = site)) +
  geom_line(aes(color=site)) + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab("survival probability") + xlab("size (cm)") + #ggtitle("Phi:site+size, p:size+dist") +
  facet_wrap(~site) +
  #facet_wrap(~site, labeller=label_context(labels, multi_line = TRUE))+
  theme_bw() 
dev.off()

# labels= c("Cabatoan", "Caridad\n Cemetery", "Caridad\n Proper", "Elementary\n School", "Gabas",
#           "Haina", "Hicgop South", "N.\n Magbangon", "Palanas", "Poroc Rose", "Poroc San\n Flower",
#           "S.\n Magbangon", "San Agustin", "Sitio\n Baybayon", "Sitio Lonas", "Sitio Tugas",
#           "Tomakin\n Dako", "Visca", "Wangag")) 

##### Recapture probability by distance and by size
# p (by dist)
p_by_dist_plot <- ggplot(data = best_fit_model_dfs$recap_dist, aes(x = dist, y = estimate_prob)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab("recapture probability") + xlab("distance (m)") + #ggtitle("Phi:site+size, p:size+dist") + 
  theme_bw() 

# p (by size)
p_by_size_plot <- ggplot(data = best_fit_model_dfs$recap_size, aes(x = size, y = estimate_prob)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab("recapture probability") + xlab("size (cm)") + #ggtitle("Phi:site+size, p:size+dist") + 
  theme_bw() 

# Put them together
pdf(file = here::here("Plots/FigureDrafts", "APP_FIG_recap_effects_Phisiteplussize_psizeplusdist.pdf"), width=7, height=4)
plot_grid(p_by_dist_plot, p_by_size_plot, labels=c("a","b"), nrow=1)
dev.off()

##### Plot of LRP and local replacement without DD accounted for
# LRP without DD 
LEP_R_plot <- ggplot(data = output_uncert_all$LEP_R_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = 'gray', fill = 'gray') +
  #geom_vline(data = LEP_R_best_est, aes(xintercept = LEP_R, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "mean offspring"))$LEP_R, color = 'black') +
  geom_vline(xintercept = LEP_R_best_est, color = "black") +
  xlab(expression("lifetime recruit production (LRP"[ DD]~")")) + #ggtitle('LRP') +
  theme_bw()

LEP_R_plot_freq <- ggplot(data = output_uncert_all$LEP_R_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color = 'gray', fill = 'gray') +
  #geom_vline(data = LEP_R_best_est, aes(xintercept = LEP_R, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "mean offspring"))$LEP_R, color = 'black') +
  geom_vline(xintercept = LEP_R_best_est, color = "black") +
  xlab(expression("lifetime recruit production (LRP"[DD]~")")) + ylab("relative frequency") + #ggtitle('LRP') +
  theme_bw()

# LRP_local without DD
LEP_R_local_plot <- ggplot(data = output_uncert_all$LEP_R_local_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = "gray", fill = "gray") +
  geom_vline(xintercept = LEP_R_local_best_est, color = "black") +
  xlab(expression("local replacement (LR"[DD]~")")) + #ggtitle("Local replacement") +
  theme_bw()

LEP_R_local_plot_freq <- ggplot(data = output_uncert_all$LEP_R_local_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color = "gray", fill = "gray") +
  geom_vline(xintercept = LEP_R_local_best_est, color = "black") +
  xlab(expression("local replacement (LR"[DD]~")")) + ylab("relative frequency") + #ggtitle("Local replacement") +
  theme_bw()

# Put them together - histograms
pdf(file = here::here('Plots/FigureDrafts', 'APP_FIG_LRP_LocalReplacement_withoutDDconsidered.pdf'), width=6, height=3)
plot_grid(LEP_R_plot, LEP_R_local_plot, labels = c("a","b"), nrow=1)
dev.off()

# Put them together - frequencies
pdf(file = here::here('Plots/FigureDrafts', 'APP_FIG_LRP_LocalReplacement_withoutDDconsidered_freq.pdf'), width=6, height=3)
plot_grid(LEP_R_plot_freq, LEP_R_local_plot_freq, labels = c("a","b"), nrow=1)
dev.off()

##### SP metrics, NP metrics, connectivity matrices without density dependence considered
# Update site name (Tomakin Dako)
best_est_metrics_mean_offspring$Cmat$org_site <- replace(best_est_metrics_mean_offspring$Cmat$org_site, 
                                                            best_est_metrics_mean_offspring$Cmat$org_site=="Tamakin Dacot", 
                                                            "Tomakin Dako")
best_est_metrics_mean_offspring$Cmat$dest_site <- replace(best_est_metrics_mean_offspring$Cmat$dest_site, 
                                                             best_est_metrics_mean_offspring$Cmat$dest_site=="Tamakin Dacot", 
                                                             "Tomakin Dako")
output_uncert_all$SP_vals_with_params$site <- replace(output_uncert_all$SP_vals_with_params$site, 
                                                         output_uncert_all$SP_vals_with_params$site=="Tamakin Dacot", 
                                                         "Tomakin Dako")
SP_best_est$site <- replace(SP_best_est$site, SP_best_est$site=="Tamakin Dacot", "Tomakin Dako")

# SP (not accounting for DD)
SP_plot <- ggplot(data = output_uncert_all$SP_vals_with_params, aes(x=reorder(site, org_geo_order), y=SP)) +
  geom_violin(fill="grey") +
  geom_point(data = SP_best_est, aes(x = site, y = SP_value), color = "black") +
  xlab("\npatch") + ylab(expression("self persistence (SP"[DD]~")")) + #ggtitle("Self-persistence") +
  #ylim(c(0,0.65)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

# NP (not accounting for DD)
NP_plot <- ggplot(data = output_uncert_all$NP_out_df, aes(x=value)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = NP_best_est, color = "black") +
  xlab(expression(lambda[DD])) + #ggtitle('Network persistence') +
  theme_bw()

NP_plot_freq <- ggplot(data = output_uncert_all$NP_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = NP_best_est, color = "black") +
  xlab(expression(lambda[DD])) + ylab("relative frequency") + #ggtitle('Network persistence') +
  theme_bw()

# realized connectivity matrix (not accounting for DD)
realized_C_plot <- ggplot(data = best_est_metrics_mean_offspring$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_R)) +
  scale_fill_gradient(high='black', low='white', name=expression('Recruits'[DD])) +
  xlab('\norigin') + ylab('destination') + #ggtitle('Realized connectivity matrix') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(legend.position = "bottom")

# put together (NP as histogram)
pdf(file = here::here('Plots/FigureDrafts', 'APP_FIG_SP_NP_connMatrixR_withoutDDcompensation.pdf'), width=11, height=5)  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(SP_plot, realized_C_plot, NP_plot, rel_widths=c(1.2,1.5,1.2), labels = c("a","b","c"), nrow=1)
dev.off()

# put together (NP as frequency)
pdf(file = here::here('Plots/FigureDrafts', 'APP_FIG_SP_NP_connMatrixR_withoutDDcompensation_freq.pdf'), width=11, height=5)  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(SP_plot, realized_C_plot, NP_plot_freq, rel_widths=c(1.2,1.5,1.2), labels = c("a","b","c"), nrow=1)
dev.off()

##### LEP by site (with DD compensation)
LEP_by_site_to_plot <- left_join(output_uncert_all_DD$LEP_by_site_out_df, site_vec_order, by = c("site"="alpha_order"))  # add in site name

pdf(file = here::here("Plots/FigureDrafts","APP_FIG_LEP_by_site.pdf"))
ggplot(data = LEP_by_site_to_plot, aes(x=value)) +
  #geom_histogram(bins=100, color = 'gray', fill = 'gray') +
  geom_histogram(aes(y=..count../sum(..count..)), bins=100, color="gray", fill="gray") + 
  #geom_vline(data = LEP_best_est, aes(xintercept = LEP, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_best_est %>% filter(recruit_size == "mean offspring"))$LEP, color='black') +
  geom_vline(xintercept = LEP_best_est, color = "black") +
  xlab('lifetime egg production (LEP)') + ylab("relative frequency") + #ggtitle('LEP by site') +
  theme_bw() +
  facet_wrap(~site_name) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### LRP by site (with DD compensation)
LRP_by_site_to_plot <- left_join(output_uncert_all_DD$LEP_R_by_site_out_df, site_vec_order, by = c("site"="alpha_order"))  # add in site name

pdf(file = here::here("Plots/FigureDrafts","APP_FIG_LRP_by_site.pdf"))
ggplot(data = LRP_by_site_to_plot, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=100, color="gray", fill="gray") + 
  #geom_histogram(bins=100, color = 'gray', fill = 'gray') +
  #geom_vline(data = LEP_best_est, aes(xintercept = LEP, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_best_est %>% filter(recruit_size == "mean offspring"))$LEP, color='black') +
  geom_vline(xintercept = LEP_R_best_est, color = "black") +
  xlab("lifetime recruit production (LRP)") + ylab("relative frequency") +
  theme_bw() +
  facet_wrap(~site_name)
dev.off()


##### Uncertainty exploration for LEP
# All together
pdf(file = here::here("Plots/FigureDrafts", "LEP_uncertainty_by_param.pdf"), width=6, height=4)
ggplot(data = LEP_uncert_DD %>% filter(uncertainty_type %in% c("start recruit size", "breeding size", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("start recruit size", "breeding size", "growth", "survival", "all"),
                               value = LEP_best_est), color = "black") +
  xlab("uncertainty type") + ylab("lifetime egg production (LEP)") + #ggtitle("Uncertainty in LEP") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  theme(text=element_text(size=12))
dev.off()

# # start recruit size
# LEP_start_recruit_uncertainty_plot <- ggplot(data = LEP_uncert %>% filter(uncertainty_type %in% c("start recruit size")), aes(x=uncertainty_type, y=value)) +
#   geom_violin(fill="grey") +
#   geom_point(data = data.frame(uncertainty_type = c("start recruit size"), value = LEP_best_est), color = "black") +
#   xlab("LEP") + ggtitle("Recruit size effects") +
#   theme_bw() +
#   theme(axis.text.x = element_blank())
# 
# # breeding size
# LEP_breeding_size_uncertainty_plot <- ggplot(data = LEP_uncert %>% filter(uncertainty_type %in% c("breeding size")), aes(x=uncertainty_type, y=value)) +
#   geom_violin(fill="grey") +
#   geom_point(data = data.frame(uncertainty_type = c("breeding size"), value = LEP_best_est), color = "black") +
#   xlab("LEP") + ggtitle("Breeding size effects") +
#   theme_bw() +
#   theme(axis.text.x = element_blank())
# 
# # growth
# LEP_growth_uncertainty_plot <- ggplot(data = LEP_uncert %>% filter(uncertainty_type %in% c("growth")), aes(x=uncertainty_type, y=value)) +
#   geom_violin(fill="grey") +
#   geom_point(data = data.frame(uncertainty_type = c("growth"), value = LEP_best_est), color = "black") +
#   xlab("LEP") + ggtitle("Growth curve effects") +
#   theme_bw() +
#   theme(axis.text.x = element_blank())
# 
# # survival
# LEP_survival_uncertainty_plot <- ggplot(data = LEP_uncert %>% filter(uncertainty_type %in% c("survival")), aes(x=uncertainty_type, y=value)) +
#   geom_violin(fill="grey") +
#   geom_point(data = data.frame(uncertainty_type = c("survival"), value = LEP_best_est), color = "black") +
#   xlab("LEP") + ggtitle("Survival effects") +
#   theme_bw() +
#   theme(axis.text.x = element_blank())
# 
# # all
# LEP_all_uncertainty_plot <- ggplot(data = LEP_uncert %>% filter(uncertainty_type %in% c("all")), aes(x=uncertainty_type, y=value)) +
#   geom_violin(fill="grey") +
#   geom_point(data = data.frame(uncertainty_type = c("all"), value = LEP_best_est), color = "black") +
#   xlab("LEP") + ggtitle("All uncertainty effects") +
#   theme_bw() +
#   theme(axis.text.x = element_blank())

# # put together
# pdf(file = here::here("Plots/FigureDrafts", "LEP_uncertainty_by_param.pdf"))
# plot_grid(LEP_uncertainty_plot, LEP_start_recruit_uncertainty_plot, LEP_breeding_size_uncertainty_plot, LEP_growth_uncertainty_plot, 
#           LEP_survival_uncertainty_plot, LEP_all_uncertainty_plot, labels=c("a","b","c","d","e","f"), nrow=3)
# dev.off()
# 
# ggplot(data = output_uncert_survival_DD$LEP_by_site_out_df, aes(x=uncertainty_type, y=value)) +
#   geom_violin(fill="grey") +
#   geom_point(data = data.frame(uncertainty_type = c("survival"), value = LEP_best_est), color = "black") +
#   xlab("LEP") + ggtitle("Survival effects") + 
#   facet_wrap(~site) +
#   theme_bw() +
#   theme(axis.text.x = element_blank())

##### Uncertainty exploration for RperE
RperE_uncert_DD_plot <- RperE_uncert_DD
RperE_uncert_DD_plot$uncertainty_type <- replace(RperE_uncert_DD_plot$uncertainty_type, 
                                                 RperE_uncert_DD_plot$uncertainty_type=="prob r", 
                                               "capture probability")
pdf(file = here::here("Plots/FigureDrafts", "RperE_uncertainty_by_param.pdf"), width=6, height=4)
ggplot(data = RperE_uncert_DD_plot %>% filter(uncertainty_type %in% c("breeding size", "assigned offspring", "capture probability", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("breeding size", "assigned offspring", "capture probability", "growth", "survival", "all"),
                               value = best_est_metrics_mean_offspring_DD$recruits_per_egg), color = "black") +
  xlab("uncertainty type") + ylab("recruits-per-egg") + #ggtitle("Uncertainty in egg-recruit survival") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(text=element_text(size=12))
dev.off()

##### Uncertainty exploration for LRP
LEP_R_uncert_DD_plot <- LEP_R_uncert_DD
LEP_R_uncert_DD_plot$uncertainty_type <- replace(LEP_R_uncert_DD_plot$uncertainty_type, 
                                                 LEP_R_uncert_DD_plot$uncertainty_type=="prob r", 
                                                 "capture probability")

pdf(file = here::here("Plots/FigureDrafts", "LRP_uncertainty_by_param.pdf"), width=6, height=4)
ggplot(data = LEP_R_uncert_DD_plot %>% filter(uncertainty_type %in% c("start recruit size", "breeding size", "assigned offspring", "capture probability", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("start recruit size", "breeding size", "assigned offspring", "capture probability", "growth", "survival", "all"),
                               value = LEP_R_best_est_DD), color = "black") +
  xlab("uncertainty type") + ylab("lifetime recruit production (LRP)") + #ggtitle("Uncertainty in LRP") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(text=element_text(size=12))
dev.off()

##### Uncertainty exploration for NP
pdf(file = here::here("Plots/FigureDrafts", "NP_uncertainty_by_param.pdf"), width=6)
ggplot(data = NP_uncert_DD %>% filter(uncertainty_type %in% c("breeding size", "assigned offspring", "prob r", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("breeding size", "assigned offspring", "prob r", "growth", "survival", "all"),
                               value = NP_best_est_DD), color = "black") +
  xlab("uncertainty type") + ylab("network persistence (NP)") + #ggtitle("Uncertainty in network persistence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
dev.off()

##### Uncertainty exploration for NP - with all included - why were start recruit size and dispersal not included before?
# rename params for easier understanding (will eventually fix in code above)
NP_uncert_DD_for_plot <- NP_uncert_DD
NP_uncert_DD_for_plot$uncertainty_type <- replace(NP_uncert_DD_for_plot$uncertainty_type, 
                                                  NP_uncert_DD_for_plot$uncertainty_type == "dispersal k",
                                                  "dispersal")
NP_uncert_DD_for_plot$uncertainty_type <- replace(NP_uncert_DD_for_plot$uncertainty_type, 
                                                  NP_uncert_DD_for_plot$uncertainty_type == "start recruit size",
                                                  "recruit census size")
NP_uncert_DD_for_plot$uncertainty_type <- replace(NP_uncert_DD_for_plot$uncertainty_type, 
                                                  NP_uncert_DD_for_plot$uncertainty_type == "prob r",
                                                  "capture probability")

pdf(file = here::here("Plots/FigureDrafts", "NP_uncertainty_by_param.pdf"), width=6, height=5)
ggplot(data = NP_uncert_DD_for_plot %>% filter(uncertainty_type %in% c("breeding size", "assigned offspring", "capture probability", "growth", "survival", "dispersal", "recruit census size", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("breeding size", "assigned offspring", "capture probability", "growth", "survival", "dispersal", "recruit census size", "all"),
                               value = NP_best_est_DD), color = "black") +
  xlab("uncertainty type") + ylab(expression(lambda)) + #ggtitle("Uncertainty in network persistence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(text = element_text(size=12))
dev.off()

##### Parameters (and their uncertainty) not shown in main text fig 
# Census (start-recruit) size 
startRecruit_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=start_recruit_size)) +
  geom_histogram(bins=40, color="gray", fill="gray") +
  geom_vline(xintercept = mean_sampled_offspring_size) +
  xlab("recruit census size (cm)") + #ggtitle("Census size") +
  theme_bw()

# Growth - Linf + k (VBL growth model) - show here too even though also in main text?
growthLinf_k_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=Linf, y=k_growth)) +
  geom_point(color = "gray", fill = "gray") +
  geom_point(x = Linf_growth_mean, y = k_growth_mean, color = "black", fill = "black") +
  xlab('Linf (cm)') + ylab("K") + #ggtitle("VBL growth model") +
  theme_bw()

# Capture probability (prob_r)
probR_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x = prob_r)) +
  geom_histogram(bins = 40, color = "gray", fill = "gray") +
  geom_vline(xintercept = prob_r_mean, color = "black") +
  xlab(bquote("P"[c])) + #ggtitle("Capture probability") +
  theme_bw()

# Assigned offspring
assignedOffspring_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x = assigned_offspring)) +
  geom_histogram(bins = 40, color = "gray", fill = "gray") +
  geom_vline(xintercept = n_offspring_matched, color = "black") +
  xlab("# assigned offspring") + #ggtitle("Assigned offspring") +
  theme_bw()

# Habitat and DD scaling
# Proportion habitat 
# P_h = cumulative proportion sites sampled across time
# P_d = proportion of the dispersal kernel from each site covered by our sampling region
# P_s = proportion of our sampling region that is habitat
# P_DD = proportion habitat available to juveniles if exclude DD

prop_scaling_vals_for_plot <- data.frame(value = c("P_h", "P_d", "P_s", "DD"),
                                         estimate = c(total_prop_hab_sampled_through_time, prop_total_disp_area_sampled_best_est,
                                                      prop_sampling_area_habitat, (perc_APCL_val+perc_UNOC_val)/perc_UNOC_val))
propHabitat_plot <- ggplot(data = prop_scaling_vals_for_plot, aes(x = value, y = estimate)) +
  geom_point() +
  xlab("recruit scaling parameter") + ylab("proportion") + #ggtitle("Habitat scaling") +
  ylim(c(0,1.8)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Put them together
pdf(file = here::here('Plots/FigureDrafts', 'APP_FIG_Parameter_inputs.pdf'))  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(startRecruit_plot, growthLinf_k_plot, probR_plot, assignedOffspring_plot, propHabitat_plot,
          labels = c("a","b","c","d","e"), nrow=2)
dev.off()

##### LEP vs. max size
pdf(file = here::here("Plots/FigureDrafts", "LEP_vs_max_size.pdf"))
ggplot(data = LEP_test_max_size , aes(x=max_size, y=LEP)) +
  geom_point() +
  geom_vline(xintercept = Linf_growth_mean, color = "red") +
  xlab("max size (cm)") + ylab("LEP") + 
  theme_bw()
dev.off()

##### Comparing min/max of grid search k and theta parameters to KC's 95% confidence calculated using Ben Bolker method
disp_95CI_surface_plot <- ggplot(data = k_theta_allyear_95CI_values, aes(x = k_eval, y = theta_eval, fill = log_like, color = log_like)) +
  geom_point(size = 4, shape = 15) +
  geom_hline(data = kernel_summary %>% filter(Year == "2012-18"), aes(yintercept = theta_CI95_upper), color = "black") +
  geom_hline(data = kernel_summary %>% filter(Year == "2012-18"), aes(yintercept = theta_CI95_lower), color = "black") +
  geom_vline(data = kernel_summary %>% filter(Year == "2012-18"), aes(xintercept = k_CI95_upper), color = "black") +
  geom_vline(data = kernel_summary %>% filter(Year == "2012-18"), aes(xintercept = k_CI95_lower), color = "black")

disp_95CI_surface_plot <- disp_95CI_surface_plot +
  geom_point(data = kernel_summary %>% filter(Year == "2012-18"), aes(x = best_k, y = best_theta), color = "white", fill = "white") +
  xlab("k") + ylab("theta")
  
pdf(file = here::here("Plots/FigureDrafts", "disp_95CI_surface_and_CI_limits.pdf"))
disp_95CI_surface_plot
dev.off()

########## Plots for WSN talk #########
### LRP (LEP_R) with DD
pdf(file = here::here("Plots/WSN_2019", "LRP_with_DD.pdf"))
ggplot(data = output_uncert_all_DD$LEP_R_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = 'gray', fill = 'gray') +
  #geom_vline(data = LEP_R_best_est, aes(xintercept = LEP_R, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "mean offspring"))$LEP_R, color = 'black') +
  geom_vline(xintercept = LEP_R_best_est_DD, color = "black") +
  xlab('LRP') + ggtitle('Lifetime recruit production') +
  theme_bw() +
  theme(text=element_text(size=25)) 
dev.off()

pdf(file = here::here("Plots/WSN_2019", "NP_with_DD.pdf"))
ggplot(data = output_uncert_all_DD$NP_out_df, aes(x=value)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  #geom_vline(data = NP_best_est, aes(xintercept = NP, color = recruit_size)) +
  #geom_vline(xintercept = (NP_best_est %>% filter(recruit_size == "mean offspring"))$NP, color='black') +
  geom_vline(xintercept = NP_best_est_DD, color = "black") +
  xlab('lambda') + ggtitle('Network persistence') +
  theme_bw() +
  theme(text=element_text(size=25))
dev.off()

# How many >= 1 of the 1000 runs?
output_uncert_all_DD$NP_out_df %>% filter(value >= 1) %>% summarize(npersist = n())

pdf(file = here::here("Plots/WSN_2019", "NP_with_DD_all_offspring.pdf"))
ggplot(data = output_uncert_all_offspring_all_DD$NP_out_df, aes(x=value)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  #geom_vline(data = NP_best_est, aes(xintercept = NP, color = recruit_size)) +
  #geom_vline(xintercept = (NP_best_est %>% filter(recruit_size == "mean offspring"))$NP, color='black') +
  geom_vline(xintercept = best_est_metrics_mean_offspring_all_offspring_DD$NP, color= "black") +
  #geom_vline(xintercept = NP_best_est_aDD, color = "black") +
  xlab('NP') + ggtitle('Network persistence') +
  theme_bw() +
  theme(text=element_text(size=25))
dev.off()

pdf(file = here::here("Plots/WSN_2019", "realized_C_DD_plot.pdf"))
ggplot(data = best_est_metrics_mean_offspring_DD$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_R)) +
  scale_fill_gradient(high='black', low='white', name='recruits') +
  xlab('origin') + ylab('destination') + #ggtitle('Realized connectivity matrix with DD') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
dev.off()

pdf(file = here::here("Plots/WSN_2019", "SP_with_DD.pdf"))
ggplot(data = output_uncert_all_DD$SP_vals_with_params, aes(x=reorder(site, org_geo_order), y=SP)) +
  geom_violin(fill="grey") +
  geom_point(data = SP_best_est_DD, aes(x = site, y = SP_value), color = "black") +
  xlab("Site") + ylab("SP") + #ggtitle("Self-persistence with DD") +
  ylim(c(0,0.65)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
dev.off()

########## Plots for figures ##########

##### Demographic parameters, relationships: survival curve, growth curve, dispersal kernel, transition size to female (?)
distance_vec <- seq(from=0, to=50, by=0.01)
connectivity_est_vec <- disp_kernel_all_years(distance_vec, k_allyears, theta_allyears)  # theta = 1, equation for p(d) in eqn. 6c in Bode et al. 2018
connectivity_lowest_k <- disp_kernel_all_years(distance_vec, min(k_connectivity_values), theta_allyears)  # minimum k in the 95% CI
connectivity_highest_k <- disp_kernel_all_years(distance_vec, max(k_connectivity_values), theta_allyears)  # maximum k in the 95% CI
dispersal_df <- data.frame(distance = distance_vec, kernel_bestfit = connectivity_est_vec, kernel_CI1 = connectivity_lowest_k, kernel_CI2 = connectivity_highest_k) %>%
  mutate(low_bound = case_when(kernel_CI1 <= kernel_CI2 ~ kernel_CI1,
                               kernel_CI1 > kernel_CI2 ~ kernel_CI2),
         upper_bound = case_when(kernel_CI1 > kernel_CI1 ~ kernel_CI1,
                                 kernel_CI1 <= kernel_CI2 ~ kernel_CI2))

# Dispersal kernel
dispersal_kernel_plot <- ggplot(data=dispersal_df, aes(x=distance, y=kernel_bestfit, ymin=kernel_CI1, ymax=kernel_CI2)) +
  geom_line(color='black') +
  geom_ribbon(alpha=0.5, color='gray') +
  xlab('distance (km)') + ylab('dispersal probability') + ggtitle('Dispersal kernel') +
  theme_bw()

# Growth curve - plot the data and linear model for fish caught about a year apart, models for only one recapture pair per fish, pairs selected randomly and models fit 1000x
growth_curve_plot <- ggplot(data = recap_pairs_year, aes(x = L1, y = L2)) +
  geom_point(shape = 1) +
  geom_abline(aes(intercept = mean(growth_info_estimate$intercept_est), slope = mean(growth_info_estimate$slope_est)), color = "black") +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", lwd = 1.5) +  #  1:1 line
  geom_ribbon(aes(x=seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1)),
                  ymin = (min(growth_info_estimate$intercept_est) + min(growth_info_estimate$slope_est)*seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1))),
                  ymax = (max(growth_info_estimate$intercept_est) + max(growth_info_estimate$slope_est)*seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1)))), fill = "light gray", alpha = 0.5) +
  xlab("length (cm)") + ylab("length (cm) next year") + ggtitle('Growth') +
  theme_bw()

# 
# growth_df <- data.frame(length1 = seq(min_size, max_size, length.out = n_bins*10)) %>%
#   mutate(length2 = VBL_growth(Linf_growth_mean, k_growth_mean, length1))
# 
# growth_curve_plot <- ggplot(data=growth_df, aes(x=length1, y=length2)) +
#   geom_point(color='black') +
#   xlab('size (cm)') + ylab('size next year') + ggtitle('b) Growth curve') +
#   ylim(0,13) +
#   theme_bw()

# Breeding size distribution
breeding_size_plot <- ggplot(data = recap_first_female, aes(x=size)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept=breeding_size_mean, color='black') +
  xlab('size (cm)') + ggtitle('Female transition size') +
  theme_bw()

# Survival plot
min_size_plot = 0
size.values <- min_size_plot+(0:30)*(max_size-min_size_plot)/30

#eall_mean.Phi.size.p.size.plus.dist.results = as.data.frame(eall_mean.Phi.size.p.size.plus.dist$results$beta)
survival_output_to_plot <- data.frame(size = size.values) %>%
  mutate(Phi_logit = survival_output$estimate[Phi_int_pos] + survival_output$estimate[Phi_size_pos]*size,
         Phi_lcl_logit = survival_output$lcl[Phi_int_pos] + survival_output$lcl[Phi_size_pos]*size,
         Phi_ucl_logit = survival_output$ucl[Phi_int_pos] + survival_output$ucl[Phi_size_pos]*size,
         Phi = logit_recip(Phi_logit),
         Phi_lcl = logit_recip(Phi_lcl_logit),
         Phi_ucl = logit_recip(Phi_ucl_logit))

# Phibysize_Phi.size.p.size.plus.dist_means <- data.frame(size = size.values) %>%
#   mutate(Phi_logit = eall_mean.Phi.size.p.size.plus.dist.results$estimate[1] + eall_mean.Phi.size.p.size.plus.dist.results$estimate[2]*size,
#          Phi_lcl_logit = eall_mean.Phi.size.p.size.plus.dist.results$lcl[1] + eall_mean.Phi.size.p.size.plus.dist.results$lcl[2]*size,
#          Phi_ucl_logit = eall_mean.Phi.size.p.size.plus.dist.results$ucl[1] + eall_mean.Phi.size.p.size.plus.dist.results$ucl[2]*size,
#          Phi = logit_recip(Phi_logit),
#          Phi_lcl = logit_recip(Phi_lcl_logit),
#          Phi_ucl = logit_recip(Phi_ucl_logit))

survival_plot <- ggplot(data = survival_output_to_plot, aes(size, Phi)) +
  geom_ribbon(aes(ymin=Phi_lcl,ymax=Phi_ucl),color="gray",fill="gray") +
  geom_line(color="black") +
  xlab("size (cm)") + ylab("probability of survival") + ggtitle("Annual survival") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

pdf(file = here('Plots/FigureDrafts','Parameter_inputs.pdf'), width=6, height=6)
plot_grid(dispersal_kernel_plot, growth_curve_plot, survival_plot, breeding_size_plot, labels = c("a","b","c","d"), nrow=2)
dev.off()


##### Metrics: LEP, recruit-egg-survival, LEP_R
# Without accounting for DD
LEP_plot <- ggplot(data = output_uncert_all$LEP_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = 'gray', fill = 'gray') +
  #geom_vline(data = LEP_best_est, aes(xintercept = LEP, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_best_est %>% filter(recruit_size == "mean offspring"))$LEP, color='black') +
  geom_vline(xintercept = LEP_best_est, color = "black") +
  xlab('LEP') + ggtitle('LEP') +
  theme_bw()

RperE_plot <- ggplot(data = output_uncert_all$RperE_out_df, aes(x=value)) +
  geom_histogram(bins=50, color='gray', fill='gray') +
  geom_vline(xintercept = recruits_per_egg_best_est, color = "black") +
  xlab('recruits-per-egg') + ggtitle('Egg-recruit surv') +
  theme_bw()

LEP_R_plot <- ggplot(data = output_uncert_all$LEP_R_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = 'gray', fill = 'gray') +
  #geom_vline(data = LEP_R_best_est, aes(xintercept = LEP_R, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "mean offspring"))$LEP_R, color = 'black') +
  geom_vline(xintercept = LEP_R_best_est, color = "black") +
  xlab('LRP') + ggtitle('LRP') +
  theme_bw()

LEP_R_local_plot <- ggplot(data = output_uncert_all$LEP_R_local_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = "gray", fill = "gray") +
  geom_vline(xintercept = LEP_R_local_best_est, color = "black") +
  xlab("LRP_local") + ggtitle("LRP_local") +
  theme_bw()

# With accounting for DD
RperE_plot_DD <- ggplot(data = output_uncert_all_DD$RperE_out_df, aes(x=value)) +
  geom_histogram(bins=50, color='gray', fill='gray') +
  geom_vline(xintercept = recruits_per_egg_best_est_DD, color = "black") +
  xlab('recruits-per-egg') + ggtitle('Egg-recruit surv,DD') +
  theme_bw()

LEP_R_plot_DD <- ggplot(data = output_uncert_all_DD$LEP_R_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = 'gray', fill = 'gray') +
  #geom_vline(data = LEP_R_best_est, aes(xintercept = LEP_R, color = recruit_size)) +
  #geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "mean offspring"))$LEP_R, color = 'black') +
  geom_vline(xintercept = LEP_R_best_est_DD, color = "black") +
  xlab('LRP') + ggtitle('LRP with DD') +
  theme_bw()

LEP_R_local_plot_DD <- ggplot(data = output_uncert_all_DD$LEP_R_local_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = "gray", fill = "gray") +
  geom_vline(xintercept = LEP_R_local_best_est_DD, color = "black") +
  xlab("LRP_local") + ggtitle("LRP_local with DD") +
  theme_bw()

pdf(file = here('Plots/FigureDrafts', 'LEP_RperE_LRP.pdf'), width=8.5, height=6)
plot_grid(LEP_plot, RperE_plot, LEP_R_plot, LEP_R_local_plot,
          NULL, RperE_plot_DD, LEP_R_plot_DD, LEP_R_local_plot_DD,
          labels = c("a","b","c","d","","e","f","g"), nrow=2)
dev.off()

##### Metrics: NP, realized connectivity matrix -- GOT A WEIRD ERROR ABOUT SCALE FOR COMPLEX? MAYBE FROM CMAT? BUT CAME WHEN I PUT THE TWO PLOTS TOGETHER SO NOT SURE?
best_est_metrics_mean_offspring$Cmat$org_site <- replace(best_est_metrics_mean_offspring$Cmat$org_site, 
                                                         best_est_metrics_mean_offspring$Cmat$org_site=="Tamakin Dacot", 
                                                      "Tomakin Dako")
best_est_metrics_mean_offspring_DD$Cmat$org_site <- replace(best_est_metrics_mean_offspring_DD$Cmat$org_site, 
                                                         best_est_metrics_mean_offspring_DD$Cmat$org_site=="Tamakin Dacot", 
                                                         "Tomakin Dako")

# Without accounting for DD
NP_plot <- ggplot(data = output_uncert_all$NP_out_df, aes(x=value)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  #geom_vline(data = NP_best_est, aes(xintercept = NP, color = recruit_size)) +
  #geom_vline(xintercept = (NP_best_est %>% filter(recruit_size == "mean offspring"))$NP, color='black') +
  geom_vline(xintercept = NP_best_est, color = "black") +
  xlab('NP') + ggtitle('Network persistence') +
  theme_bw()

realized_C_plot <- ggplot(data = best_est_metrics_mean_offspring$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_R)) +
  scale_fill_gradient(high='dark gray', low='white', name='Recruits') +
  xlab('origin') + ylab('destination') + ggtitle('Realized connectivity matrix') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 

# With accounting for DD
NP_plot_DD <- ggplot(data = output_uncert_all_DD$NP_out_df, aes(x=value)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  #geom_vline(data = NP_best_est, aes(xintercept = NP, color = recruit_size)) +
  #geom_vline(xintercept = (NP_best_est %>% filter(recruit_size == "mean offspring"))$NP, color='black') +
  geom_vline(xintercept = NP_best_est_DD, color = "black") +
  xlab('NP') + ggtitle('Network persistence with DD') +
  theme_bw()

realized_C_plot_DD <- ggplot(data = best_est_metrics_mean_offspring_DD$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_R)) +
  scale_fill_gradient(high='black', low='white', name='Recruits') +
  xlab('origin') + ylab('destination') + ggtitle('Realized connectivity matrix with DD') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 

pdf(file = here('Plots/FigureDrafts', 'NP_and_connMatrixR.pdf'), width = 10, height = 10)  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(NP_plot, realized_C_plot, NP_plot_DD, realized_C_plot_DD, labels = c("a","b","c","d"), nrow=2)
dev.off()

##### Metrics: Self persistence by site
output_uncert_all$SP_vals_with_params$site <- replace(output_uncert_all$SP_vals_with_params$site, 
                                                      output_uncert_all$SP_vals_with_params$site=="Tamakin Dacot", 
                                                      "Tomakin Dako")
output_uncert_all_DD$SP_vals_with_params$site <- replace(output_uncert_all_DD$SP_vals_with_params$site, 
                                                      output_uncert_all_DD$SP_vals_with_params$site=="Tamakin Dacot", 
                                                      "Tomakin Dako")
sites_for_total_areas_TD <- sites_for_total_areas
sites_for_total_areas_TD[13] <- "Tomakin Dako"

# Without accounting for DD
SP_plot <- ggplot(data = output_uncert_all$SP_vals_with_params %>% filter(site %in% sites_for_total_areas_TD), aes(x=reorder(site, org_geo_order), y=SP)) +
  geom_violin(fill="grey") +
  geom_point(data = SP_best_est %>% filter(site %in% sites_for_total_areas_TD), aes(x = site, y = SP_value), color = "black") +
  xlab("Site") + ylab("SP") + ggtitle("Self-persistence") +
  ylim(c(0,0.65)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

# Accounting for DD
SP_plot_DD <- ggplot(data = output_uncert_all_DD$SP_vals_with_params %>% filter(site %in% sites_for_total_areas_TD), aes(x=reorder(site, org_geo_order), y=SP)) +
    geom_violin(fill="grey") +
    geom_point(data = SP_best_est_DD %>% filter(site %in% sites_for_total_areas_TD), aes(x = site, y = SP_value), color = "black") +
    xlab("Site") + ylab("SP") + ggtitle("Self-persistence with DD") +
    ylim(c(0,0.65)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

pdf(file = here::here("Plots/FigureDrafts", "SP_hists_by_site_noSLSTCP.pdf"))
plot_grid(SP_plot, SP_plot_DD, labels = c("a","b"), nrow = 1)
dev.off()

# ggplot(data = output_uncert_all$SP_vals_with_params %>% filter(site %in% sites_for_total_areas), aes(x=SP)) +
#   geom_histogram(binwidth=0.001, color='gray', fill='gray') +
#   #geom_histogram(binwidth=0.0005, color='gray', fill='gray') +
#   #geom_vline(data=(SP_best_est %>% filter(recruit_size == "mean offspring") %>% filter(site %in% sites_for_total_areas)), 
#   #           aes(xintercept=SP_value), color='black') +   
#   geom_vline(data = SP_best_est, aes(xintercept = SP_value), color = "black") +
#   facet_wrap(~reorder(site, org_geo_order)) +
#   xlab('SP') + ggtitle('Self-persistence by site') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

########## Appendix plots: ##########

##### Inputs 
# Start-recruit size
startRecruit_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=start_recruit_size)) +
  geom_histogram(bins=40, color="gray", fill="gray") +
  geom_vline(xintercept = (start_recruit_size_options %>% filter(recruit_size == "mean offspring"))$size) +
  xlab("recruit size (cm)") + ggtitle("Census size") +
  theme_bw()

# Breeding size
breedingSize_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=breeding_size)) +
  geom_histogram(bins=40, color="gray", fill="gray") +
  geom_vline(xintercept=breeding_size_mean, color='black') +
  xlab("female breeding size (cm)") + ggtitle('Female transition') +
  theme_bw()

# Dispersal k
dispersalK_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=k_connectivity)) +
  geom_histogram(bins=40, color="gray", fill="gray") +
  geom_vline(xintercept = k_allyears, color='black') +
  xlab("k parameter") + ggtitle("Dispersal kernel") +
  theme_bw()

# Growth - Linf + k
growthLinf_k_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=Linf, y=k_growth)) +
  geom_point(color = "gray", fill = "gray") +
  geom_point(x = Linf_growth_mean, y = k_growth_mean, color = "black", fill = "black") +
  xlab('Linf (cm)') + ylab("k") + ggtitle("VBL growth model") +
  theme_bw()

# Survival
# survivalSint_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=Sint)) +
#   geom_histogram(bins=40, color='gray', fill='gray') +
#   geom_vline(xintercept = Sint_mean, color = 'black') +
#   xlab('intercept of size-survival relationship') + ggtitle('Survival') +
#   theme_bw()

run_list <- rep(NA, n_runs*length(size.values))
for(i in 1:n_runs) {
  run_list_start = (length(size.values)*(i-1))+1
  run_list_end = length(size.values)*i
  run_list[run_list_start:run_list_end] = i
}
surv_input_plot_df <- data.frame(size = rep(size.values, n_runs),
                                 run = run_list,
                                 value = NA)
for(i in 1:n_runs) {
  out_vec = logit_recip(Sint_set[i] + size.values*Sl_set[i])
  surv_input_plot_df$value[((length(size.values)*(i-1))+1):(length(size.values)*i)] = out_vec
}

survivalSintandSl_plot <- ggplot(data = surv_input_plot_df, aes(x=size, y=value, group=run)) +
  geom_line(alpha = 0.1) +
  xlab("size (cm)") + ylab("prob survival") + ggtitle("Survival") +
  theme_bw()

# Prob r
probR_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x = prob_r)) +
  geom_histogram(bins = 40, color = "gray", fill = "gray") +
  xlab("P_c") + ggtitle("Capture probability") +
  theme_bw()

# Assigned offspring
assignedOffspring_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x = assigned_offspring)) +
  geom_histogram(bins = 40, color = "gray", fill = "gray") +
  geom_vline(xintercept = n_offspring_matched, color = "black") +
  xlab("# assigned offspring") + ggtitle("Assigned offspring") +
  theme_bw()

# Proportion habitat 
# P_h = cumulative proportion sites sampled across time
# P_d = proportion of the dispersal kernel from each site covered by our sampling region
# P_s = proportion of our sampling region that is habitat
# P_DD = proportion habitat available to juveniles if exclude DD

prop_scaling_vals_for_plot <- data.frame(value = c("P_h", "P_d", "P_s", "DD"),
                                         estimate = c(total_prop_hab_sampled_through_time, prop_total_disp_area_sampled_best_est,
                                                      prop_sampling_area_habitat, (perc_APCL_val+perc_UNOC_val)/perc_UNOC_val))
# prop_hab_vals_for_plot <- data.frame(value = c("prop sampled all \n sites across time", "dispersal kernel \n area", "additional habitat \n without DD"),
#                                      estimate = c(total_prop_hab_sampled_through_time, prop_total_disp_area_sampled_best_est, (perc_APCL_val+perc_UNOC_val)/perc_UNOC_val))
propHabitat_plot <- ggplot(data = prop_scaling_vals_for_plot, aes(x = value, y = estimate)) +
  geom_point() +
  xlab("") + ylab("proportion") + ggtitle("Habitat scaling") +
  ylim(c(0,1.8)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##### Histograms (or scatters) of all inputs into uncertainty runs
pdf(file = here::here("Plots/FigureDrafts", "Uncertainty_inputs.pdf"))
plot_grid(startRecruit_plot, breedingSize_plot, dispersalK_plot, 
             growthLinf_k_plot, survivalSintandSl_plot , probR_plot, 
             assignedOffspring_plot, propHabitat_plot, labels = c("a","b","c","d","e","f","g","h"),  nrow=3)
dev.off()

##### Uncertainty in LEP
pdf(file = here::here("Plots/FigureDrafts", "LEP_uncertainty_breakdown.pdf"))
ggplot(data = LEP_uncert %>% filter(uncertainty_type %in% c("start recruit size", "breeding size", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("start recruit size", "breeding size", "growth", "survival", "all"),
                               value = LEP_best_est), color = "black") +
  xlab("LEP") + ggtitle("Uncertainty in LEP") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
dev.off()

# pdf(file = here::here("Plots/FigureDrafts", "LEP_uncertainty_breakdown.pdf"))
# ggplot(data = LEP_uncert %>% filter(uncertainty_type %in% c("start recruit size", "breeding size", "growth", "survival", "all")), aes(x = value)) +
#   geom_histogram(bins = 50, color = "gray", fill = "gray") +
#   #geom_vline(data = LEP_best_est %>% filter(recruit_size == "mean offspring"), aes(xintercept = LEP)) +
#   geom_vline(xintercept = LEP_best_est, color = "black") +
#   facet_wrap(~uncertainty_type) +
#   xlab('LEP') + ggtitle('Uncertainty in LEP') +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
# dev.off()

##### Uncertainty in LEP_R
pdf(file = here::here("Plots/FigureDrafts", "LEP_R_uncertainty_breakdown.pdf"))
ggplot(data = LEP_R_uncert %>% filter(uncertainty_type %in% c("start recruit size", "breeding size", "assigned offspring", "prob r", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("start recruit size", "breeding size", "assigned offspring", "prob r", "growth", "survival", "all"),
                               value = LEP_R_best_est), color = "black") +
  xlab("LRP") + ggtitle("Uncertainty in LRP") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
dev.off()

# LEP_R_violin_df <- rbind((LEP_R_uncert %>% 
#   filter(uncertainty_type %in% c("start recruit size", "breeding size", "assigned offspring", "prob r", "growth", "survival", "all")) %>%
#   mutate(DD_type = "DD unaccounted for")),
#   (LEP_R_uncert_DD %>%
#   filter(uncertainty_type %in% c("start recruit size", "breeding size", "assigned offspring", "prob r", "growth", "survival", "all")) %>%
#     mutate(DD_type = "DD accounted for")))
# LEP_R_violin_best_est_df <- data.frame(uncertainty_type = rep(c("start recruit size", "breeding size", "assigned offspring", "prob r", "growth", "survival", "all"),2),
#                                        value = c(rep(LEP_R_best_est, 7), rep(LEP_R_best_est_DD, 7)),
#                                        DD_type = c(rep("DD unaccounted for", 7), rep("DD accounted for", 7)))
# 
# pdf(file = here::here("Plots/FigureDrafts", "LEP_R_uncertainty_breakdown.pdf"))
# ggplot(data = LEP_R_violin_df,  aes(x=uncertainty_type, y=value, fill=DD_type, color = DD_type)) +
#   geom_violin() +
#   geom_point(data = LEP_R_violin_best_est_df, color = "black") +
#   scale_fill_manual(values = c("gray", "dark blue")) +
#   scale_color_manual(values = c("grey", "dark blue")) +
#   # geom_point(data = data.frame(uncertainty_type = c("start recruit size", "breeding size", "assigned offspring", "prob r", "growth", "survival", "all"),
#   #                              value = LEP_R_best_est), color = "black") +
#   xlab("LRP") + ggtitle("Uncertainty in LRP") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
# dev.off()

# pdf(file = here::here("Plots/FigureDrafts", "LEP_R_uncertainty_breakdown.pdf"))
# ggplot(data = LEP_R_uncert %>% filter(uncertainty_type %in% c("start recruit size", "breeding size", "assigned offspring", "prob r", "growth", "survival", "all")), 
#        aes(x = value)) +
#   geom_histogram(bins = 50, color = "gray", fill = "gray") +
#   #geom_vline(data = LEP_R_best_est %>% filter(recruit_size == "mean offspring"), aes(xintercept = LEP_R)) +
#   geom_vline(xintercept = LEP_R_best_est, color = "black") +
#   facet_wrap(~uncertainty_type) +
#   xlab('LRP') + ggtitle('Uncertainty in LRP') +
#   theme_bw()
# dev.off()

##### Uncertainty in LEP_R with DD
pdf(file = here::here("Plots/FigureDrafts", "LEP_R_uncertainty_breakdown_with_DD.pdf"))
ggplot(data = LEP_R_uncert_DD %>% filter(uncertainty_type %in% c("start recruit size", "breeding size", "assigned offspring", "prob r", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("start recruit size", "breeding size", "assigned offspring", "prob r", "growth", "survival", "all"),
                               value = LEP_R_best_est_DD), color = "black") +
  xlab("LRP") + ggtitle("Uncertainty in LRP accounting for DD") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
dev.off()

ggplot(data = LEP_R_uncert %>% filter(uncertainty_type == "assigned offspring"), aes(x=value)) +
  geom_histogram(bins = 40)


##### Uncertainty in RperE
# Without DD accounted for 
pdf(file = here::here("Plots/FigureDrafts", "RperE_uncertainty_breakdown.pdf"))
ggplot(data = RperE_uncert %>% filter(uncertainty_type %in% c("breeding size", "assigned offspring", "prob r", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("breeding size", "assigned offspring", "prob r", "growth", "survival", "all"),
                               value = recruits_per_egg_best_est), color = "black") +
  xlab("recruits-per-egg") + ggtitle("Uncertainty in egg-recruit survival") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
dev.off()

# With DD accounted for 
pdf(file = here::here("Plots/FigureDrafts", "RperE_uncertainty_breakdown_DD.pdf"))
ggplot(data = RperE_uncert_DD %>% filter(uncertainty_type %in% c("breeding size", "assigned offspring", "prob r", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("breeding size", "assigned offspring", "prob r", "growth", "survival", "all"),
                               value = recruits_per_egg_best_est_DD), color = "black") +
  xlab("recruits-per-egg") + ggtitle("Uncertainty in egg-recruit survival /n accounting for DD") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
dev.off()

# pdf(file = here::here("Plots/FigureDrafts", "RperE_uncertainty_breakdown.pdf"))
# ggplot(data = RperE_uncert %>% filter(uncertainty_type %in% c("breeding size", "assigned offspring", "prob r", "growth", "survival", "all")), 
#        aes(x = value)) +
#   geom_histogram(bins = 50, color = "gray", fill = "gray") +
#   geom_vline(xintercept = recruits_per_egg_best_est, color='black') +
#   facet_wrap(~uncertainty_type) +
#   xlab('recruits-per-egg') + ggtitle('Uncertainty in egg-recruit survival') +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
# dev.off()

##### Uncertainty in NP
# Without DD accounted for 
pdf(file = here::here("Plots/FigureDrafts", "NP_uncertainty_breakdown.pdf"))
ggplot(data = NP_uncert %>% filter(uncertainty_type != "assigned offspring and prob r"), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("breeding size", "assigned offspring", "prob r", "growth", "survival", "all", "dispersal k", "start recruit size"),
                               value = NP_best_est), color = "black") +
  xlab("NP") + ggtitle("Uncertainty in network persistence") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
dev.off()

# With DD accounted for 
pdf(file = here::here("Plots/FigureDrafts", "NP_uncertainty_breakdown_DD.pdf"))
ggplot(data = NP_uncert_DD %>% filter(uncertainty_type != "assigned offspring and prob r"), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("breeding size", "assigned offspring", "prob r", "growth", "survival", "all", "dispersal k", "start recruit size"),
                               value = NP_best_est_DD), color = "black") +
  xlab("NP") + ggtitle("Uncertainty in network persistence /n accounting for DD") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
dev.off()

# pdf(file = here::here("Plots/FigureDrafts", "NP_uncertainty_breakdown.pdf"))
# ggplot(data = NP_uncert %>% filter(uncertainty_type != "assigned offspring and prob r"), aes(x = value)) +
#   geom_histogram(bins = 50, color = "gray", fill = "gray") +
#   #geom_vline(data = NP_best_est %>% filter(recruit_size == "mean offspring"), aes(xintercept = NP)) +
#   geom_vline(xintercept = NP_best_est, color = "black") +
#   facet_wrap(~uncertainty_type) +
#   xlab('NP') + ggtitle('Uncertainty in network persistence') +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
# dev.off()

##### Proportion of total kernel area from each site covered by our sampling
all_parents_site$site <- replace(all_parents_site$site, all_parents_site$site=="Tamakin Dacot", "Tomakin Dako")

pdf(file = here('Plots/FigureDrafts', 'Prop_of_kernel_area_sampled_by_site.pdf'))
ggplot(data = all_parents_site, aes(x = reorder(site, site_geo_order), y = prop_disp_area_within_sites)) + # the geo orders are all off here...
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  xlab("site") + ylab("proportion kernel within sampled area") + ggtitle("Proportion of kernel area sampled") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### Effects of size and distance on recap prob
min_dist_plot = 0
max_dist_plot = 500
dist_values = seq(from = min_dist_plot, to = max_dist_plot, by = 1)

# Get size effect estimates ready
recap_size_output_to_plot <- data.frame(size = size.values) %>%
  mutate(p_logit = survival_output$estimate[p_int_pos] + survival_output$estimate[p_size_pos]*size,
         p_lcl_logit = survival_output$lcl[p_int_pos] + survival_output$lcl[p_size_pos]*size,
         p_ucl_logit = survival_output$ucl[p_int_pos] + survival_output$ucl[p_size_pos]*size,
         p = logit_recip(p_logit),
         p_lcl = logit_recip(p_lcl_logit),
         p_ucl = logit_recip(p_ucl_logit))

# And plot
recap_size_effect_plot <- ggplot(data = recap_size_output_to_plot, aes(size, p)) +
  geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="gray",fill="gray") +
  geom_line(color="black") +
  xlab("size (cm)") + ylab("probability of recapture") + ggtitle("a) Size effect") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

# Get distance effect estiamtes ready
recap_dist_output_to_plot <- data.frame(dist = dist_values) %>%
  mutate(p_logit = survival_output$estimate[p_int_pos] + survival_output$estimate[p_dist_pos]*dist,
         p_lcl_logit = survival_output$lcl[p_int_pos] + survival_output$lcl[p_dist_pos]*dist,
         p_ucl_logit = survival_output$ucl[p_int_pos] + survival_output$ucl[p_dist_pos]*dist,
         p = logit_recip(p_logit),
         p_lcl = logit_recip(p_lcl_logit),
         p_ucl = logit_recip(p_ucl_logit))

# And plot
recap_dist_effect_plot <- ggplot(data = recap_dist_output_to_plot, aes(dist, p)) +
  geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="gray",fill="gray") +
  geom_line(color="black") +
  xlab("distance (m)") + ylab("probability of recapture") + ggtitle("b) Distance effect") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

# Put together and save
pdf(file = here::here("Plots/FigureDrafts", "Recapture_size_distance_effects.pdf"))
grid.arrange(recap_size_effect_plot, recap_dist_effect_plot, nrow=1)
dev.off()

##### Relationships between values (taken from ########## Relationships between values - prettier plots to put in manuscript ########## below)
# Size at recruitment and LEP
recruit_size_vs_LEP_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=start_recruit_size, y=LEP)) +
  geom_point(size=2) +
  xlab("recruit size") + ylab("LEP") + ggtitle("Recruit size vs. LEP") +
  theme_bw()

# Size at female transition and egg-recruit-surv
breeding_size_vs_RperE_df <- data.frame(breeding_size = output_uncert_all$metric_vals_with_params$breeding_size,
                                        recruits_per_egg = output_uncert_all$RperE_out_df$value) 
breeding_size_vs_RperE_plot <- ggplot(data = breeding_size_vs_RperE_df, aes(x=breeding_size, y=recruits_per_egg)) +
  geom_point(size=2) +
  xlab("breeding size (cm)") + ylab("recruits per egg") + ggtitle("Female size vs. \n egg-recruit survival") +
  theme_bw()

# LEP_R and NP 
LRP_vs_NP_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=LEP_R, y=NP)) +
  geom_point(size=2) +
  xlab('LRP') + ylab('NP') + ggtitle('LRP vs. NP') +
  theme_bw()

# Prob r and NP
prob_r_vs_RperE <- data.frame(prob_r = output_uncert_all$metric_vals_with_params$prob_r,
                              recruits_per_egg = output_uncert_all$RperE_out_df$value,
                              NP = output_uncert_all$NP_out_df$value)
Prob_r_vs_NP_plot <- ggplot(data = prob_r_vs_RperE, aes(x = prob_r, y=NP)) +
  geom_point(size=2) +
  xlab('capture prob (P_c)') + ylab('NP') + ggtitle('Capture prob. vs. NP') +
  theme_bw()

# Survival intercept and LEP
Sint_vs_LEP_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=logit_recip(Sint), y=LEP)) +
  geom_point(size=2) +
  xlab("annual survival intercept") + ylab("LEP") + ggtitle("Survival vs. LEP") +
  theme_bw()

# Start recruit size and LRP
Start_recruit_vs_LRP <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=start_recruit_size, y=LEP_R)) +
  geom_point(size = 2) +
  xlab("size of recruit (cm)") + ylab("LRP") + ggtitle("Recruit size vs. LRP") +
  theme_bw()

pdf(file = here::here("Plots/FigureDrafts", "Param_metric_relationships.pdf"))
plot_grid(breeding_size_vs_RperE_plot, LRP_vs_NP_plot, 
          Prob_r_vs_NP_plot, Sint_vs_LEP_plot, Start_recruit_vs_LRP, 
          labels = c("a","b","c","d","e"), nrow=2)
dev.off()


########### What-if calcs ##########
########## What-if #1: what if all the offspring we genotype were from our sites? ##########

# LEP_R (LEP in terms of recruits) (this is plot "LEP_R_histogram_whatif_all_offspring.pdf") - where are these separate plots? Find and set to whatifs folder
LEP_R_histogram_whatif_all_offspring <- ggplot(data = output_uncert_all_offspring_all$LEP_R_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = 'gray', fill = 'gray') +
  geom_vline(xintercept = best_est_metrics_mean_offspring_all_offspring$LEP_R) +
  xlab("recruits") + ggtitle("Recruits per recruit \n with all offspring") +
  theme_bw()

# LEP_R accounting for DD
LEP_R_histogram_whatif_all_offspring_DD <- ggplot(data = output_uncert_all_offspring_all_DD$LEP_R_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = 'gray', fill = 'gray') +
  geom_vline(xintercept = best_est_metrics_mean_offspring_all_offspring_DD$LEP_R) +
  xlab("recruits") + ggtitle("Recruits per recruit \n with all offspring \n accounting for DD") +
  theme_bw()

# NP (this plot is "NP_histogram_whatif_all_offspring.pdf")
NP_histogram_whatif_all_offspring <- ggplot(data = output_uncert_all_offspring_all$NP_out_df, aes(x=value)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = best_est_metrics_mean_offspring_all_offspring$NP) +
  xlab("NP") + ggtitle("NP with all offspring") +
  theme_bw()

# NP accounting for DD
NP_histogram_whatif_all_offspring_DD <- ggplot(data = output_uncert_all_offspring_all_DD$NP_out_df, aes(x=value)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = best_est_metrics_mean_offspring_all_offspring_DD$NP) +
  xlab("NP") + ggtitle("NP with all offspring \n accounting for DD") +
  theme_bw()

pdf(file = here::here("Plots/PersistenceMetrics/Whatifs", "LEP_R_and_NP_histograms_whatif_all_offspring.pdf"), width=7, height=6)
plot_grid(LEP_R_histogram_whatif_all_offspring, NP_histogram_whatif_all_offspring, LEP_R_histogram_whatif_all_offspring_DD, NP_histogram_whatif_all_offspring_DD,
          labels = c("a","b","c","d"), nrow=2)
dev.off()

# ########## STOPPED RE-ORDERING, MOVING FIGURE AND APPENDIX PLOTS UP, HERE! ##########
# 
# ########## Plotting metrics ##########
# 
# ##### Plot the histograms of LEP, LEP_R, recruits_per_egg, and NP output
# # LEP
# pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_histogram.pdf'))
# ggplot(data = output_uncert_all$LEP_out_df, aes(x=value)) +
#   geom_histogram(bins=40, color = 'gray', fill = 'gray') +
#   geom_vline(data = LEP_best_est, aes(xintercept = LEP, color = recruit_size)) +
#   #geom_vline(xintercept = (LEP_best_est %>% filter(recruit_size == "3.5cm"))$LEP, color='black') +
#   xlab('LEP') + ggtitle('Histogram of LEP values') +
#   theme_bw()
# dev.off()
# 
# # LEP_R (LEP in terms of recruits) 
# pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_R_histogram.pdf'))
# ggplot(data = output_uncert_all$LEP_R_out_df, aes(x=value)) +
#   geom_histogram(bins=40, color = 'gray', fill = 'gray') +
#   geom_vline(data = LEP_R_best_est, aes(xintercept = LEP_R, color = recruit_size)) +
#   #geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "3.5cm"))$LEP_R, color = 'black') +
#   xlab('LEP_R') + ggtitle('Histogram of LEP_R values') +
#   theme_bw()
# dev.off()
# 
# # # Zoomed in LEP_R (since sometimes there are some really high values)
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_R_histogram_zoomed.pdf'))
# # ggplot(data = LEP_R_out_df, aes(x=value)) +
# #   geom_histogram(binwidth=1, color = 'gray', fill = 'gray') +
# #   geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "3.5cm"))$LEP_R, color = 'black') +
# #   xlim(c(0,100)) +
# #   xlab('LEP_R') + ggtitle('Histogram of LEP_R values, zoomed') +
# #   theme_bw()
# # dev.off()
# 
# # NP
# pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'NP_histogram.pdf'))
# ggplot(data = output_uncert_all$NP_out_df, aes(x=value)) +
#   geom_histogram(bins=40, color='gray', fill='gray') +
#   geom_vline(data = NP_best_est, aes(xintercept = NP, color = recruit_size)) +
#   #geom_vline(xintercept = (NP_best_est %>% filter(recruit_size == "3.5cm"))$NP, color='black') +
#   xlab('NP') + ggtitle('Histogram of NP values') +
#   theme_bw()
# dev.off()
# 
# # Recruits-per-egg
# pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'RperE_histogram.pdf'))
# ggplot(data = output_uncert_all$RperE_out_df, aes(x=value)) +
#   geom_histogram(bins=50, color='gray', fill='gray') +
#   geom_vline(xintercept = recruits_per_egg_best_est, color = "black") +
#   #geom_vline(xintercept = (NP_best_est %>% filter(recruit_size == "3.5cm"))$NP, color='black') +
#   xlab('recruits-per-egg') + ggtitle('Histogram of recruits-per-egg values') +
#   theme_bw()
# dev.off()
# 
# # Recruits-per-egg zoomed
# pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'RperE_histogram_zoomed.pdf'))
# ggplot(data = output_uncert_all$RperE_out_df, aes(x=value)) +
#   geom_histogram(bins=50, color='gray', fill='gray') +
#   geom_vline(xintercept = recruits_per_egg_best_est, color='black') +
#   xlim(0,0.002) +
#   xlab('recruits-per-egg') + ggtitle('Histogram of recruits-per-egg values zoomed') +
#   theme_bw()
# dev.off()
# 
# ##### SP at each site - this plot takes a few seconds to make
# pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty','SP_histogram.pdf'))
# #ggplot(data = SP_out_df, aes(x=value)) +
# ggplot(data = output_uncert_all$SP_vals_with_params, aes(x=SP)) +
#   #geom_histogram(binwidth=0.0005) +
#   geom_histogram(binwidth=0.001, color='gray', fill='gray') +
#   geom_vline(data=(SP_best_est %>% filter(recruit_size == "mean offspring")), aes(xintercept=SP_value), color='black') +   # now these look weird, don't seem to fit in the dists for most sites? -- not anymore...
#   #ylim(0,300) +
#   facet_wrap(~site) +
#   xlab('SP') + ggtitle('Self-persistence histograms by site') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
# dev.off()
# 
# # SP at each site - zoomed in
# pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty','SP_histogram_zoomed.pdf'))
# #ggplot(data = SP_out_df, aes(x=value)) +
# ggplot(data = output_uncert_all$SP_vals_with_params, aes(x=SP)) +
#   #geom_histogram(binwidth=0.0005) +
#   geom_histogram(binwidth=0.001, color='gray', fill='gray') +
#   geom_vline(data=(SP_best_est %>% filter(recruit_size == "mean offspring")), aes(xintercept=SP_value), color='black') +   # now these look weird, don't seem to fit in the dists for most sites? -- not anymore...
#   xlim(0,0.05) +
#   ylim(0,300) +
#   facet_wrap(~site) +
#   xlab('SP') + ggtitle('Self-persistence histograms by site') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
# dev.off()
# 
# # ########## Plotting the metrics with different kinds of uncertainty ########## STOPPED SAVING PLOTS HERE
# # 
# # # LEP
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "LEP_uncertainty.pdf"))
# # ggplot(data = LEP_uncert, aes(x = value)) +
# #   geom_histogram(bins = 50, color = "gray", fill = "gray") +
# #   geom_vline(data = LEP_best_est %>% filter(recruit_size == "mean offspring"), aes(xintercept = LEP)) +
# #   facet_wrap(~uncertainty_type) +
# #   xlab('LEP') + ggtitle('LEP with different types of uncertainty') +
# #   theme_bw()
# # dev.off()
# # 
# # # LEP zoomed 
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "LEP_uncertainty_zoomed.pdf"))
# # ggplot(data = LEP_uncert, aes(x = value)) +
# #   geom_histogram(bins = 50, color = "gray", fill = "gray") +
# #   geom_vline(data = LEP_best_est %>% filter(recruit_size == "mean offspring"), aes(xintercept = LEP)) +
# #   ylim(0,500) +
# #   facet_wrap(~uncertainty_type) +
# #   xlab('LEP') + ggtitle('LEP with different types of uncertainty') +
# #   theme_bw()
# # dev.off()
# # 
# # # LEP_R
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "LEP_R_uncertainty.pdf"))
# # ggplot(data = LEP_R_uncert, aes(x=value)) +
# #   geom_histogram(bins = 40, color = "gray", fill = "gray") +
# #   geom_vline(data = LEP_R_best_est %>% filter(recruit_size == "mean offspring"), aes(xintercept = LEP_R)) +
# #   facet_wrap(~uncertainty_type) +
# #   xlab("LEP_R") + ggtitle("LEP_R with different types of uncertainty") +
# #   theme_bw()
# # dev.off()
# # 
# # # LEP_R zoomed
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "LEP_R_uncertainty_zoomed.pdf"))
# # ggplot(data = LEP_R_uncert, aes(x=value)) +
# #   geom_histogram(bins = 40, color = "gray", fill = "gray") +
# #   geom_vline(data = LEP_R_best_est %>% filter(recruit_size == "mean offspring"), aes(xintercept = LEP_R)) +
# #   ylim(0,500) +
# #   facet_wrap(~uncertainty_type) +
# #   xlab("LEP_R") + ggtitle("LEP_R with different types of uncertainty") +
# #   theme_bw()
# # dev.off()
# # 
# # # NP
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "NP_uncertainty.pdf"))
# # ggplot(data = NP_uncert, aes(x=value)) +
# #   geom_histogram(bins = 40, color = "gray", fill = "gray") +
# #   geom_vline(data = NP_best_est %>% filter(recruit_size == "mean offspring"), aes(xintercept =NP)) +
# #   facet_wrap(~uncertainty_type) +
# #   xlab("NP") + ggtitle("NP with different types of uncertainty") +
# #   theme_bw()
# # dev.off()
# # 
# # # NP zoomed
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "NP_uncertainty_zoomed.pdf"))
# # ggplot(data = NP_uncert, aes(x=value)) +
# #   geom_histogram(bins = 40, color = "gray", fill = "gray") +
# #   geom_vline(data = NP_best_est %>% filter(recruit_size == "mean offspring"), aes(xintercept =NP)) +
# #   ylim(0,750) +
# #   facet_wrap(~uncertainty_type) +
# #   xlab("NP") + ggtitle("NP with different types of uncertainty") +
# #   theme_bw()
# # dev.off()
# # 
# # # Recruits-per-egg - why is all different than recruits-per-egg here? Was on old, fixed now?
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "RperE_uncertainty.pdf"))
# # ggplot(data = RperE_uncert %>% filter(uncertainty_type %in% c("all", "assigned offspring", "assigned offspring and prob r", "breeding size", "growth", "survival")), aes(x=value)) +
# #   geom_histogram(bins = 50, color = "gray", fill = "gray") +
# #   geom_vline(xintercept = recruits_per_egg_best_est) +
# #   facet_wrap(~uncertainty_type) +
# #   xlab("recruits-per-egg") + ggtitle("RperE with different types of uncertainty") +
# #   theme_bw() +
# #   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
# # dev.off()
# # 
# # # Recruits-per-egg  - zoomed
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "RperE_uncertainty_zoomed.pdf"))
# # ggplot(data = RperE_uncert %>% filter(uncertainty_type %in% c("all", "assigned offspring", "assigned offspring and prob r", "breeding size", "growth", "survival")), aes(x=value)) +
# #   geom_histogram(bins = 40, color = "gray", fill = "gray") +
# #   geom_vline(xintercept = recruits_per_egg_best_est) +
# #   xlim(0,0.0002) +
# #   facet_wrap(~uncertainty_type) +
# #   xlab("recruits-per-egg") + ggtitle("RperE with different types of uncertainty") +
# #   theme_bw() +
# #   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
# # dev.off()
# 
# ########## Plotting inputs (histograms of data inputs) ########## 
# 
# # start recruit size (what size we think a "recruit" is)
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'start_recruit_histogram.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=start_recruit_size)) +
#   geom_histogram(bins=40, color='gray', fill='gray') +
#   geom_vline(data = start_recruit_size_options, aes(xintercept = size, color = recruit_size)) +
#   xlab('start recruit size (cm)') + ggtitle('start recruit size values') +
#   theme_bw()
# dev.off()
# 
# # k_connectivity
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'k_connectivity_histogram.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=k_connectivity)) +
#   geom_histogram(bins=40, color='gray', fill='gray') +
#   geom_vline(xintercept = k_allyears, color='black') +
#   xlab('k_connectivity') + ggtitle('k_connectivity values') +
#   theme_bw()
# dev.off()
# 
# # Linf
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Linf_histogram.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=Linf)) +
#   geom_histogram(bins=40, color='gray', fill='gray') +
#   geom_vline(xintercept = Linf_growth_mean, color='black') +
#   xlab('Linf') + ggtitle('Linf values') +
#   theme_bw()
# dev.off()
# 
# # k growth
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'k_growth_histogram.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=k_growth)) +
#   geom_histogram(bins=40, color='gray', fill='gray') +
#   geom_vline(xintercept = k_growth_mean, color='black') +
#   xlab('k in VBL') + ggtitle('k (growth) values') +
#   theme_bw()
# dev.off()
# 
# # Sint -- #come up with a more descriptive way of titling this!
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Sint_histogram.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=Sint)) +
#   geom_histogram(bins=40, color='gray', fill='gray') +
#   geom_vline(xintercept = Sint_mean, color = 'black') +
#   xlab('Sint') + ggtitle('Sint (survival) values') +
#   theme_bw()
# dev.off()
# 
# # Breeding size
# pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Breeding_size_histogram.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=breeding_size)) +
#   geom_histogram(bins=40, color='gray', fill='gray') +
#   geom_vline(xintercept=breeding_size_mean, color='black') +
#   xlab('Breeding size') + ggtitle('Breeding size values') +
#   theme_bw()
# dev.off()
# 
# # Prob r (prob of catching a fish)
# pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "Prob_r_histogram.pdf"))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x = prob_r)) +
#   geom_histogram(bins = 40, color = "gray", fill = "gray") +
#   xlab("prob r") + ggtitle("Prob r values") +
#   theme_bw()
# dev.off()
# 
# # Number of offspring assigned to parents (using rbinom)
# pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "Assigned_offspring_histogram.pdf"))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x = assigned_offspring)) +
#   geom_histogram(bins = 40, color = "gray", fill = "gray") +
#   geom_vline(xintercept = n_offspring_parentage, color = "black") +
#   xlab("# assigned offspring") + ggtitle("Number offspring assigned to parents") +
#   theme_bw()
# dev.off()
# 
# 
# ########## Relationships between values ##########
# 
# # LEP_R and NP - do we expect these to be related perfectly linearly? I guess, when both connectivity and recruits-per-egg are static...
# pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_R_NP_scatter.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=LEP_R, y=NP)) +
#   geom_point(size=2) +
#   xlab('LEP_R') + ylab('NP') + ggtitle('Scatter of LEP_R vs NP values') +
#   theme_bw()
# dev.off()
# 
# # Breeding size and LEP
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Breedingsize_LEP_scatter.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=breeding_size, y=LEP)) +
#   geom_point(size=2) +
#   xlab('breeding size') + ylab('LEP') + ggtitle('Scatter of breeding size (female) vs LEP values') +
#   theme_bw()
# dev.off()
# 
# # Breeding size and LEP and survival
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Breedingsize_LEP_Sint_scatter.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=breeding_size, y=LEP, color = Sint)) +
#   geom_point(size=2, alpha = 0.6) +
#   xlab('breeding size') + ylab('LEP') + ggtitle('Breeding size (female) vs LEP values, Sint') +
#   theme_bw()
# dev.off()
# 
# # Breeding size and LEP and growth
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Breedingsize_LEP_Linf_scatter.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=breeding_size, y=LEP, color = Linf)) +
#   geom_point(size=2, alpha = 0.7) +
#   xlab('breeding size') + ylab('LEP') + ggtitle('Breeding size (female) vs LEP values, Linf') +
#   theme_bw()
# dev.off()
# 
# # Linf and LEP
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Linf_LEP_scatter.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=Linf, y=LEP)) +
#   geom_point(size=2) +
#   xlab('Linf') + ylab('LEP') + ggtitle('Scatter of Linf vs LEP values') +
#   theme_bw()
# dev.off()
# 
# # Linf and LEP and Sint
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Linf_LEP_Sint_scatter.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=Linf, y=LEP, color = Sint)) +
#   geom_point(size=2) +
#   xlab('Linf') + ylab('LEP') + ggtitle('Linf vs LEP vs Sint values') +
#   theme_bw()
# dev.off()
# 
# # k (connectivity) and NP
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'kConnectivity_NP_scatter.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=k_connectivity, y=NP)) +
#   geom_point(size=2) +
#   xlab('k_connectivity') + ylab('NP') + ggtitle('Scatter of k_connectivity vs NP values') +
#   theme_bw()
# dev.off()
# 
# # k (connectivity) and NP - zoomed 
# pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'kConnectivity_NP_scatter_zoomed.pdf'))
# ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=k_connectivity, y=NP)) +
#   geom_point(size=2) +
#   ylim(0, 0.075) +
#   xlab('k_connectivity') + ylab('NP') + ggtitle('Scatter of k_connectivity vs NP values') +
#   theme_bw()
# dev.off()
# 
# # Prob r and recruits-per-egg
# prob_r_vs_RperE <- data.frame(prob_r = output_uncert_all$metric_vals_with_params$prob_r,
#                               recruits_per_egg = output_uncert_all$RperE_out_df$value,
#                               NP = output_uncert_all$NP_out_df$value)
# pdf(file =  here::here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'prob_r_recruits_per_egg_scatter.pdf'))
# ggplot(data = prob_r_vs_RperE, aes(x=prob_r, y=recruits_per_egg)) +
#   geom_point(size=2) +
#   xlab('capture probability') + ylab('recruits-per-egg') + ggtitle('P_c vs. recruits-per-egg') +
#   ylim(0,0.00025) +
#   theme_bw()
# dev.off()
# 
# # Prob r and NP
# pdf(file =  here::here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'prob_r_NP_scatter.pdf'))
# #ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=prob_r, y=NP)) +
# ggplot(data = prob_r_vs_RperE, aes(x = prob_r, y=NP)) +
#   geom_point(size=2) +
#   xlab('prob r') + ylab('NP') + ggtitle('Prob r vs. NP') +
#   theme_bw()
# dev.off()
# 
# ########## Relationships between values - prettier plots to put in manuscript ##########
# # Size at recruitment and LEP
# recruit_size_vs_LEP_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=start_recruit_size, y=LEP)) +
#   geom_point(size=2) +
#   xlab("recruit size") + ylab("LEP") + ggtitle("a) Recruit size vs. LEP") +
#   theme_bw()
#   
# # Size at female transition and egg-recruit-surv
# breeding_size_vs_RperE_df <- data.frame(breeding_size = output_uncert_all$metric_vals_with_params$breeding_size,
#                                         recruits_per_egg = output_uncert_all$RperE_out_df$value) 
# breeding_size_vs_RperE_plot <- ggplot(data = breeding_size_vs_RperE_df, aes(x=breeding_size, y=recruits_per_egg)) +
#   geom_point(size=2) +
#   xlab("breeding size (cm)") + ylab("recruits per egg") + ggtitle("b) Breed size vs. egg-recruit surv") +
#   theme_bw()
# 
# # LEP_R and NP 
# LRP_vs_NP_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=LEP_R, y=NP)) +
#   geom_point(size=2) +
#   xlab('LRP') + ylab('NP') + ggtitle('c) LRP vs. NP') +
#   theme_bw()
# 
# # Prob r and NP
# Prob_r_vs_NP_plot <- ggplot(data = prob_r_vs_RperE, aes(x = prob_r, y=NP)) +
#   geom_point(size=2) +
#   xlab('capture prob (P_c)') + ylab('NP') + ggtitle('d) Capture prob. vs. NP') +
#   theme_bw()
# 
# pdf(file = here::here("Plots/FigureDrafts", "Param_metric_relationships.pdf"))
# grid.arrange(recruit_size_vs_LEP_plot, breeding_size_vs_RperE_plot, 
#              LRP_vs_NP_plot, Prob_r_vs_NP_plot, nrow=2)
# dev.off()
# 
# 
# ########## Comparing different ways of estimating different values ##########
# 
# # Prob of catching a fish, made into a distribution two ways and sampling from data
# pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "Prob_r_set_comparison_with_data.pdf"))
# ggplot(data = prob_r_comp, aes(x = values, fill = distribution)) +
#   geom_histogram(binwidth = 0.01, position = "identity", alpha = 0.6) +
#   xlab("prob r") + ggtitle("Distributions of prob r") +
#   theme_bw()
# dev.off()
# 
# # Prob of catching a fish, made into a distribution two ways, without the data in there too
# pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "Prob_r_set_comparison.pdf"))
# ggplot(data = prob_r_comp %>% filter(distribution %in% c("beta", "normal", "truncated beta")), aes(x = values, fill = distribution)) +
#   geom_histogram(binwidth = 0.01, position = "identity", alpha = 0.5) +
#   xlab("prob r") + ggtitle("Distributions of prob r") +
#   theme_bw()
# dev.off()
# 
# # Recruits-per-egg in different ways
# pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "Recruits_per_egg_3kinds_uncertainty.pdf"))
# ggplot(data = recruits_per_egg_uncertainty, aes(x=recruits_per_egg, fill=method)) +
#   geom_histogram(binwidth = 0.000001, position = "identity", alpha = 0.7) +
#   xlab("Recruits-per-egg") + ggtitle("Comparing recruits-per-egg uncertainty methods") +
#   theme_bw()
# dev.off()
# 
# # Recruits-per-egg in different ways, zoomed in - this isn't working for some reason...
# pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "Recruits_per_egg_3kinds_uncertainty_zoomed.pdf"))
# ggplot(data = recruits_per_egg_uncertainty, aes(x=recruits_per_egg, fill=method)) +
#   geom_histogram(binwidth = 0.000001, position = "identity", alpha = 0.7) +
#   xlab("Recruits-per-egg") + ggtitle("Comparing recruits-per-egg uncert. meth. zoomed") +
#   #coord_cartesian(xlim=c(0, 0.3)) +
#   theme_bw()
# dev.off()
# 
# ########## What-if #1: what if all the offspring we genotype were from our sites? ##########
# 
# # LEP_R (LEP in terms of recruits) (this is plot "LEP_R_histogram_whatif_all_offspring.pdf") - where are these separate plots? Find and set to whatifs folder
# LEP_R_histogram_whatif_all_offspring <- ggplot(data = output_uncert_all_offspring_all$LEP_R_out_df, aes(x=value)) +
#   geom_histogram(bins=40, color = 'gray', fill = 'gray') +
#   geom_vline(data = LEP_R_best_est_all_offspring %>% filter(recruit_size == "mean offspring"), aes(xintercept = LEP_R)) +
#   xlab("recruits") + ggtitle("a) Recruits per recruit with all offspring") +
#   theme_bw()
# 
# # NP (this plot is "NP_histogram_whatif_all_offspring.pdf")
# NP_histogram_whatif_all_offspring <- ggplot(data = output_uncert_all_offspring_all$NP_out_df, aes(x=value)) +
#   geom_histogram(bins=40, color='gray', fill='gray') +
#   geom_vline(data = NP_best_est_all_offspring %>% filter(recruit_size == "mean offspring"), aes(xintercept = NP)) +
#   xlab("NP") + ggtitle("b) NP with all offspring") +
#   theme_bw()
# 
# pdf(file = here::here("Plots/PersistenceMetrics/Whatifs", "LEP_R_and_NP_histograms_whatif_all_offspring.pdf"), width=7, height=3)
# grid.arrange(LEP_R_histogram_whatif_all_offspring, NP_histogram_whatif_all_offspring, nrow=1)
# dev.off()
# 
# 
# ############# STOPPED EDITING HERE 
# 
# 
# 
# ## Input distributions, parameter estimates (breeding size, growth curve, survival curve, dispersal kernel?)
# distance_vec <- seq(from=0, to=50, by=0.01)
# connectivity_est_vec <- disp_kernel_all_years(distance_vec, k_allyears, theta_allyears)  # theta = 0.5, equation for p(d) in eqn. 6c in Bode et al. 2018
# connectivity_lowest_k <- disp_kernel_all_years(distance_vec, min(k_connectivity_values), theta_allyears)  # minimum k in the 95% CI
# connectivity_highest_k <- disp_kernel_all_years(distance_vec, max(k_connectivity_values), theta_allyears)  # maximum k in the 95% CI
# dispersal_df <- data.frame(distance = distance_vec, kernel_bestfit = connectivity_est_vec, kernel_CI1 = connectivity_lowest_k, kernel_CI2 = connectivity_highest_k) %>%
#   mutate(low_bound = case_when(kernel_CI1 <= kernel_CI2 ~ kernel_CI1,
#                                kernel_CI1 > kernel_CI2 ~ kernel_CI2),
#          upper_bound = case_when(kernel_CI1 > kernel_CI1 ~ kernel_CI1,
#                                  kernel_CI1 <= kernel_CI2 ~ kernel_CI2))
# 
# # Dispersal kernel
# dispersal_kernel_plot <- ggplot(data=dispersal_df, aes(x=distance, y=kernel_bestfit, ymin=kernel_CI1, ymax=kernel_CI2)) +
#   geom_line(color='black') +
#   geom_ribbon(alpha=0.5, color='gray') +
#   xlab('distance (km)') + ylab('dispersal probability') + ggtitle('a) Dispersal kernel') +
#   theme_bw()
#   
# # Growth curve
# growth_df <- data.frame(length1 = seq(min_size, max_size, length.out = n_bins*10)) %>%
#   mutate(length2 = VBL_growth(Linf_growth_mean, k_growth_mean, length1))
# 
# growth_curve_plot <- ggplot(data=growth_df, aes(x=length1, y=length2)) +
#   geom_point(color='black') +
#   xlab('size (cm)') + ylab('size next year') + ggtitle('b) Growth curve') +
#   ylim(0,13) +
#   theme_bw()
# 
# # Breeding size distribution
# breeding_size_plot <- ggplot(data = female_sizes, aes(x=size)) +
#   geom_histogram(bins=40, color='gray', fill='gray') +
#   geom_vline(xintercept=breeding_size_mean, color='black') +
#   xlab('size (cm)') + ggtitle('d) Female size') +
#   theme_bw()
# 
# # Survival plot
# min_size_plot = 0
# size.values <- min_size_plot+(0:30)*(max_size-min_size_plot)/30
# 
# eall_mean.Phi.size.p.size.plus.dist.results = as.data.frame(eall_mean.Phi.size.p.size.plus.dist$results$beta)
# 
# Phibysize_Phi.size.p.size.plus.dist_means <- data.frame(size = size.values) %>%
#   mutate(Phi_logit = eall_mean.Phi.size.p.size.plus.dist.results$estimate[1] + eall_mean.Phi.size.p.size.plus.dist.results$estimate[2]*size,
#          Phi_lcl_logit = eall_mean.Phi.size.p.size.plus.dist.results$lcl[1] + eall_mean.Phi.size.p.size.plus.dist.results$lcl[2]*size,
#          Phi_ucl_logit = eall_mean.Phi.size.p.size.plus.dist.results$ucl[1] + eall_mean.Phi.size.p.size.plus.dist.results$ucl[2]*size,
#          Phi = logit_recip(Phi_logit),
#          Phi_lcl = logit_recip(Phi_lcl_logit),
#          Phi_ucl = logit_recip(Phi_ucl_logit))
# 
# survival_plot <- ggplot(data = Phibysize_Phi.size.p.size.plus.dist_means, aes(size, Phi)) +
#   geom_ribbon(aes(ymin=Phi_lcl,ymax=Phi_ucl),color="gray",fill="gray") +
#   geom_line(color="black") +
#   xlab("size (cm)") + ylab("probability of survival") + ggtitle("c) Annual survival") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme_bw()
# 
# pdf(file = here('Plots/FigureDrafts','Parameter_inputs.pdf'), width=6, height=6)
# grid.arrange(dispersal_kernel_plot, growth_curve_plot, survival_plot, breeding_size_plot, nrow=2)
# dev.off()
# 
# 
# 
# 
# #################### Saving things: ####################
# save(best_est_metrics, file=here('Data', 'best_est_metrics.RData'))
# save(LEP_ests, file=here('Data', 'LEP_ests.RData'))
# 
# # #################### Old code: ####################
# # 
# # disp_theta_3 <- function(d) {  # theta = 3, equation for p(d) in eqn. 6b in Bode et al. 2018
# #   z = exp(k_allyears)
# #   disp = (3*z)/(gamma(1/3))*exp(-(z*d)^3)
# #   return(disp)
# # }
# 
# # eggs_per_clutch_mean = 514.11  # need to rethink what this means with size-effect, how to do one weighted by size...
# # clutches_per_year_mean = 11.9
# # egg_intercept = -426.57
# # egg_slope = 107  # eggs/cm size of F
# 
# # # Sint_mean = eall.Phi.size.p.dist.results$estimate[1]  # survival intercept (on logit scale)
# # # Sl_mean = eall.Phi.size.p.dist.results$estimate[2]  # survival slope (on logit scale)
# # # Sint_se = eall.Phi.size.p.dist.results$se[1]  # for now using SE, should really use SD...
# # # Sint_se = eall.Phi.size.p.dist.results$se[2]  # for now using SE, should really use SD...
# # 
# # # # Breeding size (for LEP) - replacing with drawing from the actual data
# # # breeding_size_mean = (size_by_color_metrics %>% filter(color == 'YP'))$mean  # originally guessed 8, this is 8.6
# # # breeding_size_sd = (size_by_color_metrics %>% filter(color == 'YP'))$sd  # originally guessed 0.8, this is 1.6
# # 
# # 
# # # # Static connectivity matrix loaded from PersistenceMetrics.R
# # # Cmatrix <- matrix(NA,ncol=max(c_mat_allyears$org_geo_order, na.rm = TRUE), nrow=max(c_mat_allyears$org_geo_order, na.rm = TRUE))    
# # # for(i in 1:length(c_mat_allyears$org_site)) {
# # #   column = c_mat_allyears$org_geo_order[i]  # column is origin 
# # #   row = c_mat_allyears$dest_geo_order[i]  # row is destination
# # #   Cmatrix[row, column] = c_mat_allyears$prob_disp_allyears[i]
# # # }
# # 
# # # # k (connectivity) and SP
# # # pdf(file = here('Plots/PersistenceMetrics','kConnectivity_SP_scatter.pdf'))
# # # ggplot(data = SP_vals_with_params, aes(x=k_connectivity, y=SP)) +
# # #   geom_point(size=2) +
# # #   facet_wrap(~site) +
# # #   xlab('k_connectivity') + ylab('SP') + ggtitle('Scatter of k_connectivity vs SP by site') +
# # #   theme_bw() +
# # #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  #not sure why this isn't working right now...
# # # dev.off()
# # 
# # # # Find k standard deviation - this is definitely not right - check with KC about her confidence intervals and how they were estimated
# # # z_97.5 = 2.24  # z-score for 97.5% confidence interval
# # # k_sdH = (sqrt(n_runs)/z_97.5)*(k_allyears_CIh - k_allyears)
# # # k_sdL = -(sqrt(n_runs)/z_97.5)*(k_allyears_CIl - k_allyears)
# # 
# # 
# # #k_connectivity_set = runif(n_runs, min = k_allyears_CIl, max = k_allyears_CIh)  # for now, just selecting randomly from within the 97.5% confidence interval
# # #breeding_size_set = rnorm(n_runs, mean = breeding_size_mean, sd = breeding_size_sd) 
# # 
# # ### OLD CODE, moved when shifting egg-recruit-survival calc into this script
# # # # Put best-estimate parameters into one dataframe
# # # param_best_est <- data.frame(t_steps = n_tsteps) %>%
# # #   mutate(min_size = min_size, max_size = max_size, n_bins = n_bins,  # LEP-IPM matrix
# # #          eggs_per_clutch = eggs_per_clutch_mean, clutches_per_year_mean = clutches_per_year_mean,  # average fecundity info
# # #          egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect =  eyed_effect,  # size-dependent fecundity info
# # #          start_recruit_size = start_recruit_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
# # #          k_growth = k_growth_mean, s = s, Sl = Sl_mean, Linf = Linf_growth_mean, Sint = Sint_mean,
# # #          breeding_size = breeding_size_mean, recruits_per_egg = recruits_per_egg,
# # #          k_connectivity = k_allyears, theta_connectivity = theta_allyears)  # dispersal kernel parameters
# # 
# # # # Create connectivity matrix in matrix form
# # # c_mat_allyears_matrix <- c_mat_allyears %>% select(org_site, dest_site, d1_km, d2_km, org_alpha_order, org_geo_order, dest_alpha_order, dest_geo_order)
# # # for(i in 1:length(c_mat_allyears$org_site)) {
# # #   c_mat_allyears_matrix$prob_disp[i] <- integrate(disp_allyears, c_mat_allyears_matrix$d1_km[i], c_mat_allyears_matrix$d2_km[i])$value
# # # }
# # # 
# # # conn_matrix <- matrix(NA,ncol=max(Cmat$org_geo_order), nrow=max(Cmat$org_geo_order))    
# # # for(i in 1:length(Cmat$org_site)) {
# # #   column = Cmat$org_geo_order[i]  # column is origin 
# # #   row = Cmat$dest_geo_order[i]  # row is destination
# # #   conn_matrix[row, column] = Cmat$prob_disp[i]
# # # }
# # 
# # # How many potential offspring (eggs) were produced by tagged adults (parents)
# # tagged_offspring_6cm <- n_parents_parentage*LEP_6cm  # Use LEP from what size here? How to avoid double-counting if those parents mated together?
# # 
# # 
# # # Estimate survival from eggs-recruits by seeing how many "tagged" offspring we found out of eggs "tagged" parents produced
# # surv_egg_recruit <- re
# # 
# # # Estimate survival from adult-recruit
# # prop_F_M <- 0.5  # saying 50% of the "adults" we clip are males that won't make it to females -- reasonable? could check this. But LEP takes that into account, right?
# # tagged_offspring_3.5cm <- n_parents_parentage*LEP_ests$LEP_3.5cm
# # tagged_offspring_6cm <- n_parents_parentage*LEP_ests$LEP_6cm
# # recruited_offspring <- n_offspring_parentage/(total_prop_hab_sampled*mean(prob_r))  # scale up by proportion of habitat sampled and probability of catching a fish
# # surv_egg_recruit <- recruited_offspring/tagged_offspring_3.5cm
# # 
# # ##### Find 'best estimates' of some of the parameters
# # breeding_size_mean <- mean(female_sizes$size, na.rm=TRUE)
# # prob_r_mean <- mean(prob_r)  # average value of prob r from each recap dive  -- WHERE DOES THIS COME IN? OTHER SCRIPTS?
# # 
# # ##### Find 
# # 
# # ##### Generate sets of parameters
# # Linf_set = rnorm(n_runs, mean = Linf_growth_mean, sd=Linf_growth_sd)
# # Sint_set = rnorm(n_runs, mean = Sint_mean, sd= Sint_se)
# # k_connectivity_set = sample(k_connectivity_values, n_runs, replace=TRUE)  # replace should be true, right?
# # breeding_size_set = sample(female_sizes$size, n_runs, replace=TRUE)  # replace should be true, right?
# # prob_r_set = rnorm(n_runs, mean = prob_r_mean, sd=sd(prob_r))  
# # 
# # # Put static + pulled-from-distribution parameters together into one dataframe
# # param_set_full <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
# #   mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
# #          eggs_per_clutch = eggs_per_clutch_mean, clutches_per_year = clutches_per_year_mean,  # fecundity info
# #          egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
# #          start_recruit_size = start_recruit_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
# #          k_growth = k_growth_mean, s = s, Sl = Sl_mean, Linf = Linf_set, Sint = Sint_set,
# #          breeding_size = breeding_size_set, recruits_per_egg = recruits_per_egg,
# #          k_connectivity = k_connectivity_set, theta_connectivity = theta_allyears)  # dispersal kernel parameters
# # 
# # site_list <- c_mat_allyears$dest_site[1:19]
# # 
# # 
# # 
# # ##### Find the 'best-estimate' metrics
# # best_est_metrics <- calcMetrics(param_best_est, c_mat_allyears, site_list)
# # 
# # # Create connectivity matrix in matrix form
# # c_mat_allyears_matrix <- c_mat_allyears %>% select(org_site, dest_site, d1_km, d2_km, org_alpha_order, org_geo_order, dest_alpha_order, dest_geo_order)
# # for(i in 1:length(c_mat_allyears$org_site)) {
# #   c_mat_allyears_matrix$prob_disp[i] <- integrate(disp_allyears, c_mat_allyears_matrix$d1_km[i], c_mat_allyears_matrix$d2_km[i])$value
# # }
# # 
# # conn_matrix <- matrix(NA,ncol=max(Cmat$org_geo_order), nrow=max(Cmat$org_geo_order))    
# # for(i in 1:length(Cmat$org_site)) {
# #   column = Cmat$org_geo_order[i]  # column is origin 
# #   row = Cmat$dest_geo_order[i]  # row is destination
# #   conn_matrix[row, column] = Cmat$prob_disp[i]
# # }
# # 
# # 
# # # Find LEP for mean breeding size and for size 6.0
# # LEP_breeding_size_mean <- findLEP(param_best_est$min_size, param_best_est$max_size, param_best_est$n_bins, 
# #                                   param_best_est$t_steps, param_best_est$Sint, param_best_est$Sl,
# #                                   param_best_est$s, param_best_est$Linf, param_best_est$k_growth, 
# #                                   param_best_est$eggs_per_clutch, param_best_est$clutches_per_year, 
# #                                   param_best_est$breeding_size, breeding_size_mean, param_best_est$start_recruit_sd, 
# #                                   param_best_est$egg_size_slope, param_best_est$egg_size_intercept, param_best_est$eyed_effect)
# # 
# # LEP_6cm <- findLEP(param_best_est$min_size, param_best_est$max_size, param_best_est$n_bins, 
# #                    param_best_est$t_steps, param_best_est$Sint, param_best_est$Sl,
# #                    param_best_est$s, param_best_est$Linf, param_best_est$k_growth, 
# #                    param_best_est$eggs_per_clutch, param_best_est$clutches_per_year, 
# #                    param_best_est$breeding_size, 6, param_best_est$start_recruit_sd, 
# #                    param_best_est$egg_size_slope, param_best_est$egg_size_intercept, param_best_est$eyed_effect)
# # 
# # LEP_3.5cm <- findLEP(param_best_est$min_size, param_best_est$max_size, param_best_est$n_bins, 
# #                      param_best_est$t_steps, param_best_est$Sint, param_best_est$Sl,
# #                      param_best_est$s, param_best_est$Linf, param_best_est$k_growth, 
# #                      param_best_est$eggs_per_clutch, param_best_est$clutches_per_year, 
# #                      param_best_est$breeding_size, 3.5, param_best_est$start_recruit_sd, 
# #                      param_best_est$egg_size_slope, param_best_est$egg_size_intercept, param_best_est$eyed_effect)
# # 
# # LEP_ests <- list(LEP_breeding_size_mean = LEP_breeding_size_mean, LEP_6cm = LEP_6cm, LEP_3.5cm = LEP_3.5cm)
# # 
# # # Save as separate items, for plotting ease
# # LEP_best_est <- best_est_metrics$LEP
# # LEP_R_best_est <- best_est_metrics$LEP_R
# # NP_best_est <- best_est_metrics$NP
# # SP_best_est <- as.data.frame(best_est_metrics$SP)
# # 
# # ##### Run the metrics for lots of parameters
# # # Set output dataframes 
# # n_metrics = 4  # NP, LEP, LEP_R, recruits_per_egg
# # LEP_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP', inout = 'input', run = seq(1:n_runs))
# # LEP_R_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP_R', inout = 'input', run = seq(1:n_runs))
# # RperE_out_df <- data.frame(value = rep(NA, n_runs), metric = 'recruits per egg', inout = 'input', run = seq(1:n_runs))
# # NP_out_df <- data.frame(value = rep(NA, n_runs), metric = 'NP', inout = 'output', run = seq(1:n_runs))
# # 
# # metric_vals <- data.frame(run = seq(1:n_runs), LEP = NA, LEP_R = NA, recruits_per_egg = NA, NP = NA)
# # 
# # # Create vector of sites for SP dataframe      
# # runsrepped = rep(1, length(site_list))                                                         
# # for(i in 2:n_runs) {
# #   runsrepped = c(runsrepped, rep(i, length(site_list)))
# # }                                 
# # 
# # SP_out_df <- data.frame(value = rep(NA, length(site_list)*n_runs), metric = 'SP', run = NA,
# #                         site = NA, org_geo_order = NA)
# # 
# # # Calculate the metrics for each parameter set, fill into the data frames
# # for(i in 1:n_runs) {
# #   # Select parameter set
# #   params <- param_set_full[i,]
# #   
# #   # Do the run
# #   metrics_output = calcMetrics(params, c_mat_allyears, site_list)
# #   
# #   # Fill in the metrics
# #   LEP_out_df$value[i] = metrics_output$LEP
# #   LEP_R_out_df$value[i] = metrics_output$LEP_R
# #   RperE_out_df$value[i] = metrics_output$recruits_per_egg
# #   NP_out_df$value[i] = metrics_output$NP
# #   
# #   metric_vals$LEP[i] = metrics_output$LEP
# #   metric_vals$LEP_R[i] = metrics_output$LEP_R
# #   metric_vals$recruits_per_egg[i] = metrics_output$recruits_per_egg
# #   metric_vals$NP[i] = metrics_output$NP
# #   
# #   # Pull out the SP metrics
# #   start_index = (i-1)*length(site_list)+1
# #   end_index = i*length(site_list)
# #   
# #   SP_out_df$site[start_index:end_index] = metrics_output$SP$site
# #   SP_out_df$value[start_index:end_index] = metrics_output$SP$SP_value
# #   SP_out_df$run[start_index:end_index] = rep(i, length(site_list))
# #   SP_out_df$org_geo_order[start_index:end_index] = metrics_output$SP$org_geo_order
# # }
# # 
# # # Put the data frames together (easier to plot?)
# # #metrics_vals <- rbind(LEP_out_df, LEP_R_out_df, RperE_out_df, NP_out_df)
# # 
# # # Add some of the changing parameters in, so can look at in plots
# # metric_vals_with_params <- metric_vals %>%
# #   mutate(breeding_size = param_set_full$breeding_size,
# #          Linf = param_set_full$Linf,
# #          Sint = param_set_full$Sint,
# #          k_connectivity = param_set_full$k_connectivity)
# # 
# # SP_vals_with_params <- left_join(SP_out_df, metric_vals_with_params, by='run') %>%
# #   dplyr::rename(SP = value)
# # 
# # # Pull out the self-persistence values
# # 
# # ## Change size-dependent survival estimates from logit estimates in loaded file
# # #Mint = logit_recip(eall.Phi.size.p.dist.results$estimate[1])
# # #Ml = logit_recip(eall.Phi.size.p.dist.results$estimate[2])
# # 
# # #################### Plots: ####################
# # 
# # ##### Plot the histograms of LEP, LEP_R, recruits_per_egg, and NP output
# # # LEP
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_histogram.pdf'))
# # ggplot(data = LEP_out_df, aes(x=value)) +
# #   geom_histogram(bins=40, color = 'gray', fill = 'gray') +
# #   geom_vline(xintercept = LEP_best_est, color='black') +
# #   xlab('LEP') + ggtitle('Histogram of LEP values') +
# #   theme_bw()
# # dev.off()
# # 
# # # LEP_R (LEP in terms of recruits)
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_R_histogram.pdf'))
# # ggplot(data = LEP_R_out_df, aes(x=value)) +
# #   geom_histogram(bins=40, color = 'gray', fill = 'gray') +
# #   geom_vline(xintercept = LEP_R_best_est, color = 'black') +
# #   xlab('LEP_R') + ggtitle('Histogram of LEP_R values') +
# #   theme_bw()
# # dev.off()
# # 
# # # # Doing it this way looks a bit different... seems like there are more breaks, even though I also tried to ask for 30?
# # # pdf(file = here('Plots/PersistenceMetrics', 'LEP_R_histv2.pdf'))
# # # hist(LEP_R_out_df$value, breaks=30)
# # # dev.off()
# # 
# # # NP
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'NP_histogram.pdf'))
# # ggplot(data = NP_out_df, aes(x=value)) +
# #   geom_histogram(bins=25, color='gray', fill='gray') +
# #   geom_vline(xintercept = NP_best_est, color='black') +
# #   xlab('NP') + ggtitle('Histogram of NP values') +
# #   theme_bw()
# # dev.off()
# # 
# # # recruits_per_egg (but right now this is static)
# # 
# # # hist(LEP_out_df$value, breaks=30)
# # # hist(LEP_R_out$value, breaks=30)
# # # hist(NP_out_df$value, breaks=30)
# # # hist(RperE_out_df$value, breaks=30)
# # 
# # ##### SP at each site
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty','SP_histogram.pdf'))
# # ggplot(data = SP_out_df, aes(x=value)) +
# #   #geom_histogram(binwidth=0.0005) +
# #   geom_histogram(binwidth=0.001, color='gray', fill='gray') +
# #   geom_vline(data=SP_best_est, aes(xintercept=SP_value), color='black') +   # now these look weird, don't seem to fit in the dists for most sites?
# #   ylim(0,300) +
# #   facet_wrap(~site) +
# #   xlab('SP') + ggtitle('Self-persistence histograms by site') +
# #   theme_bw() +
# #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  #not sure why this isn't working right now...
# # dev.off()
# # 
# # ##### Relationships between values
# # # LEP_R and NP - do we expect these to be related perfectly linearly? I guess, when both connectivity and recruits-per-egg are static...
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_R_NP_scatter.pdf'))
# # ggplot(data = metric_vals, aes(x=LEP_R, y=NP)) +
# #   geom_point(size=2) +
# #   xlab('LEP_R') + ylab('NP') + ggtitle('Scatter of LEP_R vs NP values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Breeding size and LEP
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Breedingsize_LEP_scatter.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=breeding_size, y=LEP)) +
# #   geom_point(size=2) +
# #   xlab('breeding size') + ylab('LEP') + ggtitle('Scatter of breeding size (female) vs LEP values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Linf and LEP
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Linf_LEP_scatter.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=Linf, y=LEP)) +
# #   geom_point(size=2) +
# #   #geom_line()
# #   #geom_line(aes(x=sss)),
# #   xlab('Linf') + ylab('LEP') + ggtitle('Scatter of Linf vs LEP values') +
# #   theme_bw()
# # dev.off()
# # 
# # # k (connectivity) and NP
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'kConnectivity_NP_scatter.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=k_connectivity, y=NP)) +
# #   geom_point(size=2) +
# #   #geom_line(aes(x=sss)),
# #   xlab('k_connectivity') + ylab('NP') + ggtitle('Scatter of k_connectivity vs NP values') +
# #   theme_bw()
# # dev.off()
# # 
# # ##### Histograms of data inputs
# # # k_connectivity
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'k_connectivity_histogram.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=k_connectivity)) +
# #   geom_histogram(bins=40, color='gray', fill='gray') +
# #   geom_vline(xintercept = k_allyears, color='black') +
# #   xlab('k_connectivity') + ggtitle('k_connectivity values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Linf
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Linf_histogram.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=Linf)) +
# #   geom_histogram(bins=40, color='gray', fill='gray') +
# #   geom_vline(xintercept = Linf_growth_mean, color='black') +
# #   xlab('Linf') + ggtitle('Linf values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Sint
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Sint_histogram.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=Sint)) +
# #   geom_histogram(bins=40, color='gray', fill='gray') +
# #   geom_vline(xintercept = Sint_mean, color = 'black') +
# #   xlab('Sint') + ggtitle('Sint values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Breeding size
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Breeding_size_histogram.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=breeding_size)) +
# #   geom_histogram(bins=40, color='gray', fill='gray') +
# #   geom_vline(xintercept=breeding_size_mean, color='black') +
# #   xlab('Breeding size') + ggtitle('Breeding size values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Recruits-per-egg in different ways
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "Recruits_per_egg_3kinds_uncertainty.pdf"))
# # ggplot(data = recruits_per_egg_uncertainty, aes(x=recruits_per_egg, fill=method)) +
# #   geom_histogram(binwidth = 0.001) +
# #   xlab("Recruits-per-egg") + ggtitle("Comparing recruits-per-egg uncertainty methods") +
# #   theme_bw()
# # dev.off()
# # 
# # # Recruits-per-egg in different ways, zoomed in - this isn't working for some reason...
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "Recruits_per_egg_3kinds_uncertainty_zoomed.pdf"))
# # ggplot(data = recruits_per_egg_uncertainty, aes(x=recruits_per_egg, fill=method)) +
# #   geom_histogram(binwidth = 0.001) +
# #   xlab("Recruits-per-egg") + ggtitle("Comparing recruits-per-egg uncert. meth. zoomed") +
# #   coord_cartesian(xlim=c(0, 0.3)) +
# #   theme_bw()
# # dev.off()
# # 
# # 
# # 
# # 
# # ##### Prettier sub-figured plots for potential figures
# # ## LEP and LEP_R histograms
# # LEP_plot <- ggplot(data = LEP_out_df, aes(x=value)) +
# #   geom_histogram(bins=40, color = 'gray', fill = 'gray') +
# #   geom_vline(xintercept = LEP_best_est, color='black') +
# #   xlab('LEP') + ggtitle('a) Lifetime egg production') +
# #   theme_bw() 
# # #theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0))
# # 
# # LEP_R_plot <- ggplot(data = LEP_R_out_df, aes(x=value)) +
# #   geom_histogram(bins=40, color = 'gray', fill = 'gray') +
# #   geom_vline(xintercept = LEP_R_best_est, color = 'black') +
# #   xlab('LRP') + ggtitle('b) Lifetime recruit production') +
# #   theme_bw()
# # 
# # pdf(file = here('Plots/FigureDrafts', 'LEP_and_LRP.pdf'), width=7, height=3)
# # grid.arrange(LEP_plot, LEP_R_plot, nrow=1)
# # dev.off()
# # 
# # ## NP and realized connectivity matrix
# # NP_plot <- ggplot(data = NP_out_df, aes(x=value)) +
# #   geom_histogram(bins=25, color='gray', fill='gray') +
# #   geom_vline(xintercept = NP_best_est, color='black') +
# #   xlab('NP') + ggtitle('a) Network persistence values') +
# #   theme_bw()
# # 
# # realized_C_plot <- ggplot(data = best_est_metrics$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
# #   geom_tile(aes(fill=prob_disp_R)) +
# #   scale_fill_gradient(high='black', low='white', name='Recruits') +
# #   xlab('origin') + ylab('destination') + ggtitle('b) Realized connectivity matrix') +
# #   theme_bw() +
# #   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
# # 
# # pdf(file = here('Plots/FigureDrafts', 'NP_and_connMatrixR.pdf'), width = 9, height = 4)
# # grid.arrange(NP_plot, realized_C_plot, nrow=1)
# # dev.off()
# # 
# # ## SP by site
# # pdf(file = here('Plots/FigureDrafts','SP_hists_by_site.pdf'), width=8, height=6)
# # ggplot(data = (SP_out_df %>% filter(site != 'Sitio Lonas')), aes(x=value)) +
# #   #geom_histogram(binwidth=0.0005) +
# #   geom_histogram(binwidth=0.001, color='gray', fill='gray') +
# #   geom_vline(data=(SP_best_est %>% filter(site != 'Sitio Lonas')), aes(xintercept=SP_value), color='black') +   # now these look weird, don't seem to fit in the dists for most sites?
# #   ylim(0,150) +
# #   facet_wrap(~reorder(site, org_geo_order)) +
# #   xlab('SP') + ggtitle('Self-persistence estimates by site') +
# #   theme_bw() +
# #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  #not sure why this isn't working right now...
# # dev.off()
# # 
# # ## Input distributions, parameter estimates (breeding size, growth curve, survival curve, dispersal kernel?)
# # distance_vec <- seq(from=0, to=50, by=0.01)
# # connectivity_est_vec <- disp_kernel_all_years(distance_vec, k_allyears, theta_allyears)  # theta = 0.5, equation for p(d) in eqn. 6c in Bode et al. 2018
# # connectivity_lowest_k <- disp_kernel_all_years(distance_vec, min(k_connectivity_values), theta_allyears)  # minimum k in the 95% CI
# # connectivity_highest_k <- disp_kernel_all_years(distance_vec, max(k_connectivity_values), theta_allyears)  # maximum k in the 95% CI
# # dispersal_df <- data.frame(distance = distance_vec, kernel_bestfit = connectivity_est_vec, kernel_CI1 = connectivity_lowest_k, kernel_CI2 = connectivity_highest_k) %>%
# #   mutate(low_bound = case_when(kernel_CI1 <= kernel_CI2 ~ kernel_CI1,
# #                                kernel_CI1 > kernel_CI2 ~ kernel_CI2),
# #          upper_bound = case_when(kernel_CI1 > kernel_CI1 ~ kernel_CI1,
# #                                  kernel_CI1 <= kernel_CI2 ~ kernel_CI2))
# # 
# # # Dispersal kernel
# # dispersal_kernel_plot <- ggplot(data=dispersal_df, aes(x=distance, y=kernel_bestfit, ymin=kernel_CI1, ymax=kernel_CI2)) +
# #   geom_line(color='black') +
# #   geom_ribbon(alpha=0.5, color='gray') +
# #   xlab('distance (km)') + ylab('dispersal probability') + ggtitle('a) Dispersal kernel') +
# #   theme_bw()
# # 
# # # Growth curve
# # growth_df <- data.frame(length1 = seq(min_size, max_size, length.out = n_bins*10)) %>%
# #   mutate(length2 = VBL_growth(Linf_growth_mean, k_growth_mean, length1))
# # 
# # growth_curve_plot <- ggplot(data=growth_df, aes(x=length1, y=length2)) +
# #   geom_point(color='black') +
# #   xlab('size (cm)') + ylab('size next year') + ggtitle('b) Growth curve') +
# #   ylim(0,13) +
# #   theme_bw()
# # 
# # # Breeding size distribution
# # breeding_size_plot <- ggplot(data = female_sizes, aes(x=size)) +
# #   geom_histogram(bins=40, color='gray', fill='gray') +
# #   geom_vline(xintercept=breeding_size_mean, color='black') +
# #   xlab('size (cm)') + ggtitle('d) Female size') +
# #   theme_bw()
# # 
# # # Survival plot
# # min_size_plot = 0
# # size.values <- min_size_plot+(0:30)*(max_size-min_size_plot)/30
# # 
# # eall_mean.Phi.size.p.size.plus.dist.results = as.data.frame(eall_mean.Phi.size.p.size.plus.dist$results$beta)
# # 
# # Phibysize_Phi.size.p.size.plus.dist_means <- data.frame(size = size.values) %>%
# #   mutate(Phi_logit = eall_mean.Phi.size.p.size.plus.dist.results$estimate[1] + eall_mean.Phi.size.p.size.plus.dist.results$estimate[2]*size,
# #          Phi_lcl_logit = eall_mean.Phi.size.p.size.plus.dist.results$lcl[1] + eall_mean.Phi.size.p.size.plus.dist.results$lcl[2]*size,
# #          Phi_ucl_logit = eall_mean.Phi.size.p.size.plus.dist.results$ucl[1] + eall_mean.Phi.size.p.size.plus.dist.results$ucl[2]*size,
# #          Phi = logit_recip(Phi_logit),
# #          Phi_lcl = logit_recip(Phi_lcl_logit),
# #          Phi_ucl = logit_recip(Phi_ucl_logit))
# # 
# # survival_plot <- ggplot(data = Phibysize_Phi.size.p.size.plus.dist_means, aes(size, Phi)) +
# #   geom_ribbon(aes(ymin=Phi_lcl,ymax=Phi_ucl),color="gray",fill="gray") +
# #   geom_line(color="black") +
# #   xlab("size (cm)") + ylab("probability of survival") + ggtitle("c) Annual survival") +
# #   scale_y_continuous(limits = c(0, 1)) +
# #   theme_bw()
# # 
# # pdf(file = here('Plots/FigureDrafts','Parameter_inputs.pdf'), width=6, height=6)
# # grid.arrange(dispersal_kernel_plot, growth_curve_plot, survival_plot, breeding_size_plot, nrow=2)
# # dev.off()
# # 
# # 
# # 
# # 
# # #################### Saving things: ####################
# # save(best_est_metrics, file=here('Data', 'best_est_metrics.RData'))
# # save(LEP_ests, file=here('Data', 'LEP_ests.RData'))
# # 
# # 
# # 
# # 
# # #################### Plots: ####################
# # # 
# # 
# # 
# # # # Doing it this way looks a bit different... seems like there are more breaks, even though I also tried to ask for 30?
# # # pdf(file = here('Plots/PersistenceMetrics', 'LEP_R_histv2.pdf'))
# # # hist(LEP_R_out_df$value, breaks=30)
# # # dev.off()
# # 
# # 
# # 
# # # recruits_per_egg (but right now this is static)
# #   
# # # hist(LEP_out_df$value, breaks=30)
# # # hist(LEP_R_out$value, breaks=30)
# # # hist(NP_out_df$value, breaks=30)
# # # hist(RperE_out_df$value, breaks=30)
# # 
# # 
# # ##### Relationships between values
# # # LEP_R and NP - do we expect these to be related perfectly linearly? I guess, when both connectivity and recruits-per-egg are static...
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_R_NP_scatter.pdf'))
# # ggplot(data = metric_vals, aes(x=LEP_R, y=NP)) +
# #   geom_point(size=2) +
# #   xlab('LEP_R') + ylab('NP') + ggtitle('Scatter of LEP_R vs NP values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Breeding size and LEP
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Breedingsize_LEP_scatter.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=breeding_size, y=LEP)) +
# #   geom_point(size=2) +
# #   xlab('breeding size') + ylab('LEP') + ggtitle('Scatter of breeding size (female) vs LEP values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Linf and LEP
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Linf_LEP_scatter.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=Linf, y=LEP)) +
# #   geom_point(size=2) +
# #   #geom_line()
# #   #geom_line(aes(x=sss)),
# #   xlab('Linf') + ylab('LEP') + ggtitle('Scatter of Linf vs LEP values') +
# #   theme_bw()
# # dev.off()
# # 
# # # k (connectivity) and NP
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'kConnectivity_NP_scatter.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=k_connectivity, y=NP)) +
# #   geom_point(size=2) +
# #   #geom_line(aes(x=sss)),
# #   xlab('k_connectivity') + ylab('NP') + ggtitle('Scatter of k_connectivity vs NP values') +
# #   theme_bw()
# # dev.off()
# # 
# # ##### Histograms of data inputs
# # # k_connectivity
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'k_connectivity_histogram.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=k_connectivity)) +
# #   geom_histogram(bins=40, color='gray', fill='gray') +
# #   geom_vline(xintercept = k_allyears, color='black') +
# #   xlab('k_connectivity') + ggtitle('k_connectivity values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Linf
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Linf_histogram.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=Linf)) +
# #   geom_histogram(bins=40, color='gray', fill='gray') +
# #   geom_vline(xintercept = Linf_growth_mean, color='black') +
# #   xlab('Linf') + ggtitle('Linf values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Sint
# # pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Sint_histogram.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=Sint)) +
# #   geom_histogram(bins=40, color='gray', fill='gray') +
# #   geom_vline(xintercept = Sint_mean, color = 'black') +
# #   xlab('Sint') + ggtitle('Sint values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Breeding size
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Breeding_size_histogram.pdf'))
# # ggplot(data = metric_vals_with_params, aes(x=breeding_size)) +
# #   geom_histogram(bins=40, color='gray', fill='gray') +
# #   geom_vline(xintercept=breeding_size_mean, color='black') +
# #   xlab('Breeding size') + ggtitle('Breeding size values') +
# #   theme_bw()
# # dev.off()
# # 
# # # Recruits-per-egg in different ways
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "Recruits_per_egg_3kinds_uncertainty.pdf"))
# # ggplot(data = recruits_per_egg_uncertainty, aes(x=recruits_per_egg, fill=method)) +
# #   geom_histogram(binwidth = 0.001) +
# #   xlab("Recruits-per-egg") + ggtitle("Comparing recruits-per-egg uncertainty methods") +
# #   theme_bw()
# # dev.off()
# # 
# # # Recruits-per-egg in different ways, zoomed in - this isn't working for some reason...
# # pdf(file = here::here("Plots/PersistenceMetrics/MetricsWithUncertainty", "Recruits_per_egg_3kinds_uncertainty_zoomed.pdf"))
# # ggplot(data = recruits_per_egg_uncertainty, aes(x=recruits_per_egg, fill=method)) +
# #   geom_histogram(binwidth = 0.001) +
# #   xlab("Recruits-per-egg") + ggtitle("Comparing recruits-per-egg uncert. meth. zoomed") +
# #   coord_cartesian(xlim=c(0, 0.3)) +
# #   theme_bw()
# # dev.off()
# # 
# # # 
# # ##### Prettier sub-figured plots for potential figures - what I showed at lab meeting in early March 2019
# # ## LEP and LEP_R histograms
# # LEP_plot <- ggplot(data = LEP_out_df, aes(x=value)) +
# #   geom_histogram(bins=40, color = 'gray', fill = 'gray') +
# #   geom_vline(xintercept = LEP_best_est, color='black') +
# #   xlab('LEP') + ggtitle('a) Lifetime egg production') +
# #   theme_bw() 
# # #theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0))
# # 
# # LEP_R_plot <- ggplot(data = LEP_R_out_df, aes(x=value)) +
# #   geom_histogram(bins=40, color = 'gray', fill = 'gray') +
# #   geom_vline(xintercept = LEP_R_best_est, color = 'black') +
# #   xlab('LRP') + ggtitle('b) Lifetime recruit production') +
# #   theme_bw()
# # 
# # pdf(file = here('Plots/FigureDrafts', 'LEP_and_LRP.pdf'), width=7, height=3)
# # grid.arrange(LEP_plot, LEP_R_plot, nrow=1)
# # dev.off()
# # 
# # ## NP and realized connectivity matrix
# # NP_plot <- ggplot(data = NP_out_df, aes(x=value)) +
# #   geom_histogram(bins=25, color='gray', fill='gray') +
# #   geom_vline(xintercept = NP_best_est, color='black') +
# #   xlab('NP') + ggtitle('a) Network persistence values') +
# #   theme_bw()
# # 
# # realized_C_plot <- ggplot(data = best_est_metrics$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
# #   geom_tile(aes(fill=prob_disp_R)) +
# #   scale_fill_gradient(high='black', low='white', name='Recruits') +
# #   xlab('origin') + ylab('destination') + ggtitle('b) Realized connectivity matrix') +
# #   theme_bw() +
# #   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
# # 
# # pdf(file = here('Plots/FigureDrafts', 'NP_and_connMatrixR.pdf'), width = 9, height = 4)
# # grid.arrange(NP_plot, realized_C_plot, nrow=1)
# # dev.off()
# 
# # ## SP by site
# # pdf(file = here('Plots/FigureDrafts','SP_hists_by_site.pdf'), width=8, height=6)
# # ggplot(data = (SP_out_df %>% filter(site != 'Sitio Lonas')), aes(x=value)) +
# #   #geom_histogram(binwidth=0.0005) +
# #   geom_histogram(binwidth=0.001, color='gray', fill='gray') +
# #   geom_vline(data=(SP_best_est %>% filter(site != 'Sitio Lonas')), aes(xintercept=SP_value), color='black') +   # now these look weird, don't seem to fit in the dists for most sites?
# #   ylim(0,150) +
# #   facet_wrap(~reorder(site, org_geo_order)) +
# #   xlab('SP') + ggtitle('Self-persistence estimates by site') +
# #   theme_bw() +
# #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  #not sure why this isn't working right now...
# # dev.off()
# # 
# 
# # ########## Sub-figured plots to be able to look at things together ##########
# # 
# # ##### What if 1: all genotyped offspring considered to have come from our sites
# # 
# # # LEP_R (LEP in terms of recruits) 
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_R_histogram_whatif_all_offspring.pdf'))
# # ggplot(data = output_uncert_all_offspring_all$LEP_R_out_df, aes(x=value)) +
# #   geom_histogram(bins=40, color = 'gray', fill = 'gray') +
# #   geom_vline(data = LEP_R_best_est_all_offspring, aes(xintercept = LEP_R, color = recruit_size)) +
# #   #geom_vline(xintercept = (LEP_R_best_est %>% filter(recruit_size == "3.5cm"))$LEP_R, color = 'black') +
# #   xlab('LEP_R') + ggtitle('Histogram of LEP_R values with all offspring') +
# #   theme_bw()
# # dev.off()
# # 
# # # NP
# # pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'NP_histogram_whatif_all_offspring.pdf'))
# # ggplot(data = output_uncert_all_offspring_all$NP_out_df, aes(x=value)) +
# #   geom_histogram(bins=40, color='gray', fill='gray') +
# #   geom_vline(data = NP_best_est_all_offspring, aes(xintercept = NP, color = recruit_size)) +
# #   #geom_vline(xintercept = (NP_best_est %>% filter(recruit_size == "3.5cm"))$NP, color='black') +
# #   xlab('NP') + ggtitle('Histogram of NP values with all offspring') +
# #   theme_bw()
# # dev.off()
# # 


############### Code cut from running things section on 11/14/19 tidying of script ###############
# # with DD scaling
# param_best_est_mean_collected_offspring_DD <- data.frame(t_steps = n_tsteps) %>%
#   mutate(min_size = min_size, max_size = max_size, n_bins = n_bins,  # LEP-IPM matrix
#          eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year_mean = clutches_per_year_mean,  # average fecundity info
#          egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect =  eyed_effect,  # size-dependent fecundity info
#          start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
#          k_growth = k_growth_mean, s = s, Sl = Sl_mean, Linf = Linf_growth_mean, Sint = Sint_mean,
#          breeding_size = breeding_size_mean, recruits_per_egg = recruits_per_egg_best_est,
#          k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
#          prob_r = prob_r_mean, offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
#          total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
#          perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, DD = TRUE) 

# # Calculate the metrics for the best estimates
# best_est_metrics_mean_offspring <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info, site_vec_order$site_name, "FALSE")
# best_est_metrics_mean_offspring_DD <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info, site_vec_order$site_name, "TRUE")

# # Scale up by the proportion of site area we sampled over the time frame of finding parentage matches and prob of catching a fish
# recruited_tagged_offspring_oursites <- n_offspring_matched/((total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area*mean(prob_r))  # scale up by proportion of habitat sampled and probability of catching a fish
# 
# # Scale up by the proportion of the kernel area included within our sites (assuming patchiness of habitat outside our sites is the same as within)
# recruited_tagged_offspring_kernel <- recruited_tagged_offspring_oursites/prop_total_disp_area_sampled_best_est
# 
# # Scale up by the proportion of the sampling area (assumed to be the same as the greater region) that is habitat
# recruited_tagged_offspring_habprop <- recruited_tagged_offspring_kernel/prop_sampling_area_habitat

# # Estimate survival from eggs-recruits by seeing how many "tagged" offspring we found out of eggs "tagged" parents produced
# recruits_per_egg_best_est <- recruited_tagged_offspring_habprop/tagged_eggs_6cm
# 
# # Scale up by DD
# recruited_tagged_offspring_DDscaled <- scaleTaggedRecruitsDD(recruited_tagged_offspring_habprop, perc_UNOC_val, perc_APCL_val)
# 
# # New best est with DD
# recruits_per_egg_best_est_DD <- recruited_tagged_offspring_DDscaled/tagged_eggs_6cm
# 
# assignment_rate = n_offspring_matched/n_offspring_genotyped  # proportion of genotyped offspring that were assigned to parents in parentage analysis
# 
# # Recruits per egg just arriving to our sites (for an LRP for our sites that includes dispersal mortality in the egg-recruit surv)
# recruits_per_egg_oursites <- recruited_tagged_offspring_oursites/tagged_eggs_6cm
# recruits_per_egg_oursites_DD <- scaleTaggedRecruitsDD(recruited_tagged_offspring_oursites, perc_UNOC_val, perc_APCL_val)/tagged_eggs_6cm


# ########## Estimate egg-recruit survival, then scale up
# 
# # How many potential offspring (eggs) were produced by those potential (genotyped) parents ("tagged" adults b/c genetically marked)?
# tagged_eggs_6cm <- n_parents_genotyped*mean(LEP_different_sizes$LEP_6cm)
# tagged_eggs_3.5cm <- n_parents_genotyped*mean(LEP_different_sizes$LEP_3.5cm)
# # tagged_eggs_6cm <- n_parents_genotyped*LEP_6cm  # Use LEP from what size here? How to avoid double-counting if those parents mated together?
# # tagged_eggs_3.5cm <- n_parents_genotyped*LEP_3.5cm
# 
# 
# # Scale up by the proportion of site area we sampled over the time frame of finding parentage matches and prob of catching a fish
# recruited_tagged_offspring_oursites <- n_offspring_matched/((total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area*mean(prob_r))  # scale up by proportion of habitat sampled and probability of catching a fish
# 
# # Scale up by the proportion of the kernel area included within our sites (assuming patchiness of habitat outside our sites is the same as within)
# recruited_tagged_offspring_kernel <- recruited_tagged_offspring_oursites/prop_total_disp_area_sampled_best_est
# 
# # Scale up by the proportion of the sampling area (assumed to be the same as the greater region) that is habitat
# recruited_tagged_offspring_habprop <- recruited_tagged_offspring_kernel/prop_sampling_area_habitat
# 
# # Estimate survival from eggs-recruits by seeing how many "tagged" offspring we found out of eggs "tagged" parents produced
# recruits_per_egg_best_est <- recruited_tagged_offspring_habprop/tagged_eggs_6cm
# 
# # Scale up by DD
# recruited_tagged_offspring_DDscaled <- scaleTaggedRecruitsDD(recruited_tagged_offspring_habprop, perc_UNOC_val, perc_APCL_val)
# 
# # New best est with DD
# recruits_per_egg_best_est_DD <- recruited_tagged_offspring_DDscaled/tagged_eggs_6cm
# 
# assignment_rate = n_offspring_matched/n_offspring_genotyped  # proportion of genotyped offspring that were assigned to parents in parentage analysis
# 
# # Recruits per egg just arriving to our sites (for an LRP for our sites that includes dispersal mortality in the egg-recruit surv)
# recruits_per_egg_oursites <- recruited_tagged_offspring_oursites/tagged_eggs_6cm
# recruits_per_egg_oursites_DD <- scaleTaggedRecruitsDD(recruited_tagged_offspring_oursites, perc_UNOC_val, perc_APCL_val)/tagged_eggs_6cm


# ##### Estimate some of the "best estimates" of metrics/parameters
# breeding_size_mean <- mean(recap_first_female$size)
# prob_r_mean <- mean(prob_r)  # average value of prob r from each recap dive
# 
# ### Best-estimate survival parameters by site
# site_surv_best_est <- data.frame(site = no_space_sites_alpha, Sint = NA, Sl = best_fit_model_dfs$results$estimate[Phi_size_pos], stringsAsFactors = FALSE)
# # Cabatoan
# site_surv_best_est$Sint[1] = best_fit_model_dfs$results$estimate[1]
# # fill in rest of sites
# for(i in 2:length(no_space_sites_alpha)) {
#   site_surv_best_est$Sint[i] = best_fit_model_dfs$results$estimate[1] + best_fit_model_dfs$results$estimate[i]
# }
# 
# # LEP, starting at different sizes
# LEP_different_sizes <- site_surv_best_est %>%
#   mutate(LEP_mean_breeding_size = NA,
#          LEP_6cm = NA,
#          LEP_3.5cm = NA)
# 
# for(i in 1:length(LEP_different_sizes$site)) {
#   # Starting from mean transition to female size
#   LEP_different_sizes$LEP_mean_breeding_size[i] <- findLEP(min_size, max_size, n_bins, t_steps, LEP_different_sizes$Sint[i], LEP_different_sizes$Sl[i],
#                                                         s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts, 
#                                                         clutches_per_year_mean, breeding_size_mean, breeding_size_mean, start_recruit_sd, 
#                                                         eggs_slope_log, eggs_intercept_log, eyed_effect)
#   # Starting at tagging size (6cm)
#   LEP_different_sizes$LEP_6cm[i] <- findLEP(min_size, max_size, n_bins, t_steps, LEP_different_sizes$Sint[i], LEP_different_sizes$Sl[i],
#                      s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts, 
#                      clutches_per_year_mean, breeding_size_mean, 6, start_recruit_sd, 
#                      eggs_slope_log, eggs_intercept_log, eyed_effect)
#   
#   # Starting at fin-clip size (3.5cm)
#   LEP_different_sizes$LEP_3.5cm[i] <- findLEP(min_size, max_size, n_bins, t_steps, LEP_different_sizes$Sint[i], LEP_different_sizes$Sl[i],
#                        s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts, 
#                        clutches_per_year_mean, breeding_size_mean, 3.5, start_recruit_sd, 
#                        eggs_slope_log, eggs_intercept_log, eyed_effect)
# }

# # Test if there is LEP eviction
# max_Sint_test <- max(site_surv_best_est$Sint)
# Sl_test <- site_surv_best_est$Sl[1]
# Sl_test2 <- max(site_surv_param_sets[[9]]$Sl)
# 
# LEP_test_1 <- findLEP(min_size, max_size, n_bins, t_steps, max_Sint_test, Sl_test2, s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts,
#                       clutches_per_year_mean, breeding_size_mean, 6.0, start_recruit_sd, eggs_slope_log, eggs_intercept_log, eyed_effect)
# LEP_test_2 <- findLEP(min_size, 20, n_bins, t_steps, max_Sint_test, Sl_test2, s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts,
#                       clutches_per_year_mean, breeding_size_mean, 6.0, start_recruit_sd, eggs_slope_log, eggs_intercept_log, eyed_effect)
# LEP_test_3 <- findLEP(min_size, 25, n_bins, t_steps, max_Sint_test, Sl_test2, s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts,
#                       clutches_per_year_mean, breeding_size_mean, 6.0, start_recruit_sd, eggs_slope_log, eggs_intercept_log, eyed_effect)
# LEP_test_4 <- findLEP(min_size, 17, n_bins, t_steps, max_Sint_test, Sl_test2, s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts,
#                       clutches_per_year_mean, breeding_size_mean, 6.0, start_recruit_sd, eggs_slope_log, eggs_intercept_log, eyed_effect)


# # Starting from mean transition to female size
# LEP_breeding_size_mean <- findLEP(min_size, max_size, n_bins, t_steps, Sint_mean, Sl_mean,
#                                   s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts, 
#                                   clutches_per_year_mean, breeding_size_mean, breeding_size_mean, start_recruit_sd, 
#                                   eggs_slope_log, eggs_intercept_log, eyed_effect)
# 
# # Starting at tagging size (6cm)
# LEP_6cm <- findLEP(min_size, max_size, n_bins, t_steps, Sint_mean, Sl_mean,
#                                   s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts, 
#                                   clutches_per_year_mean, breeding_size_mean, 6, start_recruit_sd, 
#                                   eggs_slope_log, eggs_intercept_log, eyed_effect)
# 
# # Starting at fin-clip size (3.5cm)
# LEP_3.5cm <- findLEP(min_size, max_size, n_bins, t_steps, Sint_mean, Sl_mean,
#                                   s, Linf_growth_mean, k_growth_mean, mean_eggs_per_clutch_from_counts, 
#                                   clutches_per_year_mean, breeding_size_mean, 3.5, start_recruit_sd, 
#                                   eggs_slope_log, eggs_intercept_log, eyed_effect)

#################### Estimate survival from egg to recruit (method similar to Johnson et al. 2018): ####################

# ##### How much of the dispersal kernel area from each site did we actually sample? (So can scale up "tagged" recruits found to account for areas they might have gone that weren't in our sites)
# # Find the number of parents at each site (eventually, this parent file pull will go in Constants_database_common_functions). Just putting it here for now b/c going to use original 913, with same distribution as current parents
# all_parents_site <- all_parents_by_site %>%
#   group_by(site) %>%
#   summarize(nparents = n()) %>%
#   mutate(prop_parents = nparents/sum(nparents)) 
# 
# # Total dispersal kernel area (total parents*2 - total area dispersing north of site is 1 and south is 1 for each parent)
# total_parent_kernel_area = n_parents_genotyped*2
# 
# # Add in site info
# all_parents_site <- left_join(all_parents_site, site_width_info %>% select(site, site_geo_order, dist_to_N_edge_km, dist_to_S_edge_km), by = "site")
# 
# # Find area within sampling area to the north and south of each site
# all_parents_site <- all_parents_site %>%
#   mutate(disp_area_N_within_sites = NA,
#          disp_area_S_within_sites = NA)
# 
# for(i in 1:length(all_parents_site$site)) {
#   all_parents_site$disp_area_N_within_sites[i] = integrate(disp_allyears_d, 0, all_parents_site$dist_to_N_edge_km[i])$value
#   all_parents_site$disp_area_S_within_sites[i] = integrate(disp_allyears_d, 0, all_parents_site$dist_to_S_edge_km[i])$value
# }
# 
# # Find proportion of total area under dispersal kernel (where total area to INF is 2 (1 for each side)) covered within sample sites
# all_parents_site <- all_parents_site %>%
#   mutate(total_disp_area_within_sites = disp_area_N_within_sites + disp_area_S_within_sites,
#          prop_disp_area_within_sites = total_disp_area_within_sites/2,
#          total_parent_area_sampled = total_disp_area_within_sites*nparents)
# 
# all_parents_site_summarized <- all_parents_site %>%
#   summarize(total_parent_kernel_area = sum(nparents)*2,
#             sampled_parent_kernel_area = sum(total_parent_area_sampled),
#             prop_parent_kernel_area_sampled = sampled_parent_kernel_area/total_parent_kernel_area)
# 
# prop_total_disp_area_sampled_best_est <- all_parents_site_summarized$prop_parent_kernel_area_sampled
# 
# # How many potential offspring (eggs) were produced by those potential (genotyped) parents ("tagged" adults b/c genetically marked)?
# tagged_eggs_6cm <- n_parents_genotyped*mean(LEP_different_sizes$LEP_6cm)
# tagged_eggs_3.5cm <- n_parents_genotyped*mean(LEP_different_sizes$LEP_3.5cm)
# # tagged_eggs_6cm <- n_parents_genotyped*LEP_6cm  # Use LEP from what size here? How to avoid double-counting if those parents mated together?
# # tagged_eggs_3.5cm <- n_parents_genotyped*LEP_3.5cm
# 
# # Find the total prop habitat sampled over time - need to update this to have the other sites...
# #total_prop_hab_sampled_through_time <- (total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area*mean(prob_r)  # scale up by proportion of habitat sampled and probability of catching a fish
# total_prop_hab_sampled_through_time <- (total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area  # scale up by proportion of habitat sampled 
# 
# # Scale up by the proportion of site area we sampled over the time frame of finding parentage matches and prob of catching a fish
# recruited_tagged_offspring_oursites <- n_offspring_matched/((total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area*mean(prob_r))  # scale up by proportion of habitat sampled and probability of catching a fish
# 
# # Scale up by the proportion of the kernel area included within our sites (assuming patchiness of habitat outside our sites is the same as within)
# recruited_tagged_offspring_kernel <- recruited_tagged_offspring_oursites/prop_total_disp_area_sampled_best_est
# 
# # Scale up by the proportion of the sampling area (assumed to be the same as the greater region) that is habitat
# recruited_tagged_offspring_habprop <- recruited_tagged_offspring_kernel/prop_sampling_area_habitat
# 
# # Estimate survival from eggs-recruits by seeing how many "tagged" offspring we found out of eggs "tagged" parents produced
# recruits_per_egg_best_est <- recruited_tagged_offspring_habprop/tagged_eggs_6cm
# 
# # Scale up by DD
# recruited_tagged_offspring_DDscaled <- scaleTaggedRecruitsDD(recruited_tagged_offspring_habprop, perc_UNOC_val, perc_APCL_val)
# 
# # New best est with DD
# recruits_per_egg_best_est_DD <- recruited_tagged_offspring_DDscaled/tagged_eggs_6cm
# 
# assignment_rate = n_offspring_matched/n_offspring_genotyped  # proportion of genotyped offspring that were assigned to parents in parentage analysis
# 
# # Recruits per egg just arriving to our sites (for an LRP for our sites that includes dispersal mortality in the egg-recruit surv)
# recruits_per_egg_oursites <- recruited_tagged_offspring_oursites/tagged_eggs_6cm
# recruits_per_egg_oursites_DD <- scaleTaggedRecruitsDD(recruited_tagged_offspring_oursites, perc_UNOC_val, perc_APCL_val)/tagged_eggs_6cm
# 
# ##### How big are the offspring? What size should we use for recruits?
# mean_sampled_offspring_size <- mean(all_offspring$size, na.rm = TRUE)  # might be duplicate observations of fish in here, but I don't think so - should check
# #mean_sampled_offspring_size <- mean(n_offspring_genotypes_df$size, rm.na = TRUE)  # this is just one of the obs of each of these fish... not sure how many duplicates there are, should really check...
# 
# start_recruit_size = mean_sampled_offspring_size
# 
# #################### Find connectivity matrix: ####################
# 
# #################### Find metrics for "best estimate" of the various parameters: ####################
# # Put best-estimate parameters into one dataframe
# # start at mean size of actual offspring collected (about 4.45)
# param_best_est_mean_collected_offspring <- data.frame(t_steps = n_tsteps) %>%
#   mutate(min_size = min_size, max_size = max_size, n_bins = n_bins,  # LEP-IPM matrix
#          eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year_mean = clutches_per_year_mean,  # average fecundity info
#          egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect =  eyed_effect,  # size-dependent fecundity info
#          start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
#          k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean, 
#          #Sl = Sl_mean, Sint = Sint_mean,
#          breeding_size = breeding_size_mean, #recruits_per_egg = recruits_per_egg_best_est,
#          k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
#          prob_r = prob_r_mean, offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
#          total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
#          perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 
# 
# # # with DD scaling
# # param_best_est_mean_collected_offspring_DD <- data.frame(t_steps = n_tsteps) %>%
# #   mutate(min_size = min_size, max_size = max_size, n_bins = n_bins,  # LEP-IPM matrix
# #          eggs_per_clutch = mean_eggs_per_clutch_from_counts, clutches_per_year_mean = clutches_per_year_mean,  # average fecundity info
# #          egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect =  eyed_effect,  # size-dependent fecundity info
# #          start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
# #          k_growth = k_growth_mean, s = s, Sl = Sl_mean, Linf = Linf_growth_mean, Sint = Sint_mean,
# #          breeding_size = breeding_size_mean, recruits_per_egg = recruits_per_egg_best_est,
# #          k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
# #          prob_r = prob_r_mean, offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
# #          total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
# #          perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, DD = TRUE) 
# 
# # Calculate the metrics for the best estimates
# best_est_metrics_mean_offspring <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info, site_vec_order, "FALSE")
# best_est_metrics_mean_offspring_DD <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info, site_vec_order, "TRUE")
# 
# # # Calculate the metrics for the best estimates
# # best_est_metrics_mean_offspring <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info, site_vec_order$site_name, "FALSE")
# # best_est_metrics_mean_offspring_DD <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info, site_vec_order$site_name, "TRUE")
# 
# LEP_best_est <- mean(best_est_metrics_mean_offspring$LEP_by_site$LEP)  # how to deal with site variation?
# LEP_best_est_min <- min(best_est_metrics_mean_offspring$LEP_by_site$LEP)
# LEP_best_est_max <- max(best_est_metrics_mean_offspring$LEP_by_site$LEP)
# LEP_R_best_est_min <- min(best_est_metrics_mean_offspring$LEP_by_site$LEP_R)
# LEP_R_best_est_max <- max(best_est_metrics_mean_offspring$LEP_by_site$LEP_R)
# LEP_R_best_est <- best_est_metrics_mean_offspring$LEP_R_mean
# NP_best_est <- best_est_metrics_mean_offspring$NP
# SP_best_est <- best_est_metrics_mean_offspring$SP
# LEP_R_local_best_est <- best_est_metrics_mean_offspring$LEP_R_local_mean
# 
# LEP_best_est_DD <- mean(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP)  # how to deal with, show site variation?
# LEP_best_est_min_DD <- min(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP)
# LEP_best_est_max_DD <- max(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP)
# LEP_R_best_est_min_DD <- min(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP_R)
# LEP_R_best_est_max_DD <- max(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP_R)
# LEP_R_best_est_DD <- best_est_metrics_mean_offspring_DD$LEP_R_mean
# NP_best_est_DD <- best_est_metrics_mean_offspring_DD$NP
# SP_best_est_DD <- best_est_metrics_mean_offspring_DD$SP
# LEP_R_local_best_est_DD <- best_est_metrics_mean_offspring_DD$LEP_R_local_mean
# 
# LEP_R_oursites <- recruits_per_egg_oursites*LEP_best_est
# LEP_R_oursites_DD <- recruits_per_egg_oursites_DD*LEP_best_est
# 
# # Probability of capturing a fish
# prob_r_beta_params = findBetaDistParams(mean(prob_r), var(prob_r))  # find beta distribution parameters for prob r distrubtion from normal mean and variance
# prob_r_set_fodder = rbeta(n_runs, prob_r_beta_params$alpha, prob_r_beta_params$beta, 0)  # should the non-centrality parameter be 0?
# 
# prob_r_set_normal = rnorm(n_runs, mean = prob_r_mean, sd = sd(prob_r))  # not sure where this is coming in right now... would be in scaling up pops...
# prob_r_set_fromdata = sample(prob_r, n_runs, replace = TRUE)  # just sample from the 14 calculated prob_r values
# 
# # Get rid of any values lower than the lowest observed value, then re-sample from that truncated vector...
# prob_r_set_truncated = prob_r_set_fodder[prob_r_set_fodder >= min(prob_r)]  # removes about 100 obs (down to 898 in one case)
# prob_r_set = sample(prob_r_set_truncated, n_runs, replace = TRUE)
# 
# # Compare the two versions (plot down in plot section)
# prob_r_comp <- data.frame(distribution = c(rep("beta", n_runs), rep("truncated beta", n_runs), rep("normal", n_runs), rep("sample from data", n_runs)),
#                           values = c(prob_r_set_fodder, prob_r_set, prob_r_set_normal, prob_r_set_fromdata), stringsAsFactor = FALSE)
# 
# # Uncertainty in recruits-per-egg - two ways - NEED TO TIDY THIS
# # Way one: uncertainty in the number of offspring that get matched through parentage analysis
# n_offspring_parentage_set <- rbinom(n_runs, n_offspring_genotyped, assignment_rate)  # number of assigned offspring using just uncertainty in binomial (assigned/not), for recruits-per-egg est

# # I think the rest of this is now incorporated into the calcMetrics function
# scaled_tagged_recruits_set1 <- scaleTaggedRecruits(n_offspring_parentage_set, (total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area, prob_r_mean, prop_total_disp_area_sampled_best_est)
# # recruits_per_egg_set1 <- findRecruitsPerTaggedEgg(scaled_tagged_recruits_set1, LEP_6cm)  # this was the problem!! Accidentally put that those tagged recruits came from one tagged parent's egg output, not all the tagged parents
# recruits_per_egg_set1 <- findRecruitsPerTaggedEgg(scaled_tagged_recruits_set1, tagged_eggs_6cm)
# 
# # Way two: uncertainty in probability of capturing a fish (so number of offspring we find from the tagged parents) - is this worth it?
# scaled_tagged_recruits_set2 <- scaleTaggedRecruits(n_offspring_matched, (total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area, prob_r_set, prop_total_disp_area_sampled_best_est)
# # recruits_per_egg_set2 <- findRecruitsPerTaggedEgg(scaled_tagged_recruits_set2, LEP_6cm)
# recruits_per_egg_set2 <- findRecruitsPerTaggedEgg(scaled_tagged_recruits_set2, tagged_eggs_6cm)
# 
# # Way three: both together
# scaled_tagged_recruits_set3 <- scaleTaggedRecruits(n_offspring_parentage_set, (total_area_sampled_through_time %>% filter(method == "metal tags", time_frame == "2012-2018"))$total_prop_hab_sampled_area, prob_r_set, prop_total_disp_area_sampled_best_est)
# # recruits_per_egg_set3 <- findRecruitsPerTaggedEgg(scaled_tagged_recruits_set3, LEP_6cm)
# recruits_per_egg_set3 <- findRecruitsPerTaggedEgg(scaled_tagged_recruits_set3, tagged_eggs_6cm)
# 
# # Compare the uncertainty in recruits-per-egg (rbinom for assigned offspring, prob_r, both), for plot later
# recruits_per_egg_uncertainty <- data.frame(recruits_per_egg = c(recruits_per_egg_set1, recruits_per_egg_set2, recruits_per_egg_set3),
#                                            method = c(rep("n assigned offspring", n_runs), rep("prob r", n_runs), rep("both", n_runs)),
#                                            stringsAsFactors = FALSE)
# 
# # # Pick a recruits_per_egg_set to use
# # recruits_per_egg_set <- recruits_per_egg_set3  # using the set that includes both binom for genotyped offspring and draws in prob_r

# # Survival (for LEP)
# Sint_mean = survival_output$estimate[Phi_int_pos]
# Sint_se = survival_output$se[Phi_int_pos]
# Sl_mean = survival_output$estimate[Phi_size_pos]
# Sl_se = survival_output$se[Phi_size_pos]


# eall_mean.Phi.size.p.size.plus.dist.results <- as.data.frame(eall_mean.Phi.size.p.size.plus.dist$results$beta)
# Sint_mean = eall_mean.Phi.size.p.size.plus.dist.results$estimate[1]  # survival intercept (on logit scale)
# Sl_mean = eall_mean.Phi.size.p.size.plus.dist.results$estimate[2]  # survival slope (on logit scale)
# Sint_se = eall_mean.Phi.size.p.size.plus.dist.results$se[1]  # for now using SE, should really use SD...
# Sl_se = eall_mean.Phi.size.p.size.plus.dist.results$se[2]  # for now using SE, should really use SD...


