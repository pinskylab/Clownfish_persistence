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

# Load output from abundance trend (for plotting) (from Code/TimeSeriesPersistence.R)
load(file = here::here("Data/Script_outputs", "site_trends_all.RData"))
load(file = here::here("Data/Script_outputs", "site_trends_time.RData"))

##### Set parameters (for running IPM, for calculating connectivity, for uncertainty runs, etc.)

# Fecundity (for LEP) 
size_fecundity_model = length_count8llEA  # best fit model from Adam fecundity work, both length and eggs on log scale
eggs_intercept_log = size_fecundity_model$coefficients[1]  # on log-scale
eggs_slope_log = size_fecundity_model$coefficients[2]
eyed_effect = size_fecundity_model$coefficients[3]  # effect of eggs having visible eyes

# Set params for IPM structure for calculating LEP
n_bins = 100  # number of size bins
n_tsteps = 100  # number of time steps (large number, approximating infinity)
start_recruit_sd = 0.1  # standard deviation of starting recruit size (to spread initial fish across size bins)

# Parameters for incorporating growth
size_around_mean_L1_size = 0.1  # how much buffer to give for finding fish of that size (so searching for the sd of L2 of fish within 0.2 cm range)

# Parent size (for "tagged" eggs produced by parents)
parent_size = 6.0  # cm

# Parameters for what-ifs
n_sites = length(site_vec_order$site_name)  # number of sites for habitat, region width, and larval navigation simulations (same as actual number of sites)
perc_hab_vals <- seq(from=0.05,to=1.0, by=0.05)  # values of percent habitat to test for sensitivity to amount of habitat in region

# Parameters for analyzing and interpreting results
NP_decimal_points <- 2  # numbers past decimal point for NP range
quantile_lower_index = 25  # lower index for 95% quantile reporting
quantile_upper_index = 975  # upper index for 95% quantile reporting

#################### Functions: ####################
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

# Find beta distribution parameters from mean and variance (for prob_r uncertainty)
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

# Find eggs by fish size (eyed eggs) 
findEggs = function(fish_size, egg_size_intercept, egg_size_slope, eyed_effect) {
  count_logged = egg_size_intercept + egg_size_slope*log(fish_size) + eyed_effect
  raw_eggs = exp(count_logged)
  return(raw_eggs)
}

# Find scaled number of tagged recruits we would expect to have found if we sampled the whole area and caught all the fish - if we caught all recruits produced by our patches
scaleTaggedRecruits = function(offspring_assigned, total_prop_hab_sampled, prob_capture, prop_total_disp_area_sampled, prop_hab) {
  rto <- offspring_assigned/(total_prop_hab_sampled*prob_capture*prop_total_disp_area_sampled*prop_hab)
  return(rto)
}

# Just consider recruits that return to our sites (used in local replacement)
scaleTaggedRecruits_local = function(offspring_assigned, total_prop_hab_sampled, prob_capture) {
  rto <- offspring_assigned/(total_prop_hab_sampled*prob_capture)
  return(rto)
}

# Find scaled number of tagged recruits by open habitat, density dependence
scaleTaggedRecruitsDD = function(recruited_tagged_offspring, perc_UNOC, perc_APCL) {
  recruited_tagged_offspring_total_DD = ((perc_APCL+perc_UNOC)/perc_UNOC)*recruited_tagged_offspring  # scale as if currently-occupied APCL habitat is open
  return(recruited_tagged_offspring_total_DD)
}

# Find LEP - use an IPM to find LEP 
findLEP = function(min_size, max_size, n_bins, t_steps, Sint, Sl, s, Linf, k_growth, clutches_per_year,
                   breeding_size, start_recruit_size, start_recruit_sd, egg_size_slope, egg_size_intercept, eyed_effect) {

  # Create vector of lengths
  lengths_vec = seq(min_size, max_size, length.out = n_bins)
  dx = diff(lengths_vec[1:2])

  # Make matrix
  xmat = matrix(rep(lengths_vec, n_bins), nrow=n_bins, ncol=n_bins, byrow=TRUE)  # 100 rows of lengths_vec
  ymat = t(xmat)

  # Survival probability (based on size)
  S = logit_recip(Sint + lengths_vec*Sl)  # survival probability (based on size), where Sint varies by site
  S[is.na(S)] = 1  # replace the NaNs (coming from Inf/Inf in logit_recip) with 1 (shouldn't occur anymore, was sometimes an issue when survival uncertainties were wide)
  Smat = matrix(rep(S,n_bins), nrow=n_bins, ncol=n_bins, byrow = TRUE)  # survival part of kernel

  # Growth probability
  Ls = Linf - (Linf - xmat)*exp(-k_growth)  # expected length (z') in next time step based on current length (z)
  Lmat = dnorm(ymat,Ls,s)  # turn that into a matrix of probability densities

  # Make the kernel with the survival + growth inputs
  K = Lmat*Smat
  K = dx*K

  # Iterate with 1 subadult to start  # 3.5, 0.1
  N0 = dnorm(lengths_vec, start_recruit_size, start_recruit_sd)  # create a size distribution for one recruit, spread across size bins

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
  Fvec[which(lengths_vec > breeding_size)] = 1  # switch to 1 for sizes larger than size of transition to female
  
  # Vector of eyed eggs produced per clutch, dependent on female size
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

  return(LEP)
}

# Run through a metric calculation with one set of parameters
calcMetrics <- function(param_set, site_based_surv_sets, sites_and_dists, site_vec, DD) {  # param_set is set of parameters, site_based_surv_sets is survival parameters, sites_and_dists is data frame of sites and distances beween then, DD is TRUE if compensating for density dependence/FALSE if not

  sites <- site_vec$site_name

  # Define function with the right parameters
  disp_allyears <- function(d) {  # generalization of equation for p(d) in eqn. 6c in Bode et al. 2018
    z = exp(param_set$k_connectivity)
    z_front = (z*param_set$theta_connectivity)/gamma(1/param_set$theta_connectivity)
    disp = (z_front/2)*exp(-(z*d)^(param_set$theta_connectivity))
    return(disp)
  }
  
  # Create connectivity matrix, find probability of dispersing among sites
  Cmat <- sites_and_dists %>% select(org_site, dest_site, d1_km, d2_km, org_alpha_order, org_geo_order, dest_alpha_order, dest_geo_order)
  for(i in 1:length(Cmat$org_site)) {
    Cmat$prob_disp[i] <- integrate(disp_allyears, Cmat$d1_km[i], Cmat$d2_km[i])$value
    # Double probability of dispersal for self-self to include both directions (from center to patch edge to both N and S)
    if(Cmat$org_site[i] == Cmat$dest_site[i]) {
      Cmat$prob_disp[i] <- Cmat$prob_disp[i]*2
    }
  }

  # Find LEP (in terms of eggs) by site
  LEP_by_site <- data.frame(site = (site_vec %>% arrange(geo_order))$site_name, LEP = NA, LEP_parents = NA)
  
  for(i in 1:length(LEP_by_site$site)) {
    
    # Select the survival parameters for that site
    Sint_set <- (site_based_surv_sets %>% filter(site == LEP_by_site$site[i]))$Sint
    Sl_set <- (site_based_surv_sets %>% filter(site == LEP_by_site$site[i]))$Sl

    # Then calculate LEP
    LEP_by_site$LEP[i] = findLEP(param_set$min_size, param_set$max_size, param_set$n_bins, param_set$t_steps, Sint_set, Sl_set,
                                 param_set$s, param_set$Linf, param_set$k_growth, param_set$clutches_per_year,
                                 param_set$breeding_size, param_set$start_recruit_size, param_set$start_recruit_sd,
                                 param_set$egg_size_slope, param_set$egg_size_intercept, param_set$eyed_effect)
    LEP_by_site$LEP_parents[i] = findLEP(parent_size, param_set$max_size, param_set$n_bins, param_set$t_steps, Sint_set, Sl_set,
                                         param_set$s, param_set$Linf, param_set$k_growth, param_set$clutches_per_year,
                                         param_set$breeding_size, parent_size, param_set$start_recruit_sd,
                                         param_set$egg_size_slope, param_set$egg_size_intercept, param_set$eyed_effect)
  }

  # Find egg-recruit survival (recruits/egg)
  tagged_recruits_val = scaleTaggedRecruits(param_set$offspring_assigned_to_parents, param_set$total_prop_hab_sampled,
                                            param_set$prob_r, param_set$prop_total_disp_area_sampled, param_set$prop_hab)

  tagged_recruits_local = scaleTaggedRecruits_local(param_set$offspring_assigned_to_parents, param_set$total_prop_hab_sampled,
                                                    param_set$prob_r)
 
  # Find the mean of LEP starting from parent size across sites
 LEP_parents_mean <- mean(LEP_by_site$LEP_parents)

 # Use mean LEPparents to estimate the number of tagged eggs, rather than doing it by site
 tagged_eggs_val = param_set$n_parents*LEP_parents_mean  

  if(DD == "TRUE"){  # if compensating for density dependence
    tagged_recruits_val_scaled = scaleTaggedRecruitsDD(tagged_recruits_val, param_set$perc_UNOC, param_set$perc_APCL)
    tagged_recruits_local_scaled = scaleTaggedRecruitsDD(tagged_recruits_local, param_set$perc_UNOC, param_set$perc_APCL)
  } else {
    tagged_recruits_val_scaled = tagged_recruits_val
    tagged_recruits_local_scaled = tagged_recruits_local
  }
 
  recruits_per_egg = tagged_recruits_val_scaled/tagged_eggs_val

  # Find LRP (lifetime recruit production)
  LEP_by_site <- LEP_by_site %>%
    mutate(LEP_R = LEP*recruits_per_egg,
           LEP_R_local = LEP*(tagged_recruits_local_scaled/tagged_eggs_val))  # use only recruits that returned to local metapop to calculate local replacement

  # Find mean LRP across sites
  LEP_R_mean <- mean(LEP_by_site$LEP_R)

  # Find mean local replacement (LR) across sites
  LEP_R_local_mean = mean(LEP_by_site$LEP)*(tagged_recruits_local_scaled/tagged_eggs_val)

  # Put connectivity matrix into a matrix form
  conn_matrix <- matrix(NA, ncol=max(Cmat$org_geo_order), nrow=max(Cmat$org_geo_order))
  for(i in 1:length(Cmat$org_site)) {
    column = Cmat$org_geo_order[i]  # column is origin
    row = Cmat$dest_geo_order[i]  # row is destination
    conn_matrix[row, column] = Cmat$prob_disp[i]
  }

  # Make realized connectivity matrix, both in the matrix form and dataframe form
  conn_matrixR = conn_matrix  # matrix form (for eigenvalues)
  for(i in 1:length(LEP_by_site$site)) {
    conn_matrixR[,i] <- conn_matrix[,i]*LEP_by_site$LEP_R[i]
  }
  Cmat <- left_join(Cmat, LEP_by_site %>% select(site, LEP_R), by=c("org_site" = "site"))  # data frame form (for plotting)
  Cmat <- Cmat %>%
    mutate(prob_disp_R = prob_disp*LEP_R)

  # Assess network persistence
  eig_cR = eigen(conn_matrixR)

  # Pull out self-persistence values (diagonals of realized connectivity matrix)
  SP_values = data.frame(site = sites, stringsAsFactors = FALSE) %>%
    mutate(SP_value = NA, org_geo_order = NA)
  for(i in 1:length(sites)) {
    site_val = sites[i]
    SP_values$SP_value[i] = (Cmat %>% filter(org_site == site_val & dest_site == site_val))$prob_disp_R
    SP_values$org_geo_order[i] = (Cmat %>% filter(org_site == site_val & dest_site == site_val))$org_geo_order
  }

  # Put outputs together into one list
  out = list(NP = Re(eig_cR$values[1]), SP = SP_values, LEP_by_site = LEP_by_site, LEP_R_mean = LEP_R_mean, LEP_R_local_mean = LEP_R_local_mean, recruits_per_egg = recruits_per_egg,
             conn_matrix = conn_matrix, conn_matrixR = conn_matrixR, Cmat = Cmat)
  
  return(out)
}

# Calculate metrics across many runs (this is slow, if have time, should really try to write this without a for loop....)
calcMetricsAcrossRuns <- function(n_runs, param_sets, site_based_surv_sets, site_dist_info, site_vec_order, set_name, DD) {  # set name is to identify what kind of uncertainty is included in the run

  LEP_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_min_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_max_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_R_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP_R', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_R_min_out_df <-  data.frame(value = rep(NA, n_runs), metric = 'LEP_R_max', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_R_max_out_df <-  data.frame(value = rep(NA, n_runs), metric = 'LEP_R_max', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  LEP_R_local_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP_R_local', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  RperE_out_df <- data.frame(value = rep(NA, n_runs), metric = 'recruits per egg', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)
  NP_out_df <- data.frame(value = rep(NA, n_runs), metric = 'NP', run = seq(1:n_runs), uncertainty_type = set_name, stringsAsFactors = FALSE)

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

    print(i)  # Just to keep track of how many are left...

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

# Find width of sites and distance between them for different percent habitats (and sampling region lengths), sites start right at the edge of the region like in actual configuration
find_site_locations <- function(nsites, perc_hab_val, total_region) {
  
  site_info <- data.frame(org_site = seq(from=1, to=nsites, by=1))
  site_width <- (total_region*perc_hab_val)/nsites
  hab_break <- (total_region*(1-perc_hab_val))/(nsites-1)  # width of break between habitat patches
  
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
  
  # Find distances between each set of sites
  for(i in 1:length(out_df$org_site)) {
    out_df$dist_mid_to_N[i] = abs(out_df$center_org[i] - out_df$N_edge_dest[i])
    out_df$dist_mid_to_S[i] = abs(out_df$center_org[i] - out_df$S_edge_dest[i])
  }
  
  # Find d1 and d2 (min and max distances to each site)
  for(i in 1:length(out_df$org_site)) {
    out_df$d1_km[i] = min(c(out_df$dist_mid_to_S[i], out_df$dist_mid_to_N[i]))
    out_df$d2_km[i] = max(c(out_df$dist_mid_to_S[i], out_df$dist_mid_to_N[i]))
  }
  
  # Fix the self-selfs to make sure d1 is 0, d2 is half the width of the site
  out_df <- out_df %>%
    mutate(d1_km = case_when(org_site != dest_site ~ d1_km,
                             org_site == dest_site ~ 0),  # min distance is 0 for self-self distances
           d2_km = case_when(org_site != dest_site ~ d2_km,
                             org_site == dest_site ~ width_dest/2))  # max distance is half of site width for self-self distances 
  
  # Add the random other columns the calcMetrics function needs
  out_df <- out_df %>% 
    mutate(org_alpha_order = org_site,
           org_geo_order = org_site,
           dest_alpha_order = dest_site,
           dest_geo_order = dest_site)
  
  return(out_df)
}

# Find proportion habitat when including buffers for larval navigation, includes re-estimating number of scaled offspring, where scaling for loss to dispersal to non-habitat is altered for the amount of navigation
findPropHabWithLarvalNav <- function(nav_buffer, site_buffer_info_df, site_width_info_df, sampling_range) {  # nav_buffer = distance larvae can navigate in m, site_buffer_info_df = site_buffer_info, site_width_info_df = site_width_info, sampling_range = total_range_of_sampling_area
  # total sum of site widths without buffer
  sum_of_site_widths <- sum(site_width_info_df$width_m)
  buffer_to_add_N <- rep(NA, length(site_width_info_df$site))
  buffer_to_add_S <- rep(NA, length(site_width_info_df$site))
  
  # go through sites, check if buffer is less than max amount (where max is half the distance between sites, to avoid patch shadow overlap), then add either buffer or max
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
  
  total_hab <- sum_of_site_widths + sum(buffer_to_add_N) + sum(buffer_to_add_S)  # total "habitat" including patches and their buffers
  
  prop_hab <- total_hab/total_range_of_sampling_area
  
  return(prop_hab)
}

# Make new param set for new proportion hab including larval navigation buffers
makeParamSetWithSampledHab <- function(param_set, Ps) {  # param_set is param_best_est_mean_collected_offspring, Ps is prop_hab from findPropHabWithLarvalNav
  param_set_out <- param_set
  param_set_out$prop_hab <- Ps
  return(param_set_out)
}

# Re-estimate the distances data frames including navigation buffer
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
         prop_disp_area_within_sites = total_disp_area_within_sites,  # now, dispersal kernel normalized to 0.5 so total to N and S to infinity in each direction combined is 1
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

##### Survival by site (here, using median site survival for the three sites that weren't resampled and don't have survival estimates)

# Set up data frame for best-estimate survival parameters by site, for sites that were revisisted and have survival estimates
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


#################### Estimate point estimates of metrics ####################
min_size = 3.5  # start min size at size of recruit

# Put best-estimate intput parameters into one data frame (where starting recruit size is mean size of actual offspring collected ~ 4.38)
param_best_est_mean_collected_offspring <- data.frame(t_steps = n_tsteps) %>%
  mutate(min_size = min_size, max_size = max_size, n_bins = n_bins,  # for LEP-IPM matrix
         clutches_per_year_mean = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect =  eyed_effect,  # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,  # growth info
         breeding_size = breeding_size_mean,  # size of transition to female 
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time,  # Pc, Ph (for scaling up recruits)
         prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est, prop_hab = prop_sampling_area_habitat,  # Ps, Pd (for scaling up recruits)
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,  # number of offspring and parents genotyped
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val)  # proportion APCL-occupied anemones and unoccupied anemones for density dependence compensation

# Calculate the metrics for the best estimates, without and with (_DD) density-dependence compensated for
best_est_metrics_mean_offspring <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info, site_vec_order, "FALSE")
best_est_metrics_mean_offspring_DD <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est, site_dist_info, site_vec_order, "TRUE")

#################### Generate sets of parameters for uncertainty: ####################

##### Growth 
growth_set_intercept_est_mean <- mean(growth_info_estimate$intercept_est)  # mean of the intercept estimates from the fits with different growth pairs for fish with multiple recaptures
growth_set_intercept_se_mean <- mean(growth_info_estimate$intercept_se)  # mean of intercept se estimates
growth_set_slope_est_mean <- mean(growth_info_estimate$slope_est)  # mean of slope estimates from the fits with different growth pairs randomly chosen for fish with multiple recaptures
growth_set_slope_se_mean <- mean(growth_info_estimate$slope_se)  # mean of slope se estimates

# Select 1000 slopes and intercepts using the mean and se from the fits, then find k and Linf from that
growth_set_params <- data.frame(slope_est = runif(n_runs, growth_set_slope_est_mean-1.96*growth_set_slope_se_mean, growth_set_slope_est_mean+1.96*growth_set_slope_se_mean),
                                 intercept_est = runif(n_runs, growth_set_intercept_est_mean-1.96*growth_set_intercept_se_mean, growth_set_intercept_est_mean+1.96*growth_set_intercept_se_mean)) %>%
  mutate(k_est = -log(slope_est),
         Linf_est = intercept_est/(1-slope_est))

# Find mean 95% upper and lower limits 
growth_slope_est_lower = growth_set_slope_est_mean-1.96*growth_set_slope_se_mean
growth_slope_est_upper = growth_set_slope_est_mean+1.96*growth_set_slope_se_mean
growth_intercept_est_lower = growth_set_intercept_est_mean-1.96*growth_set_intercept_se_mean
growth_intercept_est_upper = growth_set_intercept_est_mean+1.96*growth_set_intercept_se_mean
k_growth_upper = -log(growth_slope_est_lower)
k_growth_lower = -log(growth_slope_est_upper)
Linf_growth_lower = growth_intercept_est_lower/(1-growth_slope_est_lower)
Linf_growth_upper = growth_intercept_est_upper/(1-growth_slope_est_upper)

##### Dispersal
# Select values from the 95% CI, weighted by log-likelihood
min_LL <- min(k_theta_allyear_95CI_values$log_like)
dispersal_sample_set <- k_theta_allyear_95CI_values %>%
  mutate(deviation = min_LL/log_like) %>%
  sample_n(size = n_runs, weight = deviation, replacement = T)

k_connectivity_set <- dispersal_sample_set$k_eval
theta_connectivity_set <- dispersal_sample_set$theta_eval

# Find 2.5th and 97.5th quantiles
k_connectivity_lower <- sort(k_connectivity_set)[quantile_lower_index]
k_connectivity_upper <- sort(k_connectivity_set)[quantile_upper_index]
theta_connectivity_lower <- sort(theta_connectivity_set)[quantile_lower_index]
theta_connectivity_upper <- sort(theta_connectivity_set)[quantile_upper_index]

##### Create survival parameter sets, sampling from within 95% confidence bounds for Cabatoan intercept, addition for each site, and size effect 
site_surv_param_sets <- list()

# Cabatoan first
site_surv_df <- data.frame(Sint_C = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[1], max=best_fit_model_dfs$results$ucl[1]),  # Cabatoan intercept
                           Sint_site = 0,  # site-specific addition to Catatoan intercept
                           Sl = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[Phi_size_pos], max=best_fit_model_dfs$results$ucl[Phi_size_pos])) %>%  # size effect
  mutate(Sint = Sint_C + Sint_site)  # add together to find site-specific intercept
site_surv_param_sets[[1]] <- site_surv_df

# then rest of sites where survival was estimated
for(i in 2:length(no_space_sites_revisited)) {
  site_surv_df <- data.frame(Sint_C = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[1], max=best_fit_model_dfs$results$ucl[1]), 
                             Sint_site = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[i], max=best_fit_model_dfs$results$ucl[i]), 
                             Sl = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[Phi_size_pos], max=best_fit_model_dfs$results$ucl[Phi_size_pos])) %>%
    mutate(Sint = Sint_C + Sint_site)
  site_surv_param_sets[[i]] <- site_surv_df
}

# then for the three that use median survival (Caridad Proper, Sitio Lonas, Sitio Tugas)
median_site_Sint_lcl = best_fit_model_dfs$results$lcl[median_site_pos]  # median site is Elementary School
median_site_Sint_ucl = best_fit_model_dfs$results$ucl[median_site_pos]

for(i in length(no_space_sites_revisited)+1:length(site_vec_order$site_name)) {
  site_surv_df <- data.frame(Sint_C = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[1], max=best_fit_model_dfs$results$ucl[1]), 
                             Sint_site = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[median_site_pos], max=best_fit_model_dfs$results$ucl[median_site_pos]), 
                             Sl = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[Phi_size_pos], max=best_fit_model_dfs$results$ucl[Phi_size_pos])) %>%
    mutate(Sint = Sint_C + Sint_site)
  site_surv_param_sets[[i]] <- site_surv_df
}

# Make a list with median Sint for all sites to use for what-if simulations
site_surv_median_Sint_param_sets <- list()
for(i in 1:length(site_vec_order$site_name)) {
  site_surv_df <- data.frame(Sint_C = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[1], max=best_fit_model_dfs$results$ucl[1]), 
                             Sint_site = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[median_site_pos], max=best_fit_model_dfs$results$ucl[median_site_pos]), 
                             Sl = runif(n=n_runs, min=best_fit_model_dfs$results$lcl[Phi_size_pos], max=best_fit_model_dfs$results$ucl[Phi_size_pos])) %>%
    mutate(Sint = Sint_C + Sint_site)
  site_surv_median_Sint_param_sets[[i]] <- site_surv_df
}

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

##### Breeding size (transition to female)
breeding_size_set = sample(recap_first_female$size, n_runs, replace = TRUE)  # transition to female size, pulled from first-observed sizes at female for fish caught as both male and female

##### Probability of capturing a fish (Pc)
prob_r_beta_params = findBetaDistParams(mean(prob_r), var(prob_r))  # find beta distribution parameters for prob r distrubtion from normal mean and variance
prob_r_set_fodder = rbeta(n_runs, prob_r_beta_params$alpha, prob_r_beta_params$beta, 0)  # generate a set of values from a beta distribution
prob_r_set_truncated = prob_r_set_fodder[prob_r_set_fodder >= min(prob_r)]  # get rid of any values lower than the lowest observed value of prob_r, then re-sample from that truncated vector... removes about 100 obs 
prob_r_set = sample(prob_r_set_truncated, n_runs, replace = TRUE)  # sample from that truncated set

##### Recruit size 
start_recruit_size_set <- runif(n_runs, min = 3.5, max = 6.0)  # pull from range of offspring sizes

###### Uncertainty in the number of offspring that get matched through parentage analysis (feeds into uncertainty in recruits-per-egg)
n_offspring_parentage_set <- rbinom(n_runs, n_offspring_genotyped, assignment_rate)  # number of assigned offspring using just uncertainty in binomial (assigned/not)

#################### Generate sets of parameters for uncertainty and do uncertainty runs: ####################
### Make parameter sets with different kinds of uncertainty included 

# Uncertainty in start recruit size only
param_set_start_recruit <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = start_recruit_size_set, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,
         breeding_size = breeding_size_mean, 
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in growth only (both Linf and k)
param_set_growth <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = growth_set_params$k_est, s=s, Linf = growth_set_params$Linf_est,
         breeding_size = breeding_size_mean, 
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in survival only (survival uncertainty enters the metrics function in a separate data frame)
param_set_survival <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean, 
         breeding_size = breeding_size_mean, 
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in breeding size only
param_set_breeding_size <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,
         breeding_size = breeding_size_set, 
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in offspring assigned to parents (affects recruits-per-egg)
param_set_offspring_assigned <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,
         breeding_size = breeding_size_mean,
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_parentage_set, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in probability catching a fish (affects recruits-per-egg)
param_set_prob_r <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean, 
         breeding_size = breeding_size_mean,
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_set, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in dispersal only
param_set_dispersal <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,
         breeding_size = breeding_size_mean, 
         k_connectivity = k_connectivity_set, theta_connectivity = theta_connectivity_set,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_matched, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Uncertainty in all parameters with uncertainty
param_set_full <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = start_recruit_size_set, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = growth_set_params$k_est, s = s, Linf = growth_set_params$Linf_est,
         breeding_size = breeding_size_set,
         k_connectivity = k_connectivity_set, theta_connectivity = theta_connectivity_set,  # dispersal kernel parameters
         prob_r = prob_r_set, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = prop_total_disp_area_sampled_best_est,
         offspring_assigned_to_parents = n_offspring_parentage_set, n_parents = n_parents_genotyped,
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat)  

### Run metrics for different types of uncertainty - no DD compensation
output_uncert_start_recruit <- calcMetricsAcrossRuns(n_runs, param_set_start_recruit, site_surv_best_est_sets, site_dist_info, site_vec_order, "start recruit size", FALSE)
output_uncert_growth <- calcMetricsAcrossRuns(n_runs, param_set_growth, site_surv_best_est_sets, site_dist_info, site_vec_order, "growth", FALSE)
output_uncert_survival <- calcMetricsAcrossRuns(n_runs, param_set_survival, site_surv_param_sets, site_dist_info, site_vec_order, "survival", FALSE)
output_uncert_breeding_size <- calcMetricsAcrossRuns(n_runs, param_set_breeding_size, site_surv_best_est_sets, site_dist_info, site_vec_order, "breeding size", FALSE)
output_uncert_offspring_assigned <- calcMetricsAcrossRuns(n_runs, param_set_offspring_assigned, site_surv_best_est_sets, site_dist_info, site_vec_order, "assigned offspring", FALSE)
output_uncert_prob_r <- calcMetricsAcrossRuns(n_runs, param_set_prob_r, site_surv_best_est_sets, site_dist_info, site_vec_order, "prob r", FALSE)
output_uncert_dispersal <- calcMetricsAcrossRuns(n_runs, param_set_dispersal, site_surv_best_est_sets, site_dist_info, site_vec_order, "dispersal k", FALSE)
output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info, site_vec_order, "all", FALSE)

### Run metrics for different types of uncertainty - with DD compensation
output_uncert_start_recruit_DD <- calcMetricsAcrossRuns(n_runs, param_set_start_recruit, site_surv_best_est_sets, site_dist_info, site_vec_order, "start recruit size", TRUE)
output_uncert_growth_DD <- calcMetricsAcrossRuns(n_runs, param_set_growth, site_surv_best_est_sets, site_dist_info, site_vec_order, "growth", TRUE)
output_uncert_survival_DD <- calcMetricsAcrossRuns(n_runs, param_set_survival, site_surv_param_sets, site_dist_info, site_vec_order, "survival", TRUE)
output_uncert_breeding_size_DD <- calcMetricsAcrossRuns(n_runs, param_set_breeding_size, site_surv_best_est_sets, site_dist_info, site_vec_order, "breeding size", TRUE)
output_uncert_offspring_assigned_DD <- calcMetricsAcrossRuns(n_runs, param_set_offspring_assigned, site_surv_best_est_sets, site_dist_info, site_vec_order, "assigned offspring", TRUE)
output_uncert_prob_r_DD <- calcMetricsAcrossRuns(n_runs, param_set_prob_r, site_surv_best_est_sets, site_dist_info, site_vec_order, "prob r", TRUE)
output_uncert_dispersal_DD <- calcMetricsAcrossRuns(n_runs, param_set_dispersal, site_surv_best_est_sets, site_dist_info, site_vec_order, "dispersal k", TRUE)
output_uncert_all_DD <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_param_sets, site_dist_info, site_vec_order, "all", TRUE)

#################### Analyze output, prepare for plotting: ####################
##### Join together runs with DD compensation with different types of uncertainty for plotting purposes 
# LEP averaged by site
LEP_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_out_df, output_uncert_growth_DD$LEP_out_df,
                    output_uncert_survival_DD$LEP_out_df, output_uncert_breeding_size_DD$LEP_out_df,
                    output_uncert_offspring_assigned_DD$LEP_out_df, output_uncert_prob_r_DD$LEP_out_df,
                    output_uncert_dispersal_DD$LEP_out_df, output_uncert_all_DD$LEP_out_df)

# LEP by site
LEP_by_site_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_by_site_out_df, output_uncert_growth_DD$LEP_by_site_out_df,
                       output_uncert_survival_DD$LEP_by_site_out_df, output_uncert_breeding_size_DD$LEP_by_site_out_df,
                       output_uncert_offspring_assigned_DD$LEP_by_site_out_df, output_uncert_prob_r_DD$LEP_by_site_out_df,
                       output_uncert_dispersal_DD$LEP_by_site_out_df, output_uncert_all_DD$LEP_by_site_out_df)

# LRP averaged across sites
LEP_R_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_R_out_df, output_uncert_growth_DD$LEP_R_out_df,
                      output_uncert_survival_DD$LEP_R_out_df, output_uncert_breeding_size_DD$LEP_R_out_df,
                      output_uncert_offspring_assigned_DD$LEP_R_out_df, output_uncert_prob_r_DD$LEP_R_out_df,
                      output_uncert_dispersal_DD$LEP_R_out_df, output_uncert_all_DD$LEP_R_out_df)

# LRP by site
LEP_R_by_site_uncert_DD <- rbind(output_uncert_start_recruit_DD$LEP_R_by_site_out_df, output_uncert_growth_DD$LEP_R_by_site_out_df,
                               output_uncert_survival_DD$LEP_R_by_site_out_df, output_uncert_breeding_size_DD$LEP_R_by_site_out_df,
                               output_uncert_offspring_assigned_DD$LEP_R_by_site_out_df, output_uncert_prob_r_DD$LEP_R_by_site_out_df,
                               output_uncert_dispersal_DD$LEP_R_by_site_out_df, output_uncert_all_DD$LEP_R_by_site_out_df)

# Egg-recruit survival
RperE_uncert_DD <- rbind(output_uncert_start_recruit_DD$RperE_out_df, output_uncert_growth_DD$RperE_out_df,
                      output_uncert_survival_DD$RperE_out_df, output_uncert_breeding_size_DD$RperE_out_df,
                      output_uncert_offspring_assigned_DD$RperE_out_df, output_uncert_prob_r_DD$RperE_out_df,
                      output_uncert_dispersal_DD$RperE_out_df, output_uncert_all_DD$RperE_out_df)

# Network persistence
NP_uncert_DD <- rbind(output_uncert_start_recruit_DD$NP_out_df, output_uncert_growth_DD$NP_out_df,
                   output_uncert_survival_DD$NP_out_df, output_uncert_breeding_size_DD$NP_out_df,
                   output_uncert_offspring_assigned_DD$NP_out_df, output_uncert_prob_r_DD$NP_out_df,
                   output_uncert_dispersal_DD$NP_out_df, output_uncert_all_DD$NP_out_df)

##### Find 95% range for estimates 

### Estimates with DD compensation
# NP with DD
NP_vec <- sort(output_uncert_all_DD$NP_out_df$value)
NP_lower <- round(NP_vec[quantile_lower_index], NP_decimal_points)
NP_upper <- round(NP_vec[quantile_upper_index], NP_decimal_points)

# LEP averaged across sites
LEP_avg_vec <- sort(output_uncert_all_DD$LEP_out_df$value)
LEP_avg_lower <- round(LEP_avg_vec[quantile_lower_index])
LEP_avg_upper <- round(LEP_avg_vec[quantile_upper_index])

# LEP at each site
LEP_site_limits <- data.frame(site = site_vec_order$site_name, lower = NA, upper = NA)
for(i in 1:n_sites) {
  LEP_site_vec <- sort((output_uncert_all_DD$LEP_by_site_out_df %>% filter(site == i))$value) 
  LEP_site_limits$lower[i] <- LEP_site_vec[quantile_lower_index]
  LEP_site_limits$upper[i] <- LEP_site_vec[quantile_upper_index]
}

# Egg-recruit survival with DD compensation
egg_recruit_survival_vec <- sort(output_uncert_all_DD$RperE_out_df$value)
egg_recruit_survival_lower <- egg_recruit_survival_vec[quantile_lower_index]
egg_recruit_survival_upper <- egg_recruit_survival_vec[quantile_upper_index]

# LRP local averaged across sites
LRP_local_average_vec <- sort(output_uncert_all_DD$LEP_R_local_out_df$value)
LRP_local_average_lower <- LRP_local_average_vec[quantile_lower_index]
LRP_local_average_upper <- LRP_local_average_vec[quantile_upper_index]

# Estimates > 1
LRP_avg_ests_above_1 <- output_uncert_all_DD$LEP_R_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)  # LRP averaged across sites
LR_est_DD_above_1 <- output_uncert_all_DD$LEP_R_local_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)  # LR 
NP_est_DD_above_1 <- output_uncert_all_DD$NP_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)

# Haina SP (site with highest SP)
Haina_SP_vec <- sort((output_uncert_all_DD$SP_vals_with_params %>% filter(site == "Haina"))$SP)
Haina_SP_lower <- Haina_SP_vec[quantile_lower_index]
Haina_SP_upper <- Haina_SP_vec[quantile_upper_index]

### Estimates without DD compensation
# NP
NP_vec_noDD <- sort(output_uncert_all$NP_out_df$value)
NP_lower_noDD <- round(NP_vec_noDD[quantile_lower_index], NP_decimal_points)
NP_upper_noDD <- round(NP_vec_noDD[quantile_upper_index], NP_decimal_points)

# Egg-recruit survival without DD compensation
egg_recruit_survival_vec_noDD <- sort(output_uncert_all$RperE_out_df$value)
egg_recruit_survival_lower_noDD <- egg_recruit_survival_vec_noDD[quantile_lower_index]
egg_recruit_survival_upper_noDD <- egg_recruit_survival_vec_noDD[quantile_upper_index]

# LEP R averaged across sites with DD compensation
LRP_average_vec <- sort(output_uncert_all_DD$LEP_R_out_df$value)
LRP_average_lower <- LRP_average_vec[quantile_lower_index]
LRP_average_upper <- LRP_average_vec[quantile_upper_index]

# LEP R averaged across sites without DD compensation
LRP_average_vec_noDD <- sort(output_uncert_all$LEP_R_out_df$value)
LRP_average_lower_noDD <- LRP_average_vec_noDD[quantile_lower_index]
LRP_average_upper_noDD <- LRP_average_vec_noDD[quantile_upper_index]

## LRP local averaged across sites without DD compensation
LRP_local_average_vec_noDD <- sort(output_uncert_all$LEP_R_local_out_df$value)
LRP_local_average_lower_noDD <- LRP_local_average_vec_noDD[quantile_lower_index]
LRP_local_average_upper_noDD <- LRP_local_average_vec_noDD[quantile_upper_index]

# Estimates > 1
LRP_avg_noDD_ests_above_1 <- output_uncert_all$LEP_R_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)  # LRP averaged across sites
NP_est_noDD_above_1 <- output_uncert_all$NP_out_df %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)

##### Find best estimates
# With density dependence compensation
LEP_best_est <- mean(best_est_metrics_mean_offspring$LEP_by_site$LEP)  # mean LEP across sites
LEP_R_best_est_DD <- best_est_metrics_mean_offspring_DD$LEP_R_mean  # mean LRP across sites
LEP_R_local_best_est_DD <- best_est_metrics_mean_offspring_DD$LEP_R_local_mean
NP_best_est_DD <- best_est_metrics_mean_offspring_DD$NP
SP_best_est_DD <- best_est_metrics_mean_offspring_DD$SP

# Without density dependence compensation
LEP_R_best_est <- best_est_metrics_mean_offspring$LEP_R_mean
LEP_R_local_best_est <- best_est_metrics_mean_offspring$LEP_R_local_mean
NP_best_est <- best_est_metrics_mean_offspring$NP
SP_best_est <- best_est_metrics_mean_offspring$SP

#################### What-if calculations, including alternative geographies and larval navigation: ####################

##### What-if calculation 1) what if all genotyped offspring came from the population? (Is the metapopulation persistent when we include immigrants)

# Parameter set to do an estimate using all offspring genotyped, rather than just those matched
param_best_est_mean_collected_offspring_all_offspring <- data.frame(t_steps = n_tsteps) %>%
  mutate(min_size = min_size, max_size = max_size, n_bins = n_bins,  # LEP-IPM matrix
         clutches_per_year_mean = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect =  eyed_effect,  # size-dependent fecundity info
         start_recruit_size = mean_sampled_offspring_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Linf = Linf_growth_mean,
         breeding_size = breeding_size_mean,  # size of transition to female
         k_connectivity = k_allyears, theta_connectivity = theta_allyears,  # dispersal kernel parameters
         prob_r = prob_r_mean, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled = 1,  # here, assuming all the offspring arriving came from these sites so don't need to scale it up for the dispersal kernel area 
         offspring_assigned_to_parents = n_offspring_genotyped, n_parents = n_parents_genotyped,  # offspring assigned are all offspring genotyped
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat) 

# Calculate the metrics for the best estimates
best_est_metrics_mean_offspring_all_offspring <- calcMetrics(param_best_est_mean_collected_offspring_all_offspring, site_surv_best_est, site_dist_info, site_vec_order, FALSE)  # without density dependence compensation
best_est_metrics_mean_offspring_all_offspring_DD <- calcMetrics(param_best_est_mean_collected_offspring_all_offspring, site_surv_best_est, site_dist_info, site_vec_order, TRUE)  # with density dependence compensation

# Do uncertainty run with all uncertainty included (except offspring assigned since using all genotyped offspring as matches in this estimate)
param_set_full_all_offspring <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = start_recruit_size_set, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = growth_set_params$k_est, s = s, Linf = growth_set_params$Linf_est,
         breeding_size = breeding_size_set,
         k_connectivity = k_connectivity_set, theta_connectivity = theta_connectivity_set,  # dispersal kernel parameters
         prob_r = prob_r_set, total_prop_hab_sampled = total_prop_hab_sampled_through_time, prop_total_disp_area_sampled_best = 1,
         offspring_assigned_to_parents = n_offspring_genotyped, n_parents = n_parents_genotyped,  # no uncertainty in offspring assigned this time
         perc_APCL = perc_APCL_val, perc_UNOC = perc_UNOC_val, prop_hab = prop_sampling_area_habitat)  

output_uncert_all_offspring_all <- calcMetricsAcrossRuns(n_runs, param_set_full_all_offspring, site_surv_param_sets, site_dist_info, site_vec_order, "all: alloff", FALSE)
output_uncert_all_offspring_all_DD <- calcMetricsAcrossRuns(n_runs, param_set_full_all_offspring, site_surv_param_sets, site_dist_info, site_vec_order, "all: alloff", TRUE)

##### What-if calculation 2) What would LRP need to be for NP to be 1? 
LRP_vec_WI2 = seq(from=0.1, to=10.0, by=0.01)  # vector of potential LRP values
NP_output_vec_WI2  <- rep(NA, length(LRP_vec_WI2))

# Find the NP value for those LRP values
for(i in 1:length(NP_output_vec_WI2)) {
  conn_matR_WI2 = best_est_metrics_mean_offspring$conn_matrix*LRP_vec_WI2[i]
  eig_CR_WI2 = eigen(conn_matR_WI2)
  NP_output_vec_WI2[i] = Re(eig_CR_WI2$values[1])
}

LRP_for_NP = LRP_vec_WI2[which(NP_output_vec_WI2 >= 1)[1]]  

##### What-if calculation 3) To get the required LRP for NP, what do egg-recruit survival (for our LEP) and LEP (for our egg-recruit-survival) need to be?
egg_recruit_survival_for_NP <- LRP_for_NP/mean(best_est_metrics_mean_offspring$LEP_by_site$LEP)

LEP_for_NP <- LRP_for_NP/best_est_metrics_mean_offspring$recruits_per_egg

LEP_for_NP_DD <- LRP_for_NP/best_est_metrics_mean_offspring_DD$recruits_per_egg

# % of estimates of LRP (with DD) that are >= LRP_for_NP
LRP_ests_above_LRP_for_NP <- output_uncert_all_DD$LEP_R_out_df %>% filter(value >= LRP_for_NP) %>% summarize(n_estimates = n()) %>% select(n_estimates)

##### Alternative geography - how much habitat would we need for the sampling region to be persistent?
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

# Find best estimate and uncertainty metrics for those habitat configurations with median-site survs, including compensation for density-dependence
perc_hab_best_ests_avgSurvs <- list()  # place to store point estimates
perc_hab_uncertainty_avgSurvs <- list()  # place to store estimates from uncertainty runs
for(i in 1:length(perc_hab_vals)) {
  perc_hab_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_med_Sint, site_dist_info_habperc[[i]], site_vec_order, TRUE)
  perc_hab_uncertainty_avgSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_median_Sint_param_sets, site_dist_info_habperc[[i]], site_vec_order, "all: hab sens, med surv", TRUE)
}

# Make dfs for plotting where each run is its own line - perc hab on x axis, NP on y axis
perc_hab_plot_df <- perc_hab_uncertainty_avgSurvs[[1]]$NP_out_df %>% select(value,run) %>% mutate(perc_hab = perc_hab_vals[1])  # one for uncertainty runs
for(i in 2:length(perc_hab_vals)) {
  perc_hab_plot_df <- rbind(perc_hab_plot_df, perc_hab_uncertainty_avgSurvs[[i]]$NP_out_df %>% select(value,run) %>% mutate(perc_hab = perc_hab_vals[i]))
}

perc_hab_best_est_NP_df <- data.frame(perc_hab = rep(NA,length(perc_hab_vals)), value = rep(NA,length(perc_hab_vals)))  # and one for best ests
for(i in 1:length(perc_hab_vals)) {
  perc_hab_best_est_NP_df$perc_hab[i] = perc_hab_vals[i]
  perc_hab_best_est_NP_df$value[i] = perc_hab_best_ests_avgSurvs[[i]]$NP
}

# Find number of runs with NP >=1 at perc_hab = 0.9 (where NP is > 1) (not sure why filtering for == 0.9 doesn't work)
perc_hab_plot_df %>% filter(0.88 < perc_hab & perc_hab < 0.92) %>% filter(value >= 1) %>% summarize(n_estimates = n()) %>% select(n_estimates)

##### Alternative geographic - same habitat density, wider region 
region_width_list <- c(region_width_km, 35, 40, 45, 50, 55)

# Go through different region widths for the real habitat percentage and find site locations, distances, etc. for equally sized and spaced sites
site_dist_info_wider_region <- list()
for(i in 1:length(region_width_list)) {
  site_dist_out_df <- make_output_with_dist(n_sites, Ps, region_width_list[i])
  site_dist_out_df$org_site <- (as.data.frame(site_vec_NS, stringsAsFactors = FALSE) %>% slice(rep(1:n(), each=n_sites)))$site_vec_NS  # replace numeric org_site with names
  site_dist_out_df$dest_site <- rep(site_vec_NS, n_sites)  # replace numeric dest_site with names
  site_dist_info_wider_region[[i]] <- site_dist_out_df
}

# Find best estimate and uncertainty metrics for those habitat configurations with median-site survs (and including compensation for density-dependence)
wider_region_best_ests_avgSurvs <- list()
wider_region_uncertainty_avgSurvs <- list()
for(i in 1:length(region_width_list)) {
  wider_region_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_med_Sint, site_dist_info_wider_region[[i]], site_vec_order, TRUE)
  wider_region_uncertainty_avgSurvs[[i]] <- calcMetricsAcrossRuns(n_runs, param_set_full, site_surv_median_Sint_param_sets, site_dist_info_wider_region[[i]], site_vec_order, "all: region sens, avg surv", TRUE)
}

# Make dfs for plotting where each run is it's own line - region width on x axis, NP on y axis
wider_region_plot_df <- wider_region_uncertainty_avgSurvs[[1]]$NP_out_df %>% select(value,run) %>% mutate(region_width = region_width_list[1])  # one for uncertainty runs
for(i in 2:length(region_width_list)) {
  wider_region_plot_df <- rbind(wider_region_plot_df, wider_region_uncertainty_avgSurvs[[i]]$NP_out_df %>% select(value,run) %>% mutate(region_width = region_width_list[i]))
}
wider_region_best_est_NP_df <- data.frame(region_width = rep(NA,length(region_width_list)), value = rep(NA,length(region_width_list)))  # and one for best ests
for(i in 1:length(region_width_list)) {
  wider_region_best_est_NP_df$region_width[i] = region_width_list[i]
  wider_region_best_est_NP_df$value[i] = wider_region_best_ests_avgSurvs[[i]]$NP
}

##### Alternative geography - 100% habitat, wider region
# Go through different region widths for 100% habitat percentage and find site locations, distances, etc. for equally sized and spaced sites
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

##### Alternative geography for % habitat and region width combined
# Find the site dist info
perc_hab_vals_long_form <- seq(from=0.01, to=1.0, by = 0.01)  # percent habitat values to use
region_width_list_long_form <- seq(from = 20, to = 55, by = 1)  # region width values to use

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

# Find best estimate of metrics for those habitat configurations with median-site survs (and including compensation for density-dependence)
wider_region_perc_hab_best_ests_avgSurvs <- list()
for(i in 1:(length(region_width_list_long_form)*length(perc_hab_vals_long_form))) {
  wider_region_perc_hab_best_ests_avgSurvs[[i]] <- calcMetrics(param_best_est_mean_collected_offspring, site_surv_best_est_med_Sint, site_dist_info_wider_region_perc_hab[[i]], site_vec_order, TRUE)
  print(i)
}

# Make dfs for plot where region width, perc hab, and NP all listed
wider_region_perc_hab_plot_df <- data.frame(NP = rep(NA, length(region_width_list_long_form)*length(perc_hab_vals_long_form)))
for(i in 1:(length(region_width_list_long_form)*length(perc_hab_vals_long_form))) {
  wider_region_perc_hab_plot_df$NP[i] = wider_region_perc_hab_best_ests_avgSurvs[[i]]$NP
}
wider_region_perc_hab_plot_df$region_width = (as.data.frame(region_width_list_long_form, stringsAsFactors = FALSE) %>% slice(rep(1:n(), each = length(perc_hab_vals_long_form))))$region_width_list_long_form
wider_region_perc_hab_plot_df$perc_hab = rep(perc_hab_vals_long_form, length(region_width_list_long_form))

# Put NP into categories: 0-.25, 0.25-0.5, 0.5-0.75, 0.75<1.0, >=1.0
wider_region_perc_hab_plot_df_categories <- wider_region_perc_hab_plot_df %>%
  mutate(NP_categories = case_when(NP < 0.25 ~ "< 0.25",
                                 (0.25 <= NP & NP < 0.5) ~ "[0.25 - 0.5)",
                                 (0.5 <= NP & NP < 0.75) ~ "[0.5 - 0.75)",
                                 (0.75 <= NP & NP < 1) ~ "[0.75 - 1)",
                                 (NP >= 1) ~ ">= 1"),
         order_5 = case_when(NP < 0.25 ~ 1,
                                  (0.25 <= NP & NP < 0.5) ~ 2,
                                  (0.5 <= NP & NP < 0.75) ~ 3,
                                  (0.75 <= NP & NP < 1) ~ 4,
                                  (NP >= 1) ~ 5)) 

##### What if larvae could navigate? Integrate kernel adding a buffer of up to 1000m on either side of each patch
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

##### Make sure upper size is sufficient that fish aren't getting evicted from the LEP matrix
max_size_test_vec <- seq(from=breeding_size_mean, to=20, by=0.1 )
LEP_test_max_size <- data.frame(max_size = max_size_test_vec, LEP = NA)

for(i in 1:length(max_size_test_vec)) {
  LEP_test_max_size$LEP[i] = findLEP(param_best_est_mean_collected_offspring$min_size, max_size_test_vec[i],
                                     param_best_est_mean_collected_offspring$n_bins, param_best_est_mean_collected_offspring$t_steps,
                                     site_surv_best_est_med_Sint$Sint[1], site_surv_best_est_med_Sint$Sl[1],
                                     param_best_est_mean_collected_offspring$s, param_best_est_mean_collected_offspring$Linf,
                                     param_best_est_mean_collected_offspring$k_growth,
                                     param_best_est_mean_collected_offspring$clutches_per_year_mean,
                                     param_best_est_mean_collected_offspring$breeding_size, param_best_est_mean_collected_offspring$start_recruit_size,
                                     param_best_est_mean_collected_offspring$start_recruit_sd, param_best_est_mean_collected_offspring$egg_size_slope,
                                     param_best_est_mean_collected_offspring$egg_size_intercept, param_best_est_mean_collected_offspring$eyed_effect)
}
  
#################### Metrics and parameters summarized for easy access: ####################
params_summary <- data.frame(param = c("k_disp","k_disp_lcl","k_disp_ucl",
                                        "theta_disp","theta_disp_lcl","theta_disp_ucl",
                                        "Linf","Linf_lcl","Linf_ucl",
                                        "k_growth","k_growth_lcl","k_growth_ucl",
                                        "growth_slope_est_lower","growth_slope_est_upper",
                                        "growth_intercept_est_lower","growth_intercept_est_upper",
                                        "n_parents_genotyped","n_offspring_genotyped",
                                        "n_offspring_matched","assignment_rate",
                                        "Ph","Pc","Ps","Pd",
                                        "p_APCL","p_UNOC"),
                              value = c(k_allyears, min(k_theta_allyear_95CI_values$k_eval), max(k_theta_allyear_95CI_values$k_eval),
                                        theta_allyears, min(k_theta_allyear_95CI_values$theta_eval), max(k_theta_allyear_95CI_values$theta_eval),
                                        Linf_growth_mean, Linf_growth_lower, Linf_growth_upper,
                                        k_growth_mean, k_growth_lower, k_growth_upper,
                                        growth_slope_est_lower, growth_slope_est_upper,
                                        growth_intercept_est_lower, growth_intercept_est_upper,
                                        n_parents_genotyped, n_offspring_genotyped, 
                                        n_offspring_matched, assignment_rate,
                                        Ph, prob_r_mean, Ps, Pd,
                                        perc_APCL_val, perc_UNOC_val), stringsAsFactors = FALSE)

metrics_summary <- data.frame(metric = c("LEP avg best est", "LEP avg lower", "LEP avg upper",
                                         "LRP_DD avg best est", "LRP_DD avg lower", "LRP_DD avg upper",
                                         "NP_DD best est", "NP_DD lower", "NP_DD upper",
                                         "egg_recruit surv best est DD", "egg recruit surv lower DD", "egg recruit surv upper DD",
                                         "LR avg best est DD", "LR avg lower DD", "LR avg upper DD",
                                         "LRP avg ests above 1 DD", "LR avg ests DD above 1", "NP ests DD above 1",
                                         "Haina SP", "Haina SP lower", "Haina SP upper",
                                         "LRP with all offspring", "LR with all offspring", "LR no DD",
                                         "LRP no DD avg best est", "LRP no DD avg lower", "LRP no DD avg upper",
                                         "NP no DD best est", "NP no DD lower", "NP no DD upper",
                                         "LR avg best est no DD", "LR avg best est lower no DD", "LR avg best est upper no DD",
                                         "LRP needed for NP >= 1", "egg recruit surv for NP >=1", "LEP for NP >= 1",
                                         "LRP ests above LRP for NP",
                                         "egg recruit surv best est no DD", "egg recruit surv lower no DD", "egg recruit surv upper no DD",
                                         "LRP avg ests above 1 no DD"
                                         ),
                              value = c(LEP_best_est, LEP_avg_lower, LEP_avg_upper,
                                        LEP_R_best_est_DD, LRP_average_lower, LRP_average_upper,
                                        NP_best_est_DD, NP_lower, NP_upper,
                                        best_est_metrics_mean_offspring_DD$recruits_per_egg, egg_recruit_survival_lower, egg_recruit_survival_upper,
                                        LEP_R_local_best_est_DD, LRP_local_average_lower, LRP_local_average_upper,
                                        LRP_avg_ests_above_1$n_estimates, LR_est_DD_above_1$n_estimates, NP_est_DD_above_1$n_estimates,
                                        (SP_best_est_DD %>% filter(site == "Haina"))$SP_value, Haina_SP_lower, Haina_SP_upper,
                                        best_est_metrics_mean_offspring_all_offspring_DD$LEP_R_mean, best_est_metrics_mean_offspring_all_offspring_DD$LEP_R_local_mean, best_est_metrics_mean_offspring$LEP_R_local_mean,
                                        LEP_R_best_est, LRP_average_lower_noDD, LRP_average_upper_noDD,
                                        NP_best_est, NP_lower_noDD, NP_upper_noDD,
                                        LEP_R_local_best_est, LRP_local_average_lower_noDD, LRP_local_average_upper_noDD,
                                        LRP_for_NP, egg_recruit_survival_for_NP, LEP_for_NP_DD,
                                        LRP_ests_above_LRP_for_NP$n_estimates,
                                        best_est_metrics_mean_offspring$recruits_per_egg, egg_recruit_survival_lower_noDD, egg_recruit_survival_upper_noDD,
                                        LRP_avg_noDD_ests_above_1$n_estimates
                                        ), stringsAsFactors = FALSE)
                              

#################### Figure prep work - vectors and values needed: ####################
# Dispersal kernel
distance_vec <- seq(from=0, to=50, by=0.01)  # distances to show
connectivity_est_vec <- disp_kernel_all_years(distance_vec, k_allyears, theta_allyears)  # vector of dispersal probs to plot

figure_vectors_and_values <- list(distance_vec = distance_vec, connectivity_est_vec = connectivity_est_vec, 
                                  region_width_km = region_width_km, region_width_list = region_width_list,
                                  region_width_list_long_form = region_width_list_long_form,
                                  perc_hab_vals = perc_hab_vals, perc_hab_vals_long_form = perc_hab_vals_long_form,
                                  nav_buffer_vec = nav_buffer_vec)

#################### Save output: ####################

### Relevant intermediate outputs
save(all_parents_site, file=here::here("Data/Script_outputs","all_parents_site.RData"))  # info on parents by site, distance to edges of sampling region

### Parameter sets
save(param_best_est_mean_collected_offspring, file=here::here("Data/Script_outputs", "param_best_est_mean_collected_offspring.RData"))  # point estimates
save(param_set_full, file=here::here("Data/Script_outputs", "param_set_full.RData"))  # parameter set with uncertainty values
save(site_surv_param_sets, file=here::here("Data/Script_outputs","site_surv_param_sets.RData"))  # survival parameter sets
save(site_surv_median_Sint_param_sets, file=here::here("Data/Script_outputs","site_surv_median_Sint_param_sets.RData"))  # survival parameter sets for median survival (for alternate geographies)

### Point estimates
save(best_est_metrics_mean_offspring_DD, file=here::here("Data/Script_outputs", "best_est_metrics_mean_offspring_DD.RData"))  # point estimates with density dependence compensation (presented in main text)
save(best_est_metrics_mean_offspring, file=here::here("Data/Script_outputs", "best_est_metrics_mean_offspring.RData"))  # point estimates without density dependence compensation (presented in appendix)

### Uncertainty estimates with density dependence compensation (presented in main text)
save(output_uncert_start_recruit_DD, file=here::here("Data/Script_outputs", "output_uncert_start_recruit_DD.RData"))
save(output_uncert_growth_DD, file=here::here("Data/Script_outputs", "output_uncert_growth_DD.RData"))
save(output_uncert_survival_DD, file=here::here("Data/Script_outputs", "output_uncert_survival_DD.RData"))
save(output_uncert_breeding_size_DD, file=here::here("Data/Script_outputs", "output_uncert_breeding_size_DD.RData"))
save(output_uncert_offspring_assigned_DD, file=here::here("Data/Script_outputs", "output_uncert_offspring_assigned_DD.RData"))
save(output_uncert_prob_r_DD, file=here::here("Data/Script_outputs", "output_uncert_prob_r_DD.RData"))
save(output_uncert_dispersal_DD, file=here::here("Data/Script_outputs", "output_uncert_dispersal_DD.RData"))
save(output_uncert_all_DD, file=here::here("Data/Script_outputs", "output_uncert_all_DD.RData"))

### Uncertainty estimates without density dependence compensation (presented in appendix)
save(output_uncert_start_recruit, file=here::here("Data/Script_outputs", "output_uncert_start_recruit.RData"))
save(output_uncert_growth, file=here::here("Data/Script_outputs", "output_uncert_growth.RData"))
save(output_uncert_survival, file=here::here("Data/Script_outputs", "output_uncert_survival.RData"))
save(output_uncert_breeding_size, file=here::here("Data/Script_outputs", "output_uncert_breeding_size.RData"))
save(output_uncert_offspring_assigned, file=here::here("Data/Script_outputs", "output_uncert_offspring_assigned.RData"))
save(output_uncert_prob_r, file=here::here("Data/Script_outputs", "output_uncert_prob_r.RData"))
save(output_uncert_dispersal, file=here::here("Data/Script_outputs", "output_uncert_dispersal.RData"))
save(output_uncert_all, file=here::here("Data/Script_outputs", "output_uncert_all.RData"))

### LEP, LRP, egg recruit survival, NP data frames for plotting contributors to uncertainty
save(LEP_uncert_DD, file=here::here("Data/Script_outputs","LEP_uncert_DD.RData"))
save(LEP_R_uncert_DD, file=here::here("Data/Script_outputs","LEP_R_uncert_DD.RData"))
save(RperE_uncert_DD, file=here::here("Data/Script_outputs","RperE_uncert_DD.RData"))
save(NP_uncert_DD, file=here::here("Data/Script_outputs","NP_uncert_DD.RData"))

### What-ifs
# All arriving recruits considered offspring
save(best_est_metrics_mean_offspring_all_offspring, file=here::here("Data/Script_outputs","best_est_metrics_mean_offspring_all_offspring.RData"))  # all arriving recruits considered offspring, without DD compensation
save(best_est_metrics_mean_offspring_all_offspring_DD, file=here::here("Data/Script_outputs","best_est_metrics_mean_offspring_all_offspring_DD.RData"))  # all arriving recruits considered offspring, with DD compensation

# Alternative geographies - sensitivity to percent habitat
save(perc_hab_best_ests_avgSurvs, file=here::here("Data/Script_outputs","perc_hab_best_ests_avgSurvs.RData"))  # point estimates
save(perc_hab_uncertainty_avgSurvs, file=here::here("Data/Script_outputs","perc_hab_uncertainty_avgSurvs.RData"))  # uncertainty
save(perc_hab_plot_df, file=here::here("Data/Script_outputs","perc_hab_plot_df.RData"))  # perc hab sensitivity in easy-to-plot form
save(perc_hab_best_est_NP_df, file=here::here("Data/Script_outputs","perc_hab_best_est_NP_df.RData"))

# Alternative geographies - sensitivity to wider region
save(wider_region_best_ests_avgSurvs, file=here::here("Data/Script_outputs","wider_region_best_ests_avgSurvs.RData"))  # point estimates
save(wider_region_uncertainty_avgSurvs, file=here::here("Data/Script_outputs","wider_region_uncertainty_avgSurvs.RData"))  # uncertainty
save(wider_region_plot_df, file=here::here("Data/Script_outputs","wider_region_plot_df.RData"))  # in an easy-to-plot form
save(wider_region_best_est_NP_df, file=here::here("Data/Script_outputs","wider_region_best_est_NP_df.RData"))
save(wider_region_perc_hab_plot_df_categories, file=here::here("Data/Script_outputs","wider_region_perc_hab_plot_df_categories.RData"))  # with NP values lumped into 5 categories

# Alternative geographies - sensitivity to wider region at 100% habitat
save(wider_region_all_hab_best_ests_avgSurvs, file=here::here("Data/Script_outputs","wider_region_all_hab_best_ests_avgSurvs.RData"))
save(wider_region_all_hab_uncertainty_avgSurvs, file=here::here("Data/Script_outputs","wider_region_all_hab_uncertainty_avgSurvs.RData"))

# Alternative geographies - sensitivity to wider region and percent habitat simultaneously
save(wider_region_perc_hab_best_ests_avgSurvs, file=here::here("Data/Script_outputs","wider_region_perc_hab_best_ests_avgSurvs.RData"))
save(wider_region_perc_hab_plot_df, file=here::here("Data/Script_outputs","wider_region_perc_hab_plot_df.RData"))

# Larval navigation sensitivity
save(larv_nav_prop_hab, file=here::here("Data/Script_outputs","larv_nav_prop_hab.RData"))  # effective proportion habitat when consider larva navigation
save(larv_nav_params_list, file=here::here("Data/Script_outputs","larv_nav_params_list.RData"))  # parameter sets that take into account effective proportion habitat
save(larv_nav_best_ests, file=here::here("Data/Script_outputs","larv_nav_best_ests.RData"))  # point estimates
save(larv_nav_uncert_all, file=here::here("Data/Script_outputs","larv_nav_uncert_all.RData"))  # uncertainty
save(larv_nav_plot_df, file=here::here("Data/Script_outputs","larv_nav_plot_df.RData"))  # in an easy-to-plot form
save(larv_nav_best_est_NP_df, file=here::here("Data/Script_outputs","larv_nav_best_est_NP_df.RData"))

# Summary of parameters and metrics and figure inputs for easy reference while writing
save(params_summary, file=here::here("Data/Script_outputs","params_summary.RData"))
save(metrics_summary, file=here::here("Data/Script_outputs","metrics_summary.RData"))
save(figure_vectors_and_values, file=here::here("Data/Script_outputs","figure_vectors_and_values.RData"))  # inputs and parameters, etc. useful for plotting figures

 
