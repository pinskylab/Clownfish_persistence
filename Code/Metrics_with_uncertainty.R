# Characterizing uncertainty in LEP, recruit survival, dispersal estimates

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))

load(file=here("Data", "eall_Phi_size_p_dist_results.RData")) #MARK output 

##### Set run info
n_runs = 100

##### Parameter info (candidates for uncertainty)
# Growth (for LEP)
s = exp(-0.0148)  # what is this?? goes into dnorm for growth part...
k_growth_mean = 0.9447194  # lowest AIC model
Linf_growth_mean = 10.50670  # lowest AIC model
Linf_growth_sd = sqrt(1.168163)  # from variance for Linf in lowest AIC model

# Eggs (for LEP)
eggs_per_clutch_mean = 514.11
clutches_per_year_mean = 11.9

# Survival (for LEP)
Sint_mean = eall.Phi.size.p.dist.results$estimate[1]  # survival intercept (on logit scale)
Sl_mean = eall.Phi.size.p.dist.results$estimate[2]  # survival slope (on logit scale)
Sint_se = eall.Phi.size.p.dist.results$se[1]  # for now using SE, should really use SD...
Sint_se = eall.Phi.size.p.dist.results$se[2]  # for now using SE, should really use SD...

# Breeding size (for LEP)
breeding_size_mean = 8  # could actually find this by finding the mean and SD of size of breeding females...
breeding_size_sd = 0.8

##### Other parameters (for running IPM, for calculating connectivity, for uncertainty runs, etc.)
# Set params for IPM structure
n_bins = 100
n_tsteps = 100

# Other size points
min_size = 0
max_size = 15 #should check this w/data...

##### Generate sets of parameters
Linf_set = rnorm(n_runs, mean = Linf_growth_mean, sd=Linf_growth_sd)
Sint_set = rnorm(n_runs, mean = Sint_mean, sd= Sint_se)
breeding_size_set = rnorm(n_runs, mean = breeding_size, sd = breeding_size_sd) 

##### Other parameters that stay static

#################### Functions: ####################
findLEP = function(min_size, max_size, n_bins, t_steps, Sint, Sl, s, Linf, k_growth, eggs_per_clutch, clutches_per_year) {
  
  # Create vector of lengths
  lengths_vec = seq(min_size, max_size, length.out = n_bins)
  dx = diff(lengths_vec[1:2])
  
  # Make matrix
  xmat = matrix(rep(lengths_vec, n_bins), nrow=n_bins, ncol=n_bins, byrow=TRUE) #100 rows of lengths_vec
  ymat = t(xmat)
  
  # Survival probability (based on size)
  S = logit_recip(Sint + lengths_vec*Sl) #is this right? Also, this should be prob of surviving, not prob of dying, right? (Which it is right now)
  Smat = matrix(rep(S,n_bins), nrow=n_bins, ncol=n_bins, byrow = TRUE) #survival part
  
  # Growth probability (this is the part that really could use updating!)
  Ls = Linf - (Linf - xmat)*exp(-k_growth)
  Lmat = dnorm(ymat,Ls,s) #growth part
  
  # Make the kernel with the survival + growth inputs
  K = Lmat*Smat
  K = dx*K
  
  # Iterate with 1 subadult to start
  N0 = dnorm(lengths_vec,3.5,0.1) # not totally clear what is happening here...
  
  # Integrate and scale to have it integrate to 1
  N0 = N0/sum(N0*dx)
  
  # Iterate through time
  N = matrix(rep(NA,n_bins*n_tsteps), nrow=n_bins, ncol=n_tsteps)
  N[,1] = N0  # initial conditions
  # iterate through time
  for (t in 2:n_tsteps) {
    N[,t] = K %*% N[, t-1]
  }
  
  # Find the number of breeding adults one recruit produces
  breeding_Adults <- colSums(N[lengths_vec>breeding_size,]*dx)
  
  # Eggs produced by one breeding adult in one year
  eggs <- eggs_per_clutch*clutches_per_year
  
  # Combine to get LEP
  LEP <- sum(breeding_Adults*eggs)
  
  return(LEP)
}

##### Run through a metric calculationg
runCalc <- function(static_param_set,) {
  
  # Find LEP
  LEP = findLEP(param_set$min_size, param_set$max_size, param_set$n_bins, param_set$t_steps, param_set)
  
  
  min_size, max_size, n_bins, t_steps, Sint, Sl, s, Linf, k_growth, eggs_per_clutch, clutches_per_year
  
}
#################### Running things: ####################











## Change size-dependent survival estimates from logit estimates in loaded file
#Mint = logit_recip(eall.Phi.size.p.dist.results$estimate[1])
#Ml = logit_recip(eall.Phi.size.p.dist.results$estimate[2])


#################### Running things: ####################

#################### Plots: ####################

#################### Saving things: ####################

#################### Old code: ####################