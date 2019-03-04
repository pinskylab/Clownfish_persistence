# Characterizing uncertainty in LEP, recruit survival, dispersal estimates

# To add: uncertainty in egg-recruit survival
# fix uncertainty in growth and uncertainty in annual survival
# uncertainty in egg-size relationship?

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))

load(file=here('Data', 'female_sizes.RData'))  # sizes of females from data
load(file=here('Data', 'eall_mean_Phi_size_p_size_plus_dist.RData'))  # MARK output (lowest AICc model)
load(file=here('Data', 'loglogFecunditySizeModel.RData'))  # size-fecundity output for best-fit model from Adam, called length_count8llEA
load(file=here('Data', 'c_mat_allyears.RData'))  # Probability of dispersing (for C matrix for now, before use kernel params to include uncertainty)
k_connectivity_values <- as.vector(readRDS(file=here('Data', 'avg_bootstrapped_k.rds')))  # values of k within the 95% confidence interval, bootstrapped - downloaded from KC parentage repository on 2/27/19

#load(file=here('Data', 'size_by_color_metrics.RData'))  # size distribution info by tail color
#load(file=here("Data", "eall_Phi_size_p_dist_results.RData")) #MARK output 

#### Set-up parameters (for running IPM, for calculating connectivity, for uncertainty runs, etc.)
# Number of runs
n_runs = 1000

# Set params for IPM structure
n_bins = 100
n_tsteps = 100
start_recruit_size = 3.5  # size of recruit that starts out the IPM for LEP
start_recruit_sd = 0.1

# Other size points
min_size = 0
max_size = 15 #should check this w/data...

##### Parameter info (candidates for uncertainty)
# # Connectivity  - estimates from Dec. 18 KC paper draft - now pulling straight from distribution of values, rather than re-creating distribution
k_allyears = -1.36  # with 2012-2015 data
theta_allyears = 0.5  # with 2012-2015 data
# k_allyears_CIh = -0.97  # upper 97.5% confidence interval of k (from KC email with screenshot)
# k_allyears_CIl = -1.94  # lower 97.5% confidence interval of k (these aren't symmetric, so not normal? Check with KC how CI were derived)

# Growth (for LEP)
s = exp(-0.0148)  # what is this?? goes into dnorm for growth part... sd around the mean size? Not sure where this estimate came from...
k_growth_mean = 0.9447194  # lowest AIC model
Linf_growth_mean = 10.50670  # lowest AIC model
Linf_growth_sd = sqrt(1.168163)  # from variance for Linf in lowest AIC model

# Eggs (for LEP)
size_fecundity_model = length_count8llEA  # assign here, in case model input from Adam changes, both length and eggs on log scale
eggs_intercept_log = size_fecundity_model$coefficients[1]  # on log-scale
eggs_slope_log = size_fecundity_model$coefficients[2]
eyed_effect = size_fecundity_model$coefficients[3]

eggs_per_clutch_mean = 514.11  # need to rethink what this means with size-effect, how to do one weighted by size...
clutches_per_year_mean = 11.9
# egg_intercept = -426.57
# egg_slope = 107  # eggs/cm size of F

# Survival (for LEP)
eall_mean.Phi.size.p.size.plus.dist.results <- as.data.frame(eall_mean.Phi.size.p.size.plus.dist$results$beta)
Sint_mean = eall_mean.Phi.size.p.size.plus.dist.results$estimate[1]  # survival intercept (on logit scale)
Sl_mean = eall_mean.Phi.size.p.size.plus.dist.results$estimate[2]  # survival slope (on logit scale)
Sint_se = eall_mean.Phi.size.p.size.plus.dist.results$se[1]  # for now using SE, should really use SD...
Sint_se = eall_mean.Phi.size.p.size.plus.dist.results$se[2]  # for now using SE, should really use SD...

# Egg-recruit survival (for getting LEP in terms of recruits)
recruits_per_egg = 8.367276e-05  # surv_egg_recruit estimating using Johnson method in PersistenceMetrics.R

##### Other parameters that stay static

#################### Functions: ####################

# Find eggs by fish size (eyed eggs) - should double check used log (ln) and not log10 (base 10 log)
findEggs = function(fish_size, egg_size_intercept, egg_size_slope, eyed_effect) {
  count_logged = egg_size_intercept + egg_size_slope*log(fish_size) + eyed_effect
  raw_eggs = exp(count_logged)
  return(raw_eggs)
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
  
  # Compare to the LEP estimate without size-dependent fecundity
  # Find the number of breeding adults one recruit produces
  breeding_Adults <- colSums(N[lengths_vec>breeding_size,]*dx)
  
  # Eggs produced by one breeding adult in one year
  eggs <- eggs_per_clutch*clutches_per_year
  
  # Combine to get LEP
  LEP_nossF <- sum(breeding_Adults*eggs)
  
  out = list(LEP=LEP, LEP_nossF = LEP_nossF)
  
  return(LEP)
}

# Run through a metric calculation with one set of parameters
calcMetrics <- function(param_set, Cmatrix, sites) {
  
  # Define function with the right parameters
  disp_allyears <- function(d) {  # theta = 0.5, equation for p(d) in eqn. 6c in Bode et al. 2018
    z = exp(param_set$k_connectivity)
    disp = (z/2)*exp(-(z*d)^(param_set$theta_connectivity))
    return(disp)
  }
  
  # Create connectivity matrix
  Cmat <- Cmatrix %>% select(org_site, dest_site, d1_km, d2_km, org_alpha_order, org_geo_order, dest_alpha_order, dest_geo_order)
  for(i in 1:length(Cmat$org_site)) {
    Cmat$prob_disp[i] <- integrate(disp_allyears, Cmat$d1_km[i], Cmat$d2_km[i])$value
  }
  
  # Find LEP (in terms of eggs) - COULD MAKE MORE OF THESE PARAMETERS PULLED FROM A DISTRIBUTION!
  LEP = findLEP(param_set$min_size, param_set$max_size, param_set$n_bins, param_set$t_steps, param_set$Sint, param_set$Sl,
                param_set$s, param_set$Linf, param_set$k_growth, param_set$eggs_per_clutch, param_set$clutches_per_year, 
                param_set$breeding_size, param_set$start_recruit_size, param_set$start_recruit_sd, 
                param_set$egg_size_slope, param_set$egg_size_intercept, param_set$eyed_effect)
  
  # Find egg-recruit survival (recruits/egg) - RIGHT NOW, USING JOHNSON-LIKE ESTIMATE BUT COULD MAKE THIS A DISTRIBUTION TOO
  recruits_per_egg = param_set$recruits_per_egg
  
  # Find LEP in terms of recruits
  LEP_R = LEP*recruits_per_egg
  
  # Find connectivity matrix - EVENTUALLY, WILL USE CONFIDENCE INTERVALS AROUND DISPERSAL KERNELS TO DO THIS - FOR ALL-YEARS ONE? NOT SURE...
  #conn_matrix = Cmatrix
  conn_matrix <- matrix(NA,ncol=max(Cmat$org_geo_order, na.rm = TRUE), nrow=max(Cmat$org_geo_order, na.rm = TRUE))    
  for(i in 1:length(Cmat$org_site)) {
    column = Cmat$org_geo_order[i]  # column is origin 
    row = Cmat$dest_geo_order[i]  # row is destination
    conn_matrix[row, column] = Cmat$prob_disp[i]
  }
  
  # Make realized connectivity matrix
  conn_matrixR = conn_matrix*LEP_R
  
  # Assess network persistence
  eig_cR = eigen(conn_matrixR)
  
  # Pull out self-persistence values (diagonals of realized connectivity matrix)
  SP_values = data.frame(site = site_list, stringsAsFactors = FALSE) %>%
    mutate(SP_value = NA)
  for(i in 1:length(site_list)) {
    SP_values$SP_value[i] = conn_matrixR[i,i]  # pull out diagonal entries of realized connectivity matrix
  }
  
  # Put outputs together into one list
  out = list(NP = eig_cR$values[1], SP = SP_values, LEP = LEP, LEP_R = LEP_R , recruits_per_egg = recruits_per_egg, 
             conn_matrix = conn_matrix, conn_matrixR = conn_matrixR)
}

#################### Running things: ####################
##### Find 'best estimates' of some of the parameters
breeding_size_mean <- mean(female_sizes$size, na.rm=TRUE)

##### Generate sets of parameters
Linf_set = rnorm(n_runs, mean = Linf_growth_mean, sd=Linf_growth_sd)
Sint_set = rnorm(n_runs, mean = Sint_mean, sd= Sint_se)
k_connectivity_set = sample(k_connectivity_values, n_runs, replace=TRUE)  # replace should be true, right?
breeding_size_set = sample(female_sizes$size, n_runs, replace=TRUE)  # replace should be true, right?

# Put static + pulled-from-distribution parameters together into one dataframe
param_set_full <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = eggs_per_clutch_mean, clutches_per_year = clutches_per_year_mean,  # fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect = eyed_effect, # size-dependent fecundity info
         start_recruit_size = start_recruit_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Sl = Sl_mean, Linf = Linf_set, Sint = Sint_set,
         breeding_size = breeding_size_set, recruits_per_egg = recruits_per_egg,
         k_connectivity = k_connectivity_set, theta_connectivity = theta_allyears)  # dispersal kernel parameters
site_list <- c_mat_allyears$dest_site[1:19]

# Put best-estimate parameters into one dataframe
param_best_est <- data.frame(t_steps = n_tsteps) %>%
  mutate(min_size = min_size, max_size = max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = eggs_per_clutch_mean, clutches_per_year_mean = clutches_per_year_mean,  # average fecundity info
         egg_size_slope = eggs_slope_log, egg_size_intercept = eggs_intercept_log, eyed_effect =  eyed_effect,  # size-dependent fecundity info
         start_recruit_size = start_recruit_size, start_recruit_sd = start_recruit_sd,  # for initializing IPM with one recruit
         k_growth = k_growth_mean, s = s, Sl = Sl_mean, Linf = Linf_growth_mean, Sint = Sint_mean,
         breeding_size = breeding_size_mean, recruits_per_egg = recruits_per_egg,
         k_connectivity = k_allyears, theta_connectivity = theta_allyears)  # dispersal kernel parameters

##### Find the 'best-estimate' metrics
best_est_metrics <- calcMetrics(param_best_est, c_mat_allyears, site_list)

# Find LEP for mean breeding size and for size 6.0
LEP_breeding_size_mean <- findLEP(param_best_est$min_size, param_best_est$max_size, param_best_est$n_bins, 
                                  param_best_est$t_steps, param_best_est$Sint, param_best_est$Sl,
                                  param_best_est$s, param_best_est$Linf, param_best_est$k_growth, 
                                  param_best_est$eggs_per_clutch, param_best_est$clutches_per_year, 
                                  param_best_est$breeding_size, param_best_est$start_recruit_size, param_best_est$start_recruit_sd, 
                                  param_best_est$egg_size_slope, param_best_est$egg_size_intercept, param_best_est$eyed_effect)

LEP_6cm <- findLEP(param_best_est$min_size, param_best_est$max_size, param_best_est$n_bins, 
                   param_best_est$t_steps, param_best_est$Sint, param_best_est$Sl,
                   param_best_est$s, param_best_est$Linf, param_best_est$k_growth, 
                   param_best_est$eggs_per_clutch, param_best_est$clutches_per_year, 
                   6, param_best_est$start_recruit_size, param_best_est$start_recruit_sd, 
                   param_best_est$egg_size_slope, param_best_est$egg_size_intercept, param_best_est$eyed_effect)

LEP_ests <- list(LEP_breeding_size_mean = LEP_breeding_size_mean, LEP_6cm = LEP_6cm)

# Save as separate items, for plotting ease
LEP_best_est <- best_est_metrics$LEP
LEP_R_best_est <- best_est_metrics$LEP_R
NP_best_est <- best_est_metrics$NP
SP_best_est <- as.data.frame(best_est_metrics$SP)


#param_set_1 <- param_set_full[1,]
#test_calcMetrics <- calcMetrics(param_set_1, c_mat_allyears, site_list)  # this doesn't look right...

##### Run the metrics for lots of parameters
# Set output dataframes 
n_metrics = 4  # NP, LEP, LEP_R, recruits_per_egg
LEP_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP', inout = 'input', run = seq(1:n_runs))
LEP_R_out_df <- data.frame(value = rep(NA, n_runs), metric = 'LEP_R', inout = 'input', run = seq(1:n_runs))
RperE_out_df <- data.frame(value = rep(NA, n_runs), metric = 'recruits per egg', inout = 'input', run = seq(1:n_runs))
NP_out_df <- data.frame(value = rep(NA, n_runs), metric = 'NP', inout = 'output', run = seq(1:n_runs))

metric_vals <- data.frame(run = seq(1:n_runs), LEP = NA, LEP_R = NA, recruits_per_egg = NA, NP = NA)

# Create vector of sites for SP dataframe      
runsrepped = rep(1, length(site_list))                                                         
for(i in 2:n_runs) {
  runsrepped = c(runsrepped, rep(i, length(site_list)))
}                                 
           
SP_out_df <- data.frame(value = rep(NA, length(site_list)*n_runs), metric = 'SP', run = NA,
                        site = NA)

# Calculate the metrics for each parameter set, fill into the data frames
for(i in 1:n_runs) {
  # Select parameter set
  params <- param_set_full[i,]
  
  # Do the run
  metrics_output = calcMetrics(params, c_mat_allyears, site_list)
  
  # Fill in the metrics
  LEP_out_df$value[i] = metrics_output$LEP
  LEP_R_out_df$value[i] = metrics_output$LEP_R
  RperE_out_df$value[i] = metrics_output$recruits_per_egg
  NP_out_df$value[i] = metrics_output$NP
  
  metric_vals$LEP[i] = metrics_output$LEP
  metric_vals$LEP_R[i] = metrics_output$LEP_R
  metric_vals$recruits_per_egg[i] = metrics_output$recruits_per_egg
  metric_vals$NP[i] = metrics_output$NP
  
  # Pull out the SP metrics
  start_index = (i-1)*length(site_list)+1
  end_index = i*length(site_list)

  SP_out_df$site[start_index:end_index] = metrics_output$SP$site
  SP_out_df$value[start_index:end_index] = metrics_output$SP$SP_value
  SP_out_df$run[start_index:end_index] = rep(i, length(site_list))
}

# Put the data frames together (easier to plot?)
#metrics_vals <- rbind(LEP_out_df, LEP_R_out_df, RperE_out_df, NP_out_df)

# Add some of the changing parameters in, so can look at in plots
metric_vals_with_params <- metric_vals %>%
  mutate(breeding_size = param_set_full$breeding_size,
         Linf = param_set_full$Linf,
         Sint = param_set_full$Sint,
         k_connectivity = param_set_full$k_connectivity)

SP_vals_with_params <- left_join(SP_out_df, metric_vals_with_params, by='run') %>%
  dplyr::rename(SP = value)

# Pull out the self-persistence values

## Change size-dependent survival estimates from logit estimates in loaded file
#Mint = logit_recip(eall.Phi.size.p.dist.results$estimate[1])
#Ml = logit_recip(eall.Phi.size.p.dist.results$estimate[2])

#################### Plots: ####################

##### Plot the histograms of LEP, LEP_R, recruits_per_egg, and NP output
# LEP
pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_histogram.pdf'))
ggplot(data = LEP_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = 'gray', fill = 'gray') +
  geom_vline(xintercept = LEP_best_est, color='black') +
  xlab('LEP') + ggtitle('Histogram of LEP values') +
  theme_bw()
dev.off()

# LEP_R (LEP in terms of recruits)
pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_R_histogram.pdf'))
ggplot(data = LEP_R_out_df, aes(x=value)) +
  geom_histogram(bins=40, color = 'gray', fill = 'gray') +
  geom_vline(xintercept = LEP_R_best_est, color = 'black') +
  xlab('LEP_R') + ggtitle('Histogram of LEP_R values') +
  theme_bw()
dev.off()

# # Doing it this way looks a bit different... seems like there are more breaks, even though I also tried to ask for 30?
# pdf(file = here('Plots/PersistenceMetrics', 'LEP_R_histv2.pdf'))
# hist(LEP_R_out_df$value, breaks=30)
# dev.off()

# NP
pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'NP_histogram.pdf'))
ggplot(data = NP_out_df, aes(x=value)) +
  geom_histogram(bins=25, color='gray', fill='gray') +
  geom_vline(xintercept = NP_best_est, color='black') +
  xlab('NP') + ggtitle('Histogram of NP values') +
  theme_bw()
dev.off()

# recruits_per_egg (but right now this is static)
  
# hist(LEP_out_df$value, breaks=30)
# hist(LEP_R_out$value, breaks=30)
# hist(NP_out_df$value, breaks=30)
# hist(RperE_out_df$value, breaks=30)

##### SP at each site
pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty','SP_histogram.pdf'))
ggplot(data = SP_out_df, aes(x=value)) +
  #geom_histogram(binwidth=0.0005) +
  geom_histogram(binwidth=0.005, color='gray', fill='gray') +
  geom_vline(data=SP_best_est, aes(xintercept=SP_value), color='black') +
  ylim(0,300) +
  facet_wrap(~site) +
  xlab('SP') + ggtitle('Self-persistence histograms by site') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  #not sure why this isn't working right now...
dev.off()

##### Relationships between values
# LEP_R and NP - do we expect these to be related perfectly linearly? I guess, when both connectivity and recruits-per-egg are static...
pdf(file = here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'LEP_R_NP_scatter.pdf'))
ggplot(data = metric_vals, aes(x=LEP_R, y=NP)) +
  geom_point(size=2) +
  xlab('LEP_R') + ylab('NP') + ggtitle('Scatter of LEP_R vs NP values') +
  theme_bw()
dev.off()

# Breeding size and LEP
pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Breedingsize_LEP_scatter.pdf'))
ggplot(data = metric_vals_with_params, aes(x=breeding_size, y=LEP)) +
  geom_point(size=2) +
  xlab('breeding size') + ylab('LEP') + ggtitle('Scatter of breeding size (female) vs LEP values') +
  theme_bw()
dev.off()

# Linf and LEP
pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Linf_LEP_scatter.pdf'))
ggplot(data = metric_vals_with_params, aes(x=Linf, y=LEP)) +
  geom_point(size=2) +
  #geom_line()
  #geom_line(aes(x=sss)),
  xlab('Linf') + ylab('LEP') + ggtitle('Scatter of Linf vs LEP values') +
  theme_bw()
dev.off()

# k (connectivity) and NP
pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'kConnectivity_NP_scatter.pdf'))
ggplot(data = metric_vals_with_params, aes(x=k_connectivity, y=NP)) +
  geom_point(size=2) +
  #geom_line(aes(x=sss)),
  xlab('k_connectivity') + ylab('NP') + ggtitle('Scatter of k_connectivity vs NP values') +
  theme_bw()
dev.off()

##### Histograms of data inputs
# k_connectivity
pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'k_connectivity_histogram.pdf'))
ggplot(data = metric_vals_with_params, aes(x=k_connectivity)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = k_allyears, color='black') +
  xlab('k_connectivity') + ggtitle('k_connectivity values') +
  theme_bw()
dev.off()

# Linf
pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Linf_histogram.pdf'))
ggplot(data = metric_vals_with_params, aes(x=Linf)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = Linf_growth_mean, color='black') +
  xlab('Linf') + ggtitle('Linf values') +
  theme_bw()
dev.off()

# Sint
pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Sint_histogram.pdf'))
ggplot(data = metric_vals_with_params, aes(x=Sint)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept = Sint_mean, color = 'black') +
  xlab('Sint') + ggtitle('Sint values') +
  theme_bw()
dev.off()

# Breeding size
pdf(file =  here('Plots/PersistenceMetrics/MetricsWithUncertainty', 'Breeding_size_histogram.pdf'))
ggplot(data = metric_vals_with_params, aes(x=breeding_size)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept=breeding_size_mean, color='black') +
  xlab('Breeding size') + ggtitle('Breeding size values') +
  theme_bw()
dev.off()

#################### Saving things: ####################
save(best_est_metrics, file=here('Data', 'best_est_metrics.RData'))
save(LEP_ests, file=here('Data', 'LEP_ests.RData'))

#################### Old code: ####################
# Sint_mean = eall.Phi.size.p.dist.results$estimate[1]  # survival intercept (on logit scale)
# Sl_mean = eall.Phi.size.p.dist.results$estimate[2]  # survival slope (on logit scale)
# Sint_se = eall.Phi.size.p.dist.results$se[1]  # for now using SE, should really use SD...
# Sint_se = eall.Phi.size.p.dist.results$se[2]  # for now using SE, should really use SD...

# # Breeding size (for LEP) - replacing with drawing from the actual data
# breeding_size_mean = (size_by_color_metrics %>% filter(color == 'YP'))$mean  # originally guessed 8, this is 8.6
# breeding_size_sd = (size_by_color_metrics %>% filter(color == 'YP'))$sd  # originally guessed 0.8, this is 1.6


# # Static connectivity matrix loaded from PersistenceMetrics.R
# Cmatrix <- matrix(NA,ncol=max(c_mat_allyears$org_geo_order, na.rm = TRUE), nrow=max(c_mat_allyears$org_geo_order, na.rm = TRUE))    
# for(i in 1:length(c_mat_allyears$org_site)) {
#   column = c_mat_allyears$org_geo_order[i]  # column is origin 
#   row = c_mat_allyears$dest_geo_order[i]  # row is destination
#   Cmatrix[row, column] = c_mat_allyears$prob_disp_allyears[i]
# }

# # k (connectivity) and SP
# pdf(file = here('Plots/PersistenceMetrics','kConnectivity_SP_scatter.pdf'))
# ggplot(data = SP_vals_with_params, aes(x=k_connectivity, y=SP)) +
#   geom_point(size=2) +
#   facet_wrap(~site) +
#   xlab('k_connectivity') + ylab('SP') + ggtitle('Scatter of k_connectivity vs SP by site') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  #not sure why this isn't working right now...
# dev.off()

# # Find k standard deviation - this is definitely not right - check with KC about her confidence intervals and how they were estimated
# z_97.5 = 2.24  # z-score for 97.5% confidence interval
# k_sdH = (sqrt(n_runs)/z_97.5)*(k_allyears_CIh - k_allyears)
# k_sdL = -(sqrt(n_runs)/z_97.5)*(k_allyears_CIl - k_allyears)


#k_connectivity_set = runif(n_runs, min = k_allyears_CIl, max = k_allyears_CIh)  # for now, just selecting randomly from within the 97.5% confidence interval
#breeding_size_set = rnorm(n_runs, mean = breeding_size_mean, sd = breeding_size_sd) 

