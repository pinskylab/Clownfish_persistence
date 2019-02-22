# Characterizing uncertainty in LEP, recruit survival, dispersal estimates

#################### Set-up: ####################
# source(here::here('Code', 'Constants_database_common_functions.R'))

load(file=here("Data", "eall_Phi_size_p_dist_results.RData")) #MARK output 
load(file=here('Data', 'c_mat_allyears.RData'))  # Probability of dispersing (for C matrix for now, before use kernel params to include uncertainty)

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

# Egg-recruit survival (for getting LEP in terms of recruits)
recruits_per_egg = 8.367276e-05  # surv_egg_recruit estimating using Johnson method in PersistenceMetrics.R


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
breeding_size_set = rnorm(n_runs, mean = breeding_size_mean, sd = breeding_size_sd) 

##### Other parameters that stay static

#################### Functions: ####################

# Find LEP (SHOULD CHECK, UPDATE THIS!)
findLEP = function(min_size, max_size, n_bins, t_steps, Sint, Sl, s, Linf, k_growth, eggs_per_clutch, clutches_per_year, breeding_size) {
  
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

# Run through a metric calculation with one set of parameters
calcMetrics <- function(param_set, Cmatrix, sites) {
  
  # Find LEP (in terms of eggs) - COULD MAKE MORE OF THESE PARAMETERS PULLED FROM A DISTRIBUTION!
  LEP = findLEP(param_set$min_size, param_set$max_size, param_set$n_bins, param_set$t_steps, param_set$Sint, param_set$Sl,
                param_set$s, param_set$Linf, param_set$k_growth, param_set$eggs_per_clutch, param_set$clutches_per_year, param_set$breeding_size)
  
  # Find egg-recruit survival (recruits/egg) - RIGHT NOW, USING JOHNSON-LIKE ESTIMATE BUT COULD MAKE THIS A DISTRIBUTION TOO
  recruits_per_egg = param_set$recruits_per_egg
  
  # Find LEP in terms of recruits
  LEP_R = LEP*recruits_per_egg
  
  # Find connectivity matrix - EVENTUALLY, WILL USE CONFIDENCE INTERVALS AROUND DISPERSAL KERNELS TO DO THIS - FOR ALL-YEARS ONE? NOT SURE...
  conn_matrix = Cmatrix
  
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

##### Put static + pulled-from-distribution parameters together into one dataframe
param_set_full <- data.frame(t_steps = rep(n_tsteps, n_runs)) %>%
  mutate(min_size = min_size, max_size=max_size, n_bins = n_bins,  # LEP-IPM matrix
         eggs_per_clutch = eggs_per_clutch_mean, clutches_per_year = clutches_per_year_mean,
         k_growth = k_growth_mean, s = s, Sl = Sl_mean, Linf = Linf_set, Sint = Sint_set,
         breeding_size = breeding_size_set, recruits_per_egg = recruits_per_egg)
site_list <- c_mat_allyears$dest_site[1:19]

# Static connectivity matrix loaded from PersistenceMetrics.R
Cmatrix <- matrix(NA,ncol=max(c_mat_allyears$org_geo_order, na.rm = TRUE), nrow=max(c_mat_allyears$org_geo_order, na.rm = TRUE))    
for(i in 1:length(c_mat_allyears$org_site)) {
  column = c_mat_allyears$org_geo_order[i]  # column is origin 
  row = c_mat_allyears$dest_geo_order[i]  # row is destination
  Cmatrix[row, column] = c_mat_allyears$prob_disp_allyears[i]
}
         
# For testing calcMetrics function        
param_set_1 <- param_set_full[1,]
test_calcMetrics <- calcMetrics(param_set_1, Cmatrix, site_list)

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
  metrics_output = calcMetrics(params, Cmatrix, site_list)
  
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
         Sint = param_set_full$Sint)

# Pull out the self-persistence values

## Change size-dependent survival estimates from logit estimates in loaded file
#Mint = logit_recip(eall.Phi.size.p.dist.results$estimate[1])
#Ml = logit_recip(eall.Phi.size.p.dist.results$estimate[2])

#################### Plots: ####################

##### Plot the histograms of LEP, LEP_R, recruits_per_egg, and NP output
# LEP
pdf(file = here('Plots/PersistenceMetrics', 'LEP_histogram.pdf'))
ggplot(data = LEP_out_df, aes(x=value)) +
  geom_histogram(bins=25) +
  xlab('LEP') + ggtitle('Histogram of LEP values') +
  theme_bw()
dev.off()

# LEP_R (LEP in terms of recruits)
pdf(file = here('Plots/PersistenceMetrics', 'LEP_R_histogram.pdf'))
ggplot(data = LEP_R_out_df, aes(x=value)) +
  geom_histogram(bins=30) +
  xlab('LEP_R') + ggtitle('Histogram of LEP_R values') +
  theme_bw()
dev.off()

# # Doing it this way looks a bit different... seems like there are more breaks, even though I also tried to ask for 30?
# pdf(file = here('Plots/PersistenceMetrics', 'LEP_R_histv2.pdf'))
# hist(LEP_R_out_df$value, breaks=30)
# dev.off()

# NP
pdf(file = here('Plots/PersistenceMetrics', 'NP_histogram.pdf'))
ggplot(data = NP_out_df, aes(x=value)) +
  geom_histogram(bins=25) +
  xlab('NP') + ggtitle('Histogram of NP values') +
  theme_bw()
dev.off()

# recruits_per_egg (but right now this is static)
  
# hist(LEP_out_df$value, breaks=30)
# hist(LEP_R_out$value, breaks=30)
# hist(NP_out_df$value, breaks=30)
# hist(RperE_out_df$value, breaks=30)

##### SP at each site
pdf(file = here('Plots/PersistenceMetrics','SP_histogram.pdf'))
ggplot(data = SP_out_df, aes(x=value)) +
  geom_histogram(binwidth=0.0005) +
  facet_wrap(~site) +
  xlab('SP') + ggtitle('Self-persistence histograms by site') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  #not sure why this isn't working right now...
dev.off()

##### Relationships between values
# LEP_R and NP - do we expect these to be related perfectly linearly? I guess, when both connectivity and recruits-per-egg are static...
pdf(file = here('Plots/PersistenceMetrics', 'LEP_R_NP_scatter.pdf'))
ggplot(data = metric_vals, aes(x=LEP_R, y=NP)) +
  geom_point(size=2) +
  xlab('LEP_R') + ylab('NP') + ggtitle('Scatter of LEP_R vs NP values') +
  theme_bw()
dev.off()

# Breeding size and LEP
pdf(file =  here('Plots/PersistenceMetrics', 'Breedingsize_LEP_scatter.pdf'))
ggplot(data = metric_vals_with_params, aes(x=breeding_size, y=LEP)) +
  geom_point(size=2) +
  xlab('breeding size') + ylab('LEP') + ggtitle('Scatter of breeding size (female) vs LEP values') +
  theme_bw()
dev.off()

# Linf and LEP
pdf(file =  here('Plots/PersistenceMetrics', 'Linf_LEP_scatter.pdf'))
ggplot(data = metric_vals_with_params, aes(x=Linf, y=LEP)) +
  geom_point(size=2) +
  geom_line(aes(x=sss),
  xlab('Linf') + ylab('LEP') + ggtitle('Scatter of Linf vs LEP values') +
  theme_bw()
dev.off()

#################### Saving things: ####################

#################### Old code: ####################