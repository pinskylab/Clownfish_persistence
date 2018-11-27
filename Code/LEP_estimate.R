# Calculate LEP from mark-recapture data using an IPM

#################### Set-up: ####################
## Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(here)

## Load input files
load(file=here("Data", "eall_Phi_size_p_dist_results.RData")) #MARK output 

## Set various constants
# Set size-dependent growth (right now from Will's prelim estimates, need to refine growth)
# Linf = 10.68
# k = 1.078
s = exp(-0.0148)

k_lowestAIC = 0.9447194
Linf_lowestAIC = 10.50670

Linf = Linf_lowestAIC
k = k_lowestAIC

# Other size points
min_size = 0
max_size = 15 #should check this w/data...
breeding_size = 8

# Egg info (mix of from Adam and from Will)
eggs_per_clutch_AY = 514.11 #from AY poster work
eggs_per_clutch_WW = 1763
clutches_per_year = 11.9 #from WW, should check with latest data
female_egg_slope = 107 #eggs/cm size of F, from AY poster work
egg_b_young = -426.57 #young eggs intercept, from AY poster work
egg_b_old = -585.93 #old eggs intercept, from AY poster work

# Size-dependent survival
Sint = eall.Phi.size.p.dist.results$estimate[1] #survival intercept (on logit scale)
Sl = eall.Phi.size.p.dist.results$estimate[2] #survival slope (on logit scale)

# Set params for IPM structure
n_bins = 100
n_tsteps = 100

#################### Functions: ####################
# Functions and constants from my GitHub function/constant collection
# script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQAXExRItr2tnepmkRr7NRt4hylZrks5aciBtwA%3D%3D", ssl.verifypeer = FALSE)
# eval(parse(text = script))
script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Functions from Michelle's GitHub helpers script
#field_helpers (similar idea to helpers, in field repository) - this might be the newer version of the helpers collection?
script <- getURL("https://raw.githubusercontent.com/pinskylab/field/master/scripts/field_helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Finds the real parameter estimate from the logit estimate
logit_recip <- function(logitval) {
  recip = (exp(logitval))/(1 + exp(logitval))
  return(recip)
}

#################### Running things: ####################
# Create vector of lengths
lengths_vec = seq(min_size, max_size, length.out=n_bins)
dx = diff(lengths_vec[1:2])

# Make matrix
xmat = matrix(rep(lengths_vec,n_bins), nrow=n_bins, ncol=n_bins, byrow=TRUE) #100 rows of lengths_vec
ymat = t(xmat)

# Survival probability (based on size)
S = logit_recip(Sint + lengths_vec*Sl) #is this right? Also, this should be prob of surviving, not prob of dying, right? (Which it is right now)
Smat = matrix(rep(S,n_bins), nrow=n_bins, ncol=n_bins, byrow = TRUE) #survival part

# Growth probability (this is the part that really could use updating!)
Ls = Linf - (Linf - xmat)*exp(-k)
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
eggs <- eggs_per_clutch_AY*clutches_per_year
eggs_Will <- eggs_per_clutch_WW*clutches_per_year

# Combine to get LEP
LEP <- sum(breeding_Adults*eggs)
LEP_W <- sum(breeding_Adults*eggs_Will)

LEP_out = data.frame(LEP = LEP, LEP_W = LEP_W)

## Change size-dependent survival estimates from logit estimates in loaded file
#Mint = logit_recip(eall.Phi.size.p.dist.results$estimate[1])
#Ml = logit_recip(eall.Phi.size.p.dist.results$estimate[2])

#################### Plots: ####################

#################### Saving output: ####################
save(LEP_out, file=here("Data","LEP_out.RData")) #save data frame with multiple versions of LEP
