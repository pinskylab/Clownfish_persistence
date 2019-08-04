# Make plots for Malin NSF proposal

#################### Set-up: ####################
## Load libraries
library(dplyr)
library(ggplot2)
library(cowplot)

## Load data/model output
# Load simple VBL growth analysis
load(file = here::here("Data/Script_outputs", "growth_info_estimate.RData"))
load(file = here::here("Data/Script_outputs", "recap_pairs_year.RData"))  # all recap pairs a year apart, for plotting purposes

# Load survival from MARK models
survival_output = readRDS(file=here::here("Data/Script_outputs", "eall_mean_Phi_size_p_size_plus_dist.RData"))  # MARK output (lowest AICc model)
Phi_int_pos = 1  # placement of intercept for survival
Phi_size_pos = 2  # placement of size effect for survival
p_int_pos = 3  # placement of intercept for recap prob
p_size_pos = 4  # placement of size effect for recap prob
p_dist_pos = 5  # placement of distance effect for recap prob

# Size transition info (Michelle analysis in genomics repo) - switch this to just females from males (is that reasonable?)
recap_first_female = readRDS(file=here::here("Data/From_other_analyses", "recap_first_female.RData"))

# Set some constants
min_size_plot = 0
max_size = 15
size.values <- min_size_plot+(0:30)*(max_size-min_size_plot)/30

t0 = 0  # age at size 0 (for vBL)

age_values <- seq(from=0, to=10, by=0.5)

#################### Functions: ####################
# Finds the real parameter estimate from the logit estimate
logit_recip <- function(logitval) {
  recip = (exp(logitval))/(1 + exp(logitval))
  return(recip)
}

# von Bertalanffy growth curve 
vBL_size_at_age <- function(Linf, k, t0, t) {
  length_out <- Linf*(1 - exp(-k*(t-t0)))
  return(length_out)
}

#################### Plots: ####################
##### Growth and VBL plot
## a) size-based survival (4c in persistence draft) (code from Metrics_with_uncertainty)

#eall_mean.Phi.size.p.size.plus.dist.results = as.data.frame(eall_mean.Phi.size.p.size.plus.dist$results$beta)
survival_output_to_plot <- data.frame(size = size.values) %>%
  mutate(Phi_logit = survival_output$estimate[Phi_int_pos] + survival_output$estimate[Phi_size_pos]*size,
         Phi_lcl_logit = survival_output$lcl[Phi_int_pos] + survival_output$lcl[Phi_size_pos]*size,
         Phi_ucl_logit = survival_output$ucl[Phi_int_pos] + survival_output$ucl[Phi_size_pos]*size,
         Phi = logit_recip(Phi_logit),
         Phi_lcl = logit_recip(Phi_lcl_logit),
         Phi_ucl = logit_recip(Phi_ucl_logit))

survival_plot <- ggplot(data = survival_output_to_plot, aes(size, Phi)) +
  geom_ribbon(aes(ymin=Phi_lcl,ymax=Phi_ucl),color="gray",fill="gray") +
  geom_line(color="black") +
  xlab("size (cm)") + ylab("probability of survival") + ggtitle("Annual survival") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

## b) size vs. age from estimated von Bertalanffy curve

Linf_mean <- mean(growth_info_estimate$Linf_est)
k_mean <- mean(growth_info_estimate$k_est)
Linf_max <- max(growth_info_estimate$Linf_est)
Linf_min <- min(growth_info_estimate$Linf_est)
k_max <- max(growth_info_estimate$k_est)
k_min <- min(growth_info_estimate$k_est)

length_at_age_best_est <- vBL_size_at_age(Linf_mean, k_mean, t0, age_values)
length_at_age_low <- vBL_size_at_age(Linf_min, k_min, t0, age_values)
length_at_age_high <- vBL_size_at_age(Linf_max, k_max, t0, age_values)

length_age_to_plot <- data.frame(age = age_values) %>%
  mutate(size = length_at_age_best_est,
         size_low = length_at_age_low,
         size_high = length_at_age_high)

vBL_plot <- ggplot(data = length_age_to_plot, aes(age, size)) +
  geom_ribbon(aes(ymin=size_low, ymax=size_high),color="gray",fill="gray") +
  geom_line(color="black") +
  xlab("age (years)") + ylab("size (cm)") + ggtitle("Size-at-age") +
  theme_bw()

## c) female transition size (plot 4d from Persistence draft)

# Breeding size distribution
breeding_size_plot <- ggplot(data = recap_first_female, aes(x=size)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept=mean(recap_first_female$size), color='black') +
  xlab('size (cm)') + ggtitle('Female transition size') +
  theme_bw()


## arrange together
pdf(file = here::here("Plots/NSF_proposal_plots", "Survival_growth_females.pdf"), width=8, height=4)
plot_grid(survival_plot, vBL_plot, breeding_size_plot, 
          labels = c("a","b","c"), nrow=1)
dev.off()


# Would you be able to put together a plot for the proposal with 
# three parts, arranged horizontally? 
#   It would be a) plot 4c from the current persistence draft; 
# b) a plot of size vs. age from the estimated von Bertalanffy curve,
# with gray confidence interval (similar style to 4c); and 
# c) plot 4d from the persistence draft. 
# 
