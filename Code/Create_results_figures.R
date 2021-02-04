# Creating figures of persistence metrics

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))  # Pull in common constants, functions, saved data files, and processed data files 

##### Load libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

##### Load input from other analyses outside this repository
# # Size transition info (Michelle analysis in genomics repo) 
# recap_first_female = readRDS(file=here::here("Data/From_other_analyses", "recap_first_female.RData"))  # sizes of fish when they were first caught as females

##### Needed functions (from Metrics_with_uncertainty_site_specific.R)
# For flexible theta and k, function of d, k, and theta
disp_kernel_all_years <- function(d, k, theta) {  # generalization of equation for p(d) in eqn. 6c in Bode et al. 2018
  z = exp(k)
  z_front = (z*theta)/gamma(1/theta)
  disp = (z_front/2)*exp(-(z*d)^(theta))
  return(disp)
}

# # Growth (in case want to show a VBL growth plot)
# VBL_growth <- function(Linf, k_growth, length) {
#   Ls = Linf - (Linf - length)*exp(-k_growth)
#   return(Ls)
# }


##### Load analyses and parameters from within this repository

## Output from abundance trend (from Code/TimeSeriesPersistence.R)
load(file = here::here("Data/Script_outputs", "site_trends_all.RData"))
load(file = here::here("Data/Script_outputs", "site_trends_time.RData"))

## Load simple VBL growth analysis (from Code/Growth_analysis.R)
load(file = here::here("Data/Script_outputs", "growth_info_estimate.RData"))
load(file = here::here("Data/Script_outputs", "recap_pairs_year.RData"))  # all recap pairs a year apart, for plotting purposes

## Load survival and recap from MARK models (from Code/ClownfishMarkRecap.R)
load(file = here::here("Data/Script_outputs", "best_fit_model_dfs.RData"))  # load best MARK model

## Relevant intermediate outputs
load(file=here::here("Data/Script_outputs","all_parents_site.RData"))  # info on parents by site, distance to edges of sampling region

## Parameter sets (from Metrics_with_uncertainty_site_specific.R)
load(file=here::here("Data/Script_outputs", "param_best_est_mean_collected_offspring.RData"))  # point estimates
load(file=here::here("Data/Script_outputs", "param_set_full.RData"))  # parameter set with uncertainty values

### Point estimates (from Metrics_with_uncertainty_site_specific.R)
load(file=here::here("Data/Script_outputs", "best_est_metrics_mean_offspring_DD.RData"))  # point estimates with density dependence compensation (presented in main text)
load(file=here::here("Data/Script_outputs", "best_est_metrics_mean_offspring.RData"))  # point estimates without density dependence compensation (presented in appendix)

### Uncertainty estimates with density dependence compensation (presented in main text) (from Metrics_with_uncertainty_site_specific.R)
load(file=here::here("Data/Script_outputs", "output_uncert_start_recruit_DD.RData"))
load(file=here::here("Data/Script_outputs", "output_uncert_growth_DD.RData"))
load(file=here::here("Data/Script_outputs", "output_uncert_survival_DD.RData"))
load(file=here::here("Data/Script_outputs", "output_uncert_breeding_size_DD.RData"))
load(file=here::here("Data/Script_outputs", "output_uncert_offspring_assigned_DD.RData"))
load(file=here::here("Data/Script_outputs", "output_uncert_prob_r_DD.RData"))
load(file=here::here("Data/Script_outputs", "output_uncert_dispersal_DD.RData"))
load(file=here::here("Data/Script_outputs", "output_uncert_all_DD.RData"))

### Uncertainty estimates without density dependence compensation (presented in appendix) (from Metrics_with_uncertainty_site_specific.R)
load(file=here::here("Data/Script_outputs", "output_uncert_all.RData"))

### Uncertainty contribution data frames for plotting (presented in appendix) (from Metrics_with_uncertainty_site_specific.R)
load(file=here::here("Data/Script_outputs","LEP_uncert_DD.RData"))
load(file=here::here("Data/Script_outputs","LEP_R_uncert_DD.RData"))
load(file=here::here("Data/Script_outputs","RperE_uncert_DD.RData"))
load(file=here::here("Data/Script_outputs","NP_uncert_DD.RData"))

### What-ifs (from Metrics_with_uncertainty_site_specific.R)
# Alternative geographies - sensitivity to percent habitat
load(file=here::here("Data/Script_outputs","perc_hab_plot_df.RData"))  # perc hab sensitivity in easy-to-plot form
load(file=here::here("Data/Script_outputs","perc_hab_best_est_NP_df.RData"))

# Alternative geographies - sensitivity to wider region
load(file=here::here("Data/Script_outputs","wider_region_plot_df.RData"))  # in an easy-to-plot form
load(file=here::here("Data/Script_outputs","wider_region_best_est_NP_df.RData"))

# Alternative geographies - sensitivity to wider region and percent habitat simultaneously
load(file=here::here("Data/Script_outputs","wider_region_perc_hab_plot_df_categories.RData"))  # with NP values lumped into 5 categories

# Larval navigation sensitivity
load(file=here::here("Data/Script_outputs","larv_nav_plot_df.RData"))  # in an easy-to-plot form
load(file=here::here("Data/Script_outputs","larv_nav_best_est_NP_df.RData"))

# Load params and metrics summaries
load(file=here::here("Data/Script_outputs","params_summary.RData"))
load(file=here::here("Data/Script_outputs","metrics_summary.RData"))

### Vectors and values useful for producing figures, like region widths for sensitivity, etc. (from Metrics_with_uncertainty_site_specific.R)
load(file=here::here("Data/Script_outputs","figure_vectors_and_values.RData"))

# Parameters for analyzing and interpreting results
NP_decimal_points <- 2  # numbers past decimal point for NP range
quantile_lower_index = 25  # lower index for 95% quantile reporting
quantile_upper_index = 975  # upper index for 95% quantile reporting

#################### Plots: ####################

########## Main text figures ##########

##### Figure 1 - life cycle and metrics schematic - made in SiteMap.R script

##### Figure 2 - map + photo - made in SiteMap.R script

##### Figure 3 - demographic and dispersal inputs: survival curve, growth curve, dispersal kernel, transition size to female)
# Prep work for dispersal kernel figure
connectivity_est_vec <- disp_kernel_all_years(figure_vectors_and_values$distance_vec, param_best_est_mean_collected_offspring$k_connectivity, param_best_est_mean_collected_offspring$theta_connectivity)  # vector of dispersal probs to plot
connectivity_sens <- matrix(ncol = length(figure_vectors_and_values$distance_vec), nrow = n_runs)
for(i in 1:n_runs) {
  connectivity_sens[i,] = disp_kernel_all_years(figure_vectors_and_values$distance_vec, param_set_full$k_connectivity[i], param_set_full$theta_connectivity[i])
}
connectivity_lowest_bound = apply(connectivity_sens, 2, min)
connectivity_highest_bound = apply(connectivity_sens, 2, max)
dispersal_df <- data.frame(distance = figure_vectors_and_values$distance_vec, kernel_bestfit = connectivity_est_vec, kernel_CI1 = connectivity_lowest_bound, kernel_CI2 = connectivity_highest_bound)

# Dispersal kernel plot
dispersal_kernel_plot <- ggplot(data=dispersal_df, aes(x=distance, y=kernel_bestfit, ymin=kernel_CI1, ymax=kernel_CI2)) +
  geom_line(color='black') +
  geom_ribbon(alpha=0.5, color='gray') +
  xlab('Distance (km)') + ylab('Dispersal probability') +
  theme_bw() 

# Growth curve - plot the data and linear model for fish caught about a year apart, models for only one recapture pair per fish, pairs selected randomly and models fit 1000x
growth_curve_plot <- ggplot(data = recap_pairs_year, aes(x = L1, y = L2)) +
  geom_point(shape = 1, color = "dark gray") +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", lwd = 0.5, linetype = "dashed") +  #  1:1 line
  geom_ribbon(aes(x=seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1)),
                  ymin = (params_summary %>% filter(param == "growth_intercept_est_lower"))$value + (params_summary %>% filter(param == "growth_slope_est_lower"))$value*seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1)),
                  ymax = (params_summary %>% filter(param == "growth_intercept_est_upper"))$value + (params_summary %>% filter(param == "growth_slope_est_upper"))$value*seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1))), fill = "dark gray", alpha = 0.75) +
                  #ymin = (min(growth_info_estimate$intercept_est) + min(growth_info_estimate$slope_est)*seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1))),
                  #ymax = (max(growth_info_estimate$intercept_est) + max(growth_info_estimate$slope_est)*seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1)))), fill = "dark gray", alpha = 0.75) +
  geom_abline(aes(intercept = mean(growth_info_estimate$intercept_est), slope = mean(growth_info_estimate$slope_est)), color = "black", lwd = 1.5) +
  xlab("Length (cm)") + ylab("Length (cm) next year") + 
  theme_bw()

# # Alternate version if want to do a VBL instead
# growth_df <- data.frame(length1 = seq(min_size, max_size, length.out = 100)) %>%
#   mutate(length2 = VBL_growth(param_best_est_mean_collected_offspring$Linf, param_best_est_mean_collected_offspring$k_growth, length1))
# 
# growth_curve_plot <- ggplot(data=growth_df, aes(x=length1, y=length2)) +
#   geom_point(color='black') +
#   xlab('Size (cm)') + ylab('Size next year') +
#   ylim(0,13) +
#   theme_bw()

# Breeding size distribution
breeding_size_plot <- ggplot(data = recap_first_female, aes(x=size)) +
  geom_histogram(bins=40, color='gray', fill='gray') +
  geom_vline(xintercept=param_best_est_mean_collected_offspring$breeding_size, color='black') +
  xlab('Transition size (cm)') + ylab("Count") +
  theme_bw()

# Survival plot - Elementary School, which has median survival, as an example site
survival_output_to_plot <- best_fit_model_dfs$surv_site_size %>% filter(site == "Elementary School")

survival_plot <- ggplot(data = survival_output_to_plot, aes(size, estimate_prob)) +
  geom_ribbon(aes(ymin=lcl_prob,ymax=ucl_prob),color="gray",fill="gray") +
  geom_line(color="black") +
  #xlab("Size (cm)") + ylab("Probability of survival \n (ex. site: Elementary School)") + 
  xlab("Size (cm)") + ylab(expression(paste("Probability of survival (", phi[i], ")"))) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()

# Put all four subplots together
pdf(file = here('Plots/FigureDrafts','Parameter_inputs.pdf'), width=6, height=6)
plot_grid(dispersal_kernel_plot, growth_curve_plot, survival_plot, breeding_size_plot, labels = c("a","b","c","d"), nrow=2)
dev.off()

##### Figure 4 - abundance + replacement metrics

# LEP - uncertainty for site-specific LEP values (LEP_i), best estimate for averaged across sites (LEP_*)
LEP_plot_freq <- ggplot(data = output_uncert_all$LEP_by_site_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=100, color = 'light gray', fill = 'light gray') +
  geom_vline(xintercept = mean(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP), color = "black") +
  xlab("Lifetime egg production (LEP"[i]~")") + ylab("Relative frequency") +
  theme_bw()

# LRP averaged across sites with DD compensation
LEP_R_95_lower <- (metrics_summary %>% filter(metric == "LRP_DD avg lower"))$value
LEP_R_95_upper <- (metrics_summary %>% filter(metric == "LRP_DD avg upper"))$value
LEP_R_vals_95 <- output_uncert_all_DD$LEP_R_out_df %>% filter(value >= LEP_R_95_lower & value <= LEP_R_95_upper)

LEP_R_plot_DD_freq <- ggplot(data = output_uncert_all_DD$LEP_R_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=50, color = 'light gray', fill = 'light gray') +
  geom_histogram(data = LEP_R_vals_95, aes(x=value, y=..count../sum(..count..)), bins=50, color = "dark gray", fill = "dark gray") +
  geom_vline(xintercept = best_est_metrics_mean_offspring_DD$LEP_R_mean, color = "black") +
  xlab(bquote("Lifetime recruit production (LRP)")) + ylab("Relative frequency") +
  theme_bw()

# LEP_R_plot_DD_freq <- ggplot(data = output_uncert_all_DD$LEP_R_out_df, aes(x=value)) +
#   geom_histogram(aes(y=..count../sum(..count..)), bins=50, color = 'gray', fill = 'gray') +
#   geom_vline(xintercept = best_est_metrics_mean_offspring_DD$LEP_R_mean, color = "black") +
#   xlab(bquote("Lifetime recruit production (LRP)")) + ylab("Relative frequency") + 
#   theme_bw()

# LR with DD
LR_95_lower <- (metrics_summary %>% filter(metric == "LR avg lower DD"))$value
LR_95_upper <- (metrics_summary %>% filter(metric == "LR avg upper DD"))$value
LR_vals_95 <- output_uncert_all_DD$LEP_R_local_out_df %>% filter(value >= LR_95_lower & value <= LR_95_upper)

LEP_R_local_plot_DD_freq <- ggplot(data = output_uncert_all_DD$LEP_R_local_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color = "light gray", fill = "light gray") +
  geom_histogram(data = LR_vals_95, aes(x=value, y=..count../sum(..count..)), bins=40, color = "dark gray", fill = "dark gray") +
  geom_vline(xintercept = best_est_metrics_mean_offspring_DD$LEP_R_local_mean, color = "black") +
  xlab(bquote("Local replacement (LR)")) + ylab("Relative frequency") +
  theme_bw()

# LEP_R_local_plot_DD_freq <- ggplot(data = output_uncert_all_DD$LEP_R_local_out_df, aes(x=value)) +
#   geom_histogram(aes(y=..count../sum(..count..)), bins=40, color = "gray", fill = "gray") +
#   geom_vline(xintercept = best_est_metrics_mean_offspring_DD$LEP_R_local_mean, color = "black") +
#   xlab(bquote("Local replacement (LR)")) + ylab("Relative frequency") +
#   theme_bw()

# Abundance trend 
Fig4_abundance_plot <- ggplot(data = site_trends_time, aes(x=year, y=mean_nF, group=site)) +
  geom_line(color="grey") +
  geom_line(data=site_trends_all, aes(x=year, y=mean_nF), color = "black", size=1.5) +
  xlab("Year") + ylab("# Females") + 
  scale_x_continuous(breaks=c(2,4,6), labels=c("2013","2015","2017")) +
  theme_bw()

# Put them together
pdf(file = here::here('Plots/FigureDrafts', 'Abundance_LEP_LRP_LocalReplacement_FreqPlots.pdf'), width=6, height=6)
plot_grid(LEP_plot_freq, LEP_R_plot_DD_freq, LEP_R_local_plot_DD_freq, Fig4_abundance_plot,
          labels = c("a","b","c","d"), nrow=2)
dev.off()

##### Figure 5 - SP metrics, NP metrics, connectivity matrices
# Find 95% quantiles for SP
SP_vals_95 <- data.frame(site = site_vec_NS, SP_lower = rep(NA, length(site_vec)), SP_upper = rep(NA, length(site_vec)),
                         SP_min = rep(NA, length(site_vec)), SP_max =  rep(NA, length(site_vec)), stringsAsFactors = FALSE)

for(i in 1:length(site_vec_NS)) {
  SP_site_vals <- (output_uncert_all_DD$SP_vals_with_params %>% filter(site == site_vec_NS[i]))$SP
  SP_site_vals_vec <- sort(SP_site_vals)
  SP_vals_95$site[i] <- site_vec_NS[i]
  SP_vals_95$SP_lower[i] <- SP_site_vals_vec[quantile_lower_index]
  SP_vals_95$SP_upper[i] <- SP_site_vals_vec[quantile_upper_index]
  SP_vals_95$SP_min[i] <- min(SP_site_vals_vec)
  SP_vals_95$SP_max[i] <- max(SP_site_vals_vec)
}

SP_vals_95_level_order <- c("Palanas","Wangag","N. Magbangon", "S. Magbangon", "Cabatoan", "Caridad Cemetery",
                            "Caridad Proper", "Hicgop South", "Sitio Tugas", "Elementary School", "Sitio Lonas",
                            "San Agustin", "Poroc San Flower", "Poroc Rose", "Visca", "Gabas",
                            "Tomakin Dako", "Haina", "Sitio Baybayon")

SP_vals_95$site <- replace(SP_vals_95$site, SP_vals_95$site == "Tamakin Dacot", "Tomakin Dako")
best_est_metrics_mean_offspring_DD$Cmat$org_site <- replace(best_est_metrics_mean_offspring_DD$Cmat$org_site, 
                                                            best_est_metrics_mean_offspring_DD$Cmat$org_site=="Tamakin Dacot", 
                                                            "Tomakin Dako")
best_est_metrics_mean_offspring_DD$Cmat$dest_site <- replace(best_est_metrics_mean_offspring_DD$Cmat$dest_site, 
                                                             best_est_metrics_mean_offspring_DD$Cmat$dest_site=="Tamakin Dacot", 
                                                             "Tomakin Dako")
output_uncert_all_DD$SP_vals_with_params$site <- replace(output_uncert_all_DD$SP_vals_with_params$site, 
                                                         output_uncert_all_DD$SP_vals_with_params$site=="Tamakin Dacot", 
                                                         "Tomakin Dako")
best_est_metrics_mean_offspring_DD$SP$site <- replace(best_est_metrics_mean_offspring_DD$SP$site, best_est_metrics_mean_offspring_DD$SP$site=="Tamakin Dacot", "Tomakin Dako")


# SP (accounting for DD)
SP_plot_DD <- ggplot(data = SP_vals_95, aes(x=factor(site, level=SP_vals_95_level_order))) +
  geom_linerange(data = SP_vals_95, aes(ymin=SP_min, ymax=SP_max), color = "grey") +
  geom_linerange(data = SP_vals_95, aes(ymin=SP_lower, ymax=SP_upper), color = "black") +
  geom_point(data = best_est_metrics_mean_offspring_DD$SP, aes(x = site, y = SP_value), color = "black") +
  xlab("\nPatch") + ylab(bquote("Self persistence (SP"[i]~")")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# SP_plot_DD <- ggplot(data = output_uncert_all_DD$SP_vals_with_params, aes(x=reorder(site, org_geo_order), y=SP)) +
#   geom_violin(fill = "grey") +
#   geom_point(data = best_est_metrics_mean_offspring_DD$SP, aes(x = site, y = SP_value), color = "black") +
#   xlab("\nPatch") + ylab(bquote("Self persistence (SP"[i]~")")) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# NP (accounting for DD)
NP_95_lower <- (metrics_summary %>% filter(metric == "NP_DD lower"))$value
NP_95_upper <- (metrics_summary %>% filter(metric == "NP_DD upper"))$value
NP_vals_95 <- output_uncert_all_DD$NP_out_df %>% filter(value >= NP_95_lower & value <= NP_95_upper)

NP_plot_DD_freq <- ggplot(data = output_uncert_all_DD$NP_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color='light gray', fill='light gray') +
  geom_histogram(data = NP_vals_95, aes(x=value, y=..count../sum(..count..)), bins=40, color = "dark gray", fill="dark gray") +
  geom_vline(xintercept = best_est_metrics_mean_offspring_DD$NP, color = "black") +
  xlab(expression(lambda[c])) + ylab("Relative frequency") +
  theme_bw() 
# 
# NP_plot_DD_freq <- ggplot(data = output_uncert_all_DD$NP_out_df, aes(x=value)) +
#   geom_histogram(aes(y=..count../sum(..count..)), bins=40, color='gray', fill='gray') +
#   geom_vline(xintercept = best_est_metrics_mean_offspring_DD$NP, color = "black") +
#   xlab(expression(lambda[c])) + ylab("Relative frequency") +
#   theme_bw() 

# realized connectivity matrix
realized_C_plot_DD <- ggplot(data = best_est_metrics_mean_offspring_DD$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_R)) +
  scale_fill_gradient(high='black', low='white', name='Recruits') +
  xlab('\nOrigin') + ylab('Destination') + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)) +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(angle=45,hjust=1)) 

# All three together - NP as frequency
pdf(file = here::here('Plots/FigureDrafts', 'SP_NP_connMatrixR_freq.pdf'), width=11, height=5)  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(SP_plot_DD, realized_C_plot_DD, NP_plot_DD_freq, rel_widths=c(1.2,1.5,1.2), labels = c("a","b","c"), nrow=1)
dev.off()

##### Figure 6 - alternative geographies and larval navigation 
# A: more habitat
perc_hab_plot <- ggplot(data = perc_hab_plot_df %>% filter(run==1), aes(x=perc_hab, y=value)) +
  geom_line(color="dark gray", alpha=0.2) +
  xlab("Percent habitat") + ylab(expression(lambda[c])) +
  theme_bw()
for(i in 2:n_runs) {
  perc_hab_plot <- perc_hab_plot +
    geom_line(data=perc_hab_plot_df %>% filter(run==i), color="dark gray", alpha=0.2)
}
perc_hab_plot <- perc_hab_plot +
  geom_line(data = perc_hab_best_est_NP_df, aes(x=perc_hab,y=value),color="black") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = param_best_est_mean_collected_offspring$prop_hab, color = "orange") 

# B: wider region, same habitat density
wider_region_plot <- ggplot(data = wider_region_plot_df %>% filter(run==1), aes(x=region_width, y=value)) +
  geom_line(color="dark gray", alpha=0.2) +
  xlab("Region width (km)") + ylab(expression(lambda[c])) +
  theme_bw()
for(i in 2:n_runs) {
  wider_region_plot <- wider_region_plot +
    geom_line(data=wider_region_plot_df %>% filter(run==i), color="dark gray", alpha=0.2)
}
wider_region_plot <- wider_region_plot +
  geom_line(data = wider_region_best_est_NP_df, aes(x=region_width,y=value),color="black") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = figure_vectors_and_values$region_width_km, color = "orange") 

# C: wider region, perc habitat density
# color schemes are from ColorBrewer: https://colorbrewer2.org/#type=sequential&scheme=YlGnBu&n=5
lambda_c_color_scheme_5 <- c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494")  # 5 colors

# 5 categories of NP values
wider_region_perc_hab_plot_5cat <- ggplot(wider_region_perc_hab_plot_df_categories, aes(x=region_width, y=perc_hab, color=reorder(NP_categories,order_5), fill=reorder(NP_categories,order_5))) +
  geom_point(size=4, shape=15) +
  scale_color_manual(values = lambda_c_color_scheme_5) +
  scale_fill_manual(values = lambda_c_color_scheme_5) +
  geom_hline(yintercept = param_best_est_mean_collected_offspring$prop_hab, color = "orange") +
  geom_vline(xintercept = figure_vectors_and_values$region_width_km, color = "orange") +
  theme(legend.title = element_text( size=2), legend.text=element_text(size=2)) +
  theme(legend.key.size = unit(1, "cm")) +
  xlab("Region width (km)") + ylab("Percent habitat") + labs(fill = expression(lambda[c]), color = expression(lambda[c])) +
  theme_bw() 

# Pull out legend from 3D plot to make spacing tighter than default
legend_3D <- get_legend(wider_region_perc_hab_plot_5cat)
wider_region_perc_hab_plot_5cat <- wider_region_perc_hab_plot_5cat + theme(legend.position = "none")
wider_region_perc_hab_plot_5cat_with_legend <- plot_grid(wider_region_perc_hab_plot_5cat, legend_3D, nrow=1, rel_widths = c(0.68, 0.32))

# D: larval navigation
larv_nav_plot <- ggplot(data = larv_nav_plot_df %>% filter(run==1), aes(x=larv_nav, y=value)) +
  geom_line(color="gray", alpha=0.2) +
  xlab("Larval navigation (km)") + ylab(expression(lambda[c])) +
  theme_bw()
for(i in 2:n_runs) {
  larv_nav_plot <- larv_nav_plot +
    geom_line(data=larv_nav_plot_df %>% filter(run==i), color="gray", alpha=0.2)
}
larv_nav_plot <- larv_nav_plot +
  geom_line(data = larv_nav_best_est_NP_df, aes(x=larv_nav,y=value),color="black") +
  geom_hline(yintercept = 1, color = "blue") 

# All together 
pdf(file=here::here("Plots/FigureDrafts","What_if_4_panels_3D.pdf"), width=6, height=6)
plot_grid(perc_hab_plot, wider_region_plot, wider_region_perc_hab_plot_5cat_with_legend, larv_nav_plot, 
          nrow=2, labels=c("a","b","c","d"))
dev.off()

########## Appendix figures ##########

##### D1 - schematic, produced outside of R
##### D2 - schematic, produced outside of R

##### D3 - parameters and their uncertainty not shown in main text Fig (recruit census size, Linf and k for growth, Pc, # assigned offspring)
# Census (start-recruit) size 
startRecruit_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=start_recruit_size)) +
  geom_histogram(bins=40, color="gray", fill="gray") +
  geom_vline(xintercept = param_best_est_mean_collected_offspring$start_recruit_size) +
  xlab("Recruit census size (cm)") + ylab("Count") +
  theme_bw()

# Growth - Linf + k (VBL growth model) 
growthLinf_k_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x=Linf, y=k_growth)) +
  geom_point(color = "gray", fill = "gray") +
  geom_point(x = param_best_est_mean_collected_offspring$Linf, y = param_best_est_mean_collected_offspring$k_growth, color = "black", fill = "black") +
  xlab('Linf (cm)') + ylab("k") + 
  theme_bw()

# Capture probability (Pc)
probR_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x = prob_r)) +
  geom_histogram(bins = 40, color = "gray", fill = "gray") +
  geom_vline(xintercept = param_best_est_mean_collected_offspring$prob_r, color = "black") +
  xlab(bquote("P"[c])) + ylab("Count") +
  theme_bw()

# Assigned offspring (Rm)
assignedOffspring_plot <- ggplot(data = output_uncert_all$metric_vals_with_params, aes(x = assigned_offspring)) +
  geom_histogram(bins = 40, color = "gray", fill = "gray") +
  geom_vline(xintercept = param_best_est_mean_collected_offspring$offspring_assigned_to_parents, color = "black") +
  xlab("# Assigned offspring") + ylab("Count") +
  theme_bw()

# Put them together
pdf(file = here::here('Plots/FigureDrafts', 'APP_FIG_Parameter_inputs.pdf'))  
plot_grid(startRecruit_plot, growthLinf_k_plot, probR_plot, assignedOffspring_plot,
          labels = c("a","b","c","d"), nrow=2)
dev.off()

##### D4 - proportion of dispersal kernel from each patch covered by sampling region
# Correct spelling of name for Tomakin Dako patch
all_parents_site$site <- replace(all_parents_site$site, all_parents_site$site=="Tamakin Dacot", "Tomakin Dako")

pdf(file = here('Plots/FigureDrafts', 'Prop_of_kernel_area_sampled_by_site.pdf'))
ggplot(data = all_parents_site, aes(x = reorder(site, site_geo_order), y = prop_disp_area_within_sites)) + # the geo orders are all off here...
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 1) +
  xlab("Patch") + ylab("Proportion kernel within sampled area (P"[d]~")") + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##### D5 - survival by patch
surv_by_site_to_plot <- best_fit_model_dfs$surv_site_size %>%  # change site from factor to character so can update spelling of Tomakin Dako
  mutate_if(is.factor, as.character) %>%
  mutate(site = if_else(site == "Tamakin Dacot", "Tomakin Dako", site))

pdf(file = here::here("Plots/FigureDrafts", "APP_FIG_surv_by_size_and_site_Phisiteplussize_psizeplusdist.pdf"))
ggplot(data = surv_by_site_to_plot, aes(x = size, y = estimate_prob)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab(expression(paste("Survival probability (",phi[i],")"))) + xlab("Size (cm)") +
  facet_wrap(~site) +
  theme_bw() 
dev.off()

##### D6 - size and distance effects on probability of recapture
# p (by dist)
p_by_dist_plot <- ggplot(data = best_fit_model_dfs$recap_dist, aes(x = dist, y = estimate_prob)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab(expression("Recapture probability (p"[r]~")")) + xlab("Distance (m)") + 
  theme_bw() 

# p (by size)
p_by_size_plot <- ggplot(data = best_fit_model_dfs$recap_size, aes(x = size, y = estimate_prob)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab(expression("Recapture probability (p"[r]~")")) + xlab("Size (cm)") + 
  theme_bw() 

# Put them together
pdf(file = here::here("Plots/FigureDrafts", "APP_FIG_recap_effects_Phisiteplussize_psizeplusdist.pdf"), width=7, height=4)
plot_grid(p_by_dist_plot, p_by_size_plot, labels=c("a","b"), nrow=1)
dev.off()

##### D7 - abundance trends by patch, produced in TimeSeriesPersistence.R

##### D8 - patch-specific LRP
# Join so have site name rather than just number
site_specific_LRP_plot_df <- left_join(output_uncert_all_DD$LEP_R_by_site_out_df, site_vec_order, by = c("site" = "alpha_order"))

pdf(file = here::here("Plots/FigureDrafts","LRP_by_site.pdf"))
ggplot(data = site_specific_LRP_plot_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=50, color='gray', fill='gray') +
  #geom_vline(xintercept = best_est_metrics_mean_offspring_DD$LEP_by_site, aes(x = site, y = LEP_R), color = "black") +
  facet_wrap(~site_name) +
  xlab(bquote("Lifetime recruit production (LRP"[i] ~")")) + ylab("Relative frequency") + 
  theme_bw() 
dev.off()

##### D9 - LRP and LR without DD compensation
# LRP without DD 
LRP_plot_freq <- ggplot(data = output_uncert_all$LEP_R_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color = 'light gray', fill = 'light gray') +
  geom_vline(xintercept = best_est_metrics_mean_offspring$LEP_R_mean, color = "black") +
  xlab(expression("Lifetime recruit production (LRP"[D]~")")) + ylab("Relative frequency") + 
  theme_bw()

# LR without DD
LR_plot_freq <- ggplot(data = output_uncert_all$LEP_R_local_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color = "light gray", fill = "light gray") +
  geom_vline(xintercept = best_est_metrics_mean_offspring$LEP_R_local_mean, color = "black") +
  xlab(expression("Local replacement (LR"[D]~")")) + ylab("Relative frequency") + 
  theme_bw()

# Put them together - frequencies
pdf(file = here::here('Plots/FigureDrafts', 'APP_FIG_LRP_LocalReplacement_withoutDDconsidered_freq.pdf'), width=6, height=3)
plot_grid(LRP_plot_freq, LR_plot_freq, labels = c("a","b"), nrow=1)
dev.off()

##### D10 - SP, NP, and connectivity matrix without density dependence compensation
# Update patch name to the correct spelling (Tomakin Dako)
best_est_metrics_mean_offspring$Cmat$org_site <- replace(best_est_metrics_mean_offspring$Cmat$org_site, 
                                                         best_est_metrics_mean_offspring$Cmat$org_site=="Tamakin Dacot", 
                                                         "Tomakin Dako")
best_est_metrics_mean_offspring$Cmat$dest_site <- replace(best_est_metrics_mean_offspring$Cmat$dest_site, 
                                                          best_est_metrics_mean_offspring$Cmat$dest_site=="Tamakin Dacot", 
                                                          "Tomakin Dako")
output_uncert_all$SP_vals_with_params$site <- replace(output_uncert_all$SP_vals_with_params$site, 
                                                      output_uncert_all$SP_vals_with_params$site=="Tamakin Dacot", 
                                                      "Tomakin Dako")
best_est_metrics_mean_offspring$SP$site <- replace(best_est_metrics_mean_offspring$SP$site, best_est_metrics_mean_offspring$SP$site=="Tamakin Dacot", "Tomakin Dako")

# SP (not accounting for DD)
SP_plot <- ggplot(data = output_uncert_all$SP_vals_with_params, aes(x=reorder(site, org_geo_order), y=SP)) +
  geom_violin(fill="light gray") +
  geom_point(data = best_est_metrics_mean_offspring$SP, aes(x = site, y = SP_value), color = "black") +
  xlab("\nPatch") + ylab(expression("Self persistence (SP"[i[D]]~")")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  

# NP (not accounting for DD)
NP_plot_freq <- ggplot(data = output_uncert_all$NP_out_df, aes(x=value)) +
  geom_histogram(aes(y=..count../sum(..count..)), bins=40, color='light gray', fill='light gray') +
  geom_vline(xintercept = best_est_metrics_mean_offspring$NP, color = "black") +
  xlab(expression(lambda[c[D]])) + ylab("Relative frequency") + 
  theme_bw()

# realized connectivity matrix (not accounting for DD)
realized_C_plot <- ggplot(data = best_est_metrics_mean_offspring$Cmat, aes(x=reorder(org_site, org_geo_order), y=reorder(dest_site, dest_geo_order))) +
  geom_tile(aes(fill=prob_disp_R)) +
  scale_fill_gradient(high='black', low='white', breaks=c(0,0.02,0.04), name=expression('Recruits'[D])) +
  xlab('\nOrigin') + ylab('Destination') + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(legend.position = "bottom")

# Put together (NP as frequency)
pdf(file = here::here('Plots/FigureDrafts', 'APP_FIG_SP_NP_connMatrixR_withoutDDcompensation_freq.pdf'), width=11, height=5)  # hacked the color scales being comparable, deal with for real if people like showing both
plot_grid(SP_plot, realized_C_plot, NP_plot_freq, rel_widths=c(1.2,1.5,1.2), labels = c("a","b","c"), nrow=1)
dev.off()

##### D11 - uncertainty in LEP (averaged across patches)
pdf(file = here::here("Plots/FigureDrafts", "LEP_uncertainty_by_param.pdf"), width=6, height=4)
ggplot(data = LEP_uncert_DD %>% filter(uncertainty_type %in% c("start recruit size", "breeding size", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("start recruit size", "breeding size", "growth", "survival", "all"),
                               value = mean(best_est_metrics_mean_offspring_DD$LEP_by_site$LEP), color = "black")) +
  xlab("Uncertainty type") + ylab(expression("Lifetime egg production (LEP"["*"]~")")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  theme(text=element_text(size=12))
dev.off()

##### D12 - uncertainty in LRP
LEP_R_uncert_DD_plot <- LEP_R_uncert_DD
LEP_R_uncert_DD_plot$uncertainty_type <- replace(LEP_R_uncert_DD_plot$uncertainty_type, 
                                                 LEP_R_uncert_DD_plot$uncertainty_type=="prob r", 
                                                 "capture probability")

pdf(file = here::here("Plots/FigureDrafts", "LRP_uncertainty_by_param.pdf"), width=6, height=4)
ggplot(data = LEP_R_uncert_DD_plot %>% filter(uncertainty_type %in% c("start recruit size", "breeding size", "assigned offspring", "capture probability", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("start recruit size", "breeding size", "assigned offspring", "capture probability", "growth", "survival", "all"),
                               value = best_est_metrics_mean_offspring_DD$LEP_R_mean, color = "black")) +
  xlab("Uncertainty type") + ylab("Lifetime recruit production (LRP)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(text=element_text(size=12))
dev.off()

##### D13 - uncertainty in egg-recruit survival
RperE_uncert_DD_plot <- RperE_uncert_DD
RperE_uncert_DD_plot$uncertainty_type <- replace(RperE_uncert_DD_plot$uncertainty_type, 
                                                 RperE_uncert_DD_plot$uncertainty_type=="prob r", 
                                                 "capture probability")
pdf(file = here::here("Plots/FigureDrafts", "RperE_uncertainty_by_param.pdf"), width=6, height=4)
ggplot(data = RperE_uncert_DD_plot %>% filter(uncertainty_type %in% c("breeding size", "assigned offspring", "capture probability", "growth", "survival", "all")), aes(x=uncertainty_type, y=value)) +
  geom_violin(fill="grey") +
  geom_point(data = data.frame(uncertainty_type = c("breeding size", "assigned offspring", "capture probability", "growth", "survival", "all"),
                               value = best_est_metrics_mean_offspring_DD$recruits_per_egg), color = "black") +
  xlab("Uncertainty type") + ylab(expression("Recruits per egg (S"[e]~")")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(text=element_text(size=12))
dev.off()

##### D14 - uncertainty in NP
# Rename params to make labels easier to understand
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
                               value = best_est_metrics_mean_offspring_DD$NP), color = "black") +
  xlab("Uncertainty type") + ylab(expression(lambda[c])) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(text = element_text(size=12))
dev.off()
