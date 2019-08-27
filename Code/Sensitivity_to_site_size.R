# For "what-ifs" where there are 19 sites, evenly spaced, of increasing proportion of area covered by habitat
# Create different habitat set-ups

# Will incoporate this into other scripts, probably Site_widths_and_distances.R, later

# Length of region (total_range_of_sampling_area is from Metrics_with_uncertainty)
region_width_km <- total_range_of_sampling_area/1000

# Proportion of habitat of actual sites - prop_sampling_area_habitat (from Metrics_with_uncertainty), say it's 20%


# Range of percent habitats:
perc_hab_vals <- seq(from=0.2,to=0.9, by=0.05)

# Find width of sites and distance between them
find_site_locations <- function(nsites, perc_hab_val, total_region) {
  
  site_info <- data.frame(org_site = seq(from=1, to=nsites, by=1))
  site_width <- (total_region*perc_hab_val)/nsites
  hab_break <- (total_region*(1-perc_hab_val))/(nsites-1)
  
  # if((site_width*nsites + hab_break*(nsites-1)) != total_region){
  #   print("Error with sum of sites and non habitat regions!")
  # }

  for(i in 1:nsites){
   site_info$S_edge_org[i] = (i-1)*(site_width + hab_break)
   site_info$N_edge_org[i] = (i-1)*(site_width + hab_break) + site_width
   site_info$center_org[i] = (i-1)*(site_width + hab_break) + (site_width/2)
   site_info$width_org[i] = site_info$N_edge_org[i] - site_info$S_edge_org[i]
  }
  
  return(site_info)
}

# Turn into a bigger data frame 
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
                             org_site == dest_site ~ width_dest))  # max distance is site width for self-self distances
  
  # Add the random other columns the calcMetrics function needs
  out_df <- out_df %>% 
    mutate(org_alpha_order = org_site,
           org_geo_order = org_site,
           dest_alpha_order = dest_site,
           dest_geo_order = dest_site)
  
 return(out_df)
  
}

# Make different sets across perc_hab_vals

###### Do runs 
# param_best_est_mean_collected_offspring
# param_set_full

site_vec_order_Sens <- data.frame(site_name = seq(from=1,to=nsites,by=1),
                                  alpha_order = seq(from=1,to=nsites,by=1),
                                  geo_order = seq(from=1,to=nsites,by=1))
# Site info
site_dist_info_0.05 <- make_output_with_dist(19, 0.05, region_width_km)
site_dist_info_0.1 <- make_output_with_dist(19, 0.1, region_width_km)
site_dist_info_0.15 <- make_output_with_dist(19, 0.15, region_width_km)
site_dist_info_0.2 <- make_output_with_dist(19, 0.2, region_width_km)
site_dist_info_0.25 <-  make_output_with_dist(19, 0.25, region_width_km)
site_dist_info_0.3 <-  make_output_with_dist(19, 0.3, region_width_km)
site_dist_info_0.35 <-  make_output_with_dist(19, 0.35, region_width_km)
site_dist_info_0.4 <-  make_output_with_dist(19, 0.4, region_width_km)
site_dist_info_0.45 <-  make_output_with_dist(19, 0.45, region_width_km)
site_dist_info_0.5 <-  make_output_with_dist(19, 0.5, region_width_km)
site_dist_info_0.55 <-  make_output_with_dist(19, 0.55, region_width_km)
site_dist_info_0.6 <-  make_output_with_dist(19, 0.6, region_width_km)
site_dist_info_0.65 <-  make_output_with_dist(19, 0.65, region_width_km)
site_dist_info_0.7 <-  make_output_with_dist(19, 0.7, region_width_km)
site_dist_info_0.75 <-  make_output_with_dist(19, 0.75, region_width_km)
site_dist_info_0.8 <-  make_output_with_dist(19, 0.8, region_width_km)
site_dist_info_0.85 <-  make_output_with_dist(19, 0.85, region_width_km)
site_dist_info_0.9 <-  make_output_with_dist(19, 0.9, region_width_km)
site_dist_info_0.95 <-  make_output_with_dist(19, 0.95, region_width_km)

# Best ests
perc_hab_0.05 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.05, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.1 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.1, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.15 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.15, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.2 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.2, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.25 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.25, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.3 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.3, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.35 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.35, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.4 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.4, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.45 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.45, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.5 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.5, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.55 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.55, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.6 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.6, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.65 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.65, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.7 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.7, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.75 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.75, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.8 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.8, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.85 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.85, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.9 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.9, site_vec_order_Sens$site_name, "TRUE")
perc_hab_0.95 <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info_0.95, site_vec_order_Sens$site_name, "TRUE")

# Do uncertainty runs
perc_hab_0.2_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.2, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.3_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.3, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.4_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.4, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.5_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.5, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.6_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.6, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.7_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.7, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.8_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.8, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.9_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.9, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.25_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.0.25, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.35_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.35, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.45_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.45, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.55_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.55, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.65_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.65, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.75_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.75, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.85_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.85, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.05_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.05, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.1_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.1, site_vec_order_Sens, "all", "TRUE")
perc_hab_0.15_output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info_0.15, site_vec_order_Sens, "all", "TRUE")



# Make data frame to plot
NP_by_perc_hab <- data.frame(perc_hab = c(0.05, 0.1, 0.15, perc_hab_vals),
                             NP = c(perc_hab_0.05$NP, perc_hab_0.1$NP, perc_hab_0.15$NP,
                                    perc_hab_0.2$NP, perc_hab_0.25$NP, perc_hab_0.3$NP,
                                    perc_hab_0.35$NP, perc_hab_0.4$NP, perc_hab_0.45$NP,
                                    perc_hab_0.5$NP, perc_hab_0.55$NP, perc_hab_0.6$NP,
                                    perc_hab_0.65$NP, perc_hab_0.7$NP, perc_hab_0.75$NP,
                                    perc_hab_0.8$NP, perc_hab_0.85$NP, perc_hab_0.9$NP),
                             NP_sd = c(sd(perc_hab_0.05_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.1_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.15_output_uncert_all$NP_out_df$value),
                                        sd(perc_hab_0.2_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.25_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.3_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.35_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.4_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.45_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.5_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.55_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.6_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.65_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.7_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.75_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.8_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.85_output_uncert_all$NP_out_df$value),
                                       sd(perc_hab_0.9_output_uncert_all$NP_out_df$value)))
NP_by_perc_hab <- NP_by_perc_hab %>% mutate(NP_min = NP-NP_sd, NP_max = NP+NP_sd)

####### Make a plot
ggplot(data = NP_by_perc_hab, aes(x = perc_hab, y = NP)) +
  geom_line() 

pdf(file = here('Plots/Poster_presentation_plots','NP_by_perc_hab.pdf'), width=11, height=11)
ggplot(data=NP_by_perc_hab, aes(x=perc_hab, y=NP, ymin=NP_min, ymax=NP_max)) +
  geom_line(color='black') +
  geom_ribbon(alpha=0.5, color='gray', fill="gray") +
  geom_hline(yintercept = 1, color = "blue") +
  geom_vline(xintercept = 0.2, color = "orange") +
  xlab('proportion habitat') + ylab(expression(lambda)) + ggtitle("Habitat for persistent metapopulation") +
  theme_bw() +
  theme(text = element_text(size=30))
dev.off()

make_space_dfs_Acr
output_uncert_all <- calcMetricsAcrossRuns(n_runs, param_set_full, site_dist_info, site_vec_order, "all", FALSE)
best_est_metrics_mean_offspring_DD <- calcMetrics(param_best_est_mean_collected_offspring, site_dist_info, site_vec_order$site_name, "TRUE")
NP_best_est_DD <- best_est_metrics_mean_offspring_DD$NP

