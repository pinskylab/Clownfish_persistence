# Find the width of each site and distances between them (for probability of dispersal)

#################### Set-up: ####################
# Load packages
library(ggplot2)
library(geosphere)

# Pull data, functions, constants 
source(here::here('Code', 'Constants_database_common_functions.R'))

#################### Functions: ####################

#################### Running things: ####################

##### Find lat/lon coordinates of anems at edges of sites
# Make a data frame for sites, boundary anems, lat lons of those anems, widths (should I just make this the same thing and put it in Constants_database_common_functions?)
site_edges_info <- site_edge_anems %>% filter(anem_loc %in% c("north", "south"))

# Pull out info on those anems so can use anem_table_ids in function to get lat/lon
site_width_anems <- anems_Processed %>%
  filter(anem_id %in% site_edges_info$anem_id) %>%
  distinct(anem_id, .keep_all = TRUE)

# Match those anem_table_ids to the anem_ids in site_width_info
site_edges_info <- left_join(site_edges_info, site_width_anems %>% select(anem_table_id, anem_id), by = "anem_id")

# Find lat/lon for each of those anems
site_edges_info <- site_edges_info %>%
  mutate(lat = NA, lon = NA)

for(i in 1:length(site_edges_info$anem_table_id)) {
  if (is.na(site_edges_info$anem_table_id[i])) {  # some sites don't have a mid-anem and anemid_latlong function doesn't work if you put NA for anem_table_id
    site_edges_info$lat[i] = NA
    site_edges_info$lon[i] = NA
  } else {
    out_anemLL <- anemid_latlong(site_edges_info$anem_table_id[i], anems_Processed, gps_Info)  
    site_edges_info$lat[i] = out_anemLL$lat
    site_edges_info$lon[i] = out_anemLL$lon
  }
}

# Remove temporary values from for loop for neatness
rm(out_anemLL)

##### Find northern-most edge and southern-most edge of sampling area
northern_edge_lat <- (site_edges_info %>% filter(site == "Palanas", anem_loc == "north"))$lat
northern_edge_lon <- (site_edges_info %>% filter(site == "Palanas", anem_loc == "north"))$lon
southern_edge_lat <- (site_edges_info %>% filter(site == "Sitio Baybayon", anem_loc == "south"))$lat
southern_edge_lon <- (site_edges_info %>% filter(site == "Sitio Baybayon", anem_loc == "south"))$lon

# Put in a data frame so can save
sampling_area_edges <- data.frame(edge = c("north", "south"),
                                  lat = c(northern_edge_lat, southern_edge_lat),
                                  lon = c(northern_edge_lon, southern_edge_lon))

##### Find the width of each site and distance to edges of sampling area
# Calculate the distance between the edge anems for each site
site_width_info <- site_edges_info %>%
  filter(anem_loc == "north") %>%
  select(site, site_geo_order) %>%
  mutate(width_m = NA,  # width of site (distance between N and S-most anems)
         dist_to_N_edge_m = NA,  # distance from middle of site to N-most edge of sampling area
         dist_to_S_edge_m = NA)  # distance from middle of site to S-most edge of sampling area

# UMMM, THESE ARE SUPER DIFFERENT THAN THE WIDTHS I HAD BEFORE - WHY??
for(i in 1:length(site_width_info$site)) {
  site_width_info$width_m[i] = distHaversine(c((site_edges_info %>% filter(anem_loc == "north", site == site_width_info$site[i]))$lon, (site_edges_info %>% filter(anem_loc == "north", site == site_width_info$site[i]))$lat), 
                                             c((site_edges_info %>% filter(anem_loc == "south", site == site_width_info$site[i]))$lon, (site_edges_info %>% filter(anem_loc == "south", site == site_width_info$site[i]))$lat))  
  site_width_info$dist_to_N_edge_m[i] = distHaversine(c((site_centers %>% filter(site == site_width_info$site[i]))$lon, (site_centers %>% filter(site == site_width_info$site[i]))$lat),
                                                      c(northern_edge_lon, northern_edge_lat))
  site_width_info$dist_to_S_edge_m[i] = distHaversine(c((site_centers %>% filter(site == site_width_info$site[i]))$lon, (site_centers %>% filter(site == site_width_info$site[i]))$lat),
                                                      c(southern_edge_lon, southern_edge_lat))
}

# Convert to km
site_width_info <- site_width_info %>%
  mutate(width_km = width_m/1000,
         dist_to_N_edge_km = dist_to_N_edge_m/1000,
         dist_to_S_edge_km = dist_to_S_edge_m/1000)


# for(i in 1:length((site_width_info %>% filter(anem_loc == "north"))$site)) {
#   site_width_info$width_m[i] = distHaversine(c((site_width_info %>% filter(anem_loc== "north"))$lon[i], (site_width_info %>% filter(anem_loc== "north"))$lat[i]),
#                                              c((site_width_info %>% filter(anem_loc== "south"))$lon[i], (site_width_info %>% filter(anem_loc== "north"))$lat[i]))
#   #site_width_info$width_m[i] = distHaversine(c(site_width_info$N_lon[i], site_width_info$N_lat[i]), c(site_width_info$S_lon[i], site_width_info$S_lat[i]))
#   site_width_info$width_m[i+length((site_width_info %>% filter(anem_loc == "north"))$site)] = site_width_info$width_m[i]   # fill in the width for the site in the south anem slot too
#   #site_width_info$width_m[i+2*length((site_width_info %>% filter(anem_loc == "north"))$site)] = site_width_info$width_m[i]   # fill in the width for the site in the south anem slot too
# }
# 
# site_width_info <- site_width_info %>% 
#   mutate(width_km = width_m/1000)


# ##### Make a data frame with distances from middle of each site to edge of sampling area
# # Set up data frame
# site_distances_to_edge <- site_vec_order %>%
#   mutate(dist_to_N_edge_m = NA, dist_to_N_edge_km = NA,
#          dist_to_S_edge_m = NA, dist_to_S_edge_km = NA)
#   
# # Go through and find distance
# for(i in 1:length(site_distances_to_edge$site_name)) {
#   site_distances_to_edge$dist_to_N_edge_m[i] = distHaversine(c((site_centers %>% filter(site == site_distances_to_edge$site_name[i]))$lon, (site_centers %>% filter(site == site_distances_to_edge$site_name[i]))$lat), c(northern_edge_lon, northern_edge_lat))
#   site_distances_to_edge$dist_to_S_edge_m[i] = distHaversine(c((site_centers %>% filter(site == site_distances_to_edge$site_name[i]))$lon, (site_centers %>% filter(site == site_distances_to_edge$site_name[i]))$lat), c(southern_edge_lon, southern_edge_lat))
#   site_distances_to_edge$dist_to_N_edge_km[i] = site_distances_to_edge$dist_to_N_edge_m[i]/1000
#   site_distances_to_edge$dist_to_S_edge_km[i] = site_distances_to_edge$dist_to_S_edge_m[i]/1000
# }

##### Find distances between sites (RIGHT NOW DOING N-N, SHOULD DO MID-MID but not all sites have a mid anem or one that has both lat/lon coordinates, at least for the observation of the anem I chose)
# Set up data frame 
org_site_list <- rep(site_vec_order$site_name[1], length(site_vec_order$site_name))  # need a list of each site repeated 19 times so can get all sites as origins and destinations for each other site

for (i in 2:length(site_vec_order$site_name)){
  org_site_list <- c(org_site_list, rep(site_vec_order$site_name[i], length(site_vec_order$site_name)))
}

site_dist_info <- data.frame(org_site = org_site_list, stringsAsFactors = FALSE) %>%
  mutate(dest_site = rep(site_vec_order$site_name, length(site_vec_order$site_name))) %>%
  mutate(dist_mid_to_S_m = NA, dist_mid_to_S_km = NA,  # distance from midpoint of origin to south side of destination in m and km
         dist_mid_to_N_m = NA, dist_mid_to_N_km = NA,  # distance from midpoint of origin to north side of destination in m and km
         dist_mid_to_mid_m = NA, dist_mid_to_mid_km = NA,  # distance from midpoint of origin to midpoint of destination in m and km
         d1_km = NA, d2_km = NA,  # d1: distance to close edge of the destination site, d2: distance to far edge of destination site
         dest_width = NA)   # width of destination site
         #dist_avg = NA, dest_width = NA,  # average distance (like mid-mid), width of destination site
         #d1a_km = NA, d2a_km = NA)  # d1 is average distance - half the width of the dest site, d2 is avg dist + half the width of the dest site (so like mid of org to close end of dest, mid of org to far end of dest)

# site_dist_info <- data.frame(org_site = org_site_list, stringsAsFactors = FALSE) %>%
#   mutate(dest_site = rep(site_vec_order$site_name, length(site_vec_order$site_name))) %>%
#   mutate(dist_N_to_S_m = NA, dist_N_to_S_km = NA,  # distance from north side of origin to south side of destination in m and km
#          dist_S_to_N_m = NA, dist_S_to_N_km = NA,  # distance from south side of origin to north side of destination in m and km
#          dist_avg = NA, dest_width = NA,  # average distance (like mid-mid), width of destination site
#          d1a_km = NA, d2a_km = NA)  # d1 is average distance - half the width of the dest site, d2 is avg dist + half the width of the dest site (so like mid of org to close end of dest, mid of org to far end of dest)

## SOMETHING IS WRONG HERE - THE TWO WAYS OF DOING DISTANCE DON'T MATCH UP!!
# Find the distances between the sites! (for now, doing N to N, should do mid to mid...)
for(i in 1:length(site_dist_info$org_site)) {
  site_org = site_dist_info$org_site[i]
  N_lat_org = (site_edges_info %>% filter(site == site_org, anem_loc == "north"))$lat  
  N_lon_org = (site_edges_info %>% filter(site == site_org, anem_loc == "north"))$lon
  S_lat_org = (site_edges_info %>% filter(site == site_org, anem_loc == "south"))$lat 
  S_lon_org = (site_edges_info %>% filter(site == site_org, anem_loc == "south"))$lon
  mid_lat_org = (site_centers %>% filter(site == site_org))$lat
  mid_lon_org = (site_centers %>% filter(site == site_org))$lon
  site_dest = site_dist_info$dest_site[i]
  N_lat_dest = (site_edges_info %>% filter(site == site_dest, anem_loc == "north"))$lat
  N_lon_dest = (site_edges_info %>% filter(site == site_dest, anem_loc == "north"))$lon
  S_lat_dest = (site_edges_info %>% filter(site == site_dest, anem_loc == "south"))$lat
  S_lon_dest = (site_edges_info %>% filter(site == site_dest, anem_loc == "south"))$lon
  mid_lat_dest = (site_centers %>% filter(site == site_dest))$lat
  mid_lon_dest = (site_centers %>% filter(site == site_dest))$lon
  
  site_dist_info$dist_mid_to_S_m[i] = distHaversine(c(mid_lon_org, mid_lat_org), c(S_lon_dest, S_lat_dest))
  site_dist_info$dist_mid_to_S_km[i] = site_dist_info$dist_mid_to_S_m[i]/1000 
  site_dist_info$dist_mid_to_N_m[i] = distHaversine(c(mid_lon_org, mid_lat_org), c(N_lon_dest, N_lat_dest))
  site_dist_info$dist_mid_to_N_km[i] = site_dist_info$dist_mid_to_N_m[i]/1000 
  site_dist_info$dist_mid_to_mid_m[i] = distHaversine(c(mid_lon_org, mid_lat_org), c(mid_lon_dest, mid_lat_dest))
  site_dist_info$dist_mid_to_mid_km[i] = site_dist_info$dist_mid_to_mid_m[i]/1000
  
  #site_dist_info$dist_avg[i] = (site_dist_info$dist_mid_to_S_km[i] + site_dist_info$dist_mid_to_N_km[i])/2
  site_dist_info$dest_width[i] = distHaversine(c(N_lon_dest, N_lat_dest), c(S_lon_dest, S_lat_dest))/1000
  #site_dist_info$dest_width2[i] <- (site_width_info %>% filter(site == site_dest))$width_km
  
  #site_dist_info$d1a_km[i] <- site_dist_info$dist_avg[i] - site_dist_info$dest_width[i]/2
  #site_dist_info$d12_km[i] <- site_dist_info$dist_mid_to_mid_km[i] - site_dist_info$dest_width[i]/2
  #site_dist_info$d2a_km[i] <- site_dist_info$dist_avg[i] + site_dist_info$dest_width[i]/2
  site_dist_info$d1_km[i] <- min(c(site_dist_info$dist_mid_to_S_km[i], site_dist_info$dist_mid_to_N_km[i]))  # d1 is smaller distance to edge of site (since will depend on N-S order of sites)
  site_dist_info$d2_km[i] <- max(c(site_dist_info$dist_mid_to_S_km[i], site_dist_info$dist_mid_to_N_km[i]))
}

# # Find the distances between the sites! (for now, doing N to N, should do mid to mid...)
# for(i in 1:length(site_dist_info$org_site)) {
#   site_org = site_dist_info$org_site[i]
#   N_lat_org = (site_width_info %>% filter(site == site_org, anem_loc == "north"))$lat  
#   N_lon_org = (site_width_info %>% filter(site == site_org, anem_loc == "north"))$lon
#   S_lat_org = (site_width_info %>% filter(site == site_org, anem_loc == "south"))$lat 
#   S_lon_org = (site_width_info %>% filter(site == site_org, anem_loc == "south"))$lon
#   site_dest = site_dist_info$dest_site[i]
#   N_lat_dest = (site_width_info %>% filter(site == site_dest, anem_loc == "north"))$lat
#   N_lon_dest = (site_width_info %>% filter(site == site_dest, anem_loc == "north"))$lon
#   S_lat_dest = (site_width_info %>% filter(site == site_dest, anem_loc == "south"))$lat
#   S_lon_dest = (site_width_info %>% filter(site == site_dest, anem_loc == "south"))$lon
#   site_dist_info$dist_N_to_S_m[i] = distHaversine(c(N_lon_org, N_lat_org), c(S_lon_dest, S_lat_dest))
#   site_dist_info$dist_N_to_S_km[i] = site_dist_info$dist_N_to_S_m[i]/1000 
#   site_dist_info$dist_S_to_N_m[i] = distHaversine(c(N_lon_dest, N_lat_dest), c(S_lon_org, S_lat_org))
#   site_dist_info$dist_S_to_N_km[i] = site_dist_info$dist_S_to_N_m[i]/1000 
#   
#   site_dist_info$dist_avg[i] = (site_dist_info$dist_N_to_S_km[i] + site_dist_info$dist_S_to_N_km[i])/2
#   site_dist_info$dest_width[i] <- (site_width_info %>% filter(site == site_dest, anem_loc == "north"))$width_km
#   
#   site_dist_info$d1a_km[i] <- site_dist_info$dist_avg[i] - site_dist_info$dest_width[i]/2
#   site_dist_info$d2a_km[i] <- site_dist_info$dist_avg[i] + site_dist_info$dest_width[i]/2
# }

# Remove temporary values from for loop for neatness
rm(site_org, N_lat_org, N_lon_org, S_lat_org, S_lon_org, site_dest, N_lat_dest, N_lon_dest, S_lat_dest, S_lon_dest)

# For self-self distance, use 0 as min and site width as max (not the only way I could do this)
site_dist_info <- site_dist_info %>%
  mutate(d1_km = case_when(org_site != dest_site ~ d1_km,
                           org_site == dest_site ~ 0),  # min distance is 0 for self-self distances
         d2_km = case_when(org_site != dest_site ~ d2_km,
                           org_site == dest_site ~ dest_width))  # max distance is site width for self-self distances
  
# # How to do self-self distance? Sum of mid to each edge? (more like what doing above but would prob need to renormalize disp kernel?) Or edge to edge?
# # Correct d1 and d2 for selfs? The difference between them is the width of the site but should I be doing integral of 0 to close edge + 0 to far edge or (what I've been doing) integral of 0 to width of site
# # Making another column here so have option to do either with distances (d1a_km and d2a_km is mid to edges for all, d1_km and d2_km are mid to edge for between sites, 0 to edge for site to self)
# site_dist_info <- site_dist_info %>%
#   mutate(d1_km = case_when(org_site != dest_site ~ d1a_km,  # same for distances between sites
#                           org_site == dest_site ~ 0),  # 0 for site to self
#          d2_km = case_when(org_site != dest_site ~ d2a_km,  # same for distances between sites
#                           org_site == dest_site ~ dest_width))  # width of site when going self-self


# Add in origin site alpha and geo order
site_dist_info <- left_join(site_dist_info, site_vec_order, by = c('org_site' = 'site_name'))
site_dist_info <- site_dist_info %>%
  dplyr::rename(org_alpha_order = alpha_order, org_geo_order = geo_order)

# Add in destination site alpha and geo order
site_dist_info <- left_join(site_dist_info, site_vec_order, by = c('dest_site' = 'site_name'))
site_dist_info <- site_dist_info %>%
  dplyr::rename(dest_alpha_order = alpha_order, dest_geo_order = geo_order)


# 
# # Find the number of parents at each site (eventually, this parent file pull will go in Constants_database_common_functions). Just putting it here for now b/c going to use original 913, with same distribution as current parents
# all_parents_site <- all_parents %>%
#   group_by(site) %>%
#   summarize(nparents = n()) %>%
#   mutate(prop_parents = nparents/sum(nparents)) %>%
#   mutate(nparents_olddata = round(prop_parents*n_parents_parentage))
# # new number of parents (b/c rounding...)
# n_parents_parentage <- sum(all_parents_site$nparents_olddata)
# n_parents_somesites <- sum((all_parents_site %>% filter(site %in% sites_for_total_areas))$nparents_olddata)   # Just for some sites...
# 
# # Total dispersal kernel area (total parents*2 - total area dispersing north of site is 1 and south is 1 for each parent)
# total_parent_kernel_area = n_parents_somesites*2
# 
# 
# # Join up anem info with parents
# all_parents_latlon <- left_join(all_parents_site, site_width_info %>% filter(anem_loc == "mid") %>% select(site, site_geo_order, lat, lon), by = "site")
# 
# # Dispersal kernel with best-fit params as a function of d - think about where to put this now that egg-recruit survival will depend on dispersal kernel
# disp_allyears_d <- function(d) {  # theta = 0.5, equation for p(d) in eqn. 6c in Bode et al. 2018
#   z = exp(k_allyears)
#   disp = (z/2)*exp(-(z*d)^(theta_allyears))
#   return(disp)
# }
# 
# # Find distance to edges (should move this into Site widths and distances script)
# # For now, proceeding just with sites with all the coords - figure out a better way to get coords for center of other sites...
# all_parents_latlon <- all_parents_latlon %>% 
#   filter(!is.na(lon) & !is.na(lat)) %>%
#   mutate(dist_to_N_edge_m = NA,
#          dist_to_S_edge_m = NA,
#          dist_to_N_edge_km = NA,
#          dist_to_S_edge_km = NA,
#          disp_area_N_within_sites = NA,
#          disp_area_S_within_sites = NA)
# 
# for(i in 1:length(all_parents_latlon$site)) {
#   all_parents_latlon$dist_to_N_edge_m[i] = distHaversine(c(all_parents_latlon$lon[i], all_parents_latlon$lat[i]),
#                                                          c(northern_edge_lon, northern_edge_lat))
#   all_parents_latlon$dist_to_S_edge_m[i] = distHaversine(c(all_parents_latlon$lon[i], all_parents_latlon$lat[i]),
#                                                          c(southern_edge_lon, southern_edge_lat))
#   all_parents_latlon$dist_to_N_edge_km[i] = all_parents_latlon$dist_to_N_edge_m[i]/1000
#   all_parents_latlon$dist_to_S_edge_km[i] = all_parents_latlon$dist_to_S_edge_m[i]/1000
#   all_parents_latlon$disp_area_N_within_sites[i] = integrate(disp_allyears_d, 0, all_parents_latlon$dist_to_N_edge_km[i])$value
#   all_parents_latlon$disp_area_S_within_sites[i] = integrate(disp_allyears_d, 0, all_parents_latlon$dist_to_S_edge_km[i])$value
# }
# 
# # Find proportion of total area under dispersal kernel (where total area to INF is 2 - 1 for each side) covered within sample sites
# all_parents_latlon <- all_parents_latlon %>%
#   mutate(total_disp_area_within_sites = disp_area_N_within_sites + disp_area_S_within_sites,
#          prop_disp_area_within_sites = total_disp_area_within_sites/2,
#          total_parent_area_sampled = total_disp_area_within_sites*nparents)
# all_parents_latlon_summarized <- all_parents_latlon %>%
#   summarize(total_parent_kernel_area = sum(nparents)*2,
#             sampled_parent_kernel_area = sum(total_parent_area_sampled),
#             prop_parent_kernel_area_sampled = sampled_parent_kernel_area/total_parent_kernel_area)
# 
# ggplot(data = all_parents_latlon, aes(x = reorder(site, site_geo_order), y = prop_disp_area_within_sites)) + # the geo orders are all off here...
#   geom_bar(position = "dodge", stat = "identity") +
#   geom_hline(yintercept = 2) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#################### Plots: ####################

#################### Saving output: ####################
save(site_width_info, file=here::here("Data/Script_outputs", "site_width_info.RData"))  # width of sites, distance to edges of sampling area
save(site_dist_info, file=here::here("Data/Script_outputs", "site_dist_info.RData"))  # distances between pairs of sites for kernel integration
save(sampling_area_edges, file=here::here("Data/Script_outputs", "sampling_area_edges.RData"))  # coordinates of northern-most anem and southern-most anem of sampling area



# ################ Code from where these calcs were originally done in PersistenceMetrics.F
# # pull out anem_table_ids for those anems so can use anemid_latlong2 function to find lat lon for boundary anems
# site_width_anems <- leyte %>%
#   tbl("anemones") %>%
#   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id, anem_obs_time) %>%
#   collect() %>%
#   filter(anem_id %in% c(site_width_info$N_anem, site_width_info$S_anem, site_width_info$M_anem)) 
# 
# # pull out the dive info associated with those anems
# site_width_dives <- leyte %>%
#   tbl("diveinfo") %>%
#   select(dive_table_id, date, site, gps) %>%
#   collect() %>%
#   filter(dive_table_id %in% site_width_anems$dive_table_id)
# 
# # join into one data frame to put into MS function (anem_latlong), get the anem obs_time in the same time zone (UTC) as the GPX data, create placeholder lat and lon columns
# site_width_anemdives <- left_join(site_width_anems, site_width_dives, by="dive_table_id") %>%
#   mutate(obs_time = force_tz(ymd_hms(str_c(date, anem_obs_time, sep = " ")), tzone = "Asia/Manila")) %>% #tell it that it is currently in Asia/Manila time zone
#   mutate(obs_time = with_tz(obs_time, tzone = "UTC")) %>% #convert to UTC so can compare with GPS data (this line and one above largely from Michelle's assign_db_gpx function)
#   mutate(month = month(obs_time), #and separate out useful components of the time (this also from Michelle's assign_db_gpx function)
#          day = day(obs_time), 
#          hour = hour(obs_time), 
#          min = minute(obs_time), 
#          sec = second(obs_time), 
#          year = year(obs_time)) %>%
#   mutate(lat = as.numeric(NA), lon = as.numeric(NA)) #put in a placeholder column for lat and lon
# 
# # find lat lon for the boundary anems -- THIS DOESN'T WORK FOR THE MID ANEMS
# for(i in 1:length(site_width_anemdives$anem_table_id)) {
#   out_anemLL <- anemid_latlong_2(site_width_anemdives$anem_table_id[i], site_width_anemdives, gps.Info)  # SHOULD COMBINE/DECIDE BETWEEN THE TWO ANEM_LATLON FUNCTIONS AT SOME POINT
#   #out_anemLL <- anemid_latlong(site_width_anemdives$anem_table_id, site_width_anemdives, gps.Info)
#   site_width_anemdives$lat[i] = out_anemLL$lat[1] #figure out why this is sometimes the wrong replacement length... (put [1] in to try to solve that issue but not sure why it's needed)
#   site_width_anemdives$lon[i] = out_anemLL$lon[1]
# }
# 
# # just pull out one row for each anem_id 
# site_width_anemdives_short <- site_width_anemdives %>%
#   distinct(anem_id, .keep_all = TRUE)
# 
# # looks like anem_id is a chr now? switch to numeric
# site_width_anemdives_short <- site_width_anemdives_short %>%
#   mutate(anem_id = as.numeric(anem_id))
# 
# # put lat and lons into info table (requires some column renaming so needs to be run in order)
# # lat/lon for N_anem
# site_width_info <- left_join(site_width_info %>% dplyr::rename(anem_id = N_anem), site_width_anemdives_short %>% select(anem_id, lat, lon), by="anem_id") #join in the lat lons for the N anem
# site_width_info <- site_width_info %>% dplyr::rename(N_anem = anem_id, N_lat = lat, N_lon = lon) #rename the anem_id, lat, lon columns to be N-specific
# # lat/lon for S_anem
# site_width_info <- left_join(site_width_info %>% dplyr::rename(anem_id = S_anem), site_width_anemdives_short %>% select(anem_id, lat, lon), by="anem_id") #join in the lat lons for the S anem
# site_width_info <- site_width_info %>% dplyr::rename(S_anem = anem_id, S_lat = lat, S_lon = lon) #rename the anem_id, lat, lon columns to be S-specific 
# # lat/lon for M_anem
# site_width_info <- left_join(site_width_info %>% dplyr::rename(anem_id = M_anem), site_width_anemdives_short %>% select(anem_id, lat, lon), by="anem_id") #join in the lat lons for the mid anem
# site_width_info <- site_width_info %>% dplyr::rename(M_anem = anem_id, M_lat = lat, M_lon = lon) #rename the anem_id, lat, lon columns to be M-specific 
# 
# #and calculate the distance!
# for(i in 1:length(site_width_info$site)) {
#   site_width_info$width_m[i] = distHaversine(c(site_width_info$N_lon[i], site_width_info$N_lat[i]), c(site_width_info$S_lon[i], site_width_info$S_lat[i]))
#   site_width_info$width_km[i] = site_width_info$width_m[i]/1000 
# }
# 
# # Find distances among sites
# # Set up data frame
# site_dist_info <- data.frame(org_site = c(rep(site_vec[1],length(site_vec)), rep(site_vec[2],length(site_vec)),
#                                           rep(site_vec[3],length(site_vec)), rep(site_vec[4],length(site_vec)),
#                                           rep(site_vec[5],length(site_vec)), rep(site_vec[6],length(site_vec)),
#                                           rep(site_vec[7],length(site_vec)), rep(site_vec[8],length(site_vec)),
#                                           rep(site_vec[9],length(site_vec)), rep(site_vec[10],length(site_vec)),
#                                           rep(site_vec[11],length(site_vec)), rep(site_vec[12],length(site_vec)),
#                                           rep(site_vec[13],length(site_vec)), rep(site_vec[14],length(site_vec)),
#                                           rep(site_vec[15],length(site_vec)), rep(site_vec[16],length(site_vec)),
#                                           rep(site_vec[17],length(site_vec)), rep(site_vec[18],length(site_vec)), 
#                                           rep(site_vec[19],length(site_vec))), stringsAsFactors = FALSE)
# site_dist_info$dest_site <- rep(site_vec, length(site_vec))
# site_dist_info$dist_N_to_S_m <- rep(NA, length(site_dist_info$org_site))
# site_dist_info$dist_N_to_S_km <- rep(NA, length(site_dist_info$org_site))
# site_dist_info$dist_S_to_N_m <- rep(NA, length(site_dist_info$org_site))
# site_dist_info$dist_S_to_N_km <- rep(NA, length(site_dist_info$org_site))
# site_dist_info$dist_avg <- rep(NA, length(site_dist_info$org_site))
# site_dist_info$dest_width <- rep(NA, length(site_dist_info$org_site))
# site_dist_info$d1_km <- rep(NA, length(site_dist_info$org_site))
# site_dist_info$d2_km <- rep(NA, length(site_dist_info$org_site))
# 
# # and find site distances! (for now, doing N to N, should do mid to mid...)
# for(i in 1:length(site_dist_info$org_site)) {
#   site_org = site_dist_info$org_site[i]
#   N_lat_org = (site_width_info %>% filter(site == site_org))$N_lat 
#   N_lon_org = (site_width_info %>% filter(site == site_org))$N_lon
#   S_lat_org = (site_width_info %>% filter(site == site_org))$S_lat 
#   S_lon_org = (site_width_info %>% filter(site == site_org))$S_lon
#   site_dest = site_dist_info$dest_site[i]
#   N_lat_dest = (site_width_info %>% filter(site == site_dest))$N_lat
#   N_lon_dest = (site_width_info %>% filter(site == site_dest))$N_lon
#   S_lat_dest = (site_width_info %>% filter(site == site_dest))$S_lat
#   S_lon_dest = (site_width_info %>% filter(site == site_dest))$S_lon
#   site_dist_info$dist_N_to_S_m[i] = distHaversine(c(N_lon_org, N_lat_org), c(S_lon_dest, S_lat_dest))
#   site_dist_info$dist_N_to_S_km[i] = site_dist_info$dist_N_to_S_m[i]/1000 
#   site_dist_info$dist_S_to_N_m[i] = distHaversine(c(N_lon_dest, N_lat_dest), c(S_lon_org, S_lat_org))
#   site_dist_info$dist_S_to_N_km[i] = site_dist_info$dist_S_to_N_m[i]/1000 
#   
#   site_dist_info$dist_avg[i] = (site_dist_info$dist_N_to_S_km[i] + site_dist_info$dist_S_to_N_km[i])/2
#   site_dist_info$dest_width[i] <- (site_width_info %>% filter(site == site_dest))$width_km
#   
#   site_dist_info$d1_km[i] <- site_dist_info$dist_avg[i] - site_dist_info$dest_width[i]/2
#   site_dist_info$d2_km[i] <- site_dist_info$dist_avg[i] + site_dist_info$dest_width[i]/2
#   # if(site_dist_info$dist_N_to_S_km[i] >= site_dist_info$dist_S_to_N_km[i]){
#   #   site_dist_info$d1_km[i] <- site_dist_info$dist_S_to_N_km[i]
#   #   site_dist_info$d2_km[i] <- site_dist_info$dist_N_to_S_km[i]
#   # } else {
#   #   site_dist_info$d1_km[i] <- site_dist_info$dist_N_to_S_km[i]
#   #   site_dist_info$d2_km[i] <- site_dist_info$dist_S_to_N_km[i]
#   # }
# }
# 
# # For sites going to self, put N-S width in
# for (i in 1:length(site_dist_info$org_site)) {
#   site_org = site_dist_info$org_site[i]
#   site_dest = site_dist_info$dest_site[i]
#   if (site_org == site_dest) {
#     N_lat_org = (site_width_info %>% filter(site == site_org))$N_lat 
#     N_lon_org = (site_width_info %>% filter(site == site_org))$N_lon
#     S_lat_org = (site_width_info %>% filter(site == site_org))$S_lat 
#     S_lon_org = (site_width_info %>% filter(site == site_org))$S_lon
#     site_dist_info$dist_N_to_S_m[i] = distHaversine(c(N_lon_org, N_lat_org), c(S_lon_org, S_lat_org))
#     site_dist_info$dist_N_to_S_km[i] = site_dist_info$dist_N_to_S_m[i]/1000 
#     site_dist_info$dist_S_to_N_m[i] = distHaversine(c(N_lon_org, N_lat_org), c(S_lon_org, S_lat_org))
#     site_dist_info$dist_S_to_N_km[i] = site_dist_info$dist_S_to_N_m[i]/1000 
#     site_dist_info$d1_km[i] = 0
#     site_dist_info$d2_km[i] = site_dist_info$dist_N_to_S_km[i]
#   }
# }
# 
# # Add in origin site alpha and geo order
# site_dist_info <- left_join(site_dist_info, site_vec_order, by = c('org_site' = 'site_name'))
# site_dist_info <- site_dist_info %>%
#   dplyr::rename(org_alpha_order = alpha_order, org_geo_order = geo_order)
# 
# # Add in destination site alpha and geo order
# site_dist_info <- left_join(site_dist_info, site_vec_order, by = c('dest_site' = 'site_name'))
# site_dist_info <- site_dist_info %>%
#   dplyr::rename(dest_alpha_order = alpha_order, dest_geo_order = geo_order)
