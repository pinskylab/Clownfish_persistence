# Find the width of each site and distances between them (for probability of dispersal), including distances when include a buffer for larval navigation

#################### Set-up: ####################
# Load packages
library(ggplot2)
library(geosphere)

# Pull data, functions, constants 
source(here::here('Code', 'Constants_database_common_functions.R'))

#################### Functions: ####################
# Find the maximum buffer distance to the north and south of each side without overlapping the next site's shadow
findBufferDistance <- function(org_site, site_buffer_df, northern_most_site_pos, southern_most_site_pos, max_buffer){  # org_site is focal site, site_buffer_df is data frame with site names and geographic position numbers and lat/lon for N and S-most anemones, nothern_most_site_pos is the geographic position number of the northern-most site sampled (1 = Palanas), southern_most_site_pos is the geographic position number of the southern-most site sampled (19 = Sitio Baybayon), max_buffer is max larval buffer distance considered
  org_site_pos <- (site_buffer_df %>% filter(site == org_site))$geo_order  # find geographic position (site order) of origin site
  # distances and buffer to the north
  if(org_site_pos == northern_most_site_pos) {  # if this is the northern-most site, no limits to the buffer to the north
    dist_to_N_m = max_buffer  
    max_buffer_N_m = max_buffer  # so assign the max larval buffer as the possible distance to the north
  } else {
    dest_site_N_pos <- org_site_pos - 1  # find the geo order of the site one to the north
    dist_to_N_m <- distHaversine(c((site_buffer_df %>% filter(site == org_site))$N_lon, (site_buffer_df %>% filter(site == org_site))$N_lat),  # find the distance between the northern edge of the origin site and the southern edge of the next site north
                                 c((site_buffer_df %>% filter(geo_order == dest_site_N_pos))$S_lon, (site_buffer_df %>% filter(geo_order == dest_site_N_pos))$S_lat))
    max_buffer_N_m = dist_to_N_m/2  # otherwise, assign the maximum buffer as half of the distance between the northern edge of this site and the southern edge of the next site north
  }
  
  # distances and buffer to the south
  if(org_site_pos == southern_most_site_pos) {  # if this is the northern-most site, no limits to the buffer to the north
    dist_to_S_m = max_buffer
    max_buffer_S_m = max_buffer
  } 
  else {
    dest_site_S_pos <- org_site_pos + 1  # find the geo order of the site one to the south
    dist_to_S_m <- distHaversine(c((site_buffer_df %>% filter(site == org_site))$S_lon, (site_buffer_df %>% filter(site == org_site))$S_lat),  # find the distance between the southern edge of the origin site and the northern edge of the next site south
                                 c((site_buffer_df %>% filter(geo_order == dest_site_S_pos))$N_lon, (site_buffer_df %>% filter(geo_order == dest_site_S_pos))$N_lat))
    max_buffer_S_m = dist_to_S_m/2 
  }
  
  out <- list(dist_to_N_m = dist_to_N_m, 
              max_buffer_N_m = max_buffer_N_m,
              dist_to_S_m = dist_to_S_m,
              max_buffer_S_m = max_buffer_S_m)
  return(out)
}


#################### Running things: ####################

##### Find lat/lon coordinates of anems at edges of sites
# Make a data frame for sites, boundary anems, lat lons of those anems, widths (should I just make this the same thing and put it in Constants_database_common_functions?)
site_edges_info <- site_edge_anems %>% filter(anem_loc %in% c("north", "south"))

# Pull out info on those anems so can use anem_table_ids in function to get lat/lon
site_width_anems <- anems_Processed %>%
  filter(anem_id %in% site_edges_info$anem_id) %>%
  distinct(anem_id, .keep_all = TRUE)  # just choose one visit to each anem (rather than trying to average GPS positions from multiple visits)

# Match those anem_table_ids to the anem_ids in site_edges_info
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

##### Find the width of each site and distance to edges of sampling area from the center of each site
# Calculate the distance between the edge anems for each site
site_width_info <- site_edges_info %>%
  filter(anem_loc == "north") %>%
  select(site, site_geo_order) %>%
  mutate(width_m = NA,  # width of site (distance between N and S-most anems)
         dist_to_N_edge_m = NA,  # distance from middle of site to N-most edge of sampling area
         dist_to_S_edge_m = NA)  # distance from middle of site to S-most edge of sampling area

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


##### Find distances between sites 
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
        
# Find the distances between the sites
for(i in 1:length(site_dist_info$org_site)) {
  site_org = site_dist_info$org_site[i]
  
  # get lat/lon coordinates for the northern edge, southern edge, and center of origin and destination site
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
  
  # find distances
  site_dist_info$dist_mid_to_S_m[i] = distHaversine(c(mid_lon_org, mid_lat_org), c(S_lon_dest, S_lat_dest))
  site_dist_info$dist_mid_to_S_km[i] = site_dist_info$dist_mid_to_S_m[i]/1000 
  site_dist_info$dist_mid_to_N_m[i] = distHaversine(c(mid_lon_org, mid_lat_org), c(N_lon_dest, N_lat_dest))
  site_dist_info$dist_mid_to_N_km[i] = site_dist_info$dist_mid_to_N_m[i]/1000 
  site_dist_info$dist_mid_to_mid_m[i] = distHaversine(c(mid_lon_org, mid_lat_org), c(mid_lon_dest, mid_lat_dest))
  site_dist_info$dist_mid_to_mid_km[i] = site_dist_info$dist_mid_to_mid_m[i]/1000
  
  site_dist_info$dest_width[i] = distHaversine(c(N_lon_dest, N_lat_dest), c(S_lon_dest, S_lat_dest))/1000
  
  site_dist_info$d1_km[i] <- min(c(site_dist_info$dist_mid_to_S_km[i], site_dist_info$dist_mid_to_N_km[i]))  # d1 is smaller distance to edge of site (since will depend on N-S order of sites)
  site_dist_info$d2_km[i] <- max(c(site_dist_info$dist_mid_to_S_km[i], site_dist_info$dist_mid_to_N_km[i]))
}

# Remove temporary values from for loop for neatness
rm(site_org, N_lat_org, N_lon_org, S_lat_org, S_lon_org, site_dest, N_lat_dest, N_lon_dest, S_lat_dest, S_lon_dest)

# For self-self distance, use 0 as min and 1/2 site width as max #### (not the only way I could do this) - now doing 1/2 site width as max because changed kernel normalization
site_dist_info <- site_dist_info %>%
  mutate(d1_km = case_when(org_site != dest_site ~ d1_km,
                           org_site == dest_site ~ 0),  # min distance is 0 for self-self distances
         d2_km = case_when(org_site != dest_site ~ d2_km,
                           org_site == dest_site ~ dest_width/2))  # max distance is 1/2 site width for self-self distances

# Add in origin site alpha and geo order
site_dist_info <- left_join(site_dist_info, site_vec_order, by = c('org_site' = 'site_name'))
site_dist_info <- site_dist_info %>%
  dplyr::rename(org_alpha_order = alpha_order, org_geo_order = geo_order)

# Add in destination site alpha and geo order
site_dist_info <- left_join(site_dist_info, site_vec_order, by = c('dest_site' = 'site_name'))
site_dist_info <- site_dist_info %>%
  dplyr::rename(dest_alpha_order = alpha_order, dest_geo_order = geo_order)


##### Find distances between edges of adjacent sites
site_buffer_info <- data.frame(site = site_vec_order$site_name, 
                               alpha_order = site_vec_order$alpha_order,
                               geo_order = site_vec_order$geo_order,
                               N_lat = NA,  # lat coordinate of north edge of site
                               N_lon = NA,  # lon coordinate of north edge of site
                               S_lat = NA,  # lat coordinate of north edge of site
                               S_lon = NA)  # lon coordinate of north edge of site

# Put in coordinates of north and south edges of sites                             
for(i in 1:length(site_buffer_info$site)) {
  site_val = site_buffer_info$site[i]
  site_buffer_info$N_lat[i] = (site_edges_info %>% filter(site == site_val, anem_loc == "north"))$lat
  site_buffer_info$N_lon[i] = (site_edges_info %>% filter(site == site_val, anem_loc == "north"))$lon
  site_buffer_info$S_lat[i] = (site_edges_info %>% filter(site == site_val, anem_loc == "south"))$lat
  site_buffer_info$S_lon[i] = (site_edges_info %>% filter(site == site_val, anem_loc == "south"))$lon
}

# Calculate max buffer distances for each site to N and S
site_buffer_info <- site_buffer_info %>%
  mutate(dist_to_N_m = NA,  # distance between this site and the next to the north
         dist_to_S_m = NA,  # distance between this site and the next to the south
         max_buffer_N_m = NA,  # maximum buffer to the north (half of the distance to nearest site to the north)
         max_buffer_S_m = NA)  # maximum buffer to the south (half of the distance to nearest site to the south))

for(i in 1:length(site_buffer_info$site)) {
  buffer_out <- findBufferDistance(site_buffer_info$site[i], site_buffer_info, northern_site_pos, southern_site_pos, max_larval_nav_buffer)
  site_buffer_info$dist_to_N_m[i] <- buffer_out$dist_to_N_m
  site_buffer_info$max_buffer_N_m[i] <- buffer_out$max_buffer_N_m
  site_buffer_info$dist_to_S_m[i] <- buffer_out$dist_to_S_m
  site_buffer_info$max_buffer_S_m[i] <- buffer_out$max_buffer_S_m
}

# Convert to km too
site_buffer_info <- site_buffer_info %>%
  mutate(dist_to_N_km = dist_to_N_m/1000,
         max_buffer_N_km = max_buffer_N_m/1000,
         dist_to_S_km = dist_to_S_m/1000,
         max_buffer_S_km = max_buffer_S_m/1000)

#################### Plots: ####################

#################### Saving output: ####################
save(site_width_info, file=here::here("Data/Script_outputs", "site_width_info.RData"))  # width of sites, distance to edges of sampling area
save(site_dist_info, file=here::here("Data/Script_outputs", "site_dist_info.RData"))  # distances between pairs of sites for kernel integration
save(sampling_area_edges, file=here::here("Data/Script_outputs", "sampling_area_edges.RData"))  # coordinates of northern-most anem and southern-most anem of sampling area
save(site_buffer_info, file=here::here("Data/Script_outputs", "site_buffer_info.RData"))  # maximum buffer from each site to the N and S


