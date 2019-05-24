# Make a site map

#################### Set-up: ####################
library(grid)
library(gridExtra)
#library(leaflet)
library(maps)
library(mapdata)
library(mapview)
library(rgdal)
#library(tmap)
library(maptools)
library(raster)
library(here)

#load(here::here('Data', 'Leyte.7z'))

# devtools::install_github('pinskylab/clownfish')  # try installing clownfish package...
# library(clownfish)

source(here::here('Code', 'Constants_database_common_functions.R'))

load(here::here('Data','AnemAllInfowLatLon.RData'))  # produces a file called 'anem.Processed' - not sure if the lat/lons are right...

#load(here::here('Data/Leyte', 'Leyte.shp'))  # 'Administrative Boundaries' data set, downloaded from http://philgis.org/province-page/leyte on 3/11/19, projection info: WGS 1984, Lat/Long

Leyte_shape <- readOGR(here::here('Data/Leyte'), 'Leyte')
Philippines_regions <- readOGR(here::here('Data/Regions'), 'Regions')  # check if it already knows the projection, otherwise figure out how to tell it

#################### Functions: ####################
# # Function to find the lat and long for an anem_id (based off of Michelle's sample_latlon function) rather than anem_table_id
# anem_id_latlong <- function(anem.id, anem.df, latlondata) {  # anem.id is one anem_id value, anem.df is dataframe with anemone info (could be straight from database anemone table), latlondata is table of GPX data from database (rather than making the function call it each time); will need to think a bit more clearly about how to handle different locations read for different visits to the same anem_id (or different with same anem_obs); for now, just letting every row in anem.Info get a lat-long
# 
#   # this is what causes the multiple entries - pulls multiple rows for a few anems (81) that have multiple entries for the same anem_table_id in the database
#   anem <- anem.df %>%
#     filter(anem_id == anem.id)  # get the relevant dive, time, site, etc. info for this anem_table_id
# 
#   # find the lat long for this anem observation
#   latloninfo <- latlondata %>%
#     filter(as.character(gps_date) %in% anem$date & unit == anem$gps)
#   # Multiple dates going on here b/c anems seen different times - for this, could just pick one, 
#   
#   %>% #filter out just the GPS unit associated with this anem observation (added since previous time)
#     filter(gps_hour == anem$anem_hour & gps_min == anem$anem_min) %>%
#     mutate(lat = as.numeric(lat)) %>%
#     mutate(lon = as.numeric(lon))
# 
#   #pull duplicates (so if sat in one place for enough time that multiple readings recorded there)
#   #(there are more digits in the lats and lons than show up on the screen so sometimes things look like duplicates but aren't)
#   dups_lat <- which(duplicated(latloninfo$lat)) #vector of positions of duplicate values
#   dups_lon <- which(duplicated(latloninfo$lon))
# 
#   #either take the mean of the lat/lon readings or the duplicated values, depending if there are duplicate points
#   if(length(dups_lat) == 0) { #if all latitude points are different
#     anem$lat <- round(mean(latloninfo$lat), digits = 5) #take the mean of the latitude values (digits = 5 b/c that is what Michelle had)
#     anem$lon <- round(mean(latloninfo$lon), digits = 5) #take the mean of the longitude values
#   }else{
#     anem$lat <- latloninfo$lat[dups_lat[1]] #if are duplicates, take the value of the first duplicated point
#     anem$lon <- latloninfo$lon[dups_lon[1]]
#     print(paste("Dups in lat lons at anem_table_id", anem$anem_table_id, "on", anem$date, "with lat", anem$lat, sep = " ")) #just have this while trouble-shooting repeat entries in the database
#   }
# 
#   return(anem)
# 
# }

#################### Running things: ####################
##### Get anem lat/lon info for edge anems at each site
# Join with anem info 
site_edge_anems_for_map <- left_join(site_edge_anems, anems_Processed %>% select(anem_table_id, anem_id, anem_obs, date, gps, anem_day, anem_hour, anem_min, anem_sec), by='anem_id')

# Just pick one of the anem_table_ids for each anem, make placeholder columns for lat/lon
site_edge_anems_for_map <- site_edge_anems_for_map %>%
  group_by(site) %>%
  distinct(anem_loc, .keep_all = TRUE) %>%
  ungroup()

####### THIS ISN'T WORKING - anemid_latlong functionin Constants_database_common_functions isn't working!!! Not sure why not!!!
# # Find lat/lons for the anems at the boundaries/middle using anem_table_ids for one of the observations of that anem
# site_edge_anems_for_map <- site_edge_anems_for_map %>%
#   mutate(lat = NA, lon = NA)  # make placeholder lat/lon columns
# 
# for(i in 1:(length(site_edge_anems_for_map$site))) {
#   out_lls <- anemid_latlong(site_edge_anems_for_map$anem_table_id[i], anems_Processed, gps_Info)
#   site_edge_anems_for_map$lat[i] = anemid_latlong(site_edge_anems_for_map$anem_table_id[i], anems_Processed, gps_Info, 'latlon')$lat
#   site_edge_anems_for_map$lon[i] = anemid_latlong(site_edge_anems_for_map$anem_table_id[i], anems_Processed, gps_Info, 'latlon')$lon
# }
# 
# #### WHY WON'T MY LAT-LON FINDING FUNCTION FIND ANY LAT LONS??????
# 
# # First, find an anem_table_id for each anem_id (so can more easily link to gps)
# site_edge_anems <- site_edge_anems %>%
#   mutate(example_anem_table_id )

# Join with lat/lons from old anem.Processed dataframe - doing this for now b/c function not working 
site_edge_anems_for_map <- left_join(site_edge_anems_for_map, anem.Processed %>% select(anem_table_id, lat, lon), by='anem_table_id')

# Just pull out north and south
site_edge_anems_for_map <- site_edge_anems_for_map %>% 
  select(site, anem_loc, lat, lon) %>%
  filter(anem_loc != 'mid')

# And map!
# Trim Philippines regions to just Leyte (or something) - this doesn't work, maybe I need Philippines_regions to be a raster or something?
ext <- extent(10,11,123,127)
cropped_region <- crop(Philippines_regions, ext)

plot(Leyte_shape)

plot(Philippines_regions)

plot(cropped_region)


points(x in lon, y in lat, color, pch)


map('world', region = c('Philippines', 'Malaysia', 'China', 'Indonesia', 'Thailand', 'Cambodia', 'Laos'))


map("world", col="grey70", fill=T, border="white", lwd=0.3, xlim=c(90,140), ylim=c(-10,35))
map.cities(x=world.cities, country = "Philippines", capitals=1)

# not working, plus way too big of a map for that for now...
points(x=site_edge_anems_for_map$lon, y=site_edge_anems_for_map$lat)



map.cities(x=world.cities, country = c("Philippines", "Thailand", "Indonesia"), capitals=1)
                                       
                                       
                                       ("cities", cities = c("Manila", "Jakarta", "Bangkok"))

# Plot world countries - from Chris Fig 6
map("world", col="grey85", fill=T, border="white", lwd=0.3,
    xlim=xlim, ylim=ylim, add=T)

# And map!
#map('world', regions = 'Philippines:Leyte')




# # And map!
# leaflet((data=site_edge_anems_for_map %>% filter(anem_loc == 'north'))) %>% 
#   addTiles() %>%
#   #addMarkers() %>%
#   addRectangles(lng1 = (site_edge_anems_for_map %>% filter(anem_loc == 'north'))$lon, 
#                 lat1 = (site_edge_anems_for_map %>% filter(anem_loc == 'north'))$lat,
#                 lng2 = (site_edge_anems_for_map %>% filter(anem_loc == 'south'))$lon,
#                 lat2 = (site_edge_anems_for_map %>% filter(anem_loc == 'south'))$lat,
#                 label = (site_edge_anems_for_map %>% filter(anem_loc == 'north'))$site)
# 
# addRectangles(
#   lng1=-118.456554, lat1=34.078039,
#   lng2=-118.436383, lat2=34.062717,
#   fillColor = "transparent"
