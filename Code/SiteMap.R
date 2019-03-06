# Make a site map

#################### Set-up: ####################
library(grid)
library(gridExtra)
library(leaflet)

devtools::install_github('pinskylab/clownfish')  # try installing clownfish package...
 
library(clownfish)


source(here::here('Code', 'Constants_database_common_functions.R'))


#################### Running things: ####################
# Get anem lat/lon info for edge anems at each site
site_edge_anems <- site_edge_anems %>%
  mutate(lat = get)
leaflet() %>%
  addTiles()

get_anem() %>% filter(anem_id == 2001)
