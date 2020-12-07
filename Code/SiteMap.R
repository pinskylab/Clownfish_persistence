# Make figure with metrics schematics and life cycle schematic (Fig. 1), map and photo figure (Fig. 2), and map for recruit scaling schematic (Fig. D2)

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))  # pulls in processed clownfish and anemone data (common to all scripts)

# Load relevant libraries
library(mapdata)  # provides Philippines outline
library(rgdal)
library(cowplot)  # for arranging plots into multi-paneled figure
library(raster)
library(magick)
library(ggplot2)
library(RColorBrewer)

#################### Figure 1 (metrics schematics and life cycle schematic): ####################

# Life cycle schematic
life_cycle_plot <- ggdraw() +
  draw_image(here::here("Plots/LifeCycleSchematic", "lifecycle3.png"))

# Metrics schematic parts
LRP_schematic_plot <- ggdraw() +
  draw_image(here::here("Plots/LifeCycleSchematic", "metapopulation_diagram_LRP.png"))
SP_schematic_plot <- ggdraw() +
  draw_image(here::here("Plots/LifeCycleSchematic", "metapopulation_diagram_SP.png"))
NP_schematic_plot <- ggdraw() +
  draw_image(here::here("Plots/LifeCycleSchematic", "metapopulation_diagram_NP.png"))
LR_schematic_plot <- ggdraw() +
  draw_image(here::here("Plots/LifeCycleSchematic", "metapopulation_diagram_LR.png"))

metrics_schematic <- plot_grid(LRP_schematic_plot, NULL, SP_schematic_plot, NULL,
                               NP_schematic_plot, NULL, LR_schematic_plot,
                               rel_heights = c(1, 0.1, 1, 0.1, 1, 0.1, 1),
                               labels = c("a","","b","","c","","d"), nrow=7)

# Put them together
pdf(file=here::here("Plots/LifeCycleSchematic","metrics_life_cycle_schematics.pdf"), height=6)
plot_grid(metrics_schematic, NULL, life_cycle_plot, labels=c("","","e"), ncol=3, rel_widths = c(1,0.02,1))
dev.off()

#################### Figure 2 (map + photo): ####################
############### Figure 2 set-up: ###############

##### Load shapefiles for hulls and coast (downloaded from amphiprion, hulls created by Mario in 2016), missing Sitio Hicgop and Hicgop South
# Site hulls
patchfiles = list.files(here::here("Data/Map_data/Site_hulls/"), pattern = ".shp")  # get list of files with site hulls
patches = vector("list", length(patchfiles))
for(i in 1:length(patches)){
  temp = readOGR(paste(here::here("Data/Map_data/Site_hulls/"), patchfiles[i], sep=''))
  patches[[i]] = spTransform(temp, CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'))
}

# Coastline
coast = readOGR(here::here("Data/Map_data/", "PHL_adm3_leyte_studyarea_trimcoastline.shp"))  # downloaded from Amphiprion, need to find source

# Philippines
country = map("world", fill = TRUE, col = "grey")

##### Set constants
kmlen = 1/111.12  # length of a km, in °N

# Limits for maps
xlims = c(124.64, 124.80) # site map area
ylims = c(10.63, 10.87)
xlims_Palanas = c(124.710, 124.713)  # limits for zoomed-in map of Palanas patch
ylims_Palanas = c(10.872, 10.875)
xlims_PHL = c(118,127)  # extent of Philippines map
ylims_PHL = c(6,18)
ylims_schematic = c(10.58, 10.92)  # for filling in missing parts of the land on the map
xlims_schematic = c(124.64, 124.85)

# Set coordinates for boxes, lines, etc. (group numbers are unimportant, just so it plots easier)
corner_coords = data.frame(long = c(124.76, 124.81, 124.81), lat = c(10.6, 10.65, 10.6), group = 5)  # missing piece of coast to add back in
km_line_coords = data.frame(long = rep(xlims[1],2), lat = c(ylims[1], ylims[1]+5*kmlen), group = 6)  # 5-km scale line
red_box_coords <- data.frame(long = c(xlims_Palanas, rev(xlims_Palanas)), lat = c(rep(ylims_Palanas[1],2), rep(ylims_Palanas[2],2)), group = 8)
zoomed_area_coords <- data.frame(long = c(xlims_Palanas[1]-0.2, xlims_Palanas[2]+0.2, xlims_Palanas[2]+0.2, xlims_Palanas[1]-0.2),
                                   lat = c(ylims_Palanas[1]-0.2, ylims_Palanas[1]-0.2, ylims_Palanas[2]+0.2, ylims_Palanas[2]+0.2), group = 9)
upper_corner_coords = data.frame(long = c(124.81, 124.69, 124.81), lat = c(ylims[2], ylims_schematic[2]+0.2, ylims_schematic[2]+0.2), group=10)  # fill in missing chunk of land outside of sampling area in map for scaling schematic
inland_box_coords = data.frame(long = c(124.80, 124.80, 124.85, 124.85), lat = c(ylims_schematic[1], ylims_schematic[2], ylims_schematic[2], ylims_schematic[1]), group = 11)  # extend inland a bit for scaling up schematic

# Set colors
land_color = "grey90"
patch_color = "blue3"
colslist = brewer.pal(3, name="Dark2")  # checked on colorbrewer2 website, is color-blind safe
colslist[3] = "grey35"  # probably can't see blue patches on a blue background

# Sites in zoomed-in box
zoom_site = "Palanas"

# Other clownfish species (though decided not to show these in anemone map)
other_clownfish_spp = c("APFR","APOC","APPE","PRBI")

# Set which set of anems to show
anem_fish_species_to_show = c("empty","A. clarkii")

# Shift right for site labels from site center
long_shift = 0.008

############### Process anems and patch ids: ###############
##### Process anems - find just those in zoomed sites (Palanas) in the anem sample year (winter 2015), attach lat-lon coords
data_anems <- anems_Processed_all %>% filter(site %in% zoom_site, year == 2015, month %in% winter_months) %>%  # filter out just anems from 2015 anem survey
  mutate(lat = NA, lon = NA)  

# Find the lat lon coordinates of each of those anems
for(i in 1:length(data_anems$lon)) {
  anem_out <- anemid_latlong(data_anems$anem_table_id[i], anems_Processed_all, gps_Info)
  data_anems$lat[i] <- anem_out$lat
  data_anems$lon[i] <- anem_out$lon
}

# Match with clownfish so know which species was on each anem
fish_to_join <- rbind(fish_db %>% dplyr::select(fish_spp, anem_table_id), fish_seen_db %>% dplyr::select(fish_spp, anem_table_id)) %>%
  distinct(anem_table_id, .keep_all = TRUE)
data_anems <- left_join(data_anems, fish_to_join, by = "anem_table_id")

##### Make list of patch ids from N-S and coordinates to plot id on plot so can label in caption
patch_id_coords <- data.frame(site = site_vec_order$site_name,
                             patch_id = site_vec_order$geo_order)

# Add in coords of center of patch
patch_id_coords <- left_join(patch_id_coords, site_centers, by = "site")

# Shift coords to the right slightly from the center of the patch, then make fine-scale adjustments
patch_id_coords <- patch_id_coords %>% 
  mutate(id_lat = lat,
         id_lon = lon + long_shift)
patch_id_coords$id_lon[1] = patch_id_coords$id_lon[1]+0.003  # move Cabatoan label right slightly
patch_id_coords$id_lat[8] = patch_id_coords$id_lat[8]+0.001  # move N. Magbangon label up slightly
patch_id_coords$id_lon[9] = patch_id_coords$id_lon[9]+0.0001  # move Palanas label right slightly
patch_id_coords$id_lat[11] = patch_id_coords$id_lat[11]-0.0015  # move Poroc San Flower label down slightly
patch_id_coords$id_lat[14] = patch_id_coords$id_lat[14]-0.0015  # move Sitio Baybayon label down slightly
patch_id_coords$id_lon[6] = patch_id_coords$id_lon[6]+0.005  # move Haina label right slightly

############### Make sub-plots and overall plot: ###############
##### Main plot of sites
site_area <- ggplot(data = coast, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims, ylim = ylims, 1) +  # had 1.3 earlier, depends where you are on the globe?
  geom_polygon(colour = land_color, fill = land_color) +
  geom_polygon(data = corner_coords, aes(x = long, y = lat, group = group), fill = land_color, colour = land_color) +
  geom_polygon(data = patches[[1]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[2]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[3]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[4]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[5]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[6]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[7]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[8]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[9]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[10]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[11]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[12]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[13]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[14]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[15]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[16]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_line(data = km_line_coords, aes(x = long, y = lat, group = group), color = "black", lwd = 2) +
  annotate(geom = 'text', x = xlims[1]+0.003, y = 10.65, label = '5 km', size = 5, hjust = 0) +  # 5 km scale bar
  xlab("Longitude (°E)") + ylab('Latitude (°N)') +
  geom_polygon(data = red_box_coords, aes(x = long, y = lat, group = group), fill = NA, color = "red", lwd = 0.5)  # red box around Palanas

# Add in site id numbers
for(i in 1:length(patch_id_coords$site)) {
  site_area <- site_area +
    annotate("text", x=patch_id_coords$id_lon[i], y=patch_id_coords$id_lat[i], label=patch_id_coords$patch_id[i], size=2.8)
}

##### Inset map of Philippines
inset_map <- ggplot(data =  country, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims_PHL, ylim = ylims_PHL, 1) +  # had 1.3 for ratio earlier, depends where you are on the globe?
  geom_polygon(colour = "dark grey", fill = "light grey") +
  panel_border(colour = "black", size = 1, linetype = 1, remove = FALSE) +  # this is from cowplot, not sure why it didn't work with theme commands...
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank()) +
  geom_polygon(data = zoomed_area_coords, aes(x = long, y = lat, group = group), fill = NA, color = "black", lwd = 0.75) +
  annotate(geom = 'text', x = 121.8, y = 8.3, label = 'Philippines', cex=5) 

##### Example habitat in a site (Palanas)
# Identify fish species on each anem in example site
anems_to_plot <- data_anems %>%
  mutate(species = case_when(is.na(fish_spp) ~ "empty",
                             fish_spp == "APCL" ~ "A. clarkii",
                             fish_spp %in% other_clownfish_spp ~ "other clownfish")) %>%
  filter(species %in% anem_fish_species_to_show) %>%
  mutate(group = 11)

# Example site to show is Palanas (7 in current patchlist) 
Palanas_patch <- 7

# Make plot
Albuera_patch <- ggplot(data = coast, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims_Palanas, ylim = ylims_Palanas, 1.3) +
  geom_polygon(color = land_color, fill = land_color) +
  geom_polygon(data = patches[[Palanas_patch]], aes(x = long, y = lat, group = group), color = patch_color, fill = patch_color, alpha = 0.2) +  # Palanas
  geom_point(data = anems_to_plot, aes(x = lon, y = lat, group = group, color = species), alpha=0.7, size=0.5) +
  scale_color_manual(values = colslist[1:3]) +
  theme(legend.position = c(0.1, 0.25)) + 
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=2, alpha=1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=8), axis.text.y = element_text(size=8)) +
  xlab("") + ylab("") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

##### Clownfish photo
photo_plot <- ggdraw() +
  draw_image(here::here("Data", "Fish.jpg"))

##### Arrange all together
# Add inset to site map 
sites_with_inset <- ggdraw(site_area) + 
  draw_plot(inset_map + theme(legend.justification = "top"), 0.27, 0.27, 0.5, 0.5)

# Put the zoomed in patches and photo together
right_side <- plot_grid(Albuera_patch, photo_plot, labels = c("b","c"), nrow=2)

# Add map with sites in
pdf(file=here::here("Plots/FigureDrafts","Map_and_photo_2.pdf"), height=6)
plot_grid(sites_with_inset, NULL, right_side, labels=c("a","",""), ncol=3, rel_widths = c(1.5,0.02,1))
dev.off()

##### Map for connectivity paper - shrink to one column width, adjust text size
width_1_column = 80  # mm

# Map of sites
site_area_connectivity <- ggplot(data = coast, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims, ylim = ylims, 1) +  # had 1.3 earlier, depends where you are on the globe?
  geom_polygon(colour = land_color, fill = land_color) +
  geom_polygon(data = corner_coords, aes(x = long, y = lat, group = group), fill = land_color, colour = land_color) +
  geom_polygon(data = patches[[1]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[2]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[3]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[4]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[5]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[6]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[7]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[8]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[9]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[10]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[11]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[12]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[13]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[14]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[15]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_polygon(data = patches[[16]], aes(x = long, y = lat, group = group), fill = patch_color, colour = patch_color) +
  geom_line(data = km_line_coords, aes(x = long, y = lat, group = group), color = "black", lwd = 2) +
  annotate(geom = 'text', x = xlims[1]+0.003, y = 10.65, label = '5 km', size = 3, hjust = 0) +  # 5 km scale bar
  xlab("Longitude (°E)") + ylab('Latitude (°N)') +
  theme(axis.text=element_text(size=rel(0.5)),axis.title=element_text(size=rel(0.5)))
 
# Add in site id numbers
for(i in 1:length(patch_id_coords$site)) {
  site_area_connectivity <- site_area_connectivity +
    annotate("text", x=patch_id_coords$id_lon[i], y=patch_id_coords$id_lat[i], label=patch_id_coords$patch_id[i], size=2.8)
}

# Inset map of Philippines
inset_map_connectivity <- ggplot(data =  country, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims_PHL, ylim = ylims_PHL, 1) +  # had 1.3 for ratio earlier, depends where you are on the globe?
  geom_polygon(colour = "dark grey", fill = "light grey") +
  panel_border(colour = "black", size = 1, linetype = 1, remove = FALSE) +  # this is from cowplot, not sure why it didn't work with theme commands...
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank()) +
  geom_polygon(data = zoomed_area_coords, aes(x = long, y = lat, group = group), fill = NA, color = "black", lwd = 0.75) +
  annotate(geom = 'text', x = 120.8, y = 8.7, label = 'Philippines', cex=3.5) 

# Add inset to edited site map
sites_with_inset_connectivity <- ggdraw(site_area_connectivity) +
  draw_plot(inset_map_connectivity + theme(legend.justification = "top"), 0.27, 0.27, 0.5, 0.5)

# Save plot as png, then pdf
ggsave(sites_with_inset_connectivity, filename=here::here("Plots","Site_map_connectivity.pdf"), dpi=300, device=cairo_pdf, width=width_1_column, units = "mm")


#################### Figure D2 (map for scaling up recruits schematic): ####################

##### Map of just sites for scaling-up-recruits figure
pdf(file=here::here("Plots/FigureDrafts", "Coastline_patch_map.pdf"), width=5)
ggplot(data = coast, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims_schematic, ylim = ylims_schematic, 1) +  # had 1.3 earlier, depends where you are on the globe?
  geom_polygon(colour = "grey", fill = "grey") +
  geom_polygon(data = corner_coords, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey") +
  geom_polygon(data = upper_corner_coords, aes(x = long, y = lat, group = group), fill = "grey", color = "grey") +
  geom_polygon(data = patches[[1]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[2]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[3]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[4]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[5]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[6]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[7]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[8]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[9]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[10]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[11]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[12]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[13]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[14]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[15]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  geom_polygon(data = patches[[16]], aes(x = long, y = lat, group = group), fill = colslist[3], colour = colslist[3]) +
  theme_nothing()
dev.off()

