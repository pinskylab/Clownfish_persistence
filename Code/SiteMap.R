# Make a site map

#################### Set-up: ####################
# Load relevant libraries
library(mapdata)  # provides Philippines outline
library(rgdal)
library(cowplot)  # arranging plots into multi-paneled figure
library(raster)
library(magick)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(here)

##### Load anem lat/lons - REDO SO GET THE ANEM INFO FROM DATABASE + A SCRIPT?
data_anems = as.data.frame(read.csv(here::here("Data", "GPSSurvey.anemlatlong2015-12-16.csv"), row.names=1, stringsAsFactors = FALSE))

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
xlims = c(124.64, 124.80) # site areas
ylims = c(10.63, 10.87)
xlims_Alb = c(124.708, 124.7174) # Albuera patch 
ylims_Alb = c(10.867, 10.875)
xlims_PHL = c(118,127)  # extent of Philippines map
ylims_PHL = c(6,18)
ylims_schematic = c(10.58, 10.92)
xlims_schematic = c(124.64, 124.85)

# xlims_schematic_panel = c(0,2)
# ylims_schematic_panel = c(0,4)
# xlims_coast_panel = c()

# Set coordinates for boxes, lines, etc. (group numbers are unimportant, just so it plots easier)
corner_coords = data.frame(long = c(124.76, 124.81, 124.81), lat = c(10.6, 10.65, 10.6), group = 5)  # missing piece of coast to add back in
km_line_coords = data.frame(long = rep(xlims[1],2), lat = c(ylims[1], ylims[1]+5*kmlen), group = 6)  # 5-km scale line
hab_line_coords = data.frame(long = rep(xlims[1],2), lat = c(10.85, 10.90), group = 7)  # line for habitat patch legend
red_box_coords <- data.frame(long = c(xlims_Alb, rev(xlims_Alb)), lat = c(rep(ylims_Alb[1],2), rep(ylims_Alb[2],2)), group = 8)
zoomed_area_coords <- data.frame(long = c(xlims_Alb[1]-0.2, xlims_Alb[2]+0.2, xlims_Alb[2]+0.2, xlims_Alb[1]-0.2),
                                   lat = c(ylims_Alb[1]-0.2, ylims_Alb[1]-0.2, ylims_Alb[2]+0.2, ylims_Alb[2]+0.2), group = 9)
#upper_corner_coords = data.frame(long = c(xlims[2], 124.79, xlims[2]), lat = c(ylims[2], ylims_schematic[2], ylims_schematic[2]), group=10)  # fill in missing chunk of land outside of sampling area in map for scaling schematic
upper_corner_coords = data.frame(long = c(124.81, 124.69, 124.81), lat = c(ylims[2], ylims_schematic[2]+0.2, ylims_schematic[2]+0.2), group=10)  # fill in missing chunk of land outside of sampling area in map for scaling schematic
inland_box_coords = data.frame(long = c(124.80, 124.80, 124.85, 124.85), lat = c(ylims_schematic[1], ylims_schematic[2], ylims_schematic[2], ylims_schematic[1]), group = 11)  # extend inland a bit for scaling up schematic

# Set colors
colslist = brewer.pal(3, name="Dark2")

#################### Make sub-plots and overall plot: ####################
##### Main plot of sites
site_area <- ggplot(data = coast, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims, ylim = ylims, 1) +  # had 1.3 earlier, depends where you are on the globe?
  geom_polygon(colour = "grey", fill = "grey") +
  geom_polygon(data = corner_coords, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey") +
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
  geom_line(data = km_line_coords, aes(x = long, y = lat, group = group), color = "black", lwd = 2) +
  #annotate(geom = 'text', x = xlims[1]+0.025, y = 10.65, label = '5 km', size = 4) +
  annotate(geom = 'text', x = xlims[1]+0.012, y = 10.65, label = '5 km', size = 5) +
  xlab("Longitude (°E)") + ylab('Latitude (°N)') +
  #geom_polygon(data = red_box_coords, aes(x = long, y = lat, group = group), fill = NA, color = "red", lwd = 1) +
  geom_line(data = hab_line_coords, aes(x = long, y = lat, group = group), color = colslist[3], lwd = 2) +
  #annotate(geom = "text", x = xlims[1]+0.033, y = 10.865, label = "Habitat \n patches", size = 4)
  annotate(geom = "text", x = xlims[1]+0.012, y = 10.865, label = "Habitat \n patches", size = 5)

##### Inset map of Philippines
inset_map <- ggplot(data =  country, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims_PHL, ylim = ylims_PHL, 1) +  # had 1.3 for ratio earlier, depends where you are on the globe?
  geom_polygon(colour = "dark grey", fill = "light grey") +
  panel_border(colour = "black", size = 1, linetype = 1, remove = FALSE) +  # this is from cowplot, not sure why it didn't work with theme commands...
  #theme(panel.border = element_rect(colour = "black", fill = NA)) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank()) +
  geom_polygon(data = zoomed_area_coords, aes(x = long, y = lat, group = group), fill = NA, color = "black", lwd = 1) +
  #annotate(geom = 'text', x = 120.8, y = 8.3, label = 'Philippines', cex=4) 
  annotate(geom = 'text', x = 120.8, y = 8.3, label = 'Philippines', cex=5) 

##### Example habitat in a site
# Filter out anems in the example sites, 
anems_to_plot <- data_anems %>%
  filter(Name %in% c("Palanas", "Wangag")) %>%
  filter(Spp %in% c("","APCL")) %>%
  mutate(species = case_when(Spp == "" ~ "empty",
                             Spp == "APCL" ~ "A. clarkii")) %>%
  mutate(group = 11)

# Example sites to show are Palanas (7 in current patchlist) and Wangag (16 in current patchlist) 
Palanas_patch <- 7
Wangag_patch <- 16

Albuera_patch <- ggplot(data = coast, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims_Alb, ylim = ylims_Alb, 1.3) +
  geom_polygon(color = "grey", fill = "grey") +
  geom_polygon(data = patches[[Palanas_patch]], aes(x = long, y = lat, group = group), color = "blue", fill = "blue", alpha = 0.2) +  # Palanas
  geom_polygon(data = patches[[Wangag_patch]], aes(x = long, y = lat, group = group), color = "blue", fill = "blue", alpha = 0.2) +  # Wangag
  geom_polygon(data = red_box_coords, aes(x = long, y = lat, group = group), fill = NA, color = "red", lwd = 1) +
  geom_point(data = anems_to_plot, aes(x = lon, y = lat, group = group, color = species), alpha=0.7, size=0.6) +
  scale_color_manual(values = colslist[1:2]) +
  theme(legend.position = c(0.2, 0.3)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") + ylab("")

##### Clownfish photo
photo_plot <- ggdraw() +
  draw_image(here::here("Data", "Fish.jpg"))

##### Arrange all together
# Add inset to site map 
sites_with_inset <- ggdraw(site_area) + 
  draw_plot(inset_map + theme(legend.justification = "top"), 0.3, 0.3, 0.5, 0.5)

# Put the two maps together
top_row <- plot_grid(sites_with_inset, Albuera_patch, labels = c("a", "b", ncol = 1, align = "h"))

pdf(file=here::here("Plots/FigureDrafts", "Map_and_photo.pdf"))
#plot_grid(sites_with_inset, right_column, labels = c("a", ""), ncol = 2)
plot_grid(top_row, photo_plot, labels = c("","c"), nrow=2, rel_heights = c(1.7,1))
dev.off()


##### Map with inset (but no red box) for poster
sites_with_inset_for_poster <- ggdraw(site_area) + 
  draw_plot(inset_map + theme(legend.justification = "top"), 0.25, 0.2, 0.5, 0.5) +
  theme(text = element_text(size=40))

pdf(file=here::here("Plots/Poster_presentation_plots","Map.pdf"), width = 8, height=11)
sites_with_inset_for_poster
dev.off()

##### Map of just sites for scaling-up-recruits figure
pdf(file=here::here("Plots/FigureDrafts", "Coastline_patch_map.pdf"), width=5)
ggplot(data = coast, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims_schematic, ylim = ylims_schematic, 1) +  # had 1.3 earlier, depends where you are on the globe?
  geom_polygon(colour = "grey", fill = "grey") +
  geom_polygon(data = corner_coords, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey") +
  geom_polygon(data = upper_corner_coords, aes(x = long, y = lat, group = group), fill = "grey", color = "grey") +
  #geom_polygon(data = inland_box_coords, aes(x = long, y = lat, group = group), fill = "grey", color = "grey") +
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

# Put example patch and photo together - right side of figure
#right_column <- plot_grid(Albuera_patch, photo_plot, labels = c("b", "c"), ncol = 1, align = "v", rel_heights = c(1.5, 1))
#plot_grid(sites_with_inset, Albuera_patch, photo_plot, labels = c("a","b","c"), ncol = 3, align = "h")


#################### Old code: ####################
# ######### Altered Malin code that runs (except for inset map part)
# # get google map
# #g <- gmap(insetext, type='roadmap', lonlat=TRUE, scale=2)
# 
# # map	
# quartz(height=5, width=6.5)
# # pdf(height=5, width=6, file='Figures/anem_map2.pdf')
# par(mai=c(0.4,0.5,0.1, 0.2), las=1, cex.axis=0.8, mgp=c(2.5,0.5,0), tcl=-0.3)
# layout(matrix(c(1,1,1,2,3,4), nrow=3), widths=c(3,2))
# plot(0,0, xlim=xlims, ylim=ylims, xlab='', ylab='')
# #map('worldHires', 'Philippines', add=TRUE, fill=TRUE, col='grey', border=NA)
# plot(coast, col='grey', border=NA, add=TRUE)
# polygon(x=c(124.76, 124.81, 124.81), y=c(10.6, 10.65, 10.6), col='grey', border=NA) # fill in a missing corner
# #plot(patches[[1]], add = TRUE)
# for(i in 1:length(patches)) plot(patches[[i]], add=TRUE, col=colslist[3], border=colslist[3], lwd=2)
# #for(i in 1:5) plot(patches[[i]], add=TRUE, col=colslist[3], border=colslist[3], lwd=2)
# #points(data$lon, data$lat, pch=16, cex=0.4, col=cols)
# lines(x = rep(xlims[1],2), y=c(ylims[1], ylims[1]+5*kmlen), lwd=5, lend=3)
# text(labels="5 km", x=xlims[1]+0.015, y=10.65, cex=1.5)
# polygon(c(xlims2, rev(xlims2)), rep(ylims2, c(2,2)), border='red', lwd=1.5)
# polygon(c(xlims3, rev(xlims3)), rep(ylims3, c(2,2)), border='red', lwd=1.5)
# polygon(c(xlims4, rev(xlims4)), rep(ylims4, c(2,2)), border='red', lwd=1.5)
# mtext(side=1, text='Longitude (°E)', line=1.8, cex=0.8)
# mtext(side=2, text='Latitude (°N)', line=2.5, cex=0.8, las=0)
# legend('topleft', legend='Habitat patch', lwd=3, col=colslist[3], bty='n')
# 
# # Albuera
# par(mai=c(0.3,0.4,0.1, 0.1), las=1)
# plot(0,0, xlim=xlims4, ylim=ylims4, xlab='', ylab='')
# for(i in 1:length(patches)) plot(patches[[i]], add=TRUE, col=patchfadecol, border=patchfadecol, lwd=7)
# points(data$lon[inds], data$lat[inds], pch=16, cex=0.8, col=cols[inds])
# #map('worldHires', 'Philippines', add=TRUE, fill=TRUE, col='grey',border=NA)
# plot(coast, col='grey', border=NA, add=TRUE)
# lines(x = rep(xlims4[1],2), y=rep(10.85,2)+c(0,0.5*kmlen), lwd=3, lend=3)
# text(labels="500 m", x=124.708, y=10.852, cex=1)
# 
# # Marcos
# plot(0,0, xlim=xlims3, ylim=ylims3, xlab='', ylab='')
# for(i in 1:length(patches)) plot(patches[[i]], add=TRUE, col=patchfadecol, border=NA, lwd=2)
# points(data$lon[inds], data$lat[inds], pch=16, cex=0.8, col=cols[inds])
# #map('worldHires', 'Philippines', add=TRUE, fill=TRUE, col='grey',border=NA)
# plot(coast, col='grey', border=NA, add=TRUE)
# lines(x = rep(xlims3[1],2), y=rep(ylims3[1],2)+c(0,0.1*kmlen), lwd=3, lend=3)
# text(labels="100 m", x=xlims3[1]+0.001, y=ylims3[1]+0.0005, cex=1)
# legend('topleft', legend=c(colsnames, 'Habitat patch'), col=c(colslist[c(1,2)], patchfadecol), pch=c(16, 16, NA), lwd = c(NA, NA, 4), cex=1, text.font=1, bty='o', bg='white', box.col='black')
# 
# # Punta
# plot(0,0, xlim=xlims2, ylim=ylims2, xlab='', ylab='')
# for(i in 1:length(patches)) plot(patches[[i]], add=TRUE, col=patchfadecol, border=NA)
# points(data$lon[inds], data$lat[inds], pch=16, cex=0.8, col=cols[inds])
# #map('worldHires', 'Philippines', add=TRUE, fill=TRUE, col='grey', border=NA)
# plot(coast, col='grey', border=NA, add=TRUE)
# lines(x = rep(124.775,2), y=rep(ylims2[1],2)+c(0,0.1*kmlen), lwd=3, lend=3)
# text(labels="100 m", x=124.7753, y=ylims2[1]+0.0005, cex=1)
# 
# # inset map
# par(fig = c(0.1, 0.35, 0.35, 0.75), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE) # set inset figure extent
# plot(0,0, xlim=insetext[c(1,2)], ylim=insetext[c(3,4)], xaxt='n', yaxt='n', xlab='', ylab='')
# #plot(g, add=TRUE, interpolate=TRUE)
# #plot(g, add = TRUE)
# plot(inset_map, add=TRUE)
# polygon(x=rep(xlims, c(2,2)), y=c(ylims, rev(ylims)), col=NA, border='red', lwd=2)

######### Other old code
#library(grid)
#library(gridExtra)
#library(leaflet)
#library(maps)
#library(mapview)
#library(tmap)
#library(maptools)

# country_border = readOGR(here::here("Data/Map_data/phl_admbnda_adm0_psa_namria_itos_20180130", "phl_admbnda_adm0_psa_namria_itos_20180130.shp"))
#country = readOGR(here::here("Data/Regions/", "Regions.shp"))

#inset_point_coords <- data.frame(lat = ylims[1], long=xlims[1], group = 10)  # extent of inset map of Philippines 

