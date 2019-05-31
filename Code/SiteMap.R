# Make a site map

#################### Set-up: ####################
# Load relevant libraries
#library(grid)
#library(gridExtra)
#library(leaflet)
#library(maps)
library(mapdata)
#library(mapview)
library(rgdal)
#library(tmap)
#library(maptools)
library(cowplot)
library(raster)
library(magick)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(here)

##### Load anem lat/lons - REDO SO GET THE ANEM INFO FROM DATABASE + A SCRIPT?
#load(here::here("Data", "AnemAllInfowLatLon2.RData"))
#load(here::here("Data/Script_outputs", "anems_Processed.RData"))
# source(here::here('Code', 'Constants_database_common_functions.R'))
#data = anems.Processed2
data_anems = as.data.frame(read.csv(here::here("Data", "GPSSurvey.anemlatlong2015-12-16.csv"), row.names=1, stringsAsFactors = FALSE))

##### Load shapefiles for hulls and coast (downloaded from amphiprion, hulls created by Mario in 2016)
# Site hulls
patchfiles = list.files(here::here("Data/Map_data/Site_hulls/"), pattern = ".shp")  # get list of files with site hulls
patches = vector("list", length(patchfiles))
for(i in 1:length(patches)){
  temp = readOGR(paste(here::here("Data/Map_data/Site_hulls/"), patchfiles[i], sep=''))
  patches[[i]] = spTransform(temp, CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'))
  #patches[[i]] = temp
}

# Coastline
coast = readOGR(here::here("Data/Map_data/", "PHL_adm3_leyte_studyarea_trimcoastline.shp"))

# Philippines
# country_border = readOGR(here::here("Data/Map_data/phl_admbnda_adm0_psa_namria_itos_20180130", "phl_admbnda_adm0_psa_namria_itos_20180130.shp"))
country = readOGR(here::here("Data/Regions/", "Regions.shp"))
g = map("world", fill = TRUE, col = "grey")


##### Set constants
kmlen = 1/111.12 # length of a km, in °N

#################### Make sub-plots and overall plot: ####################
# Limits for maps
xlims = c(124.64, 124.80) # main map
ylims = c(10.63, 10.87)
xlims_Alb = c(124.708, 124.718) # Albuera patch box on main map
ylims_Alb = c(10.867, 10.875)
#xlims_Mar = c(124.78, 124.79) # Marcos patch box on main map
#ylims_Mar = c(10.752, 10.762)

# Set colors 
paste(sort(unique(data$Spp)), collapse=', ')
colslist = brewer.pal(3, name="Dark2")
colsnames = c('Open anemone', 'Anemone with A. clarkii')
cols = colslist[c(1,2,3,3,3,3,3,3)][as.numeric(data$Spp)] # treat all other species as 'other'
patchfadecol = rgb(t(col2rgb(colslist[3])), maxColorValue=255, alpha=100)
inds = data$Spp %in% c('', 'APCL')

insetext = c(118,127,6,18) # extent of inset map


# Philippines map inset subplot
g <- plot(country, col = "grey", border = NA)

# Clownfish photo subplot
photo_plot <- ggdraw() +
  draw_image(here::here("Data", "Fish.jpg"), scale=0.9)

g <- ggplot(data = country, aes(x = long, y = lat, group = group)) +
  geom_polygon(colour = "grey", fill = "grey")


test <- ggplot(data = coast, aes(x = long, y = lat, group = group)) +
  geom_polygon(colour = "grey", fill = "grey") +
  geom_polygon(data = patches[[1]], aes(x = long, y = lat, group = group), colour = "blue", fill = "blue")


geom_polygon(data = borders, aes(x = long, y = lat, group = group), colour = "black", fill = "white")
# Map of anemones
xlims = c(124.64, 124.80) # main map
ylims = c(10.63, 10.87)
xlims2 = c(124.775, 124.778) # Punta Baybayon
ylims2 = c(10.633, 10.63555)
xlims3 = c(124.78, 124.79) # Marcos
ylims3 = c(10.752, 10.762)
xlims4 = c(124.705, 124.73) # Albuera
ylims4 = c(10.85, 10.876)



# set up colors
paste(sort(unique(data$Spp)), collapse=', ')
colslist = brewer.pal(3, name="Dark2")
colsnames = c('Open anemone', 'Anemone with A. clarkii')
cols = colslist[c(1,2,3,3,3,3,3,3)][as.numeric(data$Spp)] # treat all other species as 'other'
patchfadecol = rgb(t(col2rgb(colslist[3])), maxColorValue=255, alpha=100)
inds = data$Spp %in% c('', 'APCL')

data_anems <- data_anems %>%
  filter(Name %in% c("Palanas", "Wangag")) %>%
  filter(Spp %in% c("","APCL")) %>%
  mutate(species = case_when(Spp == "" ~ "empty",
                               Spp == "APCL" ~ "APCL")) %>%
  mutate(group = 9)

# get google map
g <- gmap(insetext, type='roadmap', lonlat=TRUE, scale=2)

# map	
quartz(height=5, width=6.5)
# pdf(height=5, width=6, file='Figures/anem_map2.pdf')
par(mai=c(0.4,0.5,0.1, 0.2), las=1, cex.axis=0.8, mgp=c(2.5,0.5,0), tcl=-0.3)
layout(matrix(c(1,1,1,2,3,4), nrow=3), widths=c(3,2))
plot(0,0, xlim=xlims, ylim=ylims, xlab='', ylab='')
#map('worldHires', 'Philippines', add=TRUE, fill=TRUE, col='grey', border=NA)
plot(coast, col='grey', border=NA, add=TRUE)
polygon(x=c(124.76, 124.81, 124.81), y=c(10.6, 10.65, 10.6), col='grey', border=NA) # fill in a missing corner
#plot(patches[[1]], add = TRUE)
for(i in 1:length(patches)) plot(patches[[i]], add=TRUE, col=colslist[3], border=colslist[3], lwd=2)
#for(i in 1:5) plot(patches[[i]], add=TRUE, col=colslist[3], border=colslist[3], lwd=2)
#points(data$lon, data$lat, pch=16, cex=0.4, col=cols)
lines(x = rep(xlims[1],2), y=c(ylims[1], ylims[1]+5*kmlen), lwd=5, lend=3)
text(labels="5 km", x=xlims[1]+0.015, y=10.65, cex=1.5)
polygon(c(xlims2, rev(xlims2)), rep(ylims2, c(2,2)), border='red', lwd=1.5)
polygon(c(xlims3, rev(xlims3)), rep(ylims3, c(2,2)), border='red', lwd=1.5)
polygon(c(xlims4, rev(xlims4)), rep(ylims4, c(2,2)), border='red', lwd=1.5)
mtext(side=1, text='Longitude (°E)', line=1.8, cex=0.8)
mtext(side=2, text='Latitude (°N)', line=2.5, cex=0.8, las=0)
legend('topleft', legend='Habitat patch', lwd=3, col=colslist[3], bty='n')

# Albuera
par(mai=c(0.3,0.4,0.1, 0.1), las=1)
plot(0,0, xlim=xlims4, ylim=ylims4, xlab='', ylab='')
for(i in 1:length(patches)) plot(patches[[i]], add=TRUE, col=patchfadecol, border=patchfadecol, lwd=7)
points(data$lon[inds], data$lat[inds], pch=16, cex=0.8, col=cols[inds])
#map('worldHires', 'Philippines', add=TRUE, fill=TRUE, col='grey',border=NA)
plot(coast, col='grey', border=NA, add=TRUE)
lines(x = rep(xlims4[1],2), y=rep(10.85,2)+c(0,0.5*kmlen), lwd=3, lend=3)
text(labels="500 m", x=124.708, y=10.852, cex=1)

# Marcos
plot(0,0, xlim=xlims3, ylim=ylims3, xlab='', ylab='')
for(i in 1:length(patches)) plot(patches[[i]], add=TRUE, col=patchfadecol, border=NA, lwd=2)
points(data$lon[inds], data$lat[inds], pch=16, cex=0.8, col=cols[inds])
#map('worldHires', 'Philippines', add=TRUE, fill=TRUE, col='grey',border=NA)
plot(coast, col='grey', border=NA, add=TRUE)
lines(x = rep(xlims3[1],2), y=rep(ylims3[1],2)+c(0,0.1*kmlen), lwd=3, lend=3)
text(labels="100 m", x=xlims3[1]+0.001, y=ylims3[1]+0.0005, cex=1)
legend('topleft', legend=c(colsnames, 'Habitat patch'), col=c(colslist[c(1,2)], patchfadecol), pch=c(16, 16, NA), lwd = c(NA, NA, 4), cex=1, text.font=1, bty='o', bg='white', box.col='black')

# Punta
plot(0,0, xlim=xlims2, ylim=ylims2, xlab='', ylab='')
for(i in 1:length(patches)) plot(patches[[i]], add=TRUE, col=patchfadecol, border=NA)
points(data$lon[inds], data$lat[inds], pch=16, cex=0.8, col=cols[inds])
#map('worldHires', 'Philippines', add=TRUE, fill=TRUE, col='grey', border=NA)
plot(coast, col='grey', border=NA, add=TRUE)
lines(x = rep(124.775,2), y=rep(ylims2[1],2)+c(0,0.1*kmlen), lwd=3, lend=3)
text(labels="100 m", x=124.7753, y=ylims2[1]+0.0005, cex=1)

# inset map
par(fig = c(0.1, 0.35, 0.35, 0.75), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE) # set inset figure extent
plot(0,0, xlim=insetext[c(1,2)], ylim=insetext[c(3,4)], xaxt='n', yaxt='n', xlab='', ylab='')
#plot(g, add=TRUE, interpolate=TRUE)
#plot(g, add = TRUE)
plot(inset_map, add=TRUE)
polygon(x=rep(xlims, c(2,2)), y=c(ylims, rev(ylims)), col=NA, border='red', lwd=2)

plot(patches[[1]], col ="blue")

#################### Trying one with ggplot ####################
data = read.csv(here::here("Data", "GPSSurvey.anemlatlong2015-12-16.csv"), row.names=1)

##### Load shapefiles for hulls and coast (downloaded from amphiprion, hulls created by Mario in 2016)
# Site hulls
patchfiles = list.files(here::here("Data/Map_data/Site_hulls/"), pattern = ".shp")  # get list of files with site hulls
patches = vector("list", length(patchfiles))
for(i in 1:length(patches)){
  temp = readOGR(paste(here::here("Data/Map_data/Site_hulls/"), patchfiles[i], sep=''))
  patches[[i]] = spTransform(temp, CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'))
  #patches[[i]] = temp
}

# Coastline
coast = readOGR(here::here("Data/Map_data/", "PHL_adm3_leyte_studyarea_trimcoastline.shp"))

# Philippines
# country_border = readOGR(here::here("Data/Map_data/phl_admbnda_adm0_psa_namria_itos_20180130", "phl_admbnda_adm0_psa_namria_itos_20180130.shp"))
country = readOGR(here::here("Data/Regions/", "Regions.shp"))


##### Set coordinates
inset_point_coords <- data.frame(lat=ylims[1], long=xlims[1], group =10)

# Inset map of Philippines
inset_map <- ggplot(data =  country, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = c(118,127), ylim = c(6,18), 1.3) +
  geom_polygon(colour = "black", fill = "light grey") +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        axis.line=element_blank()) +
  geom_point(data = inset_point_coords, aes(x = long, y = lat, group=group), colour="red", cex = 3) +
  annotate(geom = 'text', x = 120.5, y = 8.3, label = 'Philippines', cex=3) 

# Site area map
corner_coords = data.frame(long = c(124.76, 124.81, 124.81), lat = c(10.6, 10.65, 10.6), group = 5)
km_line_coords = data.frame(long = rep(xlims[1],2), lat=c(ylims[1], ylims[1]+5*kmlen), group = 6)
red_box_coords = data.frame(long = c(xlims_Alb[1], xlims_Alb[2], xlims_Alb[2], xlims_Alb[1]), 
                            lat= c(ylims_Alb[1], ylims_Alb[1], ylims_Alb[2], ylims_Alb[2]), group = 7)
hab_line_coords = data.frame(long = rep(xlims[1], 2), lat = c(10.85, 10.90), group = 8)
  
site_area <- ggplot(data = coast, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims, ylim = ylims, 1.3) +
  geom_polygon(colour = "grey", fill = "grey") +
  geom_polygon(data = corner_coords, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey") +
  geom_polygon(data = patches[[1]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[2]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[3]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[4]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[5]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[6]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[7]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[8]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[9]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[10]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[11]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[12]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[13]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[14]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[15]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_polygon(data = patches[[16]], aes(x = long, y = lat, group = group), fill = "blue", colour = "blue") +
  geom_line(data = km_line_coords, aes(x = long, y = lat, group = group), color = "black", lwd = 2) +
  annotate(geom = 'text', x = xlims[1]+0.025, y = 10.65, label = '5 km', size = 4) +
  xlab("Longitude (°E)") + ylab('Latitude (°N)') +
  geom_polygon(data = red_box_coords, aes(x = long, y = lat, group = group), fill = NA, color = "red", lwd = 1) +
  geom_line(data = hab_line_coords, aes(x = long, y = lat, group = group), color = "blue", lwd = 2) +
  annotate(geom = "text", x = xlims[1]+0.033, y = 10.865, label = "Habitat \n patches", size = 4)

# Example habitat patch - patch 7 is Palanas, patch 16 is Wangag
#Palanas_patch <- patchfiles[which("Palanas Hull.shp")]
Palanas_patch <- 7
Wangag_patch <- 16
red_box_coords_zoomed <- data.frame(long = c(124.7085, 124.7174, 124.7174, 124.7085),
                                    lat = c(10.867, 10.867, 10.875, 10.875),
                                    group = 7)
Albuera_patch <- ggplot(data = coast, aes(x = long, y = lat, group = group)) +
  coord_fixed(xlim = xlims_Alb, ylim = ylims_Alb, 1.3) +
  geom_polygon(color = "grey", fill = "grey") +
  geom_polygon(data = patches[[Palanas_patch]], aes(x = long, y = lat, group = group), color = "blue", fill = "blue", alpha = 0.2) +  # Palanas
  geom_polygon(data = patches[[Wangag_patch]], aes(x = long, y = lat, group = group), color = "blue", fill = "blue", alpha = 0.2) +  # Wangag
  geom_polygon(data = red_box_coords_zoomed, aes(x = long, y = lat, group = group), fill = NA, color = "red", lwd = 1) +
  geom_point(data = data_anems, aes(x = lon, y = lat, group = group, color = species), alpha=0.7) +
  scale_color_manual(values = colslist[1:2]) +
  theme(legend.position = c(0.2, 0.3)) +
  #theme(panel.background = element_rect(colour = "red", fill = NA)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") + ylab("")
  
# Clownfish picture
photo_plot <- ggdraw() +
  draw_image(here::here("Data", "Fish.jpg"))
  

##### Arrange all together
site_with_inset <- ggdraw() +
  draw_plot(site_area, 0,0,0.5,1) +
  draw_plot(inset_map, 0.02,0.3,0.55,0.5)

site_with_inset <- ggdraw() +
  draw_plot(site_area + theme(legend.justification = "bottom"), 0,0,1,1) +
  draw_plot(inset_map + theme(legend.justification = "top"), 0.5, 0.52, 0.5, 0.5)

site_with_inset <- ggdraw(site_area) + 
  draw_plot(inset_map + theme(legend.justification = "top"), 0.3, 0.3, 0.5, 0.5)


ggdraw() +
  draw_plot(plot.diamonds + theme(legend.justification = "bottom"), 0, 0, 1, 1) +
  draw_plot(plot.mpg + scale_color_viridis(discrete = TRUE) + 
              theme(legend.justification = "top"), 0.5, 0.52, 0.5, 0.4) +
  draw_plot_label(c("A", "B"), c(0, 0.5), c(1, 0.92), size = 15)

right_column <- plot_grid(Albuera_patch, photo_plot, labels = c("b", "c"), ncol = 1, align = "v", rel_heights = c(1.5, 1))

pdf(file=here::here("Plots/FigureDrafts", "Map_and_photo.pdf"))
plot_grid(site_with_inset, right_column, labels = c("a", ""), ncol = 2)
dev.off()

# bottom_row <- plot_grid(plot.mpg, plot.diamonds, labels = c('B', 'C'), align = 'h', rel_widths = c(1, 1.3))
# plot_grid(plot.iris, bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1, 1.2))



p <- qplot(1:10, 1:10)
# draw into the top-right corner of a larger plot area
ggdraw() + draw_plot(p, .6, .6, .4, .4)

# Albuera
par(mai=c(0.3,0.4,0.1, 0.1), las=1)
plot(0,0, xlim=xlims4, ylim=ylims4, xlab='', ylab='')
for(i in 1:length(patches)) plot(patches[[i]], add=TRUE, col=patchfadecol, border=patchfadecol, lwd=7)
points(data$lon[inds], data$lat[inds], pch=16, cex=0.8, col=cols[inds])
#map('worldHires', 'Philippines', add=TRUE, fill=TRUE, col='grey',border=NA)
plot(coast, col='grey', border=NA, add=TRUE)
lines(x = rep(xlims4[1],2), y=rep(10.85,2)+c(0,0.5*kmlen), lwd=3, lend=3)
text(labels="500 m", x=124.708, y=10.852, cex=1)

lines(x = rep(xlims[1],2), y=c(ylims[1], ylims[1]+5*kmlen), lwd=5, lend=3)
text(labels="5 km", x=xlims[1]+0.015, y=10.65, cex=1.5) 
  
ggplot(data = patches[[1]], aes(x=long, y=lat, group=group)) +
  geom_polygon(colour="blue", fill="blue")





par(mai=c(0.4,0.5,0.1, 0.2), las=1, cex.axis=0.8, mgp=c(2.5,0.5,0), tcl=-0.3)
layout(matrix(c(1,1,1,2,3,4), nrow=3), widths=c(3,2))
plot(0,0, xlim=xlims, ylim=ylims, xlab='', ylab='')
#map('worldHires', 'Philippines', add=TRUE, fill=TRUE, col='grey', border=NA)
plot(coast, col='grey', border=NA, add=TRUE)
polygon(x=c(124.76, 124.81, 124.81), y=c(10.6, 10.65, 10.6), col='grey', border=NA) # fill in a missing corner
#plot(patches[[1]], add = TRUE)
for(i in 1:length(patches)) plot(patches[[i]], add=TRUE, col=colslist[3], border=colslist[3], lwd=2)
#for(i in 1:5) plot(patches[[i]], add=TRUE, col=colslist[3], border=colslist[3], lwd=2)
#points(data$lon, data$lat, pch=16, cex=0.4, col=cols)
lines(x = rep(xlims[1],2), y=c(ylims[1], ylims[1]+5*kmlen), lwd=5, lend=3)
text(labels="5 km", x=xlims[1]+0.015, y=10.65, cex=1.5)
polygon(c(xlims2, rev(xlims2)), rep(ylims2, c(2,2)), border='red', lwd=1.5)
polygon(c(xlims3, rev(xlims3)), rep(ylims3, c(2,2)), border='red', lwd=1.5)
polygon(c(xlims4, rev(xlims4)), rep(ylims4, c(2,2)), border='red', lwd=1.5)
mtext(side=1, text='Longitude (°E)', line=1.8, cex=0.8)
mtext(side=2, text='Latitude (°N)', line=2.5, cex=0.8, las=0)
legend('topleft', legend='Habitat patch', lwd=3, col=colslist[3], bty='n')

inset_coords <- data.frame(lat = ylims, long = xlims, group = 10)

inset_coords <- data.frame(lat = xlims[1], x2 = xlims[2], y1 = ylims[1], y2 = ylims[2], group=1)
xlims = c(124.64, 124.80) # main map
ylims = c(10.63, 10.87)

eom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.5)

pdf(file=here::here('Manuscript','Fig1.pdf')) 
ggplot(data = land, aes(x = long, y = lat, group = group)) + 
  coord_fixed(xlim = c(-125,-115), ylim = c(33.5, 42), 1.3) +
  geom_polygon(colour = "black", fill = "white") +
  #geom_polygon(data = land, aes(x = long, y = lat, group = group), colour = "black", fill = "white") +
  geom_polygon(data = borders, aes(x = long, y = lat, group = group), colour = "black", fill = "white") +
  geom_line(data = rivers_3, aes(x = long, y = lat, group = group), colour = "gray") +
  geom_line(data = rivers_2, aes(x = long, y = lat, group = group), colour = "gray") +
  geom_line(data = rivers, aes(x = long, y = lat, group = group), colour = "gray") +
  geom_point(data = hatch_map_df %>% filter(group %in% c(11,12,13,14,15,16,17)), mapping = aes(color = Hatchery, shape = Hatchery), stroke = 2.5, size = 4) +
  scale_shape_manual(values = hatchery_shapes) +
  scale_color_manual(values = hatchery_palette) +
  geom_point(data = hatch_map_df %>% filter(group %in% c(21, 22)), fill = 'black', size = 3) +
  annotate(geom = 'text', x = -123.7034, y = 37.6989, label = 'San \n Francisco', size = 4) +
  annotate(geom = 'text', x = -116.6010, y = 33.9734, label = 'Los Angeles', size = 4) +
  xlab("longitude") + ylab("latitude") +
  #geom_polygon(data = ca_df, aes(x = long, y = lat, group = group), colour = "black", fill = "white") +
  theme_bw() +
  theme(text = element_text(size=15))
dev.off()


#load(here::here('Data', 'Leyte.7z'))

# devtools::install_github('pinskylab/clownfish')  # try installing clownfish package...
# library(clownfish)



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
