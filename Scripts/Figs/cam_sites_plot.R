#plot of sites we have data for Aug 22, 2018



rm(list = ls())



# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman', repos = "http://cran.case.edu")
}

pacman::p_load(dplyr, rgdal, rgeos, ggplot2, raster, plotrix)




# cam sites and lat/lons --------------------------------------------------

if (Sys.info()[[1]] == 'Windows')
{
  setwd('C:/Users/Lynch Lab 7/Google_Drive/R/penguin_watch_model/Data/')
} else {
  setwd('~/Google_Drive/R/penguin_watch_model/Data/')
}


cams_Aug_22_2018 <- read.csv('test.csv')
SLL <- cams_Aug_22_2018[,1:2]
colnames(SLL) <- c('LAT', 'LON')


#remove name column
p_SLL <- cbind(SLL$LON, SLL$LAT)

#points are in 4326 (uses lat/lon)
col_points <- SpatialPoints(p_SLL, proj4string = CRS('+init=epsg:4326'))

#load Antarctic polygon
setwd('Coastline_medium_res_polygon/')
Ant <- rgdal::readOGR('Coastline_medium_res_polygon.shp')

#AP
setwd('../peninsula/')
AP_p <- rgdal::readOGR('GADM_peninsula.shp')
AP <- spTransform(AP_p, CRS(proj4string(Ant)))

#CCAMLR small scale management units (SSMU)
setwd('../ssmu-shapefile-WGS84/')
sm <- rgdal::readOGR('ssmu-shapefile-WGS84.shp')
SSMU <- spTransform(sm, CRS(proj4string(Ant)))
SSMU_names <- as.character(unique(SSMU@data$ShortLabel))

#exclude non-AP SSMUs
exl <- c('APPA', 'SOPA', 'SOW', 'SONE', 'SOSE', 'SGPA', 
         'SGW','SGE', 'SSPA', 'SSI')
SSMU_names <- SSMU_names[!SSMU_names %in% exl]

for (i in 1:length(SSMU_names))
{
  #i <- 1
  assign(SSMU_names[i], SSMU[which(SSMU@data$ShortLabel == SSMU_names[i]),])
}


#convert colony points to 3031 (rgeos expects projected spatial object)
t_col_points <- spTransform(col_points, CRS(proj4string(Ant)))

#create buffers around sites - 150km
all_site_buffers_150 <- rgeos::gBuffer(t_col_points, width = 150000)
all_site_buffers_100 <- rgeos::gBuffer(t_col_points, width = 100000)
all_site_buffers_50 <- rgeos::gBuffer(t_col_points, width = 50000)
all_site_buffers_25 <- rgeos::gBuffer(t_col_points, width = 25000)


#without KRILL
colonies <- data.frame(long = t_col_points@coords[,1],
                       lat = t_col_points@coords[,2])





ggplot(Ant, aes(long, lat, group = group),
       inherit.aes = FALSE,
       fill = 'grey') +
  geom_polygon() +
  # geom_path(data = Ant, aes(long, lat, group = group), 
  #           inherit.aes = FALSE,
  #           color = 'black') +
  #plot settings
  geom_tile() +
  theme_void() +
  coord_equal(xlim = c(-3800000, -2050000), 
              ylim = c(750000, 3500000)) +
  geom_point(data = colonies,
             inherit.aes = FALSE,
             size = 5,
             alpha = 0.9,
             aes(long, lat, color = 'red')) +
  #krill point outlines
  geom_point(data = colonies,
             inherit.aes = FALSE,
             size = 5,
             shape = 21,
             alpha = 0.6,
             stroke = 1,
             color = 'black',
             aes(long, lat)) +
  theme(legend.position='none')
