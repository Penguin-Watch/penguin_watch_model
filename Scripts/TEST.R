

# Clear environment -------------------------------------------------------

rm(list = ls())



# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman', repos = "http://cran.case.edu")
}

pacman::p_load(dplyr, rgdal, rgeos, ggplot2)





# create site buffers ----------------------------------------------------------


#created buffers around actual lat/lons - don't want low-res land mask to interfere with krill trawls

#determine which sites we have PW data for
setwd('~/Google_Drive/R/penguin_watch_model/Data/PW_data/')
PW_data <- read.csv('Markrecap_data_15.05.18.csv', stringsAsFactors = FALSE)

un_sites <- unique(PW_data$site)

site_year <- c()
for(i in 1:length(un_sites))
{
  #i <- 1
  tdat <- filter(PW_data, site == un_sites[i])
  tsy <- unique(tdat$season_year)
  tout <- data.frame(SITE = rep(un_sites[i], length(tsy)), YEAR = tsy)
  site_year <- rbind(site_year, tout)
}

#site/years we have data for:
#site_year


#BOOT is now PCHA in MAPPPD database
pos <- which(un_sites == 'BOOT')
cam_sites_p <- un_sites[-pos]
cam_sites <- c(cam_sites_p, 'PCHA')

PW_data$site[which(PW_data$site == 'BOOT')] <- 'PCHA'


#determine lat/lons for all site we have data for
site_ll <- data.frame(SITE = cam_sites, LON = rep(NA, length(cam_sites)), LAT = rep(NA, length(cam_sites)))
for (i in 1:length(cam_sites))
{
  #i <- 1
  temp <- filter(PW_data, site == cam_sites[i])[1,]
  site_ll$LON[i] <- temp$col_lon
  site_ll$LAT[i] <- temp$col_lat
}


#remove name column
p_site_ll <- site_ll[,-1]

#points are in 4326 (uses lat/lon)
col_points <- SpatialPoints(p_site_ll, proj4string = CRS('+init=epsg:4326'))

#load Antarctic polygon
setwd('../Coastline_medium_res_polygon/')
Ant <- rgdal::readOGR('Coastline_medium_res_polygon.shp')

#AP
setwd('../peninsula/')
AP_p <- rgdal::readOGR('GADM_peninsula.shp')
AP <- spTransform(AP_p, CRS(proj4string(Ant)))

#CCAMLR zones
setwd('../asd-shapefile-WGS84/')
mz <- rgdal::readOGR('asd-shapefile-WGS84.shp')
CCAMLR_zones <- spTransform(mz, CRS(proj4string(Ant)))
sub_481 <- CCAMLR_zones[which(CCAMLR_zones@data$Name == 'Subarea 48.1'),]
#small scale management units (SSMU)
setwd('../ssmu-shapefile-WGS84/')
sm <- rgdal::readOGR('ssmu-shapefile-WGS84.shp')
SSMU <- spTransform(sm, CRS(proj4string(Ant)))
SSMU_names <- as.character(unique(SSMU@data$ShortLabel))
#use names from 'CCAMLR/CCAMLR_2017_Statistical_Bulletin_Volume_29_Data_Files/ReferenceDataGeographicArea.csv

# SSMU_names <- c('APPA', 'APW', 'APDPW',
#                 'APDPE', 'APBSW', 'APBSE',
#                 'APEI', 'SOPA', 'SOW',
#                 'SONE', 'SOSE', 'SGPA',
#                 'SGW', 'SGE', 'SSPA',
#                 'SS', 'APE')
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
all_site_buffers_25 <- rgeos::gBuffer(t_col_points, width = 1)