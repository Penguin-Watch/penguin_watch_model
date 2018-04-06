#################
# Penguin Watch Model - 1 - Process data
#
# 1 - process data (krill data)
# 1 - simulate data
# 2 - run model
# 3 - analyze output

#
#Created: Apr 5, 2018
#################



# Clear environment -------------------------------------------------------

rm(list = ls())


# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman', repos = "http://cran.case.edu")
}
pacman::p_load(dplyr, rgdal, rgeos)



# setwd -------------------------------------------------------------------

setwd('Data/Site_data/')


# create buffers ----------------------------------------------------------

#load in sites (from MAPPPD/SiteCovariates/Locations)
sites <- rgdal::readOGR('SitesEPSG3031.shp')



#----------------------------#
#check to make sure everything works properly

#units for 3031 projection are in m
#proj4string(sites)

#create 150km buffer (150,000m) around all sites
site_buffer <- rgeos::gBuffer(sites, width = 150000)
#plot(site_buffer)

#transform buffer to 4326 (uses lat/lon)
#nsb <- spTransform(site_buffer, CRS("+init=epsg:4326"))
#plot(nsb)

#test points - Admiralty Bay, Cape Bird
points <- data.frame(Longitude = c(-58.42, 166.43), 
                     Latitude = c(-62.21, -77.22))

#points are in 4326 (uses lat/lon)
np <- SpatialPoints(points, proj4string = CRS("+init=epsg:4326"))
#transform to 3031 (to match site buffer)
tnp <- spTransform(np, CRS(proj4string(site_buffer)))

#check
plot(site_buffer)
points(tnp, col = 'red', pch = 19)
#----------------------------#


setwd('../Krill_data')


#sites with processed penguin watch data:
cam_sites <- c('AITC', 'BAIL', 'BOOT', 'CUVE', 'DAMO', 'DANC', 'GEOR', 
               'HALF', 'MAIV', 'NEKO', 'PETE', 'SPIG', 'SSIS')


test <- read.csv('test.csv')

#lat/lons
pts <- data.frame(Longitude = test$startLon, 
                     Latitude = test$startLat)

#points are in 4326 (uses lat/lon)
spatial_pts <- SpatialPoints(pts, proj4string = CRS("+init=epsg:4326"))
#transform to 3031 (to match everything else)
tr_spatial_pts <- spTransform(spatial_pts, CRS(proj4string(sites)))


#to select individual sites

master_krill <- data.frame()

for (i in 1:length(cam_sites))
{
  i <- 1
  temp_site <- subset(sites, site_id == cam_sites[i])
  #transform to degrees
  #spTransform(temp_site, CRS("+init=epsg:4326"))
  
  #150km (150,000m)
  temp_buffer <- rgeos::gBuffer(temp_site, width = 150000)
  
  #which trawls started within colony i buffer
  temp_in_buff <- which(as.numeric(over(tr_spatial_pts, temp_buffer)) == 1)
  temp_trawls <- mutate(test[temp_in_buff,], col_id = cam_sites[i])
  master_krill <- rbind(master_krill, temp_trawls)

  #---------------#  
  #plot everything
  #plot(site_buffer)
  #points(temp_site, pch = 19, col = 'red')
  #plot(temp_buffer, add = TRUE)
  #points(tr_spatial_pts, pch = '+', col = 'green')
  #---------------#
}







setwd('../Krill_data/A')



