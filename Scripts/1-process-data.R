#################
# Penguin Watch Model - 1 - Process data
#
# 1 - process data
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


# penguin data -----------------------------------------------------------




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





#sites with processed penguin watch data:
cam_sites <- c('AITC', 'BAIL', 'BOOT', 'CUVE', 'DAMO', 'DANC', 'GEOR', 
               'HALF', 'MAIV', 'NEKO', 'PETE', 'SPIG', 'SSIS')


#to select individual sites

for (i in 1:length(cam_sites))
{
  i <- 1
  temp_site <- subset(sites, site_id == cam_sites[i])
  
  #150km (150,000m)
  temp_buffer <- rgeos::gBuffer(temp_site, width = 150000)
  #plot(site_buffer)
  #points(nsite, pch = 19, col = 'red')
  #plot(temp_buffer, add = TRUE)
  
  
}







setwd('../Krill_data/A')



