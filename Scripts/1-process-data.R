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

setwd('~/Google_Drive/R/penguin_watch_model/Data/Site_data/')


# create site buffers ----------------------------------------------------------

#load in sites (from MAPPPD/SiteCovariates/Locations)
sites <- rgdal::readOGR('SitesEPSG3031.shp')
data.frame(sites$site_id, sites$site_name)


#----------------------------#
#check to make sure everything works properly

#units for 3031 projection are in m
#proj4string(sites)

#create 150km buffer (150,000m) around all sites
#site_buffer <- rgeos::gBuffer(sites, width = 150000)
#plot(site_buffer)

#transform buffer to 4326 (uses lat/lon)
#nsb <- spTransform(site_buffer, CRS("+init=epsg:4326"))
#plot(nsb)

#test points - Admiralty Bay, Cape Bird
#points <- data.frame(Longitude = c(-58.42, 166.43), 
#                     Latitude = c(-62.21, -77.22))

#points are in 4326 (uses lat/lon)
#np <- SpatialPoints(points, proj4string = CRS("+init=epsg:4326"))
#transform to 3031 (to match site buffer)
#tnp <- spTransform(np, CRS(proj4string(site_buffer)))

#check
#plot(site_buffer)
#points(tnp, col = 'red', pch = 19)
#----------------------------#

#create buffers around sites
site_buffer <- rgeos::gBuffer(sites, width = 150000)


#determine which sites we have PW data for
setwd('../PW_data/RAW_Fiona_Apr_15_2018/')
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


#sites with PW data that Fiona sent over:
#BOOT is now PCHA in MAPPPD database
pos <- which(un_sites == 'BOOT')
cam_sites_p <- un_sites[-pos]
cam_sites <- c(cam_sites_p, 'PCHA')

#sites Fiona initally said we have PW data for:
#cam_sites <- c('AITC', 'BAIL', 'BOOT', 'CUVE', 'DAMO', 'DANC', 'GEOR', 
#'HALF', 'MAIV', 'NEKO', 'PETE', 'SPIG', 'SSIS')

#site/years we have data for (only 1 or two years of data for each site as of Apri 16, 2018)
#site_year




# intersection of krill data with PW data sites ---------------------------

#load krill data - determine which entires lie within sites we have data for

setwd('../../Krill_data')

krill_data <- read.csv('krill_data_CLEAN.csv')
points <- data.frame(longitude = krill_data$lon_st, 
                     latitude = krill_data$lat_st,
                     krill = krill_data$krill_green_weight,
                     id = 1:NROW(krill_data))

#plot points on map
#points are in 4326 (uses lat/lon)
np <- SpatialPoints(points, proj4string = CRS("+init=epsg:4326"))
#transform to 3031 (to match everything else)
tnp <- spTransform(np, CRS(proj4string(site_buffer)))

#plot(site_buffer)
#points(tnp, col = 'red', pch = '.')




#deterine which krill trawls fall within site we have PW data for
#to select individual sites
master_krill <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 3
  temp_site <- subset(sites, site_id == cam_sites[i])
  #transform to degrees
  #spTransform(temp_site, CRS("+init=epsg:4326"))
  
  #150km (150,000m) buffer
  temp_buffer <- rgeos::gBuffer(temp_site, width = 150000)
  
  #which trawls started within colony i buffer
  temp_in_buff <- which(as.numeric(over(tnp, temp_buffer)) == 1)
  temp_trawls <- mutate(krill_data[temp_in_buff,], col_id = cam_sites[i])
  master_krill <- rbind(master_krill, temp_trawls)
  
  #---------------#  
  #plot everything
  #plot(site_buffer)
  #points(temp_site, pch = 19, col = 'red')
  #plot(temp_buffer, add = TRUE)
  #points(tnp, pch = '+', col = 'green')
  #---------------#
}

#master_krill is all trawls that occured within the site buffers for the sites we have PW data for - all years (not filtered for years we have PW data for)



#determine which years match for both krill and PW data (need to determine how to make the cutoff between years [and remember that Fiona used 2007 for 2006/2007 season - double check this])
str(master_krill)

for (i in 1:length(cam_sites))
{
  i <- 1
  temp <- filter(master_krill, col_id == cam_sites[i])
  str(temp)
  
  site_year  
}





# thoughts ----------------------------------------------------------------


#-might want to look at over trawl density (average of where most krill is caught spatially) - could be covariate in model
#-additionally/alternatively might want to look at temporal component of krill catch







# plot all krill data -----------------------------------------------------

#need to find better way to visualize krill catch

require(ggplot2)
#credit to C. Che-Castaldo for mapping code
world <- map_data('world')
worldmap <- ggplot(world, aes(x = long, y = lat, group = group)) + 
  geom_polygon(color = "black", fill = "white") + 
  geom_path() + 
  scale_y_continuous(breaks = (-2:2) * 30) + 
  scale_x_continuous(breaks = (-4:4) * 45) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank())

#90,0,0 for south pole
(antarctica <- worldmap + coord_map("ortho", orientation = c(-58, -50, 0), ylim = c(-60, -50)))


antarctica + 
  geom_point(data = points,
             aes(x = longitude,
                 y = latitude,
                 group = id,
                 size = krill),
             color = 'red',
             alpha = 0.3)

antarctica + 
  stat_density2d(data = points,
             aes(x = longitude,
                 y = latitude,
                 group = id),
                 geom = 'polygon')


require(ggmap)

#AP
map <- get_map(location = c(lon = -59, lat = -63), zoom = 6)
#AP and SG
map <- get_map(location = c(lon = -52, lat = -60), zoom = 4)
ggmap(map) + 
  geom_density2d(data = points, 
                 aes(x = longitude,
                     y = latitude)
                 ) +
  stat_density2d(data = points,
                 aes(x = longitude,
                     y = latitude,
                     fill = ..level.., alpha = ..level..),
                  bins = 20, geom = 'polygon') +
  scale_fill_gradient(low = 'green', high = 'red') +
  scale_alpha(range(0, 0.3), guide = FALSE)




#plot colony locations on map



setwd('../Krill_data/A')



