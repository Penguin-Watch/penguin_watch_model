#################
# Penguin Watch Model - 1 - Process data
#
# 1 - process data (krill data)
# 2 - run model
# 3 - analyze output
#################



# Clear environment -------------------------------------------------------

rm(list = ls())




# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman', repos = "http://cran.case.edu")
}
pacman::p_load(dplyr, rgdal, rgeos)





# create site buffers ----------------------------------------------------------


#created buffers around actual lat/lons - don't want low-res land mask to interfere with krill trawls

#determine which sites we have PW data for
setwd('Data/PW_data/RAW_Fiona_Apr_15_2018/')
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

#convert colony points to 3031 (rgeos expects projected spatial object)
t_col_points <- spTransform(np, CRS(proj4string(Ant)))

#create buffers around sites - 150km
all_site_buffers_150 <- rgeos::gBuffer(t_col_points, width = 150000)
all_site_buffers_100 <- rgeos::gBuffer(t_col_points, width = 100000)
all_site_buffers_50 <- rgeos::gBuffer(t_col_points, width = 50000)
all_site_buffers_25 <- rgeos::gBuffer(t_col_points, width = 25000)





# intersection of krill data with site buffers ---------------------------

#load krill data - determine which entires lie within sites we have data for

setwd('../Krill_data')

krill_data <- read.csv('krill_data_CLEAN.csv', stringsAsFactors = FALSE)
krill_df <- data.frame(LON = krill_data$lon_st, 
                     LAT = krill_data$lat_st,
                     KRILL = krill_data$krill_green_weight,
                     ID = 1:NROW(krill_data))

#transform krill data to spatial points
#points are in 4326 (uses lat/lon)
krill_points <- SpatialPoints(krill_df, proj4string = CRS("+init=epsg:4326"))
#transform to 3031 (to match everything else)
t_krill_points <- spTransform(krill_points, CRS(proj4string(Ant)))


#--------------#
#plot check

#continent
# plot(Ant)
# #site buffers
# plot(all_site_buffers, col = 'lightblue', add = TRUE)
# #krill
# points(t_krill_points, col = 'purple', pch = '.')
# #PW sites
# points(t_col_points, col = 'red', pch = '*')
#--------------#



#deterine which krill trawls fall within site we have PW data for
#to select individual sites
master_krill_150 <- data.frame()
master_krill_100 <- data.frame()
master_krill_50 <- data.frame()
master_krill_25 <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 1
  temp_site <- filter(site_ll, SITE == cam_sites[i])
  
  #convert to spatial points
  temp_site_sp <- SpatialPoints(temp_site[-1], proj4string = CRS('+init=epsg:4326'))
  
  #transform to 3031
  temp_site <- spTransform(temp_site_sp, CRS(proj4string(Ant)))
  
  #150km (150,000m) buffer
  temp_buffer_150 <- rgeos::gBuffer(temp_site, width = 150000)
  #100km
  temp_buffer_100 <- rgeos::gBuffer(temp_site, width = 100000)
  #50km
  temp_buffer_50 <- rgeos::gBuffer(temp_site, width = 50000)
  #25km
  temp_buffer_25 <- rgeos::gBuffer(temp_site, width = 25000)
  
  #which trawls started within colony i buffer
  temp_in_buff_150 <- which(as.numeric(over(t_krill_points, temp_buffer_150)) == 1)
  temp_trawls_150 <- mutate(krill_data[temp_in_buff_150,], col_id = cam_sites[i])
  master_krill_150 <- rbind(master_krill_150, temp_trawls_150)
  
  temp_in_buff_100 <- which(as.numeric(over(t_krill_points, temp_buffer_100)) == 1)
  temp_trawls_100 <- mutate(krill_data[temp_in_buff_100,], col_id = cam_sites[i])
  master_krill_100 <- rbind(master_krill_100, temp_trawls_100)
  
  temp_in_buff_50 <- which(as.numeric(over(t_krill_points, temp_buffer_50)) == 1)
  temp_trawls_50 <- mutate(krill_data[temp_in_buff_50,], col_id = cam_sites[i])
  master_krill_50 <- rbind(master_krill_50, temp_trawls_50)
  
  temp_in_buff_25 <- which(as.numeric(over(t_krill_points, temp_buffer_25)) == 1)
  temp_trawls_25 <- mutate(krill_data[temp_in_buff_25,], col_id = cam_sites[i])
  master_krill_25 <- rbind(master_krill_25, temp_trawls_25)
  
  #---------------#
  #plot everything for the single site
  
  # t_krill <- data.frame(LON = temp_trawls$lon_st, 
  #                       LAT = temp_trawls$lat_st)
  # temp_kp <- SpatialPoints(t_krill, proj4string = CRS("+init=epsg:4326"))
  # t_kp <- spTransform(temp_kp, CRS(proj4string(Ant)))
  # 
  # plot(all_site_buffers_150, col = 'lightblue')
  # plot(temp_buffer_150, add = TRUE, col = 'pink')
  # points(temp_points, pch = 19, col = 'red')
  # points(t_kp, pch = '+', col = 'green')
  #---------------#
}


#master_krill is all trawls that occured within the site buffers for the sites we have PW data for - all years (not filtered for years we have PW data for)





# thoughts ----------------------------------------------------------------

#2 ways to look at impact of krill
#-for each season (june - mid Feb [or wherever creche point is set])
#-do years with more krill fishing lead to lower breeding success at that site
#-total krill caught (average of year totals) within that site buffer
#-does this impact the colony intercept for breeding success




# Total krill catch for each site -----------------------------------------


total_krill <- c()
for (i in 1:length(cam_sites))
{
  #i <- 1
  #krill data
  temp <- filter(master_krill, col_id == cam_sites[i])
  kr_total <- sum(temp$krill_green_weight)
  tt <- data.frame(site = cam_sites[i], kr_total = kr_total)
  total_krill <- rbind(total_krill, tt)
}





# Krill catch in each season ----------------------------------------------

#determine which years match for both krill and PW data (remember that Fiona used 2007 for 2006/2007 season)


#Convert dates to POSIX - extract year, month, day, julian 
pos_dates <- as.POSIXct(master_krill$date_st, format = "%d-%B-%y")
#as.numeric(format(pos_dates, format = "%Y"))
#as.numeric(format(pos_dates, format = "%m"))
#as.numeric(format(pos_dates, format = "%d"))
#as.numeric(format(pos_dates, format = "%j"))



