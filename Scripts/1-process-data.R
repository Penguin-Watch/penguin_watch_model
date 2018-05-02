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

pacman::p_load(dplyr, rgdal, rgeos, ggplot2)





# create site buffers ----------------------------------------------------------


#created buffers around actual lat/lons - don't want low-res land mask to interfere with krill trawls

#determine which sites we have PW data for
setwd('~/Google_Drive/R/penguin_watch_model/Data/PW_data/RAW_Fiona_Apr_15_2018/')
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
setwd('../../Coastline_medium_res_polygon/')
Ant <- rgdal::readOGR('Coastline_medium_res_polygon.shp')

setwd('../peninsula/')
AP_p <- rgdal::readOGR('GADM_peninsula.shp')
AP <- spTransform(AP_p, CRS(proj4string(Ant)))

#convert colony points to 3031 (rgeos expects projected spatial object)
t_col_points <- spTransform(col_points, CRS(proj4string(Ant)))

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
# #plot(Ant)
# plot(AP)
# #site buffers
# plot(all_site_buffers_150, col = rgb(0,0,1,0.1), add = TRUE)
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
  #i <- 2
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
  #plot everything for site i
  
  # t_krill <- data.frame(LON = temp_trawls_150$lon_st,
  #                       LAT = temp_trawls_150$lat_st)
  # temp_kp <- SpatialPoints(t_krill, proj4string = CRS("+init=epsg:4326"))
  # t_kp <- spTransform(temp_kp, CRS(proj4string(Ant)))
  # 
  # plot(AP)
  # plot(all_site_buffers_150, col = rgb(0,0,1,0.1), add = TRUE)
  # plot(temp_buffer_150, col = rgb(1,0,0,0.1), add = TRUE)
  # points(temp_site, pch = "*", col = 'red')
  # points(t_kp, pch = '.', col = 'purple')
  #---------------#
}


#master_krill is all trawls that occured within the site buffers for the sites we have PW data for - all years (not filtered for years we have PW data for)





# Effect of krill fishing during each breeding season ---------------------

#25km radius for Dec - Feb in each year

#NO FISHING EFFORT WITHIN 25KM AT SITE WHERE WE CURRENTLY HAVE DATA

kbs <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 1
  temp_krill <- filter(master_krill_25, col_id == cam_sites[i])
  temp_PW <- filter(PW_data, site == cam_sites[i])
  
  #PW year (1999/2000 season is PW year 2000)
  yrs <- unique(temp_PW$season_year)
  pos_dates <- as.Date(temp_krill$date_st, format = "%d-%B-%y")
  
  for (j in 1:length(yrs))
  {
    #j <- 1
    #used Dec 1 - Feb 1 (Feb 1 was perscribed end of PW data)
    FIRST <- as.Date(paste0(yrs[j] - 1, '-12-01'))
    LAST <- as.Date(paste0(yrs[j], '-02-01'))
    dates <- which(pos_dates > FIRST & pos_dates < LAST)
    
    if (length(dates) > 0)
    {
      #total krill caught over this period
      t_krill <- sum(temp_krill[dates,]$krill_green_weight)
      #number of trawls (effort)
      n_trawls <- length(dates)
      #krill/trawl (CPUE)
      cpue_krill <- t_krill/n_trawls
      #output
      t_out <- data.frame(SITE = cam_sites[i], 
                          YEAR = yrs[j], 
                          T_KRILL = t_krill, 
                          N_TRAWLS = n_trawls, 
                          CPUE = cpue_krill)
      #merge with final output
      kbs <- rbind(kbs, t_out)
    } else {
      t_out <- data.frame(SITE = cam_sites[i], 
                          YEAR = yrs[j], 
                          T_KRILL = 0, 
                          N_TRAWLS = 0, 
                          CPUE = NA)
      kbs <- rbind(kbs, t_out)
    }
  }
}




# Effect of krill fishing during whole season ----------------------------


#150km radius for March - Feb (March 1999 - Feb 2000 for 1999/2000 breeding season)

kws <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 2
  temp_krill <- filter(master_krill_150, col_id == cam_sites[i])
  temp_PW <- filter(PW_data, site == cam_sites[i])
  
  #PW year (1999/2000 season is PW year 2000)
  yrs <- unique(temp_PW$season_year)
  pos_dates <- as.Date(temp_krill$date_st, format = "%d-%B-%y")
  
  for (j in 1:length(yrs))
  {
    #j <- 1
    #used March 1 - Feb 1 (Feb 1 was perscribed end of PW data)
    FIRST <- as.Date(paste0(yrs[j] - 1, '-3-01'))
    LAST <- as.Date(paste0(yrs[j], '-02-01'))
    dates <- which(pos_dates > FIRST & pos_dates < LAST)
    
    if (length(dates) > 0)
    {
      #total krill caught over this period
      t_krill <- sum(temp_krill[dates,]$krill_green_weight)
      #number of trawls (effort)
      n_trawls <- length(dates)
      #krill/trawl (CPUE)
      cpue_krill <- t_krill/n_trawls
      #output
      t_out <- data.frame(SITE = cam_sites[i], 
                          YEAR = yrs[j], 
                          T_KRILL = t_krill, 
                          N_TRAWLS = n_trawls, 
                          CPUE = cpue_krill)
      #merge with final output
      kws <- rbind(kws, t_out)
    } else {
      t_out <- data.frame(SITE = cam_sites[i], 
                          YEAR = yrs[j], 
                          T_KRILL = 0, 
                          N_TRAWLS = 0, 
                          CPUE = NA)
      kws <- rbind(kbs, t_out)
    }
  }
}





# Effect of krill fishing across years ------------------------------------

#150km radius average (or total) across years (PW seasons 2011-2017) at each site

yrs <- 2011:2017

kay <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 8
  temp_krill <- filter(master_krill_150, col_id == cam_sites[i])
  temp_PW <- filter(PW_data, site == cam_sites[i])


  #PW year (1999/2000 season is PW year 2000)
  pos_dates <- as.Date(temp_krill$date_st, format = "%d-%B-%y")
  
  yr_krill <- data.frame()
  for (j in 1:length(yrs))
  {
    #j <- 3
    #used March 1 - Feb 1 (Feb 1 was perscribed end of PW data)
    FIRST <- as.Date(paste0(yrs[j] - 1, '-3-01'))
    LAST <- as.Date(paste0(yrs[j], '-02-01'))
    dates <- which(pos_dates > FIRST & pos_dates < LAST)
    
    if (length(dates) > 0)
    {
      #total krill caught over this period
      t_krill <- sum(temp_krill[dates,]$krill_green_weight)
      #number of trawls (effort)
      n_trawls <- length(dates)
      #krill/trawl (CPUE)
      cpue_krill <- t_krill/n_trawls
      
      tyk <- data.frame(SITE = cam_sites[i], 
                          T_KRILL = t_krill,
                          N_TRAWLS = n_trawls,
                          CPUE = cpue_krill)
      
      #merge with final output
      yr_krill <- rbind(yr_krill, tyk)
    } else {
      tyk <- data.frame(SITE = cam_sites[i], 
                        T_KRILL = 0,
                        N_TRAWLS = 0,
                        CPUE = NA)
      yr_krill <- rbind(yr_krill, tyk)
    }
  }
  
  OT_KRILL <- sum(yr_krill$T_KRILL)
  MN_KRILL <- mean(yr_krill$T_KRILL)
  MN_CPUE <- mean(yr_krill$CPUE, na.rm = TRUE)

  tk <- data.frame(SITE = cam_sites[i],
                   TOTAL_KRILL = OT_KRILL,
                   MN_KRILL = MN_KRILL,
                   MN_CPUE = MN_CPUE)

  kay <- rbind(kay, tk)
}





# which sites have the most fishing at them? -------------------------------


ggplot(kay, aes(SITE, TOTAL_KRILL)) +
  geom_col() +
  ylab('Total krill catch') +
  theme_bw() +
  ggtitle('Total krill catch across sites (2010-2017)')

ggplot(kay, aes(SITE, MN_KRILL)) +
  geom_col() +
  ylab('Mean krill catch') +
  theme_bw() +
  ggtitle('Mean krill catch across sites (2010-2017)')

ggplot(kay, aes(SITE, MN_CPUE)) +
  geom_col() +
  ylab('Mean Catch Per Unit Effort') +
  theme_bw() +
  ggtitle('CPUE across sites (2010-2017)')




# When is krill fishing most intense in this region? ----------------------

#AT THE SITES WE HAVE DATA FOR: March-May of each year

kdates <- as.Date(master_krill_150$date_st, format = "%d-%B-%y")
j_kdates <- as.numeric(format(kdates, '%j'))

hist(j_kdates, col = 'grey50')
rect(xleft = 0, ybottom = 0, xright = 30, ytop = 12000, col = rgb(0.8,0,0,0.2))
rect(xleft = 60, ybottom = 0, xright = 90, ytop = 12000, col = rgb(0.8,0,0,0.2))
rect(xleft = 120, ybottom = 0, xright = 150, ytop = 12000, col = rgb(0.8,0,0,0.2))
rect(xleft = 180, ybottom = 0, xright = 210, ytop = 12000, col = rgb(0.8,0,0,0.2))
rect(xleft = 240, ybottom = 0, xright = 270, ytop = 12000, col = rgb(0.8,0,0,0.2))
rect(xleft = 300, ybottom = 0, xright = 330, ytop = 12000, col = rgb(0.8,0,0,0.2))




