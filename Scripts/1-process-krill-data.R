#################
# Penguin Watch Model - 1 - Process krill data
#
# 1-process-krill-data.R | process krill data
# 2-process-SIC-data.R | process SIC data
# 3-model.R | penguin model
# 3-run-model.pbs | pbs script to run penguin model on HPC resources
# 4-analyze-output.R | analyze model output
#
# Author: Casey Youngflesh
#################



# Krill covariates -------------------------------------------------------------------

#Several ways in which krill fishing may impact BS:

# * Effect of krill fishing during breeding season
# * Effect of krill fishing during entire season
# * Effect of krill fishing across years



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
SSMU_names_p <- unique(SSMU@data$Name)
SSMU_names <- make.names(SSMU_names_p)
#use names in Table A2.1 from 'CCAMLR_krill_report.pdf' - no data for 48.4 is that table
SSMU_names <- c('APPA_481', 'APW_481', 'APDPW_481',
                'APDPE_481', 'APBSW_481', 'APBSE_481',
                'APEI_481', 'SOPA_482', 'SOW_482',
                'SONE_482', 'SOSE_482', 'SGPA_483',
                'SGW_483', 'SGE_483', 'SSPA_484',
                'SS_484', 'APE_481')
for (i in 1:length(SSMU_names))
{
  assign(SSMU_names[i], SSMU[which(SSMU@data$Name == SSMU_names_p[i]),])
}



#convert colony points to 3031 (rgeos expects projected spatial object)
t_col_points <- spTransform(col_points, CRS(proj4string(Ant)))

#create buffers around sites - 150km
all_site_buffers_150 <- rgeos::gBuffer(t_col_points, width = 150000)
all_site_buffers_100 <- rgeos::gBuffer(t_col_points, width = 100000)
all_site_buffers_50 <- rgeos::gBuffer(t_col_points, width = 50000)
all_site_buffers_25 <- rgeos::gBuffer(t_col_points, width = 25000)






# intersection of krill data with site buffers ---------------------------

#load krill data - determine which entires lie within sites we have data for

setwd('../Krill_data/Processed/')

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


#function takes spatial object (3031) as input and determine which krill trawl points are located within that shp file
#if TYPE = 'points' -> spatial points that fall within spatial object are returned
#if TYPE = 'data' -> rows of krill data that fall within spatial object are returned
k_sp_fun <- function(INPUT, LABEL, TYPE)
{
  #INPUT <- get(SSMU_names[15])
  #LABEL <- SSMU_names[15]
  tind <- which(as.numeric(over(t_krill_points, geometry(INPUT))) == 1)
  k_temp <- mutate(krill_data[tind,], CCAMLR_Region = paste0(LABEL))
  k_ll <- data.frame(LON = k_temp$lon_st,
                     LAT = k_temp$lat_st)
  if (NROW(k_ll) > 0)
  {
    k_pre <- SpatialPoints(k_ll, proj4string = CRS("+init=epsg:4326"))
    k_pts <- spTransform(k_pre, CRS(proj4string(Ant)))
  } else {
    k_pre <- NA
    k_pts <- NA
  }
  
  if (TYPE == 'points')
  {
    return(k_pts)
  }
  if (TYPE == 'data')
  {
    return(k_temp)
  }
}



#determine which krill trawls fall within:
#CCAMLR region 48.1
k481_pts <- k_sp_fun(sub_481, '48.1', TYPE = 'points')

#CCAMLE SSMUs
for (i in 1:length(SSMU_names))
{
  #i <- 15
  assign(paste0(SSMU_names[i], '_pts'), k_sp_fun(get(SSMU_names[i]), LABEL = SSMU_names[i], TYPE = 'points'))
  assign(paste0(SSMU_names[i], '_krill'), k_sp_fun(get(SSMU_names[i]), LABEL = SSMU_names[i], TYPE = 'data'))
}


#--------------#
#plot check

#plot(Ant)
# plot(AP)
# #plot CCAMLR SSMUs
# gg_color_hue <- function(n, ALPHA = 1)
# {
#   hues = seq(15, 375, length=n+1)
#   hcl(h=hues, l=65, c=100, alpha = ALPHA)[1:n]
# }
# cols <- gg_color_hue(length(SSMU_names), ALPHA = 0.2)
# 
# for (i in 1:length(SSMU_names))
# {
#   #i <- 1
#   plot(get(SSMU_names[i]), add = TRUE, col = cols[i])
# }
# #site buffers
# #plot(all_site_buffers_150, col = rgb(0,0,1,0.5), add = TRUE)
# #PW sites
# points(t_col_points, col = rgb(0,0,1,0.5), pch = 19)
# #krill
# #points(t_krill_points, col = rgb(1,0,0,0.05), pch = '.')
# #CCAMLR zone 48.1
# #plot(sub_481, add = TRUE)
# 
# #overlay krill trawls:
# #just within 48.1
# #points(k481_pts, col = rgb(0,1,0,0.05), pch = '.')
# #SSMUs
# cols <- gg_color_hue(length(SSMU_names), ALPHA = 0.3)
# for (i in 1:length(SSMU_names))
# {
#   #i <- 3
#   points(get(paste0(SSMU_names[i], '_pts')), col = cols[i], pch = '.')
# }
#--------------#


#PW year is year t-1/t season (year 2000 is 1999/2000 season)
#CCCAMLR year is t/t+1 Dec 1 - Nov 30 (year 1999 is Dec 1, 1999 - Nov 30, 2000)


CCAMLR_krill <- read.csv('krill_table_A2_1.csv')
cn_CCAMLR <- colnames(CCAMLR_krill)
yrs <- 2011:2016
SSMU_out <- data.frame()
for (i in 1:length(SSMU_names))
{
  #i <- 1
  temp_data <- get(paste0(SSMU_names[i], '_krill'))
  pos_dates <- as.Date(temp_data$date_st, format = "%d-%B-%y")
  
  cn <- which(cn_CCAMLR == SSMU_names[i])
  zone_CCAMLR <- CCAMLR_krill[,c(1,cn)]
  
  for (j in 1:length(yrs))
  {
    #j <- 2
    if (length(cn) > 0)
    {
      zcy <- filter(zone_CCAMLR, Season == yrs[j])[,2]
    } else {
      zcy <- NA
    }

    FIRST <- as.Date(paste0(yrs[j], '-12-01'))
    LAST <- as.Date(paste0(yrs[j]+1, '-11-30'))
    dates <- which(pos_dates > FIRST & pos_dates < LAST)
    
    if (length(dates) > 0)
    {
      #total krill caught over this period - IN TONNES (AKER IN KG)
      dd <- temp_data[dates,]
      k_ll <- data.frame(LON = dd$lon_st,
                         LAT = dd$lat_st)
      k_pre <- SpatialPoints(k_ll, proj4string = CRS("+init=epsg:4326"))
      k_pts <- spTransform(k_pre, CRS(proj4string(Ant)))
      
      #---------------------#
      #plot(Ant)
      plot(AP)
      #plot CCAMLR SSMUs
      plot(get(SSMU_names[i]), add = TRUE, col = rgb(1,0,0,0.1))
      #krill trawls
      points(k_pts, col = rgb(1,0,0,0.6), pch = '.')
      #---------------------#
      
      sum(unique(dd$krill_green_weight)/1000)
      t_krill <- sum(dd$krill_green_weight)/1000
      #output
      t_out <- data.frame(ZONE = SSMU_names[i], 
                          YEAR = yrs[j], 
                          T_KRILL_TONNES = t_krill,
                          FRAC_AKER = round(t_krill/zcy, digits = 3))
      #merge with final output
      SSMU_out <- rbind(SSMU_out, t_out)
    } else {
      t_out <- data.frame(ZONE = SSMU_names[i], 
                          YEAR = yrs[j],
                          T_KRILL_TONNES = 0,
                          FRAC_AKER = NA)
      SSMU_out <- rbind(SSMU_out, t_out)
    }
  }
}






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
  # 
  # t_krill <- data.frame(LON = temp_trawls_150$lon_st,
  #                       LAT = temp_trawls_150$lat_st)
  # temp_kp <- SpatialPoints(t_krill, proj4string = CRS("+init=epsg:4326"))
  # t_kp <- spTransform(temp_kp, CRS(proj4string(Ant)))
  # 
  # plot(AP)
  # plot(all_site_buffers_150, col = rgb(0,0,1,0.1), add = TRUE)
  # points(temp_site, col = rgb(0,0,1,0.5), pch = 19)
  # plot(temp_buffer_150, col = rgb(1,0,0,0.1), add = TRUE)
  # points(t_kp, col = rgb(1,0,0,0.05), pch = '.')
  #---------------#
}


#master_krill is all trawls that occured within the site buffers for the sites we have PW data for - all years (not filtered for years we have PW data for)



# Time frame for krill processing -----------------------------------------

#years included in krill data output

yrs <- 2012:2017



# Effect of krill fishing during each breeding season ---------------------

#25km radius for Dec - Feb in each year
#YEAR is PW year

#NOT MUCH FISHING EFFORT AT SITES/TIMES WE HAVE DATA FOR

kbs <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 1
  temp_krill <- filter(master_krill_25, col_id == cam_sites[i])
  temp_PW <- filter(PW_data, site == cam_sites[i])
  
  #PW year (1999/2000 season is PW year 2000)
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


# write.csv(kbs, 'krill_breeding_season.csv', row.names = FALSE)




# Effect of krill fishing during whole season ----------------------------


#150km radius for March - Feb (e.g., March 1999 - Feb 2000 for 1999/2000 breeding season)
#YEAR is PW year

kws <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 1
  temp_krill <- filter(master_krill_150, col_id == cam_sites[i])
  temp_PW <- filter(PW_data, site == cam_sites[i])
  
  #PW year (1999/2000 season is PW year 2000)
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
      kws <- rbind(kws, t_out)
    }
  }
}


# write.csv(kws, 'krill_entire_season.csv', row.names = FALSE)




# Effect of krill fishing across years ------------------------------------

#150km radius average (or total) across all years at each site
#YEAR is PW year

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


# write.csv(kay, 'krill_average.csv', row.names = FALSE)





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

