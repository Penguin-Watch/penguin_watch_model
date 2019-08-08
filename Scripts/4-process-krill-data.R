#################
# Process krill data
#
# CCAMLR krill data: https://www.ccamlr.org/en/data/statistical-bulletin
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



# top level dir -----------------------------------------------------------


dir <- '~/Google_Drive/R/penguin_watch_model/'

OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2019-07-17'


# Load packages -----------------------------------------------------------

library(dplyr)
library(rgdal)
library(rgeos)
library(ggplot2)



# create site buffers ----------------------------------------------------------


#created buffers around actual lat/lons - don't want low-res land mask to interfere with krill trawls
#determine which sites we have PW data for

setwd(paste0(dir, 'Data'))

#read in lat/lon
site_md <- read.csv('site.csv', stringsAsFactors = FALSE)
SLL_p <- unique(site_md[,c('site_id', 'longitude', 'latitude')])
colnames(SLL_p) <- c('SITE', 'col_lon', 'col_lat')

#add SG sites
sg_add <- data.frame(SITE = c('COOP', 'GODH', 'MAIV', 'OCEA'), 
                     col_lon = c(-35.83, -36.26, -36.50, -36.27),
                     col_lat = c(-54.78, -54.29, -54.24, -54.34))
                              
SLL <- rbind(SLL_p, sg_add)

# #AP
# setwd('../peninsula/')
# AP_p <- rgdal::readOGR('GADM_peninsula.shp')
# AP <- sp::spTransform(AP_p, CRS(proj4string(Ant)))

#CCAMLR zones
# setwd('../asd-shapefile-WGS84/')
# mz <- rgdal::readOGR('asd-shapefile-WGS84.shp')
# CCAMLR_zones <- sp::spTransform(mz, CRS(proj4string(Ant)))
# sub_481 <- CCAMLR_zones[which(CCAMLR_zones@data$Name == 'Subarea 48.1'),]

#small scale management units (SSMU)
setwd('../ssmu-shapefile-WGS84/')
sm <- rgdal::readOGR('ssmu-shapefile-WGS84.shp')
SSMU <- sp::spTransform(sm, CRS(proj4string(Ant)))
SSMU_names <- as.character(unique(SSMU@data$ShortLabel))


for (i in 1:length(SSMU_names))
{
  #i <- 1
  assign(SSMU_names[i], SSMU[which(SSMU@data$ShortLabel == SSMU_names[i]),])
}


#buffers for each colony
# p_SLL <- cbind(SLL$col_lon, SLL$col_lat)

#points are in 4326 (uses lat/lon)
# col_points <- sp::SpatialPoints(p_SLL, proj4string = CRS('+init=epsg:4326'))

#convert colony points to 3031 (rgeos expects projected spatial object)
# t_col_points <- sp::spTransform(col_points, CRS(proj4string(Ant)))

# #create buffers around sites - 150km
# all_site_buffers_150 <- rgeos::gBuffer(t_col_points, width = 150000)
# all_site_buffers_100 <- rgeos::gBuffer(t_col_points, width = 100000)
# all_site_buffers_50 <- rgeos::gBuffer(t_col_points, width = 50000)
# all_site_buffers_25 <- rgeos::gBuffer(t_col_points, width = 25000)



# Spatial intersection of SSMU and site buffers - Weight krill catch ---------------------------

#CCAMLR
#PW year is year t-1/t season (year 2000 is 1999/2000 season)
#CCCAMLR year is t/t+1 Dec 1 - Nov 30 (year 1999 is Dec 1, 1999 - Nov 30, 2000)
#weight is in TONNES
setwd(paste0(dir, 'Data/Krill_data/CCAMLR/CCAMLR_2019_Statistical_Bulletin_Volume_31_Data_Files/'))

#1985-2018 calendar years (through 2017/2018 season)
CCAMLR_krill <- read.csv('AggregatedKrillCatch.csv')

#sometimes there are 0s in CCAMLR datasheet, sometimes no entries at all. Assume that no entries means that there was no fishing done in that sector in that month/year, while zero means there was fishing done, but nothing was caught.

#function to calculate weight krill catch based on percent overlap between buffer zone and SSMU
#WEIGHT IS IN TONNES
#Buffer size in km

weight_krill_fun <- function(BUFFER_SIZE = 150)
{
  #BUFFER_SIZE = 25
  
  #create buffers for each site
  #buffer size in KM
  for (i in 1:NROW(SLL))
  {
    #i <- 738
    
    tt_site <- cbind(SLL$col_lon[i], SLL$col_lat[i])
    #convert to spatial points
    temp_site_sp <- sp::SpatialPoints(tt_site, proj4string = CRS('+init=epsg:4326'))
  
    #transform to 3031
    temp_site <- sp::spTransform(temp_site_sp, CRS(proj4string(Ant)))
    #width is in m so multiple by 1k
    assign(SLL$SITE[i], rgeos::gBuffer(temp_site, width = BUFFER_SIZE*1000))
  }


  #create data.frame that shows which SSMU each buffer intersects (more than 10% buffer area)
  #remove APPA (AP pelagic area)
  SSMU_n <- SSMU[which(!SSMU@data$ShortLabel %in% c('APPA')),]
  
  #empty data.frame with SSMU as colnames
  zone_ovl <- data.frame(matrix(vector(), 
                                NROW(SLL), 
                                length(SSMU_n@data$ShortLabel)+1))
  colnames(zone_ovl) <- c('SITE', as.character(SSMU_n@data$ShortLabel))
  
  for (i in 1:NROW(SLL))
  {
    #i <- 738
    t_data <- get(SLL$SITE[i])
    vals <- which(rgeos::gIntersects(t_data, SSMU_n, byid = TRUE) == TRUE)
  
    int <- rgeos::gIntersection(t_data, SSMU_n, byid = TRUE)
  
    #-----------#
    #plot check
    # plot(AP)
    # plot(SSMU_n, add = TRUE)
    # plot(t_data, add = TRUE, col = rgb(1,0,0,0.2))
    # plot(int, add = T, col = rgb(0,1,0,0.5))
    #-----------#
  
    if (!is.null(int))
    {
      int_area <- raster::area(int)
      buff_area <- raster::area(t_data)
      per_area <- round(int_area/buff_area, digits = 3)
      #n_vals <- vals[which(per_area > 0.1)]
      zones <- SSMU_n@data$ShortLabel[vals]
  
      zone_ovl[i,c(1,(vals+1))] <- c(SLL$SITE[i], per_area)
    } else {
      zone_ovl[i,c(1,(vals+1))] <- c(SLL$SITE[i])
    }
  }
  
  #weight krill catch based on percent overlap between buffer zone and SSMU
  #WEIGHT IS IN TONNES
  yrs <- range(CCAMLR_krill$Calendar_Year)[1]:range(CCAMLR_krill$Calendar_Year)[2]
  cn <- colnames(zone_ovl)[-1]
  krill_weighted <- data.frame()
  
  #for each cam site
  for (i in 1:NROW(zone_ovl))
  {
    #i <- 738
    #which cols are not NA
    tsite <- zone_ovl[i, 1]
    tzone <- zone_ovl[i, ]
    tind <- which(!is.na(tzone[-1]))
    tper <- tzone[which(!is.na(zone_ovl[i, -1]))+1]
  
    #what does the area add up to
    total_area <- sum(as.numeric(tper))
    #fraction of that area that is made up by each of the SSMU
    final_per <- as.numeric(tper)/total_area
  
    #years
    for (k in 1:length(yrs))
    {
      #k <- 1
      temp_yr <- dplyr::filter(CCAMLR_krill, Calendar_Year == yrs[k])
    
      #each month
      for (m in 1:12)
      {
        #m <- 12
        temp_mn <- dplyr::filter(temp_yr, Month == m)
      
        t2_k <- c()
        if (length(tind) > 0)
        {
          for (j in 1:length(tind))
          {
            #j <- 1
            temp_z <- dplyr::filter(temp_mn, SSMU_Code == cn[tind[j]])
        
            #krill catch weighted by spatial overlap - final_per
            t_k <- temp_z$Krill_Green_Weight * as.numeric(final_per[j])
      
            if (length(t_k) > 0)
            {
              t2_k <- c(t2_k, t_k)
            } else {
              t2_k <- c(t2_k, 0)
            }
          }
        }
        
        #sum the weighted krill across all SSMU
        temp_mn_weight_kr <- round(sum(t2_k), digits = 2)
        
        t_mn_df <- data.frame(SITE = tsite, YEAR = yrs[k], MONTH = m, 
                              WEIGHTED_KRILL = temp_mn_weight_kr)
        krill_weighted <- rbind(krill_weighted, t_mn_df)
      }
    }
  }
  return(krill_weighted)
}

#krill weighted is total krill catch within site buffer, weighted by percent overlap with SSMU
krill_weighted_25 <- weight_krill_fun(BUFFER_SIZE = 25)
krill_weighted_150 <- weight_krill_fun(BUFFER_SIZE = 150)




# Time frame for krill covariate processing -----------------------------------------

#PW years included in krill data output (1999/2000 season is PW year 2000)
yrs <- 2000:2018



# Effect of krill fishing during each breeding season ---------------------

#25km radius for Dec - Feb in each year
#YEAR is PW year

#CCAMLR data
CCAMLR_kr_BS <- data.frame()
for (i in 1:NROW(SLL))
{
  #i <- 1
  temp_krill <- dplyr::filter(krill_weighted_25, SITE == SLL$SITE[i])

  #PW year (1999/2000 season is PW year 2000)
  for (j in 1:length(yrs))
  {
    #j <- 1
    #used Dec 1 - Feb 1 (Feb 1 was perscribed end of PW data)

    yr_one <- dplyr::filter(temp_krill, YEAR == yrs[j]-1)
    DEC <- dplyr::filter(yr_one, MONTH == 12)
    yr_two <- dplyr::filter(temp_krill, YEAR == yrs[j])
    JAN <- dplyr::filter(yr_two, MONTH == 1)

    #total krill caught over this period
    t_krill <- sum(c(DEC$WEIGHTED_KRILL, JAN$WEIGHTED_KRILL), na.rm = TRUE)

    if (!is.na(t_krill))
    {
      #output
      t_out <- data.frame(SITE = SLL$SITE[i],
                          YEAR = yrs[j],
                          T_KRILL = t_krill)
      #merge with final output
      CCAMLR_kr_BS <- rbind(CCAMLR_kr_BS, t_out)
    } else {
      t_out <- data.frame(SITE = SLL$SITE[i],
                          YEAR = yrs[j],
                          T_KRILL = 0)
      CCAMLR_kr_BS <- rbind(CCAMLR_kr_BS, t_out)
    }
  }
}

setwd('../Processed_CCAMLR/')
write.csv(CCAMLR_kr_BS, 'CCAMLR_krill_breeding_season.csv', row.names = FALSE)



# Effect of krill fishing during whole season ----------------------------


#150km radius for March - Jan (e.g., March 1999 - Feb 2000 for 1999/2000 breeding season)
#YEAR is PW year

CCAMLR_kr_WS <- data.frame()
for (i in 1:NROW(SLL))
{
  #i <- 1
  temp_krill <- dplyr::filter(krill_weighted_150, SITE == SLL$SITE[i])
  
  #PW year (1999/2000 season is PW year 2000)
  for (j in 1:length(yrs))
  {
    #j <- 1
    #used March 1 - Feb 1 (Chicks appear to be creche by Feb 1 - catch after this period would be irrelevant for survival)
    yr_one <- dplyr::filter(temp_krill, YEAR == yrs[j]-1)
    MAR_DEC <- dplyr::filter(yr_one, MONTH > 2)
    yr_two <- dplyr::filter(temp_krill, YEAR == yrs[j])
    JAN <- dplyr::filter(yr_two, MONTH == 1)
    
    #total krill caught over this period
    t_krill <- sum(c(MAR_DEC$WEIGHTED_KRILL, JAN$WEIGHTED_KRILL), na.rm = TRUE)
    
    if (!is.na(t_krill))
    {
      #output
      t_out <- data.frame(SITE = SLL$SITE[i],
                          YEAR = yrs[j],
                          T_KRILL = t_krill)
      #merge with final output
      CCAMLR_kr_WS <- rbind(CCAMLR_kr_WS, t_out)
    } else {
      t_out <- data.frame(SITE = SLL$SITE[i],
                          YEAR = yrs[j],
                          T_KRILL = 0)
      CCAMLR_kr_WS <- rbind(CCAMLR_kr_WS, t_out)
    }
  }
}


setwd('../Processed_CCAMLR/')
write.csv(CCAMLR_kr_WS, 'CCAMLR_krill_entire_season.csv', row.names = FALSE)



# Effect of krill fishing across years ------------------------------------

#MARCH - FEB
#150km radius average (or total) across all years at each site
#YEAR is PW year

CCAMLR_kr_AY <- data.frame()
for (i in 1:NROW(SLL))
{
  #i <- 2
  temp_krill <- dplyr::filter(krill_weighted_150, SITE == SLL$SITE[i])
  
  #PW year (1999/2000 season is PW year 2000)
  t_yr_krill <- c()
  for (j in 1:length(yrs))
  {
    #j <- 1
    #used March 1 - Feb 1 (Feb 1 was perscribed end of PW data)
    yr_one <- dplyr::filter(temp_krill, YEAR == yrs[j]-1)
    MAR_DEC <- dplyr::filter(yr_one, MONTH > 2)
    yr_two <- dplyr::filter(temp_krill, YEAR == yrs[j])
    JAN_FEB <- dplyr::filter(yr_two, MONTH <= 2)
    
    #total krill caught over this period
    t_krill <- sum(c(MAR_DEC$WEIGHTED_KRILL, JAN_FEB$WEIGHTED_KRILL), na.rm = TRUE)
    
    t_yr_krill <- c(t_yr_krill, t_krill)
  }
  
  mn_yr_krill <- mean(t_yr_krill, na.rm = TRUE)
  
  if (!is.na(mn_yr_krill))
  {
    #output
    t_out <- data.frame(SITE = SLL$SITE[i],
                        T_KRILL = mn_yr_krill)
    #merge with final output
    CCAMLR_kr_AY <- rbind(CCAMLR_kr_AY, t_out)
  } else {
    t_out <- data.frame(SITE = SLL$SITE[i],
                        T_KRILL = 0)
    CCAMLR_kr_AY <- rbind(CCAMLR_kr_AY, t_out)
  }
}

write.csv(CCAMLR_kr_AY, 'CCAMLR_krill_average.csv', row.names = FALSE)



# which sites have the most fishing at them? -------------------------------

#CCAMLR
ggplot(CCAMLR_kr_AY, aes(SITE, T_KRILL)) +
  geom_col() +
  ylab('Mean krill catch (tonnes)') +
  theme_bw() +
  ggtitle('CCAMLR - Mean krill catch March - Jan (2000-2018)')



# When is krill fishing most intense in this region? ----------------------

#CCAMLR
krill_time <- data.frame(MONTH = as.integer(1:12), KRILL = rep(NA, 12))
for (i in 1:12)
{
  #i <- 1
  temp <- dplyr::filter(krill_weighted_150, MONTH == i)
  swk <- sum(temp$WEIGHTED_KRILL, na.rm = TRUE)
  krill_time[i,2] <- swk
}

ggplot(krill_time, aes(x = MONTH, y = KRILL)) +
  geom_col() + 
  theme_bw() +
  ggtitle('CCAMLR - Krill catch by month')

