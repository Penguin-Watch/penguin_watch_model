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
OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2019-08-23'


# Load packages -----------------------------------------------------------

library(dplyr)
library(rgdal)
library(rgeos)
library(raster)
library(sp)


# read in merged data -----------------------------------------------------

setwd(OUTPUT)

#read in RDS from 3-analyze-output.R
bs_precip_mrg <- readRDS('bs_precip_mrg.rds')



# create site buffers ----------------------------------------------------------

setwd(paste0(dir, 'Data'))

#read in lat/lon for Antarctic continent sites
site_md <- read.csv('site.csv', stringsAsFactors = FALSE)
SLL_p <- unique(site_md[,c('site_id', 'longitude', 'latitude')])
colnames(SLL_p) <- c('SITE', 'col_lon', 'col_lat')

#add SG sites
sg_add <- data.frame(SITE = c('COOP', 'GODH', 'MAIV', 'OCEA'), 
                     col_lon = c(-35.83, -36.26, -36.50, -36.27),
                     col_lat = c(-54.78, -54.29, -54.24, -54.34))
                              
SLL <- rbind(SLL_p, sg_add)


# #AP
setwd(paste0(dir, 'Data/peninsula'))
AP_p <- rgdal::readOGR('GADM_peninsula.shp')
AP <- sp::spTransform(AP_p, CRS("+init=epsg:3031"))
#Sub antarctic
setwd(paste0(dir, 'Data/Sub-antarctic_coastline_low_res_polygon'))
SA_p <- rgdal::readOGR('Sub-antarctic_coastline_low_res_polygon.shp')
SA <- sp::spTransform(SA_p, CRS("+init=epsg:3031"))
AP_SA <- raster::union(AP, SA)



#CCAMLR zones
# setwd(paste0(dir, 'Data/asd-shapefile-WGS84'))
# mz <- rgdal::readOGR('asd-shapefile-WGS84.shp')
# CCAMLR_zones <- sp::spTransform(mz, CRS("+init=epsg:3031"))
# sub_481 <- CCAMLR_zones[which(CCAMLR_zones@data$Name == 'Subarea 48.1'),]

#CCAMLR small scale management units (SSMU)
setwd(paste0(dir, 'Data/ssmu-shapefile-WGS84/'))
sm <- rgdal::readOGR('ssmu-shapefile-WGS84.shp')
SSMU <- sp::spTransform(sm, CRS("+init=epsg:3031"))
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
# t_col_points <- sp::spTransform(col_points, CRS("+init=epsg:3031"))

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

weight_krill_fun <- function(BUFFER_SIZE = 150, SITES = SLL, PLOT = FALSE)
{
  #BUFFER_SIZE = 100
  #SITES = sel_sites
  
  #create buffers for each site
  #buffer size in KM
  for (i in 1:NROW(SITES))
  {
    #i <- 1
    
    tt_site <- cbind(SITES$col_lon[i], SITES$col_lat[i])
    #convert to spatial points
    temp_site_sp <- sp::SpatialPoints(tt_site, proj4string = CRS('+init=epsg:4326'))
  
    #transform to 3031
    temp_site <- sp::spTransform(temp_site_sp, CRS('+init=epsg:3031'))
    #width is in m so multiple by 1k
    assign(SITES$SITE[i], rgeos::gBuffer(temp_site, width = BUFFER_SIZE*1000))
  }


  #create data.frame that shows which SSMU each buffer intersects (more than 10% buffer area)
  #remove APPA (AP pelagic area), SGPA (SG pelagic area), ad SOPA (SO pelagic area)
  SSMU_n <- SSMU[which(!SSMU@data$ShortLabel %in% c('APPA', 'SGPA', 'SOPA',
                                                    'SSI', 'SSPA')),]
  
  if (PLOT == TRUE)
  {
    #-----------#
    #plot check
    plot(AP_SA, xlim = c(-1598157, -1048731), ylim = c(706494, 3310191))
    plot(SSMU_n, add = TRUE)
  }
  
  #empty data.frame with SSMU as colnames
  zone_ovl <- data.frame(matrix(vector(), 
                                NROW(SITES), 
                                length(SSMU_n@data$ShortLabel)+1))
  colnames(zone_ovl) <- c('SITE', as.character(SSMU_n@data$ShortLabel))
  
  for (i in 1:NROW(SITES))
  {
    #i <- 4
    #polygon for site
    t_data <- get(SITES$SITE[i])
    
    vals <- which(rgeos::gIntersects(t_data, SSMU_n, byid = TRUE) == TRUE)
    #intersection of SSMU and site
    int <- rgeos::gIntersection(t_data, SSMU_n, byid = TRUE)
  
    #see if polygon is segmented - if so, exclude segmented part (happens when Weddell Sea is included with sites on Western AP)
    if (!is.null(int))
    {
      #if there is more than one feature in polygon (more than one SSMU)
      if (length(int) > 0)
      {
        #see which segments are touching
        touching <- rgeos::gTouches(int, byid = TRUE)
        #which segments are not touching any others
        zz <- which(apply(touching, 2, sum) == 0)
        
        if (length(zz) > 0 & NROW(touching) > 1)
        {
          #if only two segments, take larger
          if (length(int) == 2)
          {
            mint <- which.max(area(int))
            nint <- int[mint,]
            vals2 <- vals[mint]
          } else {
            nint <- int[-zz,]
            vals2 <- vals[-zz]
          }
        } else {
          nint <- int
          vals2 <- vals
        }
      } else {
        nint <- int
        vals2 <- vals
      }
      
      #assign(paste0(SITES$SITE[i], '_nint'), nint)
      
      if (PLOT == TRUE)
      {
        #-----------#
        #plot check
        # plot(t_data, add = TRUE, col = rgb(1,0,0,0.2))
        plot(nint, add = T, col = rgb(1,0,0,0.5))
        #-----------#
      }
      
      #area of SSMU zones
      int_area <- raster::area(nint)
      
      zones <- SSMU_n@data$ShortLabel[vals2]
  
      #area of each SSMU
      zone_ovl[i,c(1,(vals2+1))] <- c(SITES$SITE[i], int_area)
    } else {
      zone_ovl[i,c(1,(vals2+1))] <- c(SITES$SITE[i])
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
    #i <- 427
    #which cols are not NA
    tsite <- zone_ovl[i, 1]
    tzone <- zone_ovl[i, ]
    tind <- which(!is.na(tzone[-1]))
    tarea <- tzone[which(!is.na(zone_ovl[i, -1])) + 1]
  
    #areas
    final_areas <- as.numeric(tarea)
  
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
            #j <- 4
            temp_z <- dplyr::filter(temp_mn, SSMU_Code == cn[tind[j]])
        
            #want mean krill
            
            #krill catch
            t_k <- temp_z$Krill_Green_Weight
      
            if (length(t_k) > 0)
            {
              t2_k <- c(t2_k, t_k)
            } else {
              t2_k <- c(t2_k, 0)
            }
          }
        }
        
        #proportion of total area in each SSMU
        per <- final_areas / sum(final_areas)
        
        #sum the weighted mean krill across all SSMU
        temp_weight_kr <- round(weighted.mean(t2_k, per), 3)
        
        t_mn_df <- data.frame(SITE = tsite, YEAR = yrs[k], MONTH = m, 
                              WEIGHTED_KRILL = temp_weight_kr)
        krill_weighted <- rbind(krill_weighted, t_mn_df)
      }
    }
  }
  return(krill_weighted)
}

#krill weighted is average krill catch within site buffer, weighted by percent overlap with SSMU
#just study sties
sel_sites <- unique(bs_precip_mrg[,c('SITE', 'col_lon', 'col_lat')])

krill_weighted_25 <- weight_krill_fun(BUFFER_SIZE = 25,
                                      SITES = sel_sites,
                                      PLOT = FALSE)
krill_weighted_100 <- weight_krill_fun(BUFFER_SIZE = 100,
                                       SITES = sel_sites, 
                                       PLOT = FALSE)
krill_weighted_150 <- weight_krill_fun(BUFFER_SIZE = 150,
                                       SITES = sel_sites, 
                                       PLOT = FALSE)

setwd('../Processed_CCAMLR/')
saveRDS(krill_weighted_25, 'krill_weighted_25.rds')
saveRDS(krill_weighted_100, 'krill_weighted_100.rds')
saveRDS(krill_weighted_150, 'krill_weighted_150.rds')



# Time frame and sites for krill covariate processing -----------------------------

#PW years included in krill data output (1999/2000 season is PW year 2000)
yrs <- 2000:2018
SITES <- sel_sites



# Effect of krill fishing during each breeding season ---------------------

#25km radius for Dec - Feb in each year
#YEAR is PW year

#CCAMLR data
CCAMLR_kr_BS <- data.frame()
for (i in 1:NROW(SITES))
{
  #i <- 1
  temp_krill <- dplyr::filter(krill_weighted_25, SITE == SITES$SITE[i])

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
      t_out <- data.frame(SITE = SITES$SITE[i],
                          YEAR = yrs[j],
                          T_KRILL = t_krill)
      #merge with final output
      CCAMLR_kr_BS <- rbind(CCAMLR_kr_BS, t_out)
    } else {
      t_out <- data.frame(SITE = SITES$SITE[i],
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
for (i in 1:NROW(SITES))
{
  #i <- 1
  temp_krill <- dplyr::filter(krill_weighted_150, SITE == SITES$SITE[i])
  
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
      t_out <- data.frame(SITE = SITES$SITE[i],
                          YEAR = yrs[j],
                          T_KRILL = t_krill)
      #merge with final output
      CCAMLR_kr_WS <- rbind(CCAMLR_kr_WS, t_out)
    } else {
      t_out <- data.frame(SITE = SITES$SITE[i],
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
for (i in 1:NROW(SITES))
{
  #i <- 2
  temp_krill <- dplyr::filter(krill_weighted_150, SITE == SITES$SITE[i])
  
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
    t_out <- data.frame(SITE = SITES$SITE[i],
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




# merge krill with bs_precip data -----------------------------------------------

setwd(OUTPUT)

#krill catch breeding season
mrg_BS <- dplyr::left_join(bs_precip_mrg, CCAMLR_kr_BS, by = c('SITE', 'YEAR'))
#krill catch entire year
mrg_BS_WS <- dplyr::left_join(mrg_BS, CCAMLR_kr_WS, by = c('SITE', 'YEAR'))
#average krill catch entire year
mrg_BS_WS_AY <- dplyr::left_join(mrg_BS_WS, CCAMLR_kr_AY, by = 'SITE')

#new colnames
coln <- colnames(mrg_BS_WS_AY)
cidx<- grep('T_KRILL', coln)
colnames(mrg_BS_WS_AY)[cidx] <- c('krill_BR', 'krill_WS', 'krill_AY')

#save as rds
saveRDS(mrg_BS_WS_AY, 'bs_precip_krill_mrg.rds')

