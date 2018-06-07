######################
#Script to plot map with all sites and buffers
#
#
######################


# Clear environment -------------------------------------------------------


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


#cam sites is taken from 'Camera List' tab of google doc
#removed a number of sites that were not in ASI sites expanded (Mostly SG, Falklands, SSI, I think)
cam_sites <- as.character(read.csv('cam_sites.csv', header = FALSE)[,1])
ASI <- read.csv('ASI_sites_expanded.csv')


SLL <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 4
  idx <- ASI$Site.code %in% cam_sites[i]
  lat <- ASI$Lat[idx]
  lon <- ASI$Lon[idx]
  temp <- data.frame(SITE = cam_sites[i], 
                     NAME = ASI$Hotspot.Name[idx],
                     LAT = lat, 
                     LON = lon)
  SLL <- rbind(SLL, temp)
}


# create site buffers ----------------------------------------------------------


#created buffers around actual lat/lons - don't want low-res land mask to interfere with krill trawls


#if using actual data

# #determine which sites we have PW data for
# setwd('~/Google_Drive/R/penguin_watch_model/Data/PW_data/')
# PW_data <- read.csv('Markrecap_data_15.05.18.csv', stringsAsFactors = FALSE)
# 
# un_sites <- unique(PW_data$site)
# 
# site_year <- c()
# for(i in 1:length(un_sites))
# {
#   #i <- 1
#   tdat <- filter(PW_data, site == un_sites[i])
#   tsy <- unique(tdat$season_year)
#   tout <- data.frame(SITE = rep(un_sites[i], length(tsy)), YEAR = tsy)
#   site_year <- rbind(site_year, tout)
# }
# 
# #site/years we have data for:
# #site_year
# 
# #BOOT is now PCHA in MAPPPD database
# pos <- which(un_sites == 'BOOT')
# cam_sites_p <- un_sites[-pos]
# cam_sites <- c(cam_sites_p, 'PCHA')
# 
# PW_data$site[which(PW_data$site == 'BOOT')] <- 'PCHA'
# 
# 
# #determine lat/lons for all site we have data for
# SLL <- data.frame(SITE = cam_sites, LON = rep(NA, length(cam_sites)), LAT = rep(NA, length(cam_sites)))
# for (i in 1:length(cam_sites))
# {
#   #i <- 1
#   temp <- filter(PW_data, site == cam_sites[i])[1,]
#   SLL$LON[i] <- temp$col_lon
#   SLL$LAT[i] <- temp$col_lat
# }




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




# calculate CCAMLR SSMU krill catch for all years, all sites --------------------------


setwd('../Krill_data/CCAMLR/CCAMLR_2017_Statistical_Bulletin_Volume_29_Data_Files/')

#1985-2016 calendar years
CCAMLR_krill <- read.csv('AggregatedKrillCatch.csv')



#function to calculate weight krill catch based on percent overlap between buffer zone and SSMU
#WEIGHT IS IN TONNES
#Buffer size in km

weight_krill_fun <- function(BUFFER_SIZE = 150)
{
  #BUFFER_SIZE = 150
  #create buffers for each site
  #buffer size in KM
  for (i in 1:length(cam_sites))
  {
    #i <- 1
    temp_site <- filter(SLL, SITE == cam_sites[i])
    
    tt <- cbind(temp_site$LON, temp_site$LAT)
    
    #convert to spatial points
    temp_site_sp <- SpatialPoints(tt, proj4string = CRS('+init=epsg:4326'))
    
    #transform to 3031
    temp_site <- spTransform(temp_site_sp, CRS(proj4string(Ant)))
    #width is in m so multiple by 1k
    assign(cam_sites[i], rgeos::gBuffer(temp_site, width = BUFFER_SIZE*1000))
  }
  
  
  #create data.frame that shows which SSMU each buffer intersects (more than 10% buffer area)
  #remove APPA (AP pelagic area) and APE (AP East)
  SSMU_n <- SSMU[which(!SSMU@data$ShortLabel %in% c('APPA', 'APE')),]
  #empty data.frame with SSMU as colnames
  zone_ovl <- data.frame(matrix(vector(), 
                                length(cam_sites), 
                                length(SSMU_n@data$ShortLabel)+1))
  colnames(zone_ovl) <- c('SITE', as.character(SSMU_n@data$ShortLabel))
  
  for (i in 1:length(cam_sites))
  {
    #i <- 1
    t_data <- get(cam_sites[i])
    vals <- which(gIntersects(t_data, SSMU_n, byid = TRUE) == TRUE)
    
    int <- gIntersection(t_data, SSMU_n, byid = TRUE)
    
    #-----------#
    #plot check
    # plot(AP)
    # plot(SSMU_n, add = TRUE)
    # plot(t_data, add = TRUE, col = rgb(1,0,0,0.2))
    # plot(int, add = T, col = rgb(0,1,0,0.5))
    #-----------#
    
    int_area <- raster::area(int)
    buff_area <- raster::area(t_data)
    per_area <- round(int_area/buff_area, digits = 3)
    #n_vals <- vals[which(per_area > 0.1)]
    zones <- SSMU_n@data$ShortLabel[vals]
    
    zone_ovl[i,c(1,(vals+1))] <- c(cam_sites[i], per_area)
  }
  
  #weight krill catch based on percent overlap between buffer zone and SSMU
  #WEIGHT IS IN TONNES
  yrs <- range(CCAMLR_krill$Calendar_Year)[1]:range(CCAMLR_krill$Calendar_Year)[2]
  cn <- colnames(zone_ovl)[-1]
  krill_weighted <- data.frame()
  for (i in 1:NROW(zone_ovl))
  {
    #i <- 1
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
      temp_yr <- filter(CCAMLR_krill, Calendar_Year == yrs[k])
      
      for (m in 1:12)
      {
        #m <- 1
        temp_mn <- filter(temp_yr, Month == m)
        
        #if data exists for that month
        if (NROW(temp_mn) > 0)
        {
          t2_k <- c()
          for (j in 1:length(tind))
          {
            #j <- 1
            temp_z <- filter(temp_mn, SSMU_Code == cn[tind[j]])
            
            #krill catch weighted by spatial overlap - final_per
            t_k <- temp_z$Krill_Green_Weight * as.numeric(final_per[j])
            
            if (length(t_k) > 0)
            {
              t2_k <- c(t2_k, t_k)
            } else {
              t2_k <- c(t2_k, 0)
            }
          }
          
          #sum the weighted krill across all SSMU
          temp_mn_weight_kr <- round(sum(t2_k), digits = 2)
          
          t_mn_df <- data.frame(SITE = tsite, YEAR = yrs[k], MONTH = m, WEIGHTED_KRILL = temp_mn_weight_kr)
          krill_weighted <- rbind(krill_weighted, t_mn_df)
        } else {
          temp_mn_weight_kr <- NA
          t_mn_df <- data.frame(SITE = tsite, YEAR = yrs[k], MONTH = m, WEIGHTED_KRILL = temp_mn_weight_kr)
          krill_weighted <- rbind(krill_weighted, t_mn_df)
        }
      }
    }
  }
  return(krill_weighted)
}

#krill weighted is total krill catch within site buffer, weighted by percent overlap with SSMU
#krill_weighted_25 <- weight_krill_fun(BUFFER_SIZE = 25)
krill_weighted_150 <- weight_krill_fun(BUFFER_SIZE = 150)


#CHECK NA vals
#krill_weighted_150[which(is.na(krill_weighted_150$WEIGHTED_KRILL)),]



# Time frame for krill covariate processing -----------------------------------------

#PW years included in krill data output

#NO KRILL DATA FOR END OF 2016
yrs <- 2012:2016



CCAMLR_kr_WS <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 1
  temp_krill <- filter(krill_weighted_150, SITE == cam_sites[i])
  
  #PW year (1999/2000 season is PW year 2000)
  for (j in 1:length(yrs))
  {
    #j <- 1
    #used March 1 - Feb 1 (Feb 1 was perscribed end of PW data)
    yr_one <- filter(temp_krill, YEAR == yrs[j]-1)
    MAR_DEC <- filter(yr_one, MONTH > 2)
    yr_two <- filter(temp_krill, YEAR == yrs[j])
    JAN <- filter(yr_two, MONTH == 1)
    
    #total krill caught over this period
    t_krill <- sum(c(MAR_DEC$WEIGHTED_KRILL, JAN$WEIGHTED_KRILL), na.rm = TRUE)
    
    if (!is.na(t_krill))
    {
      #output
      t_out <- data.frame(SITE = cam_sites[i],
                          YEAR = yrs[j],
                          T_KRILL = t_krill)
      #merge with final output
      CCAMLR_kr_WS <- rbind(CCAMLR_kr_WS, t_out)
    } else {
      t_out <- data.frame(SITE = cam_sites[i],
                          YEAR = yrs[j],
                          T_KRILL = 0)
      CCAMLR_kr_WS <- rbind(CCAMLR_kr_WS, t_out)
    }
  }
}





# Effect of krill fishing across years ------------------------------------

#150km radius average (or total) across all years at each site
#YEAR is PW year

CCAMLR_kr_AY <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 1
  temp_krill <- filter(krill_weighted_150, SITE == cam_sites[i])
  
  #PW year (1999/2000 season is PW year 2000)
  t_yr_krill <- c()
  for (j in 1:length(yrs))
  {
    #j <- 1
    #used March 1 - Feb 1 (Feb 1 was perscribed end of PW data)
    yr_one <- filter(temp_krill, YEAR == yrs[j]-1)
    MAR_DEC <- filter(yr_one, MONTH > 2)
    yr_two <- filter(temp_krill, YEAR == yrs[j])
    JAN <- filter(yr_two, MONTH == 1)
    
    #total krill caught over this period
    t_krill <- sum(c(MAR_DEC$WEIGHTED_KRILL, JAN$WEIGHTED_KRILL), na.rm = TRUE)
    
    t_yr_krill <- c(t_yr_krill, t_krill)
  }
  
  mn_yr_krill <- mean(t_yr_krill, na.rm = TRUE)
  
  if (!is.na(mn_yr_krill))
  {
    #output
    t_out <- data.frame(SITE = cam_sites[i],
                        T_KRILL = mn_yr_krill)
    #merge with final output
    CCAMLR_kr_AY <- rbind(CCAMLR_kr_AY, t_out)
  } else {
    t_out <- data.frame(SITE = cam_sites[i],
                        T_KRILL = 0)
    CCAMLR_kr_AY <- rbind(CCAMLR_kr_AY, t_out)
  }
}



# which sites have the most fishing at them? -------------------------------

#CCAMLR
ggplot(CCAMLR_kr_AY, aes(SITE, T_KRILL)) +
  geom_col() +
  ylab('Mean krill catch (tonnes)') +
  theme_bw() +
  ggtitle('CCAMLR - Mean krill catch (2012-2017)')


# #AKER
# ggplot(aker_kr_ay, aes(SITE, TOTAL_KRILL)) +
#   geom_col() +
#   ylab('Total krill catch') +
#   theme_bw() +
#   ggtitle('Aker - Total krill catch across sites (2010-2017)')
# 
# ggplot(aker_kr_ay, aes(SITE, MN_KRILL)) +
#   geom_col() +
#   ylab('Mean krill catch') +
#   theme_bw() +
#   ggtitle('Aker - Mean krill catch across sites (2010-2017)')
# 
# ggplot(aker_kr_ay, aes(SITE, MN_CPUE)) +
#   geom_col() +
#   ylab('Mean Catch Per Unit Effort') +
#   theme_bw() +
#   ggtitle('Aker - CPUE across sites (2010-2017)')




# When is krill fishing most intense in this region? ----------------------


#CCAMLR
krill_time <- data.frame(MONTH = as.integer(1:12), KRILL = rep(NA, 12))
for (i in 1:12)
{
  #i <- 1
  temp <- filter(krill_weighted_150, MONTH == i)
  swk <- sum(temp$WEIGHTED_KRILL, na.rm = TRUE)
  krill_time[i,2] <- swk
}

ggplot(krill_time, aes(x = MONTH, y = KRILL)) +
  geom_col() + 
  theme_bw() +
  ggtitle('CCAMLR - Krill catch by month')






# Average winter SIC 2012-2017 --------------------------------------------


if (Sys.info()[[1]] == 'Windows')
{
  setwd('C:/Users/Lynch Lab 7/Google_Drive/R/Project_archive/MAPPPD/SiteCovariates/PassiveMicrowaveSIC')
} else {
  setwd("~/Google_Drive/R/Project_archive/MAPPPD/SiteCovariates/PassiveMicrowaveSIC")
}


#calculate mean winter SIC in each cell for each year - years PW years (1999/2000 season is 2000)
yrs <- 2012:2017
for (i in 1:length(yrs))
{
  #i <- 2
  june <- raster::raster(paste0('nt_', yrs[i] - 1, '06_f17_v1.1_s.tif'))
  july <- raster::raster(paste0('nt_', yrs[i] - 1, '07_f17_v1.1_s.tif'))
  aug <- raster::raster(paste0('nt_', yrs[i] - 1, '08_f17_v1.1_s.tif'))
  sep <- raster::raster(paste0('nt_', yrs[i] - 1, '09_f17_v1.1_s.tif'))
  
  win <- raster::stack(june, july, aug, sep)
  mn_win <- mean(win)
  #plot(mn_win, main = yrs[i])
  
  assign(paste0('mn_win_', yrs[i]), mn_win)
}


#calculate mean across years for each cell
ay <- raster::stack(mn_win_2012, mn_win_2013, mn_win_2014, 
                    mn_win_2015, mn_win_2016, mn_win_2017)
mn_ay <- mean(ay)



# Fluctuations in krill and SIC over time at each site -------------------------------


#Order sites from N -> S latitude

site_lats <- SLL[,c(1,3)]
k_lat_t <- left_join(CCAMLR_kr_WS, site_lats, by = 'SITE')
k_lat_t2 <- k_lat_t[order(k_lat_t$LAT, decreasing = TRUE),]
k_lat <- mutate(k_lat_t2, idx = rep(1:length(unique(k_lat_t2$LAT)), each = 5))


#krill - 150km buffer - years PW years (1999/2000 season is 2000)
#color denotes N -> S latitude

ggplot(k_lat, aes(YEAR, T_KRILL, group = SITE, color = idx)) +
  geom_line(size = 1.2) + 
  theme_bw() +
  scale_color_gradient(low = 'grey', high = 'black') + 
  ggtitle('KRILL - 150 km buffer - grey = high lat; black = low lat')


# #mean changes in krill over time
# yrs <- 2012:2016
# avg_kr <- data.frame()
# for (i in 1:length(yrs))
# {
#   #i <- 1
#   temp <- filter(CCAMLR_kr_WS, YEAR == yrs[i])
#   tk <- mean(temp$T_KRILL)
#   tout <- data.frame(yrs[i], tk)
#   avg_kr <- rbind(avg_kr, tout)
# }
# 
# colnames(avg_kr) <- c('YEAR', 'T_KRILL')
# 
# ggplot(avg_kr, aes(YEAR, T_KRILL)) +
#   geom_line(size = 1.2) + 
#   theme_bw()



#SIC - 150km buffer

if (Sys.info()[[1]] == 'Windows')
{
  setwd('C:/Users/Lynch Lab 7/Google_Drive/R/penguin_watch_model/Data/SIC_data/RAW/')
} else {
  setwd('~/Google_Drive/R/penguin_watch_model/Data/SIC_data/RAW')
}

data_150 <- read.csv('ALLSITES_SIC_150_MEAN.csv')

SIC_150_W <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 1
  temp_SIC <- filter(data_150, site_id == cam_sites[i])
  
  T_OUT <- data.frame()
  for (j in 1979:2016)
  {
    #j <- 1979
    temp <- filter(temp_SIC, year == j)
    jun <- temp$SIC_MONTH_6
    jul <- temp$SIC_MONTH_7
    aug <- temp$SIC_MONTH_8
    sep <- temp$SIC_MONTH_9
    mn_jjas <- mean(c(jun, jul, aug, sep))
    temp2 <- data.frame(YEAR = j+1,
                        SITE = cam_sites[i],
                        JUN = jun,
                        JUL = jul,
                        AUG = aug,
                        SEP = sep,
                        W_MN = mn_jjas)
    T_OUT <- rbind(T_OUT, temp2)
  }
  SIC_150_W <- rbind(SIC_150_W, T_OUT)
}

n_SIC_1 <- filter(SIC_150_W, YEAR >= 2012, YEAR <= 2016)
n_SIC_2 <- left_join(n_SIC_1, site_lats, by = 'SITE')
n_SIC_3 <- n_SIC_2[order(n_SIC_2$LAT, decreasing = TRUE),]
SIC_lat <- mutate(n_SIC_3, idx = rep(1:length(unique(n_SIC_3$LAT)), each = 5))


ggplot(SIC_lat, aes(YEAR, W_MN, group = SITE, color = idx)) +
  geom_line(size = 1.2) +
  theme_bw() + 
  scale_color_gradient(low = 'grey', high = 'black') + 
  ggtitle('SIC - 150 km buffer - grey = high lat; black = low lat')



# create map plots ------------------------------------------------------------

#colors
gg_color_hue <- function(n, ALPHA = 1)
{
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100, alpha = ALPHA)[1:n]
}
cols <- gg_color_hue(length(SSMU_names), ALPHA = 1)




#base plot
# #AP shapefile
# plot(AP)
# plot(mn_ay)
# plot(AP, add = TRUE)
# 
# #plot CCAMLR SSMUs
# for (i in 1:length(SSMU_names))
# {
#   #i <- 2
#   SSMU_names[i]
#   plot(get(SSMU_names[i]), add = TRUE, col = cols[i])
# }
# 
# #plot PW sites
# points(t_col_points, col = rgb(0,0,1,0.5), pch = 19)
# #plot site buffers
# plot(all_site_buffers_150, col = rgb(0,0,1,0.2), add = TRUE)


#ggplot
colonies <- data.frame(SITE = SLL$SITE,
                   long = t_col_points@coords[,1],
                   lat = t_col_points@coords[,2],
                   KRILL = CCAMLR_kr_AY$T_KRILL)

# #without KRILL
# colonies <- data.frame(SITE = SLL$SITE,
#                        long = t_col_points@coords[,1],
#                        lat = t_col_points@coords[,2])



#zoom in on AP sites

#convert raster to dataframe
mn_ay_spdf <- as(mn_ay, 'SpatialPixelsDataFrame')
mn_ay_df <- as.data.frame(mn_ay_spdf)
colnames(mn_ay_df) <- c('value', 'x', 'y')

#NA vals for 251 (cicular mask), 252 (unused), 253 (coastlines), 254 (landmask), and 255 (missing data)
mn_ay_df$value[which(mn_ay_df$value > 250)] <- NA
#scale 0 - 100
mn_ay_df$value <- (mn_ay_df$value/250)*100


# #just SIC
# 
# ggplot(mn_ay_df, aes(x=x, y=y, fill = value)) +
#   geom_tile() +
#   theme_void() +
#   coord_equal(xlim = c(-2750000, -2200000),
#               ylim = c(1000000, 2000000)) +
#   scale_fill_gradient2('Mean SIC - June - Sep',
#                        limits = c(0,
#                                 100),
#                        na.value = 'white',
#                        low = 'white',
#                        high = 'light blue') +
#   geom_polygon(data = AP, aes(long, lat, group = group),
#                inherit.aes = FALSE,
#                fill = 'grey') +
#   geom_path(data = AP, aes(long, lat, group = group),
#             inherit.aes = FALSE,
#             color = 'black')

# 
# #just krill catch
# ggplot(AP, aes(long, lat, group = group)) +
#   geom_polygon(fill = 'grey') + 
#   geom_path(color = 'black') +
#   coord_equal(xlim = c(-2750000, -2200000), 
#             ylim = c(1000000, 2000000)) +
#   theme_void() +
#   # geom_polygon(data = APW,
#   #              aes(long, lat),
#   #              fill = NA,
#   #              color = 'black',
#   #              #color = cols[1],
#   #              inherit.aes = FALSE) +
#   # geom_polygon(data = APDPW,
#   #              aes(long, lat),
#   #              fill = NA,
#   #              color = 'black',
#   #              #color = cols[2],
#   #              inherit.aes = FALSE) +
#   # geom_polygon(data = APDPE,
#   #              aes(long, lat),
#   #              fill = NA,
#   #              color = 'black',
#   #              #color = cols[3],
#   #              inherit.aes = FALSE) +
#   # geom_polygon(data = APBSW,
#   #              aes(long, lat),
#   #              fill = NA,
#   #              color = 'black',
#   #              #color = cols[4],
#   #              inherit.aes = FALSE) +
#   # geom_polygon(data = APBSE,
#   #              aes(long, lat),
#   #              fill = NA,
#   #              color = 'black',
#   #              #color = cols[5],
#   #              inherit.aes = FALSE) +
#   # geom_polygon(data = APEI,
#   #              aes(long, lat),
#   #              fill = NA,
#   #              color = 'black',
#   #              #color = cols[6],
#   #              inherit.aes = FALSE) +
#   # geom_polygon(data = APE,
#   #              aes(long, lat),
#   #              fill = NA,
#   #              color = 'black',
#   #              #color = cols[7],
#   #              inherit.aes = FALSE) +
#   geom_point(data = colonies,
#              inherit.aes = FALSE,
#              size = 8,
#              alpha = 0.6,
#              aes(long, lat, color = KRILL)) +
#   scale_color_gradient('Commercial Krill Catch (tons)',
#                        limits = c(min(colonies$KRILL),
#                                   max(colonies$KRILL)),
#                        low = 'white',
#                        high = 'red') +
#   geom_point(data = colonies,
#              inherit.aes = FALSE,
#              size = 8,
#              shape = 21,
#              alpha = 0.6,
#              stroke = 1.2,
#              color = 'black',
#              aes(long, lat))
#   ggtitle('Average yearly Commercial Krill Catch - 2012 - 2017') +
#   geom_tile()


#SIC and krill

ggplot(mn_ay_df, aes(x=x, y=y, fill = value)) + 
  #plot settings
  geom_tile() +
  theme_void() +
  coord_equal(xlim = c(-2850000, -2050000), 
              ylim = c(1000000, 2000000)) +
  #SIC
  scale_fill_gradient2('Mean SIC - June - Sep',
                       limits = c(0,
                                  100),
                       na.value = 'white',
                       low = 'white',
                       high = 'light blue') +
  #AP shp file
  geom_polygon(data = AP, aes(long, lat, group = group),
               inherit.aes = FALSE,
               fill = 'grey') + 
  geom_path(data = AP, aes(long, lat, group = group), 
            inherit.aes = FALSE,
            color = 'black') +
  #krill
  geom_point(data = colonies,
             inherit.aes = FALSE,
             size = 8,
             alpha = 0.6,
             aes(long, lat, color = KRILL)) +
  scale_color_gradient('Commercial Krill Catch (tons)',
                       limits = c(min(colonies$KRILL),
                                  max(colonies$KRILL)),
                       low = 'white',
                       high = 'red') +
  #krill point outlines
  geom_point(data = colonies,
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.6,
             stroke = 1,
             color = 'black',
             aes(long, lat)) +
  theme(legend.position='none')





#legend KRILL
#pdf(file = 'Legend.pdf', width = 5, height = 5, useDingbats = FALSE)
plot(1:100, 1:100, type = "n", 
     ann = TRUE, xaxt = 'n', yaxt = "n", 
     bty = "n", ylab = NA, xlab = NA,
     main = paste0('Krill'))

colfunc <- colorRampPalette(c('white', 'red'))
cols <- colfunc(201)

plotrix::gradient.rect(30, 0, 70, 100, 
              col = cols, gradient = "v", 
              border = 'black')
#dev.off()


#legend SIC
#pdf(file = 'Legend.pdf', width = 5, height = 5, useDingbats = FALSE)
plot(1:100, 1:100, type = "n", 
     ann = TRUE, xaxt = 'n', yaxt = "n", 
     bty = "n", ylab = NA, xlab = NA,
     main = paste0('SIC'))

colfunc <- colorRampPalette(c('white','light blue'))
cols <- colfunc(201)

plotrix::gradient.rect(30, 0, 70, 100, 
                       col = cols, gradient = "v", 
                       border = 'black')
#dev.off()


#legend latitude
#pdf(file = 'Legend.pdf', width = 5, height = 5, useDingbats = FALSE)
plot(1:100, 1:100, type = "n", 
     ann = TRUE, xaxt = 'n', yaxt = "n", 
     bty = "n", ylab = NA, xlab = NA,
     main = paste0('Site Latitude'))

colfunc <- colorRampPalette(c('grey', 'black'))
cols <- colfunc(201)

plotrix::gradient.rect(30, 0, 70, 100, 
                       col = cols, gradient = "v", 
                       border = 'black')
#dev.off()







# #entirity of Antarctica
# ggplot(Ant, aes(long, lat, group = group)) +
#   geom_polygon(fill = 'grey') + 
#   geom_path(color = 'black') +
#   coord_equal() +
#   theme_void()


  
