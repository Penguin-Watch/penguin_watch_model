#################
# Analyze results
#
# Author: Casey Youngflesh
#################


# Clear environment -------------------------------------------------------

rm(list = ls())



# dir ---------------------------------------------------------------------


dir <- '~/Google_Drive/R/penguin_watch_model/'
OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2019-07-17'


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(boot)
library(dplyr)
library(sp)
library(rgdal)
library(ggplot2)


# Load data -------------------------------------------------------

setwd(OUTPUT)

#read in master data object
master_output <- readRDS('master_output.rds')



# plot BS ------------------------------------------

#just penguin watch BS estimates
PW_output <- dplyr::filter(master_output, SOURCE == 'PW')

setwd(OUTPUT)

pdf('site_bs.pdf')
ggplot(PW_output, aes(YEAR, mn_bs, color = SITE)) + 
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = PW_output,
                aes(ymin = LCI_bs, ymax = UCI_bs,
                    color = SITE), width = 1.1,
                position = position_dodge(width = 0.7)) +
  theme_bw() +
  ylim(c(0, 2)) +
  ylab('Breeding Success') +
  xlab('Year')
dev.off()




# BS map ------------------------------------------------------------------

#map of breeding success

#read in maps of AP and sub Antarctic
# #AP
setwd(paste0(dir, 'Data/peninsula'))
AP_p <- rgdal::readOGR('GADM_peninsula.shp')
AP <- sp::spTransform(AP_p, CRS("+init=epsg:3031"))
#Sub Antarctic
setwd(paste0(dir, 'Data/Sub-antarctic_coastline_low_res_polygon'))
SA_p <- rgdal::readOGR('Sub-antarctic_coastline_low_res_polygon.shp')
SA <- sp::spTransform(SA_p, CRS("+init=epsg:3031"))
AP_SA <- raster::union(AP, SA)
#transform to 4326
AP_SA_ll <- sp::spTransform(AP_SA, CRS("+init=epsg:4326"))
#Antarctic continent
setwd(paste0(dir, 'Data/Coastline_low_res_polygon'))
Ant <- rgdal::readOGR('Coastline_low_res_polygon.shp')



#mean BS at each site
mrg_agg <- aggregate(mn_bs ~ SITE + col_lat + col_lon, 
                     data = master_output, mean)

#merge with source
mrg_agg2 <- dplyr::left_join(mrg_agg, unique(master_output[,c('SITE', 'SOURCE')]), 
                             by = c('SITE'))

#change PETE to all
PETE_idx <- which(mrg_agg2$SITE == 'PETE')
mrg_agg3 <- mrg_agg2[-PETE_idx[2:3],]
mrg_agg3[which(mrg_agg3$SITE == 'PETE'),'SOURCE'] <- 'ALL'


setwd(OUTPUT)

#BLACK circles = Mean breeding success for site comes from this study
#RED circles = Mean breeding success for site comes from Hinke et al. 2018
#PURPE cirlces = Mean breeding success for site comes from this study, Hinke et al. 2018, and Lynch et al. 2009 (three studies)

#how transparent points are when plotted
ALPHA_PT <- 0.9


#all sites
# pdf('BS_map.pdf')
# #AP shp file
# ggplot(data = AP_SA_ll, aes(long, lat, group = group)) +
#   geom_polygon(fill = 'grey') + 
#   geom_path(data = AP_SA_ll, aes(long, lat, group = group), 
#             inherit.aes = FALSE,
#             color = 'black') +
#   #lat/lon limits
#   coord_map(xlim = c(-70, -30),
#             ylim = c(-67, -51)) +
#   #theme_void() +
#   theme_bw() +
#   #BS
#   geom_point(data = mrg_agg3,
#              inherit.aes = FALSE,
#              size = 8,
#              alpha = ALPHA_PT,
#              aes(col_lon, col_lat, color = mn_bs)) +
#   scale_color_gradient('Chicks per pair',
#                        limits = c(min(mrg_agg3$mn_bs),
#                                   max(mrg_agg3$mn_bs)),
#                        low = '#2c7fb8',
#                        high = '#edf8b1') +
#   # #point outlines
#   # geom_point(data = mrg_agg3,
#   #            inherit.aes = FALSE,
#   #            size = 8,
#   #            shape = 21,
#   #            alpha = 0.8,
#   #            stroke = 1,
#   #            color = 'black',
#   #            aes(col_lon, col_lat)) +
#   #point outlines - PW
#   geom_point(data = mrg_agg3[which(mrg_agg3$SOURCE == 'PW'),],
#              inherit.aes = FALSE,
#              size = 8,
#              shape = 21,
#              alpha = 0.8,
#              stroke = 1,
#              color = 'black',
#              aes(col_lon, col_lat)) +
#   #point outlines - Hinke
#   geom_point(data = mrg_agg3[which(mrg_agg3$SOURCE == 'Hinke'),],
#              inherit.aes = FALSE,
#              size = 8,
#              shape = 21,
#              alpha = 0.8,
#              stroke = 1,
#              color = 'red',
#              aes(col_lon, col_lat)) +
#   #point outlines - ALL THREE SOURCES
#   geom_point(data = mrg_agg3[which(mrg_agg3$SOURCE == 'ALL'),],
#              inherit.aes = FALSE,
#              size = 8,
#              shape = 21,
#              alpha = 0.8,
#              stroke = 1,
#              color = 'purple',
#              aes(col_lon, col_lat))
# #theme(legend.position='none') +
# dev.off()



#just AP sites
pdf('BS_map_AP.pdf')
#AP shp file
ggplot(data = AP_SA_ll, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey') + 
  geom_path(data = AP_SA_ll, aes(long, lat, group = group), 
            inherit.aes = FALSE,
            color = 'black') +
  #Just AP lat/lon
  coord_map(xlim = c(-68, -46),
            ylim = c(-66.5, -59)) +
  #theme_void() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #BS
  geom_point(data = mrg_agg3,
             inherit.aes = FALSE,
             size = 8,
             alpha = ALPHA_PT,
             aes(col_lon, col_lat, color = mn_bs)) +
  scale_color_gradient('Chicks per pair',
                       limits = c(min(mrg_agg3$mn_bs),
                                  max(mrg_agg3$mn_bs)),
                       low = '#2c7fb8',
                       high = '#edf8b1') +
  # #point outlines
  # geom_point(data = mrg_agg3,
  #            inherit.aes = FALSE,
  #            size = 8,
  #            shape = 21,
  #            alpha = 0.8,
  #            stroke = 1,
  #            color = 'black',
  #            aes(col_lon, col_lat)) +
  #point outlines - PW
  geom_point(data = mrg_agg3[which(mrg_agg3$SOURCE == 'PW'),],
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'black',
             aes(col_lon, col_lat)) +
  #point outlines - Hinke
  geom_point(data = mrg_agg3[which(mrg_agg3$SOURCE == 'Hinke'),],
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'red',
             #color = 'black',
             aes(col_lon, col_lat)) + 
  #point outlines - ALL THREE SOURCES
  geom_point(data = mrg_agg3[which(mrg_agg3$SOURCE == 'ALL'),],
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'purple',
             #color = 'black',
             aes(col_lon, col_lat))
#theme(legend.position='none') +
dev.off()


PT_SZ <- 40
#just SG sites
pdf('BS_map_SG.pdf')
#SG shp file
ggplot(data = AP_SA_ll, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey') + 
  geom_path(data = AP_SA_ll, aes(long, lat, group = group), 
            inherit.aes = FALSE,
            color = 'black') +
  #just SG lat/lon
  coord_map(xlim = c(-39.5, -34.5),
            ylim = c(-53, -55.5)) +
  #theme_void() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(-39, -37, -35)) +
  scale_y_continuous(breaks = c(-53, -54, -55)) +
  #BS
  geom_point(data = mrg_agg3,
             inherit.aes = FALSE,
             size = PT_SZ,
             alpha = ALPHA_PT,
             aes(col_lon, col_lat, color = mn_bs)) +
  scale_color_gradient('Chicks per pair',
                       limits = c(min(mrg_agg3$mn_bs),
                                  max(mrg_agg3$mn_bs)),
                       low = '#2c7fb8',
                       high = '#edf8b1') +
  #point outlines - PW
  geom_point(data = mrg_agg3[which(mrg_agg3$SOURCE == 'PW'),],
             inherit.aes = FALSE,
             size = PT_SZ,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'black',
             aes(col_lon, col_lat)) +
  #point outlines - Hinke
  geom_point(data = mrg_agg3[which(mrg_agg3$SOURCE == 'Hinke'),],
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'red',
             aes(col_lon, col_lat)) + 
  #point outlines - ALL THREE SOURCES
  geom_point(data = mrg_agg3[which(mrg_agg3$SOURCE == 'ALL'),],
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'purple',
             aes(col_lon, col_lat))
#theme(legend.position='none') +
dev.off()


#inset
pdf('Ant_continent.pdf')
#SG shp file
ggplot(data = Ant, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey') + 
  geom_path(data = Ant, aes(long, lat, group = group), 
            inherit.aes = FALSE,
            color = 'black') +
  theme_void()
  #theme_bw()
#theme(legend.position='none') +
dev.off()




# just PW sites, no BS -------------------------------------------------------

uPW <- unique(PW_output[,c('SITE', 'col_lat', 'col_lon')])
# uH <- dplyr::filter(master_output, SOURCE == 'Hinke')
# uL <- dplyr::filter(master_output, SOURCE == 'Lynch')

#all Penguin Watch sites
pdf('PW_site_map.pdf')
#AP shp file
ggplot(data = AP_SA_ll, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey') +
  geom_path(data = AP_SA_ll, aes(long, lat, group = group),
            inherit.aes = FALSE,
            color = 'black') +
  #lat/lon limits
  coord_map(xlim = c(-70, -30),
            ylim = c(-67, -51)) +
  theme_void() +
  #theme_bw() +
  #BS
  geom_point(data = uPW,
             inherit.aes = FALSE,
             size = 6,
             alpha = 1,
             aes(col_lon, col_lat, color = 'red')) +
  #point outlines
  geom_point(data = uPW,
             inherit.aes = FALSE,
             size = 6,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'black',
             aes(col_lon, col_lat)) +
  theme(legend.position='none')
dev.off()





# BS ~ precip -------------------------------------------------------------

plot(master_output$tsnow, master_output$mn_bs, 
     xlab = 'Number of large snow events at site',
     ylab = 'Breeding success')
plot(master_output$train, master_output$mn_bs,
     xlab = 'Number of large rain events at site',
     ylab = 'Breeding success')




# BS ~ krill --------------------------------------------------------------

#krill caught during breeding season
plot(master_output$krill_BR, master_output$mn_bs)
#krill caught during entire previous year (most krill is caught during Austral winter)
plot(master_output$krill_WS, master_output$mn_bs)
#average krill catch across years
plot(master_output$krill_AY, master_output$mn_bs)



# BS ~ tourism ------------------------------------------------------------

#need data
#
#




