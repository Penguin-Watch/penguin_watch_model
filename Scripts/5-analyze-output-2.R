#################
# Analyze covariates
#
# Author: Casey Youngflesh
#################


# Clear environment -------------------------------------------------------

rm(list = ls())



# dir ---------------------------------------------------------------------

OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2019-07-17'



# Load packages -----------------------------------------------------------

library(MCMCvis)
library(boot)
library(dplyr)



# Load data -------------------------------------------------------

setwd(OUTPUT)

#read in master data object
master_output <- readRDS('master_output.rds')



# plot BS ------------------------------------------


library(ggplot2)
p <- ggplot(mrg2, aes(YEAR, mn_mu_phi, color = SITE)) + 
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = mrg2,
                aes(ymin = LCI_mu_phi, ymax = UCI_mu_phi,
                    color = SITE), width = 1.1,
                position = position_dodge(width = 0.7)) +
  theme_bw() +
  ylim(c(0, 2)) +
  ylab('Breeding Success') +
  xlab('Year')

ggsave('site_bs.pdf', p)



# BS map ------------------------------------------------------------------

#map of breeding success
require(raster)
setwd('~/Google_Drive/R/penguin_watch_model/Data/peninsula/')
AP <- rgdal::readOGR('GADM_peninsula.shp')
setwd('~/Google_Drive/R/penguin_watch_model/Data/Sub-antarctic_coastline_low_res_polygon/')
SG <- rgdal::readOGR('Sub-antarctic_coastline_low_res_polygon.shp')
SG2 <- sp::spTransform(SG, CRS(proj4string(AP)))


mrg_agg <- aggregate(mn_mu_phi ~ SITE + col_lat + col_lon, 
                     data = mrg9, mean)

#merge with source
mrg_agg2 <- dplyr::left_join(mrg_agg, unique(mrg8[,c('SITE', 'SOURCE')]), 
                             by = c('SITE'))

#change PETE to all
PETE_idx <- which(mrg_agg2$SITE == 'PETE')
mrg_agg3 <- mrg_agg2[-PETE_idx[2:3],]
mrg_agg3[which(mrg_agg3$SITE == 'PETE'),'SOURCE'] <- 'ALL'


setwd(OUTPUT)

#all sites
pdf('BS_map.pdf')
#AP shp file
ggplot(data = AP, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey') + 
  geom_path(data = AP, aes(long, lat, group = group), 
            inherit.aes = FALSE,
            color = 'black') +
  #all coord
  coord_map(xlim = c(-68, -33),
            ylim = c(-67, -51)) +
  #theme_void() +
  theme_bw() +
  #BS
  geom_point(data = mrg_agg3,
             inherit.aes = FALSE,
             size = 8,
             alpha = 0.9,
             aes(col_lon, col_lat, color = mn_mu_phi)) +
  scale_color_gradient('Chicks per pair',
                       limits = c(min(mrg_agg3$mn_mu_phi),
                                  max(mrg_agg3$mn_mu_phi)),
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


#just AP sites
pdf('BS_map_AP.pdf')
#AP shp file
ggplot(data = AP, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey') + 
  geom_path(data = AP, aes(long, lat, group = group), 
            inherit.aes = FALSE,
            color = 'black') +
  #AP coord
  coord_map(xlim = c(-68, -46),
            ylim = c(-66.5, -59)) +
  # #all coord
  # coord_map(xlim = c(-68, -33),
  #           ylim = c(-67, -51)) +
  #theme_void() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #BS
  geom_point(data = mrg_agg3,
             inherit.aes = FALSE,
             size = 8,
             alpha = 0.9,
             aes(col_lon, col_lat, color = mn_mu_phi)) +
  scale_color_gradient('Chicks per pair',
                       limits = c(min(mrg_agg3$mn_mu_phi),
                                  max(mrg_agg3$mn_mu_phi)),
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
             #color = 'red',
             color = 'black',
             aes(col_lon, col_lat)) + 
  #point outlines - ALL THREE SOURCES
  geom_point(data = mrg_agg3[which(mrg_agg3$SOURCE == 'ALL'),],
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             #color = 'purple',
             color = 'black',
             aes(col_lon, col_lat))
#theme(legend.position='none') +
dev.off()


PT_SZ <- 40
#just SG sites
pdf('BS_map_SG.pdf')
#SG shp file
ggplot(data = SG2, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey') + 
  geom_path(data = SG2, aes(long, lat, group = group), 
            inherit.aes = FALSE,
            color = 'black') +
  # #AP coord
  # coord_map(xlim = c(-68, -53),
  #           ylim = c(-66.5, -61)) +
  #SG coord
  coord_map(xlim = c(-39.5, -34.5),
            ylim = c(-53, -55.5)) +
  # #all coord
  # coord_map(xlim = c(-68, -33),
  #           ylim = c(-67, -51)) +
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
             alpha = 0.9,
             aes(col_lon, col_lat, color = mn_mu_phi)) +
  scale_color_gradient('Chicks per pair',
                       limits = c(min(mrg_agg3$mn_mu_phi),
                                  max(mrg_agg3$mn_mu_phi)),
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
pdf('continent.pdf')
#SG shp file
ggplot(data = SG2, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey') + 
  geom_path(data = SG2, aes(long, lat, group = group), 
            inherit.aes = FALSE,
            color = 'black') +
  # #AP coord
  # coord_map(xlim = c(-68, -53),
  #           ylim = c(-66.5, -61)) +
  #SG coord
  coord_map(xlim = c(-39.5, -34.5),
            ylim = c(-53, -55.5)) +
  # #all coord
  # coord_map(xlim = c(-68, -33),
  #           ylim = c(-67, -51)) +
  #theme_void() +
  theme_bw() +
  #BS
  geom_point(data = mrg_agg3,
             inherit.aes = FALSE,
             size = 8,
             alpha = 0.9,
             aes(col_lon, col_lat, color = mn_mu_phi)) +
  scale_color_gradient('Chicks per pair',
                       limits = c(min(mrg_agg3$mn_mu_phi),
                                  max(mrg_agg3$mn_mu_phi)),
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




# BS ~ precip -------------------------------------------------------------

plot(master_output$tsnow, master_output$mn_mu_phi, 
     xlab = 'Number of large snow events at site',
     ylab = 'Breeding success')
plot(master_output$train, master_output$mn_mu_phi)



# BS ~ krill --------------------------------------------------------------

plot(master_output$krill_BR, master_output$mn_mu_phi)
plot(master_output$krill_WS, master_output$mn_mu_phi)
plot(master_output$krill_AY, master_output$mn_mu_phi)



# BS ~ tourism ------------------------------------------------------------

#
#




