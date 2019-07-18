#################
# Analyze covariates
#
# Author: Casey Youngflesh
#################


# Clear environment -------------------------------------------------------

rm(list = ls())



# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman')
}

pacman::p_load(MCMCvis, boot, dplyr)




# Load data -------------------------------------------------------

#read in RDS from 3-analyze-output.R
mrg6 <- readRDS('mrg6.rds')

#merge with krill data
setwd('~/Google_Drive/R/penguin_watch_model/Data/Krill_data/CCAMLR/Processed_CCAMLR/')

#krill catch from entire year
krill <- read.csv('CCAMLR_krill_entire_season.csv')
colnames(krill)[3] <- 'YR_KRILL'
mrg7 <- dplyr::left_join(mrg6, krill, by = c('SITE', 'YEAR'))

#krill catch from just breeding season
krill2 <- read.csv('CCAMLR_krill_breeding_season.csv')
colnames(krill2)[3] <- 'BR_KRILL'
mrg8 <- dplyr::left_join(mrg7, krill2, by = c('SITE', 'YEAR'))



# covariates/plots --------------------------------------------------------------

#compare BS to total precip - no relationship
plot(mrg8$tsnow, mrg8$mn_mu_phi)
plot(mrg8$train, mrg8$mn_mu_phi)


#compare BS to phenology
#convert creche date to julian day
jd_1 <- as.numeric(strftime(mrg8$creche_date, format = "%j"))
idx <- which(jd_1 > 300)
jd_2 <- jd_1
jd_2[idx] <- jd_1[idx] - 300
jd_2[-idx] <- jd_1[-idx] + 65

mrg8$c_jd <- jd_2

plot(mrg8$c_jd, mrg8$mn_mu_phi)
plot(mrg8$col_lat, mrg8$c_jd)


#compare BS to lat and substrate type
# to.rm <- which(mrg4$col_lat > -60)
# mrg5 <- mrg4[-to.rm,]
# plot(mrg5$col_lat, mrg5$mn_mu_phi)
# f1 <- lm(mn_mu_phi ~ col_lat, data = mrg5)
# abline(f1, col = 'red')

#compare BS to krill
library(quantreg)
plot(mrg8$YR_KRILL, mrg8$mn_mu_phi)



# BS map ------------------------------------------------------------------

#map of breeding success
require(raster)
setwd('~/Google_Drive/R/penguin_watch_model/Data/peninsula/')
AP <- rgdal::readOGR('GADM_peninsula.shp')

mrg_agg <- aggregate(mn_mu_phi ~ SITE + col_lat + col_lon + SOURCE, data = mrg8, mean)

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
  theme_void() +
  #theme_bw()
  #BS
  geom_point(data = mrg_agg,
             inherit.aes = FALSE,
             size = 8,
             alpha = 0.9,
             aes(col_lon, col_lat, color = mn_mu_phi)) +
  scale_color_gradient('Chicks per pair',
                       limits = c(min(mrg_agg$mn_mu_phi),
                                  max(mrg_agg$mn_mu_phi)),
                       low = '#2c7fb8',
                       high = '#edf8b1') +
  # #point outlines
  # geom_point(data = mrg_agg,
  #            inherit.aes = FALSE,
  #            size = 8,
  #            shape = 21,
  #            alpha = 0.8,
  #            stroke = 1,
  #            color = 'black',
  #            aes(col_lon, col_lat)) +
  #point outlines - PW
  geom_point(data = mrg6_agg[which(mrg_agg$SOURCE == 'PW'),],
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'black',
             aes(col_lon, col_lat)) +
  #point outlines - Hinke
  geom_point(data = mrg_agg[which(mrg_agg$SOURCE == 'Hinke'),],
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'red',
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
  coord_map(xlim = c(-68, -53),
            ylim = c(-66.5, -61)) +
  # #all coord
  # coord_map(xlim = c(-68, -33),
  #           ylim = c(-67, -51)) +
  theme_void() +
  #theme_bw()
  #BS
  geom_point(data = mrg_agg,
             inherit.aes = FALSE,
             size = 8,
             alpha = 0.9,
             aes(col_lon, col_lat, color = mn_mu_phi)) +
  scale_color_gradient('Chicks per pair',
                       limits = c(min(mrg_agg$mn_mu_phi),
                                  max(mrg_agg$mn_mu_phi)),
                       low = '#2c7fb8',
                       high = '#edf8b1') +
  # #point outlines
  # geom_point(data = mrg_agg,
  #            inherit.aes = FALSE,
  #            size = 8,
  #            shape = 21,
  #            alpha = 0.8,
  #            stroke = 1,
  #            color = 'black',
  #            aes(col_lon, col_lat)) +
  #point outlines - PW
  geom_point(data = mrg_agg[which(mrg_agg$SOURCE == 'PW'),],
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'black',
             aes(col_lon, col_lat)) +
  #point outlines - Hinke
  geom_point(data = mrg_agg[which(mrg_agg$SOURCE == 'Hinke'),],
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'red',
             aes(col_lon, col_lat))
#theme(legend.position='none') +
dev.off()


