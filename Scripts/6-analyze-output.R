#################
# Analyze results
#
# Author: Casey Youngflesh
#################


# Clear environment -------------------------------------------------------

rm(list = ls())



# dir ---------------------------------------------------------------------


dir <- '~/Google_Drive/R/penguin_watch_model/'
OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2019-09-08'


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(dplyr)
library(ggplot2)
library(sp)
library(rgdal)
library(boot)
library(rstanarm)


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

master2 <- dplyr::filter(master_output, SOURCE == 'PW')
sum_precip <- master2$train + master2$tsnow

rs_fun <- function(y, x, XLAB, YLAB, obj)
{
fit <- rstanarm::stan_glm(y ~ x, 
                          chains = 4)

alpha_ch <- c()
beta_ch <- c()
for (i in 1:4)
{
  a1 <- fit$stanfit@sim$samples[[i]]$alpha
  b1 <- fit$stanfit@sim$samples[[i]]$beta
  alpha_ch <- c(alpha_ch, a1)
  beta_ch <- c(beta_ch, b1)
}

sim_x <- seq(min(x), max(x), length = 100)
mf <- matrix(nrow = length(alpha_ch), ncol = 100)
for (i in 1:length(sim_x))
{
  mf[,i] <- alpha_ch + beta_ch * sim_x[i]
}

med_mf <- apply(mf, 2, median)
LCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
UCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.975))

FIT_PLOT <- data.frame(MN = med_mf,
                       MN_X = sim_x,
                       LCI = LCI_mf,
                       UCI = UCI_mf)

DATA_PLOT2 <- data.frame(x = x,
                         y = y)

p <- ggplot(data = DATA_PLOT2, aes(x, y)) +
  #model fit
  geom_ribbon(data = FIT_PLOT,
              aes(x = MN_X, ymin = LCI, ymax = UCI),
              fill = 'grey', alpha = 0.7,
              inherit.aes = FALSE) +
  geom_line(data = FIT_PLOT, aes(MN_X, MN), color = 'red',
            alpha = 0.9,
            inherit.aes = FALSE,
            size = 1.4) +
  #latent state
  geom_point(data = DATA_PLOT2, aes(x, y), color = 'black',
             inherit.aes = FALSE, size = 3, alpha = 0.3) +
  theme_bw() +
  #scale_x_discrete(limits = c(seq(18,30, by = 2))) +
  ylab(YLAB) +
  xlab(XLAB) +
  theme(
    plot.title = element_text(size = 22),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
    axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
ggsave(p, filename = paste0(obj), width = 5, height = 5)

return(fit)
}

fit1 <- rs_fun(y = master2$mn_bs, 
               x = sum_precip, 
               XLAB = 'Total # of precipitation events', 
               YLAB = 'Breeding Success (chicks / pair)',
               obj = 'bs_precip.jpg')

MCMCvis::MCMCsummary(fit1, params = 'x')


# BS ~ krill --------------------------------------------------------------

#krill caught in previous year

fit2 <- rs_fun(y = master2$mn_bs, 
               x = log(master2$krill_WS), XLAB = 'log(Krill catch) (tons)', 
               YLAB = 'Breeding Success (chicks / pair)',
               obj = 'bs_krill.jpg')

MCMCvis::MCMCsummary(fit2, params = 'x')



# BS ~ tourism ------------------------------------------------------------

#filter SG (don't have data for sites there yet)
master3 <- dplyr::filter(master2, t_visitors > 1000)

fit3 <- rs_fun(y = master3$mn_bs,
               x = (master3$t_visitors/1000), XLAB = 'Thousands of visitors',
               YLAB = 'Breeding Success (chicks / pair)',
               obj = 'bs_tourism.jpg')

MCMCvis::MCMCsummary(fit3, params = 'x')




# BS ~ precip + krill + tourism ------------------------------------------------------------

# z_idx <- which(master2$t_visitors < 1000)
# fit4 <- rstanarm::stan_glm(master3$mn_bs ~ log(master3$t_visitors) + 
#                             sum_precip[-z_idx] + 
#                              log(master3$krill_WS), 
#                           chains = 4)
# 
# MCMCsummary(fit4, round = 7)

