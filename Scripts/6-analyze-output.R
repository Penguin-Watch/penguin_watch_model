#################
# Analyze results
#
# Author: Casey Youngflesh
#################


# Clear environment -------------------------------------------------------

rm(list = ls())


# dir ---------------------------------------------------------------------


dir <- '~/Google_Drive/R/penguin_watch_model/'
OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2020-03-16'


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(dplyr)
library(ggplot2)
library(sp)
library(rgdal)
library(rjags)


# Load data -------------------------------------------------------

setwd(OUTPUT)

#read in master data object
master_output <- readRDS('master_output.rds')
sy_chicks_rep <- readRDS('sy_chicks_rep.rds')
data <- readRDS('jagsData.rds')


# posterior predictive check ----------------------------------------------

#data generated from model
sy_chicks_rep_ch <- MCMCvis::MCMCchains(sy_chicks_rep)

#observed data
sy <- sum(data$y, na.rm = TRUE)


#create dir for PPC figs if doesn't exist
ifelse(!dir.exists(paste0(OUTPUT, '/PPC_plots')), 
       dir.create(paste0(OUTPUT, '/PPC_plots')), 
       FALSE)

setwd(paste0(OUTPUT, '/PPC_plots'))

cn <- colnames(sy_chicks_rep_ch)
ppp <- rep(NA, length(cn))
for (i in 1:length(cn))
{
  #i <- 8
  #get row/col
  sp1 <- strsplit(cn[i], split = '[', fixed = TRUE)[[1]][2]
  sp2 <- strsplit(sp1, split = ']', fixed = TRUE)[[1]][1]
  sp3 <- as.numeric(strsplit(sp2, split = ',', fixed = TRUE)[[1]])
  
  #total chicks for y for that site/year
  tsy <- sum(data$y[,,sp3[1],sp3[2]], na.rm = TRUE)
  #total for yrep for that site/year
  tsyr <- sy_chicks_rep_ch[,i]
  
  YEAR <- data$yrs_array[sp3[1], sp3[2]]
  SITE <- data$unsites[sp3[2]]
  
  #ppcheck
  ppp[i] <- sum(tsyr > tsy) / length(tsyr)
  txt <- paste0('Bayes p: ', round(ppp[i], 2))
  
  jpeg(paste0(SITE, '-', YEAR, '-2020-03-13.jpg'), 
       width = 5, height = 5, units = 'in', res = 300)
  hist(tsyr, main = paste0(SITE, ' - ', 
                           YEAR-1, '-', YEAR),
       xlab = 'Simulated number of chicks') 
  abline(v = tsy, col = 'red', lwd = 3)
  legend('topleft',legend = txt, bty ="n", pch=NA)
  dev.off()
}


# plot BS ------------------------------------------

#just penguin watch BS estimates
PW_output <- dplyr::filter(master_output, SOURCE == 'PW')

min(PW_output$mn_bs)
max(PW_output$mn_bs)
mean(PW_output$mn_bs)

setwd(OUTPUT)

pdf('site_bs.pdf')
ggplot(PW_output, aes(YEAR, mn_bs, color = SITE)) + 
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = PW_output,
                aes(ymin = (mn_bs - sd_bs), ymax = (mn_bs + sd_bs),
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
Ant_3031 <- sp::spTransform(Ant, CRS("+init=epsg:3031"))
Ant_SA <- raster::union(Ant_3031, SA)


#mean BS at each site
mrg_agg <- aggregate(mn_bs ~ SITE + col_lat + col_lon, 
                     data = PW_output, mean)

#merge with source
mrg_agg2 <- dplyr::left_join(mrg_agg, unique(PW_output[,c('SITE', 'SOURCE')]), 
                             by = c('SITE'))

#change PETE to all
# PETE_idx <- which(mrg_agg2$SITE == 'PETE')
# mrg_agg3 <- mrg_agg2[-PETE_idx[2:3],]
# mrg_agg3[which(mrg_agg3$SITE == 'PETE'),'SOURCE'] <- 'ALL'


setwd(OUTPUT)

#BLACK circles = Mean breeding success for site comes from this study

#how transparent points are when plotted
ALPHA_PT <- 0.9

#just AP sites
pdf('BS_map_AP.pdf', height = 6, width = 6)
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
  geom_point(data = mrg_agg2,
             inherit.aes = FALSE,
             size = 8,
             alpha = ALPHA_PT,
             aes(col_lon, col_lat, color = mn_bs)) +
  scale_color_gradient('Chicks per pair',
                       limits = c(min(mrg_agg2$mn_bs),
                                  max(mrg_agg2$mn_bs)),
                       low = '#2c7fb8',
                       high = '#edf8b1') +
  # #point outlines
  geom_point(data = mrg_agg2,
             inherit.aes = FALSE,
             size = 8,
             shape = 21,
             alpha = 0.8,
             stroke = 1,
             color = 'black',
             aes(col_lon, col_lat)) #+
  #theme(legend.position='none')
dev.off()


PT_SZ <- 30
#just SG sites
pdf('BS_map_SG.pdf', width = 6, height = 6)
#SG shp file
ggplot(data = AP_SA_ll, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey') + 
  geom_path(data = AP_SA_ll, aes(long, lat, group = group), 
            inherit.aes = FALSE,
            color = 'black',
            size = 1.2) +
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
  geom_point(data = mrg_agg2,
             inherit.aes = FALSE,
             size = PT_SZ,
             alpha = ALPHA_PT,
             aes(col_lon, col_lat, color = mn_bs)) +
  scale_color_gradient('Chicks per pair',
                       limits = c(min(mrg_agg2$mn_bs),
                                  max(mrg_agg2$mn_bs)),
                       low = '#2c7fb8',
                       high = '#edf8b1') +
  #point outlines - PW
  geom_point(data = mrg_agg2,
             inherit.aes = FALSE,
             size = PT_SZ,
             shape = 21,
             alpha = 0.8,
             stroke = 3.5,
             color = 'black',
             aes(col_lon, col_lat)) +
  theme(legend.position='none')
dev.off()


#inset
pdf('Ant_SA.pdf')
#SG shp file
ggplot(data = Ant_SA, aes(long, lat, group = group)) +
  geom_polygon(fill = 'grey') + 
  geom_path(data = Ant_SA, aes(long, lat, group = group), 
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


# functions to fit models and plot results --------------------------------

fit_env_fun <- function(ENV, n_adapt, n_burn, n_draw, n_thin)
{
  DATA <- list(y_obs = PW_output$mn_bs,
               sigma_y = PW_output$sd_bs + 0.000001,
               x = scale(ENV, scale = FALSE)[,1],
               N = length(PW_output$mn_bs))
  
  {
    sink('bs_env.jags')
    
    cat("
        model {
      
        #site
        for (i in 1:N)
        {
          y_obs[i] ~ dnorm(y_true[i], pow(sigma_y[i], -2))
          mu[i] = alpha + beta * x[i]
          y_true[i] ~ dnorm(mu[i], tau)
        } 
  
        alpha ~ dnorm(0, 10)
        beta ~ dnorm(0, 1)
        tau <- pow(sigma, -2) 
        sigma ~ dunif(0, 5)
        
        }",fill = TRUE)
    
    sink()
  }
  
  jm <- rjags::jags.model(data = DATA,
                          file = 'bs_env.jags',
                          n.chains = 4,
                          n.adapt = n_adapt)
  
  stats::update(jm, n.iter = n_burn)
  
  samples <- rjags::coda.samples(jm,
                                 n.iter = n_draw,
                                 variable.names = c('alpha',
                                                    'beta',
                                                    'sigma',
                                                    'y_true'),
                                 thin = n_thin)
  
  file.remove('bs_env.jags')
  
  return(samples)
}


plot_env_fun <- function(ENV, samples, XLAB, YLAB, obj)
{
  ENVsc <- scale(ENV, scale = FALSE)
  sim_x <- seq(min(ENVsc), max(ENVsc), length = 100)
  
  alpha_ch <- MCMCvis::MCMCchains(samples, params = 'alpha')
  beta_ch <- MCMCvis::MCMCchains(samples, params = 'beta')
  
  mf <- matrix(nrow = length(alpha_ch), ncol = 100)
  for (i in 1:length(sim_x))
  {
    mf[,i] <- alpha_ch + beta_ch * sim_x[i]
  }
  
  med_mf <- apply(mf, 2, median)
  LCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
  UCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.975))
  
  FIT_PLOT <- data.frame(MN = med_mf,
                         MN_X = sim_x + attr(ENVsc, 'scaled:center'),
                         LCI = LCI_mf,
                         UCI = UCI_mf)
  
  yt_mn <- MCMCvis::MCMCpstr(samples, params = 'y_true', func = mean)[[1]]
  yt_sd <- MCMCvis::MCMCpstr(samples, params = 'y_true', func = sd)[[1]]
  
  DATA_PLOT2 <- data.frame(x = ENVsc + attr(ENVsc, 'scaled:center'),
                           y = yt_mn)
  
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
}


# BS ~ precip -------------------------------------------------------------

#total # precip events
sum_precip <- PW_output$train + PW_output$tsnow

#fit model
fit_precip <- fit_env_fun(ENV = sum_precip, 
                          n_adapt = 5000, 
                          n_burn = 10000, 
                          n_draw = 10000, 
                          n_thin = 1)

#analyze output
MCMCvis::MCMCsummary(fit_precip, round = 3)

#plot results
plot_env_fun(sum_precip, fit_precip, 
             XLAB = 'Total # of precipitation events', 
             YLAB = 'Breeding Success (chicks / pair)',
             obj = 'bs_precip.jpg')


# BS ~ krill --------------------------------------------------------------

#fit model
fit_krill <- fit_env_fun(log(PW_output$krill_WS), 
                         n_adapt = 5000, 
                         n_burn = 10000, 
                         n_draw = 10000, 
                         n_thin = 1)

#analyze output
MCMCvis::MCMCsummary(fit_krill, round = 3)

#plot results
plot_env_fun(log(PW_output$krill_WS), fit_krill, 
             XLAB = 'log(Krill catch) (tons)', 
             YLAB = 'Breeding Success (chicks / pair)',
             obj = 'bs_krill.jpg')


# BS ~ tourism ------------------------------------------------------------

#fit model
fit_tourism <- fit_env_fun((PW_output$t_visitors/1000), 
                           n_adapt = 5000, 
                           n_burn = 10000, 
                           n_draw = 10000, 
                           n_thin = 1)

#analyze output
MCMCvis::MCMCsummary(fit_tourism, round = 3)

#plot results
plot_env_fun((PW_output$t_visitors/1000), fit_tourism, 
             XLAB = 'Thousands of visitors',
             YLAB = 'Breeding Success (chicks / pair)',
             obj = 'bs_tourism.jpg')

