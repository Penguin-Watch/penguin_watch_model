#################
# Analyze krill
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

setwd(paste0(dir, 'Data/Krill_data/CCAMLR/Processed_CCAMLR'))

#read in weighted krill catch
krill_weighted_150 <- readRDS('krill_weighted_150.rds')

#read in average (across years) kril catch over entire season
CCAMLR_kr_AY <- read.csv('CCAMLR_krill_average.csv')


# krill figs --------------------------------------------------------------

# which sites have the most fishing at them?

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
