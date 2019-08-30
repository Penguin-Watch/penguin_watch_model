#################
# Process tourism data
#
# From IAATO
#
# Author: Casey Youngflesh
#################



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

mrg <- readRDS('bs_precip_krill_mrg.rds')

#read in model input data
data <- readRDS('jagsData.rds')

#start date of season (creche - 60 days)
start <- as.numeric(format(as.Date(as.numeric(julian(mrg$creche_date)) - 59, 
                          origin = '1970-01-01'), '%j'))


# read in tourism data ----------------------------------------------------

setwd(paste0(dir, 'Data/tourism_data/IAATO_site_visits_reports'))

fls <- list.files()

#get dates for each site/year (lay to fledge)
#sum total visitors for each site/year based on those dates
#plot # of visitors each day over course of season
for (i in 1:length(fls))
{
  #i <- 3
  temp <- read.csv(fls[i])
  temp$total_visiting
}


bs_precip_mrg <- readRDS('bs_precip_mrg.rds')
