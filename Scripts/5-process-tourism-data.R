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

mrg_p <- readRDS('bs_precip_krill_mrg.rds')

mrg <- dplyr::filter(mrg_p, SOURCE == 'PW')

#read in model input data
data <- readRDS('jagsData.rds')

#start date of season (creche - 60 days)
start <- as.numeric(format(as.Date(as.numeric(julian(mrg$creche_date)) - 59, 
                          origin = '1970-01-01'), '%j'))


# read in tourism data ----------------------------------------------------

setwd(paste0(dir, 'Data/tourism_data/IAATO_site_visits_reports-2019-09-02/csv'))

fls <- list.files()

#get dates for each site/year (lay to fledge)
#sum total visitors for each site/year based on those dates
#plot # of visitors each day over course of season
t_df <- data.frame(SITE = rep(NA, length(fls)), YEAR = NA, t_visitors = NA)
for (i in 1:length(fls))
{
  #i <- 4
  temp <- read.csv(fls[i])
  
  nm <- strsplit(fls[i], '_')[[1]]
  t_df$SITE[i] <- nm[1]
  t_df$YEAR[i] <- as.numeric(strsplit(nm[2], split = '.', fixed = TRUE)[[1]][1])
  t_df$t_visitors[i] <- sum(temp$total_visiting, na.rm = TRUE)
}



# merge with rest of data -------------------------------------------------

setwd(OUTPUT)

mrg2 <- left_join(mrg, t_df, by = c('SITE', 'YEAR'))

#save as rds
saveRDS(mrg2, 'master_output.rds')

