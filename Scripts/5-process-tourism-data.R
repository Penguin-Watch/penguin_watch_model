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
OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2020-03-16'


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


# read in tourism data ----------------------------------------------------

setwd(paste0(dir, 'Data/tourism_data/IAATO_site_visits_reports-2019-09-02/csv'))

fls <- list.files()

#get dates for each site/year (lay to fledge)
#sum total visitors for each site/year based on those dates
#plot # of visitors each day over course of season
t_df <- data.frame(SITE = rep(NA, length(fls)), YEAR = NA, t_visitors = NA)
for (i in 1:length(fls))
{
  #i <- 5
  temp <- read.csv(fls[i])

  nm <- strsplit(fls[i], '_')[[1]]
  t_df$SITE[i] <- nm[1]
  t_df$YEAR[i] <- as.numeric(strsplit(nm[2], split = '.', fixed = TRUE)[[1]][1])
  
  tmrg <- dplyr::filter(mrg, SITE == t_df$SITE[i], YEAR == t_df$YEAR[i])
  
  #visit dates
  v_dates <- as.numeric(format(as.Date(as.numeric(julian(as.Date(temp$visit_date))), 
                                       origin = '1970-01-01'), '%j'))

  #start/end of season (creche - 60 days)
  start <- as.numeric(format(as.Date(as.numeric(julian(tmrg$creche_date)) - 60, 
                                     origin = '1970-01-01'), '%j'))
  end <- as.numeric(format(as.Date(as.numeric(julian(tmrg$creche_date)), 
                                   origin = '1970-01-01'), '%j'))
  
  #which visits are in valid period
  v_ind <- which(v_dates >= start | v_dates <= end)
  
  t_df$t_visitors[i] <- sum(temp$total_visiting[v_ind], na.rm = TRUE)
}


# add SG tourism data -----------------------------------------------------

setwd(paste0(dir, 'Data/tourism_data/SG_tourism_2019-09-29/csv'))

fls <- list.files()

#get dates for each site/year (lay to fledge)
#sum total visitors for each site/year based on those dates
#plot # of visitors each day over course of season
t_df_SG <- data.frame(SITE = rep(NA, length(fls)), YEAR = NA, t_visitors = NA)
for (i in 1:length(fls))
{
  #i <- 1
  temp <- read.csv(fls[i])
  
  nm <- strsplit(fls[i], '_')[[1]]
  t_df_SG$SITE[i] <- nm[1]
  t_df_SG$YEAR[i] <- as.numeric(strsplit(nm[2], split = '.', fixed = TRUE)[[1]][1])
  
  tmrg <- dplyr::filter(mrg, SITE == t_df_SG$SITE[i], YEAR == t_df_SG$YEAR[i])
  
  #visit dates
  v_dates <- as.numeric(format(as.Date(as.numeric(julian(as.Date(temp$DATE))), 
                                       origin = '1970-01-01'), '%j'))
  
  #start/end of season (creche - 60 days)
  start <- as.numeric(format(as.Date(as.numeric(julian(tmrg$creche_date)) - 60, 
                                     origin = '1970-01-01'), '%j'))
  end <- as.numeric(format(as.Date(as.numeric(julian(tmrg$creche_date)), 
                                   origin = '1970-01-01'), '%j'))
  
  #which visits are in valid period
  v_ind <- which(v_dates >= start | v_dates <= end)
  
  t_df_SG$t_visitors[i] <- sum(c(temp$PASSENGERS[v_ind],
                              temp$STAFF[v_ind],
                              temp$CREW[v_ind],
                              temp$OTHERS[v_ind]),
                              na.rm = TRUE)
}

#combine AP and SG tourism data
t_df_comb <- rbind(t_df, t_df_SG)


# merge with rest of data -------------------------------------------------

setwd(OUTPUT)

mrg2 <- left_join(mrg, t_df_comb, by = c('SITE', 'YEAR'))

#save as rds
saveRDS(mrg2, 'master_output.rds')

