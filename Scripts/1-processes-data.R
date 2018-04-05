#################
# Penguin Watch Model - 1 - Process data
#
# 1 - process data
# 1 - simulate data
# 2 - run model
# 3 - analyze output

#
#Created: Apr 5, 2018
#################



# Clear environment -------------------------------------------------------

rm(list = ls())


# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman', repos = "http://cran.case.edu")
}
pacman::p_load(dplyr, rgdal, rgeos)



# setwd -------------------------------------------------------------------

setwd('Data/Site_data/')


# penguin data -----------------------------------------------------------



# krill data --------------------------------------------------------------

#150 km buffer around each site
sites <- rgdal::readOGR('SitesEPSG3031.shp')

#units for projection are in m
#proj4string(sites)

#create 150km buffer (150,000m) around each site
site_buffer <- rgeos::gBuffer(sites, width = 150000)
#plot(site_buffer)





