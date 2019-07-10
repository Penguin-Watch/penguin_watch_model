#MODIFIED FROM MAPPPD REPO - written by Chris Che-Castaldo (https://github.com/CCheCastaldo)

# ________________________________________________________________________________
# load packages

if('pacman' %in% rownames(installed.packages()) == FALSE) {
  
  install.packages('pacman', repos = "http://cran.case.edu")
  
}

setwd('~/Google_Drive/R/penguin_watch_model/Scripts/2-process-SIC-data/')

pacman::p_load(rgdal, sp, gdata)

work1 <- read.csv('../../Data/site_ll.csv', stringsAsFactors = FALSE, na.strings = NA)

coordinates(work1) <- ~longitude + latitude

work1.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

proj4string(work1) <- work1.geo  # define projection system of our data

work2 <- sp::spTransform(work1, CRS("+init=epsg:3031"))

rgdal::writeOGR(work2, dsn = 'Locations', layer = 'SitesEPSG3031', driver = "ESRI Shapefile")

