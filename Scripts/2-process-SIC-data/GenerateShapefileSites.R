# ________________________________________________________________________________
# load packages

if('pacman' %in% rownames(installed.packages()) == FALSE) {
  
  install.packages('pacman', repos = "http://cran.case.edu")
  
}

pacman::p_load(rgdal, sp, gdata)

work1 <- read.xls(paste0(wd, "/SiteLocations.xlsx"), sheet = 1, stringsAsFactors = FALSE, na.strings = NA)

coordinates(work1) <- ~longitude + latitude

work1.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

proj4string(work1) <- work1.geo  # define projection system of our data

work2 <- spTransform(work1, CRS("+init=epsg:3031"))

writeOGR(work2, paste0(wd2, "/SiteCovariates/Locations"), "SitesEPSG3031", driver = "ESRI Shapefile")
