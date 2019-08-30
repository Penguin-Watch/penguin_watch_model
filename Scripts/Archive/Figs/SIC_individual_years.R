######################
#Script to plot SIC rasters for individaul years/months
#
#
######################


if (Sys.info()[[1]] == 'Windows')
{
  setwd('C:/Users/Lynch Lab 7/Google_Drive/R/Project_archive/MAPPPD/SiteCovariates/PassiveMicrowaveSIC')
} else {
  setwd("~/Google_Drive/R/Project_archive/MAPPPD/SiteCovariates/PassiveMicrowaveSIC")
}

mn_ay <- raster::raster(paste0('nt_201409_f17_v1.1_s.tif'))


mn_ay_spdf <- as(mn_ay, 'SpatialPixelsDataFrame')
mn_ay_df <- as.data.frame(mn_ay_spdf)
colnames(mn_ay_df) <- c('value', 'x', 'y')

#NA vals for 251 (cicular mask), 252 (unused), 253 (coastlines), 254 (landmask), and 255 (missing data)
mn_ay_df$value[which(mn_ay_df$value > 250)] <- NA
#scale 0 - 100
mn_ay_df$value <- (mn_ay_df$value/250)*100



ggplot(mn_ay_df, aes(x=x, y=y, fill = value)) +
  geom_tile() +
  theme_void() +
  coord_equal(xlim = c(-2750000, -2200000),
              ylim = c(1000000, 2000000)) +
  # scale_fill_gradient2('Mean SIC - June - Sep',
  #                      limits = c(0,
  #                               100),
  #                      na.value = 'white',
  #                      low = 'white',
  #                      high = 'light blue') +
  geom_polygon(data = AP, aes(long, lat, group = group),
               inherit.aes = FALSE,
               fill = 'grey') +
  geom_path(data = AP, aes(long, lat, group = group),
            inherit.aes = FALSE,
            color = 'black')



