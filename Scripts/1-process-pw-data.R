#################
# 1 - Process PW data
#
# Author: Casey Youngflesh
#################

#*create images with nest polygons for each site/year
#*classify images in penguin watch pro
#*create QC images to check user classifications
#*manually QC images - place in 'Manual_QC_data_files/
#*retain only relevant portions of post-QC classifications for each site (post-first chick sighting, pre-creche) - place in 'Model_input/'
#*run this script
#*run model script
#*analyze model output


# Clear environment -------------------------------------------------------

rm(list = ls())


# DIR ---------------------------------------------------------------------


DATE <- '2019-07-10'

setwd(paste0('~/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/Model_input_', DATE, '/'))




# Load packages -----------------------------------------------------------

library(dplyr)



# load QC data ------------------------------------------------------------

files_p <- list.files()[grep('.csv', list.files())]
#don't yet have 2018 krill data
#to.rm <- grep('2018', files)
#files <- files_p[-to.rm]
files <- files_p

full_df <- data.frame()
for (i in 1:length(files))
{
  #i <- 24
  temp <- read.csv(paste0(files[i]), stringsAsFactors = FALSE, header = TRUE)
  full_df <- dplyr::bind_rows(full_df, temp)
}


# write to csv ------------------------------------------------------------

setwd('~/Google_Drive/R/penguin_watch_model/Data/PW_data/')

write.csv(full_df, file = paste0('PW_data_', DATE, '.csv'), row.names = FALSE)
