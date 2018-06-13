#################
# Penguin Watch Model - 3 - Process PW data
#
# 0-detect-params.R | recovering generating values using one detection param vs many detection params
# 00-recover-params.R | recovering generating values using realistic data/model
# 1-process-krill-data.R | process krill data
# 2-process-SIC-data.R | process SIC data
# 3-process-pw-data.R | process PW Pro data
# 4-model.R | penguin model
# 4-run-model.pbs | pbs script to run penguin model on HPC resources
# 5-analyze-output.R | analyze model output
#
# Author: Casey Youngflesh
#################



# Clear environment -------------------------------------------------------

rm(list = ls())


# DIR ---------------------------------------------------------------------

#laptop

if (Sys.info()[[1]] == 'Windows')
{
  setwd('C:/Users/Lynch Lab 7/Research/Projects/Penguin_watch/PW_surv_model_data/Manual_QC_data_files/')
} else {
  setwd('~/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/Manual_QC_data_files/')
}



# Load packages -----------------------------------------------------------

library(dplyr)



# load QC data ------------------------------------------------------------


files <- list.files()[grep('.csv', list.files())]

full_df <- data.frame()
for (i in 1:length(files))
{
  temp <- read.csv(paste0(files[i]), stringsAsFactors = FALSE)
  full_df <- dplyr::bind_rows(full_df, temp)
}


# write to csv ------------------------------------------------------------

if (Sys.info()[[1]] == 'Windows')
{
  setwd('C:/Users/Lynch Lab 7/Google_Drive/R/penguin_watch_model/Data/PW_data/')
} else {
  setwd('~/Google_Drive/R/penguin_watch_model/Data/PW_data/')
}


write.csv(full_df, file = 'PW_data_June_12_2018.csv', row.names = FALSE)
