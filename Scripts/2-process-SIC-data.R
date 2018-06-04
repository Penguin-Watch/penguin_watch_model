#################
# Penguin Watch Model - 2 - Process SIC data
#
# 1-process-krill-data.R | process krill data
# 2-process-SIC-data.R | process SIC data
# 3-model.R | penguin model
# 3-run-model.pbs | pbs script to run penguin model on HPC resources
# 4-analyze-output.R | analyze model output
#
# Author: Casey Youngflesh
#################



# Sea ice covariates -----------------------------------------------------------------

#Several ways in which SIC may impact BS:

# * Effect of SIC during breeding season (access to food resources)
# * Effect of SIC during winter previous to breeding season (body condtion from previous winter)
# * Effect of SIC during previous 5 years (conditions that may or may not be favorable to krill)



# Clear environment -------------------------------------------------------

rm(list = ls())



# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman')
}

pacman::p_load(dplyr)





# Effect of SIC during breeding season -------------------------------------

#50km radius for Dec - Feb in each year

setwd('~/Google_Drive/R/penguin_watch_model/Data/SIC_data/RAW')

#year for OUT is PW year (e.g., year 2000 in OUT is 1999/2000 breeding season)


#NEED TO RUN MAPPPD CODE FOR 50KM

# data_50 <- read.csv('ALLSITES_SIC_50_MEAN.csv')
# 
# SIC_50_BS <- data.frame()
# for (i in 1:length(cam_sites))
# {
#   #i <- 1
#   temp_SIC <- filter(data_50, site_id == cam_sites[i])
# 
#   T_OUT <- data.frame()
#   for (j in 1979:2016)
#   {
#     #j <- 1979
#     temp <- filter(temp_SIC, year == j)
#     dec <- temp$SIC_MONTH_12
#     jan <- temp$SIC_MONTH_1
#     feb <- temp$SIC_MONTH_2
#     mn_djf <- mean(c(dec, jan, feb))
#     temp2 <- data.frame(YEAR = j+1,
#                         SITE = cam_sites[i],
#                         DEC = dec, 
#                         JAN = jan,
#                         FEB = feb,
#                         S_MN = mn_djf)
#     T_OUT <- rbind(T_OUT, temp2)
#   }
#   SIC_50_BS <- rbind(SIC_50_BS, T_OUT)
# }
# 
# setwd('../Processed')
# write.csv(SIC_50_BS, 'SIC_50_BS.csv', row.names = FALSE)





# Effect of SIC during winter previous to breeding season -----------------

#500km radius for June - Sep

setwd('../RAW')
data_500 <- read.csv('ALLSITES_SIC_500_MEAN.csv')

SIC_500_W <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 1
  temp_SIC <- filter(data_500, site_id == cam_sites[i])
  
  T_OUT <- data.frame()
  for (j in 1979:2016)
  {
    #j <- 1979
    temp <- filter(temp_SIC, year == j)
    jun <- temp$SIC_MONTH_6
    jul <- temp$SIC_MONTH_7
    aug <- temp$SIC_MONTH_8
    sep <- temp$SIC_MONTH_9
    mn_jjas <- mean(c(jun, jul, aug, sep))
    temp2 <- data.frame(YEAR = j+1,
                        SITE = cam_sites[i],
                        JUN = jun,
                        JUL = jul,
                        AUG = aug,
                        SEP = sep,
                        W_MN = mn_jjas)
    T_OUT <- rbind(T_OUT, temp2)
  }
  SIC_500_W <- rbind(SIC_500_W, T_OUT)
}

# setwd('../Processed/')
# 
# write.csv(SIC_500_W, 'SIC_500_W.csv', row.names = FALSE)




#150km radius for June - Sep


data_150 <- read.csv('ALLSITES_SIC_150_MEAN.csv')

SIC_150_W <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 1
  temp_SIC <- filter(data_150, site_id == cam_sites[i])
  
  T_OUT <- data.frame()
  for (j in 1979:2016)
  {
    #j <- 1979
    temp <- filter(temp_SIC, year == j)
    jun <- temp$SIC_MONTH_6
    jul <- temp$SIC_MONTH_7
    aug <- temp$SIC_MONTH_8
    sep <- temp$SIC_MONTH_9
    mn_jjas <- mean(c(jun, jul, aug, sep))
    temp2 <- data.frame(YEAR = j+1,
                        SITE = cam_sites[i],
                        JUN = jun,
                        JUL = jul,
                        AUG = aug,
                        SEP = sep,
                        W_MN = mn_jjas)
    T_OUT <- rbind(T_OUT, temp2)
  }
  SIC_150_W <- rbind(SIC_150_W, T_OUT)
}

# setwd('../Processed/')
# 
# write.csv(SIC_150_W, 'SIC_150_W.csv', row.names = FALSE)



# Effect of SIC during previous 5 years -----------------------------------

#500km radius for Jun - Sep for previous 5 winters

#Justification (From Che-Castaldo et al. 2017 Nat Comms): We chose covariates hypothesized to influence Adélie breeding success that were available for all sites and years included in this study. Krill abundance relies on winter sea ice conditions, with krill populations requiring at least one year of heavy winter sea ice every 4-5 years to replenish standing stocks2. Based on this, we defined peak winter sea ice concentration (wsic) in the Adélie model as the maximum monthly concentration (June through September) across the previous 5 winters in a 500 km radius around each Adélie site. This covariate was standardized and has a mean and standard deviation of 81 and 12, respectively.


SIC_500_MAX <- data.frame()
for (i in 1:length(cam_sites))
{
  #i <- 1
  #input from SIC_500_W uses PW year
  temp_SIC <- filter(SIC_500_W, SITE == cam_sites[i])
  
  T_OUT <- data.frame()
  for (j in 1984:2017)
  {
    #j <- 1984
    #5 years prior to 1984 PW season
    t2_SIC <- filter(temp_SIC, YEAR >= j - 4, YEAR <= j)
    mx_SIC<- max(t2_SIC$W_MN)
    
    temp2 <- data.frame(YEAR = j,
                        SITE = cam_sites[i],
                        MAX_SIC = mx_SIC)
    
    T_OUT <- rbind(T_OUT, temp2)
  }
  SIC_500_MAX <- rbind(SIC_500_MAX, T_OUT)
}

# setwd('../Processed/')
# 
# write.csv(SIC_500_MAX, 'SIC_500_MAX.csv', row.names = FALSE)





