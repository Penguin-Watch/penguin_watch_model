#################
# Penguin Watch Model - 4 - Penguin model (run with PBS script; results saved to Results)
#
# model to look at when is season nests fail
#
# Author: Casey Youngflesh
#################

#include data starting at first chick sighting, ending at creche point. No NA needed for night.
#can control number of days starting data from first chick sighting with DAYS_BUFF_FIRST_CHICK


#include data starting at first chick sighting, ending at creche point. No NA needed for night.
#can control number of days starting data from first chick sighting with DAYS_BUFF_FIRST_CHICK


# Clear environment -------------------------------------------------------


rm(list = ls())


# DIR ---------------------------------------------------------------------


#laptop
dir <- c('~/Google_Drive/R/penguin_watch_model/Data/PW_data/',
         '../Krill_data/CCAMLR/Processed_CCAMLR/',
         '../../../SIC_data/Processed/',
         '~/Google_Drive/R/penguin_watch_model/Results/')

#HPC
#dir <- c('../Data', '../Data', '../Data', '../Results')




# Load packages -----------------------------------------------------------

#devtools::install_github('caseyyoungflesh/jagsRun')

library(dplyr)
library(jagsRun)


# determine PW dates to use -----------------------------------------------

setwd(dir[1])

#make sure only periods of data that have been QCed are read in here (NA vals will be added to fill the rest of the period)
#unused nests should be marked with all NAs

#PW_data <- read.csv('PW_data_2019-03-23.csv', stringsAsFactors = FALSE)
PW_data <- read.csv('PW_data_2019-04-06.csv', stringsAsFactors = FALSE)

#remove specified colonies
un_sites_p <- sort(unique(PW_data$site))
to_rm <- c('AITC')
un_sites <- un_sites_p[-which(un_sites_p %in% to_rm)]
#un_sites <- un_sites_p


#determine years
yrs <- c()
for (k in 1:length(un_sites))
{
  #k <- 1
  temp <- dplyr::filter(PW_data, site == un_sites[k])
  un_yrs <- unique(temp$season_year)
  
  yrs <- c(yrs, un_yrs)
}



# Create nests_array -------------------------------------------

#array with number of chicks seen at each nest in each time step

#find just nest time series columns
#all colnames
cols <- colnames(PW_data)
#just columns with 'nest'
ind <- grep('nest', cols)
#which columns have x.coord and y.coord
to.rm1 <- grep('x.coord', cols)
to.rm2 <- grep('y.coord', cols)
to.rm <- c(to.rm1, to.rm2)
#which columns just have 'nest'
tog <- c(ind, to.rm)
tog2 <- tog[!(duplicated(tog) | duplicated(tog, fromLast = TRUE))]


#set dates to use for model - days before first data point in input .csv (first chick sighting)
DAYS_BEFORE <- 30
DAYS_BUFF_FIRST_CHICK <- 0
DAYS_AFTER <- 29 - DAYS_BUFF_FIRST_CHICK

#number of time steps (rows) in response data
n_ts <- DAYS_BEFORE + DAYS_BUFF_FIRST_CHICK + DAYS_AFTER + 1

#number of nests (columns) in response data - i
n_nests <- length(tog2)

#number of years (3rd dim) in response data - j
d_yrs <- sort(unique(yrs))
n_yrs <- length(d_yrs)

#number of sites (4th dim) in response data - k
n_sites <- length(un_sites)


#aggregate data by day (max number of chicks for each day)

#create blank array
nests_array <- array(NA, dim = c(n_ts, n_nests, n_yrs, n_sites))

#FILL RESPONSE DATA ARRAY
#adds NA buffer to beginning and end of data
#nests with NAs are simply removed (e.g., if there are 4 nests, and nest 3 is all NAs, nest 4 becomes nest 3)

yrs_array <- array(NA, dim = c(n_yrs, n_sites))
for (k in 1:n_sites)
{
  #k <- 5
  temp <- dplyr::filter(PW_data, site == un_sites[k])
  
  j_idx <- 1
  for (j in 1:n_yrs)
  {
    #j <- 1
    temp2 <- dplyr::filter(temp, season_year == d_yrs[j])
    
    if (NROW(temp2) > 0)
    {
      temp_dates <- as.Date(temp2$datetime, format = "%Y:%m:%d %H:%M:%S")
      
      #date range (includes days that might be missing for some reason)
      date_rng <- seq(temp_dates[1], temp_dates[length(temp_dates)], by = 'day')
      
      temp_agg <- data.frame()
      for (t in 1:length(date_rng))
      {
        #t <- 1
        td_filt <- which(temp_dates == date_rng[t])
        temp3 <- temp2[td_filt,]
        temp_max <- suppressWarnings(apply(temp3[,tog2], 2, function(x) max(x, na.rm = TRUE)))
        temp4 <- data.frame(datetime = date_rng[t], t(temp_max))
        temp_agg <- rbind(temp_agg , temp4)
      }
      
      #replace -Inf (from max) with NA
      temp_agg[which(temp_agg == -Inf, arr.ind = TRUE)] <- NA
      
      
      #Specify FIRST and LAST days of season - number of days before first data point (first chick sighting) - days after specified data start day
      FIRST <- min(date_rng) - DAYS_BEFORE
      LAST <- min(date_rng) + DAYS_BUFF_FIRST_CHICK + DAYS_AFTER
      
      # #Use same dates across all sites
      # first_date <- format(as.Date('12-15', format = '%m-%d'), '%m-%d')
      # last_date <- format(as.Date('02-15', format = '%m-%d'), '%m-%d')
      # FIRST <- as.Date(paste0((d_yrs[j]-1), '-', first_date), format = "%Y-%m-%d")
      # LAST <- as.Date(paste0(d_yrs[j], '-', last_date), format = "%Y-%m-%d")
      
      #which dates are within the designated period
      valid_dates <- which(temp_agg$datetime >= min(date_rng) + 
                             DAYS_BUFF_FIRST_CHICK & temp_agg$datetime <= LAST)
      sel_dates <- temp_agg$datetime[valid_dates]
      
      if (min(sel_dates) > FIRST)
      {
        #add NA vals to front
        #num_first <- length(which(sel_dates == min(sel_dates))) #used when not aggregating days
        
        lna <- min(sel_dates) - FIRST
        na_first <- matrix(NA, ncol = length(tog2), nrow = (lna)) # + (num_first))
      } else {
        na_first <- NULL
      }
      
      if (max(sel_dates) < LAST)
      {
        #add NA vals to end
        #num_last <- length(which(sel_dates == max(sel_dates))) #used when not aggregating days
        
        lna <- LAST - max(sel_dates)
        na_last <- matrix(NA, ncol = length(tog2), nrow = (lna))# + (num_last))
      } else {
        na_last <- NULL
      }
      
      #add buffers to front and back (if needed)
      vals <- as.matrix(temp_agg[, -1])
      n_vals <- rbind(na_first, vals[valid_dates, ], na_last)
      
      #determines if there are any nests with NA values for the entire column (removed during the QC step)
      #first time step with value at nest 1 (unless all NA, then move to next nest)
      counter <- 1
      while (sum(!is.na(n_vals[,counter])) < 1)
      {
        counter <- counter + 1
      }
      ft <- n_vals[min(which(!is.na(n_vals[,counter]))),]
      #last nest with value at this time step
      lnv <- max(which(!is.na(ft)))
      #are there any NA vals in this row?
      nst_na <- which(is.na(ft[1:lnv]))
      if (length(nst_na) > 0)
      {
        fill_mat <- matrix(NA, nrow = NROW(n_vals), ncol = length(nst_na))
        f_n_vals <- cbind(n_vals[,-nst_na], fill_mat)
      } else {
        f_n_vals <- n_vals
      }
      
      #appropriate date range and appropriate columns for nests
      nests_array[,,j_idx,k] <- f_n_vals
      yrs_array[j_idx,k] <- d_yrs[j]
      j_idx <- j_idx + 1
    }
  }
}


# create real_nests matrix ------------------------------------------------

#create matrix that has number of nests at each site/year

#rows are years, columns are sites
real_nests <- matrix(NA, nrow = n_yrs, ncol = n_sites)
for (k in 1:dim(nests_array)[4])
{
  #k <- 2
  for (j in 1:dim(nests_array)[3])
  {
    #j <- 1
    #just nest 1 - which positions are not NA
    idx_nna <- which(!is.na(nests_array[,1,j,k]))
    
    if (length(idx_nna) > 0)
    {
      #values at that time point
      ft <- nests_array[min(idx_nna),,j,k]
      #last nest with value at this time step
      lnv <- max(which(!is.na(ft)))
    } else {
      lnv <- 0
    }
    
    #max is going to be number of nests in image
    real_nests[j,k] <- lnv
  }
}



# observations with > 2 chicks --------------------------------------------

# check to see if any nests have observations with more than 2 chicks
# which(nests_array > 2, arr.ind = TRUE)
# num_o2 <- length(nests_array[which(nests_array > 2, arr.ind = TRUE)])
# total <- length(nests_array[which(!is.na(nests_array), arr.ind = TRUE)])
# num_o2/total

#determine which observation have more than two chicks observed and change them to 2
#ind.g2 <- which(nests_array > 2, arr.ind = TRUE)
#nests_array[ind.g2] <- 2




# create z_array ----------------------------------------------------------

#z_array - array of known true values

z_array <- nests_array

for (k in 1:dim(nests_array)[4])
{
  #k <- 2
  for (j in 1:dim(nests_array)[3])
  {
    #j <- 1
    if (real_nests[j,k] > 0)
    {
      for (i in 1:real_nests[j,k])
      {
        #two chicks in first position (alive at time step one)
        z_array[1,i,j,k] <- 2
        
        if (sum(z_array[,i,j,k] == 2, na.rm = TRUE) > 1)
        {
          #last sight with two chicks
          n2 <- max(which(z_array[,i,j,k] == 2))
          #fill 2 for all between first val and last sight of 2
          z_array[1:n2,i,j,k] <- 2
        }
      }
    }
  }
}

#after there aren't two, don't know if there are actually 2, 1, or 0 so NA
zeros <- which(z_array == 0, arr.ind = TRUE)
ones <- which(z_array == 1, arr.ind = TRUE)
z_array[zeros] <- NA
z_array[ones] <- NA


# Create Data for JAGS ---------------------------------------------------------

#nests_array:
#dim1 (rows) [t] = time steps
#dim2 (cols) [i] = nests
#dim3 [j] = years (d_yrs)
#dim4 [k] = sites (un_sites)


#number of years for each site
NJ <- rep(NA, NCOL(yrs_array))
for (i in 1:length(NJ))
{
  #i <- 1
  NJ[i] <- max(which(!is.na(yrs_array[,i])))
}



#just LOCK, just one year
s_ind <- which(un_sites == 'LOCK')
n_array <- nests_array[,,2,s_ind]
n_NJ <- NJ[s_ind]
n_real_nests <- real_nests[2,s_ind]
n_z_array <- z_array[,,2,s_ind]

DATA <- list(
  y = n_array, #response
  NK = 1, #dim(n_array)[4], #number of sites
  NJ = n_NJ, #number of years covered for each site
  NI = n_real_nests, #number of nests j,k [year, site]
  NT = dim(n_array)[1], #number of time steps
  z = n_z_array, #known points of bird being alive
  x = scale(as.numeric(1:dim(n_array)[1]), scale = FALSE)[,1],
  unsites = un_sites,
  yrs_array = yrs_array)


setwd(dir[4])

# Model -------------------------------------------------------------------

{
  sink("pwatch_surv.jags")
  
  cat("
      
      model {
      
      #nests
      for (i in 1:NI)
      {
      #both chicks alive at time step 1 (z[1,i,j] = 2)
      
      #time step
      for (t in 2:NT)
      {
      #state model
      z[t,i] ~ dbinom(p_alive[t,i], z[t-1,i])
      p_alive[t,i] <- ifelse(z[t-1,i] < 2, 
      phi[t,i] * z[t-1,i],
      phi[t,i])
      
      #observation model
      y[t,i] ~ dbinom(p_sight[t,i], z[t,i])
      p_sight[t,i] <- ifelse(z[t,i] < 2,
      p[t,i] * z[t,i],
      p[t,i])
      
      }
      }
      
      #transforms
      #nests      
      for (i in 1:NI)
      {
      #time
      for (t in 1:NT)
      {
    
      logit(phi[t,i]) <- mu_phi[t]
      logit(p[t,i]) <- mu_p + beta_p*x[t]
      
      } #t
      } #i
    
      for (t in 1:NT)
      {
      mu_phi[t] ~ dnorm(3, 0.1)
      }

  
      #priors - p and phi
      beta_p ~ dnorm(0.1, 10) T(0, 0.5)
      mu_p ~ dnorm(2, 0.1)
      
      
      }",fill = TRUE)

  sink()
}



# Starting values ---------------------------------------------------------


Inits_1 <- list(mu_phi = rep(4, dim(DATA$y)[1]),
                mu_p = 2,
                beta_p = 0.1,
                .RNG.name = "base::Mersenne-Twister", 
                .RNG.seed = 1)

Inits_2 <- list(mu_phi = rep(4, dim(DATA$y)[1]),
                mu_p = 2,
                beta_p = 0.1,
                .RNG.name = "base::Wichmann-Hill", 
                .RNG.seed = 2)

Inits_3 <- list(mu_phi = rep(4, dim(DATA$y)[1]),
                mu_p = 2,
                beta_p = 0.1,
                .RNG.name = "base::Marsaglia-Multicarry", 
                .RNG.seed = 3)

# Inits_4 <- list(mu_phi = matrix(4, nrow = dim(DATA$y)[3], ncol = dim(DATA$y)[4]),
#                 mu_p = 2,
#                 beta_p = 0.1,
#                 sigma_mu_phi = 1,
#                 sigma_nu_p = 1,
#                 alpha_theta = 0,
#                 pi_theta = 0,
#                 rho_theta = 0,
#                 .RNG.name = "base::Mersenne-Twister", 
#                 .RNG.seed = 4)
# 
# Inits_5 <- list(mu_phi = matrix(4, nrow = dim(DATA$y)[3], ncol = dim(DATA$y)[4]),
#                 mu_p = 2,
#                 beta_p = 0.1,
#                 sigma_mu_phi = 1,
#                 sigma_nu_p = 1,
#                 alpha_theta = 0,
#                 pi_theta = 0,
#                 rho_theta = 0,
#                 .RNG.name = "base::Wichmann-Hill",
#                 .RNG.seed = 5)
# 
# Inits_6 <- list(mu_phi = matrix(4, nrow = dim(DATA$y)[3], ncol = dim(DATA$y)[4]),
#                 mu_p = 2,
#                 beta_p = 0.1,
#                 sigma_mu_phi = 1,
#                 sigma_nu_p = 1,
#                 alpha_theta = 0,
#                 pi_theta = 0,
#                 rho_theta = 0,
#                 .RNG.name = "base::Wichmann-Hill", 
#                 .RNG.seed = 6)

F_Inits <- list(Inits_1, Inits_2, Inits_3)#, Inits_4, Inits_5, Inits_6)



# Parameters to track -----------------------------------------------------

Pars <- c('mu_phi',
          'mu_p',
          'beta_p'
)


# Run model ---------------------------------------------------------------

#make sure model compiles
# jagsRun(jagsData = DATA,
#         jagsModel = 'pwatch_surv.jags',
#         jagsInits = F_Inits,
#         DEBUG = TRUE)


jagsRun(jagsData = DATA, 
        jagsModel = 'pwatch_surv.jags',
        jagsInits = F_Inits,
        params = Pars,
        jagsID = 'PW_30k_2019-04-05_LOCK_time_v',
        jagsDsc = '1 sites, 1 years, timing
        logit(phi) <- mu_phi_t;
        logit(p) <- mu_p + beta*x[t]',
        db_hash = 'PW_data_2019-04-06.csv',
        n_chain = 3,
        n_adapt = 5000,
        n_burn = 30000,
        n_draw = 30000,
        n_thin = 20,
        EXTRA = FALSE,
        Rhat_max = 1.1,
        n_max = 100000,
        save_data = TRUE)


