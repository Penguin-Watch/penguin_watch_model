#################
# Penguin Watch Model - 4 - Penguin model (run with PBS script; results saved to Results)
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

#include data starting at first chick sighting, ending at creche point. No NA needed for night.
#can control number of days starting data from first chick sighting with DAYS_BUFF_FIRST_CHICK


# Clear environment -------------------------------------------------------


rm(list = ls())


# DIR ---------------------------------------------------------------------


#laptop
# dir <- c('~/Google_Drive/R/penguin_watch_model/Data/PW_data/',
#          '../Krill_data/CCAMLR/Processed_CCAMLR/',
#          '../../../SIC_data/Processed/',
#          '~/Google_Drive/R/penguin_watch_model/Results/')

#HPC
dir <- c('../Data', '../Data', '../Data', '../Results')




# Load packages -----------------------------------------------------------

#devtools::install_github('caseyyoungflesh/jagsRun')

library(dplyr)
library(jagsRun)



# determine PW dates to use -----------------------------------------------

setwd(dir[1])

#make sure only periods of data that have been QCed are read in here (NA vals will be added to fill the rest of the period)
#unused nests should be marked with all NAs

PW_data <- read.csv('PW_data_June_12_2018.csv', stringsAsFactors = FALSE)

boots <- which(PW_data$site == 'BOOT')
PW_data$site[boots] <- 'PCHA'

un_sites_p <- unique(PW_data$site)

#remove colonies to make dataset smaller (faster model run)
#to_rm <- c('AITC')
#un_sites <- un_sites_p[-which(un_sites_p %in% to_rm)]

un_sites <- un_sites_p


#determine years
yrs <- c()
for (k in 1:length(un_sites))
{
  #k <- 1
  temp <- filter(PW_data, site == un_sites[k])
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
d_yrs <- sort(yrs[!duplicated(yrs)])
n_yrs <- length(d_yrs)

#number of sites (4th dim) in response data - k
n_sites <- length(un_sites)


#aggregate data by day (max number of chicks for each day)

#create blank array
nests_array <- array(NA, dim = c(n_ts, n_nests, n_yrs, n_sites))

#FILL RESPONSE DATA ARRAY
#adds NA buffer to beginning and end of data
#nests with NAs are simply removed (e.g., if there are 4 nests, and nest 3 is all NAs, nest 4 becomes nest 3)
for (k in 1:n_sites)
{
  #k <- 1
  temp <- filter(PW_data, site == un_sites[k])
  
  for (j in 1:n_yrs)
  {
    #j <- 2
    temp2 <- filter(temp, season_year == d_yrs[j])
    
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
      valid_dates <- which(temp_agg$datetime >= min(date_rng) + DAYS_BUFF_FIRST_CHICK & temp_agg$datetime <= LAST)
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
      #first time step with value at nest 1
      ft <- n_vals[min(which(!is.na(n_vals[,1]))),]
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
      nests_array[,,j,k] <- f_n_vals
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




# Krill covariate ---------------------------------------------------------------

#total krill caught over the previous winter and current breeding season
#150km radius for March - Feb

#dim1 [j] = years (d_yrs)
#dim2 [k] = sites (un_sites)

un_sites_ncc <- substr(un_sites, start = 1, stop = 4)

setwd(dir[2])

krill <- read.csv('CCAMLR_krill_entire_season.csv')

i_KRILL <- matrix(nrow = length(d_yrs), ncol = length(un_sites_ncc))
for (k in 1:length(un_sites_ncc))
{
  #k <- 5
  temp <- filter(krill, SITE == un_sites_ncc[k])
  
  for (j in 1:length(d_yrs))
  {
    #j <- 1
    temp2 <- filter(temp, YEAR == d_yrs[j])
    i_KRILL[j,k] <- temp2$T_KRILL
  }
}

#standardize krill catch
t1 <- as.vector(i_KRILL)
t2 <- scale(t1)[,1]
KRILL <- matrix(t2, nrow = NROW(i_KRILL))



# SIC covariate -----------------------------------------------------------

#SIC for the previous winter
#150km radius for June - Sep
#COULD ALSO USE 500 KM
#COULD ALSO USE MAX OVER LAST 5 YEARS


setwd(dir[3])

sea_ice <- read.csv('SIC_150_W.csv')

i_SIC <- matrix(nrow = length(d_yrs), ncol = length(un_sites_ncc))
for (k in 1:length(un_sites_ncc))
{
  #k <- 1
  temp <- filter(sea_ice, SITE == un_sites_ncc[k])
  
  for (j in 1:length(d_yrs))
  {
    #j <- 1
    temp2 <- filter(temp, YEAR == d_yrs[j])
    i_SIC[j,k] <- temp2$W_MN
  }
}


#standardize krill catch
t1 <- as.vector(i_SIC)
t2 <- scale(t1)[,1]
SIC <- matrix(t2, nrow = NROW(i_SIC))



# Create Data for JAGS ---------------------------------------------------------

#nests_array:
#dim1 (rows) [t] = time steps
#dim2 (cols) [i] = nests
#dim3 [j] = years (d_yrs)
#dim4 [k] = sites (un_sites)

DATA <- list(
  y = nests_array, #response
  NK = dim(nests_array)[4], #number of sites
  NJ = dim(nests_array)[3], #number of years covered for all sites
  NI = real_nests, #matrix with number of nests for each site/year NI[j,k]
  NT = dim(nests_array)[1], #number of time steps
  z = z_array, #known points of bird being alive
  x = scale(as.numeric(1:dim(nests_array)[1]), scale = FALSE)[,1], #time steps for increase in surv/detection over time 
  KRILL = KRILL, #standardized krill catch data (total krill caught over the previous winter and current breeding season)
  SIC = SIC) #standardized SIC for previous winter


#checks:
# nests_array[1:50, 1:10, 2, 1]
# z_array[1:10, 1:11, 2, 1]
# DATA$NI[2,1]



setwd(dir[4])

# Model -------------------------------------------------------------------

{
  sink("pwatch_surv.jags")
  
  cat("

      model {
      
      #sites
      for (k in 1:NK)
      {
      #years
      for (j in 1:NJ)
      {
      #nests
      for (i in 1:NI[j,k])
      {
      #both chicks alive at time step 1 (z[1,i,j,k] = 2)
      
      #time step
      for (t in 2:NT)
      {
      #state model
      z[t,i,j,k] ~ dbinom(p_alive[t,i,j,k], z[t-1,i,j,k])
      p_alive[t,i,j,k] <- ifelse(z[t-1,i,j,k] < 2, 
      phi[t,i,j,k] * z[t-1,i,j,k],
      phi[t,i,j,k])
      
      #observation model
      y[t,i,j,k] ~ dbinom(p_sight[t,i,j,k], z[t,i,j,k])
      p_sight[t,i,j,k] <- ifelse(z[t,i,j,k] < 2,
      p[t,i,j,k] * z[t,i,j,k],
      p[t,i,j,k])
      
      }
      }
      }
      }
      
      
      #transforms
      #sites
      for (k in 1:NK)
      {
      #years
      for (j in 1:NJ)
      {
      #nests
      for (i in 1:NI[j,k])
      {
      #time
      for (t in 1:NT)
      {
      #phi = survival prob
      #mu_phi = grand mean for all sites/years
      #eta_phi = effect of site
      #gamma_phi = effect of year
      #pi_phi = effect of SIC on survival
      #rho_phi = effect of KRILL on survival
      
      logit(phi[t,i,j,k]) <- mu_phi + gamma_phi[j] + eta_phi[k] + pi_phi * SIC[j,k] + rho_phi * KRILL[j,k]

      #p = detection prob
      #mu_p = grand mean for all sites/years
      #beta_phi = slope for increasing detection over time (older chicks have higher detection p)
      #nu_p = detection probability random effect for site/year
      
      logit(p[t,i,j,k]) <- mu_p + beta_p*x[t] + nu_p[j,k]

      } #t
      } #i
      } #j
      } #k
      

      #priors - covariates
      pi_phi ~ dnorm(0, 0.386)
      rho_phi ~ dnorm(0, 0.386)

      #priors - phi
      mu_phi ~ dnorm(3, 0.1)
      
      for (k in 1:NK)
      {
      #eta_phi[k] ~ dnorm(0, tau_eta_phi)
      eta_phi[k] ~ dnorm(0, 0.386)    
      }
      
      for (j in 1:NJ)
      {
      #gamma_phi[j] ~ dnorm(0, tau_gamma_phi)
      gamma_phi[j] ~ dnorm(0, 0.386)
      }
      
      #tau_eta_phi <- pow(sigma_eta_phi, -2)
      #sigma_eta_phi ~ dunif(0, 5)
      
      #tau_gamma_phi <- pow(sigma_gamma_phi, -2)
      #sigma_gamma_phi ~ dunif(0, 5)
      
      #priors - p
      mu_p ~ dnorm(2, 0.1)
      beta_p ~ dnorm(0.1, 10) T(0, 0.5)

      for (k in 1:NK)
      {
       for (j in 1:NJ)
       {
         nu_p[j,k] ~ dnorm(0, tau_nu_p)
       }
      }
      
      tau_nu_p ~ dunif(0, 25)

      
      }",fill = TRUE)

  sink()
}



# Starting values ---------------------------------------------------------


Inits_1 <- list(mu_phi = 4,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                pi_phi = 0,
                rho_phi = 0,
                mu_p = 2,
                beta_p = 0.1,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                tau_nu_p = 15,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mu_phi = 4,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                pi_phi = 0,
                rho_phi = 0,
                mu_p = 2,
                beta_p = 0.1,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                tau_nu_p = 15,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mu_phi = 4,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                pi_phi = 0,
                rho_phi = 0,
                mu_p = 2,
                beta_p = 0.1,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                tau_nu_p = 15,
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

Inits_4 <- list(mu_phi = 4,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                pi_phi = 0,
                rho_phi = 0,
                mu_p = 2,
                beta_p = 0.1,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                tau_nu_p = 15,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 4)

Inits_5 <- list(mu_phi = 4,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                pi_phi = 0,
                rho_phi = 0,
                mu_p = 2,
                beta_p = 0.1,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                tau_nu_p = 15,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 5)

Inits_6 <- list(mu_phi = 4,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                pi_phi = 0,
                rho_phi = 0,
                mu_p = 2,
                beta_p = 0.1,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                tau_nu_p = 15,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 6)

F_Inits <- list(Inits_1, Inits_2, Inits_3)#, Inits_4, Inits_5, Inits_6)



# Parameters to track -----------------------------------------------------

Pars <- c('mu_phi',
          'eta_phi',
          'gamma_phi',
          'pi_phi',
          'rho_phi',
          #'sigma_eta_phi',
          #'sigma_gamma_phi',
          'mu_p',
          'beta_p',
          'nu_p',
          'tau_nu_p'
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
               jagsID = 'June_12_2018_100k',
               jagsDsc = '5 sites, 3 years
            100k iter
            Long queue
            First QC data - before POLAR
            non-hierarchical gamma, etc, hierarchical nu
            logit(phi) <- mu + gamma + eta + pi + rho; 
            logit(p) <- mu + beta*x + nu',
               db_hash = 'PW_data_June_12_2018.csv',
               n_chain = 6,
               n_adapt = 5000,
               n_burn = 100000,
               n_draw = 100000,
               n_thin = 20,
               EXTRA = FALSE,
               Rhat_max = 1.1,
               n_max = 100000,
               save_data = TRUE)


