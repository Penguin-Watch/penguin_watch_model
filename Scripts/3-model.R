#################
# Penguin Watch Model - 3 - Penguin model (run with PBS script; results saved to Results)
#
# 1-process-krill-data.R | process krill data
# 2-process-SIC-data.R | process SIC data
# 3-model.R | penguin model
# 3-run-model.pbs | pbs script to run penguin model on HPC resources
# 4-analyze-output.R | analyze model output
#
# Author: Casey Youngflesh
#################


# Clear environment -------------------------------------------------------


rm(list = ls())


# DIR ---------------------------------------------------------------------


#laptop
# dir <- c('~/Google_Drive/R/penguin_watch_model/Data/PW_data/RAW_Fiona_Apr_15_2018/',
#          '../../Krill_data/CCAMLR/Processed_CCAMLR/',
#          '../../../SIC_data/Processed/',
#          '~/Google_Drive/R/penguin_watch_model/Results/')

#HPC
dir <- c('../Data', '../Data', '../Data', '../Results')




# Load packages -----------------------------------------------------------

#devtools::install_github('caseyyoungflesh/jagsRun')

library(abind)
library(dplyr)
library(jagsRun)



# determine PW dates to use -----------------------------------------------

#ensures that row dimension (time steps within season) will have the same dimension

setwd(dir[1])

PW_data <- read.csv('Markrecap_data_15.05.18.csv', stringsAsFactors = FALSE)

boots <- which(PW_data$site == 'BOOT')
PW_data$site[boots] <- 'PCHA'

#remove colonies to make dataset smaller (faster model run)
to_rm <- c('BAIL', 'HALF', 'NEKO', 'ORNE')

un_sites_p <- unique(PW_data$site)
un_sites <- un_sites_p[-which(un_sites_p %in% to_rm)]


#first date of each season
yrs <- c()
min_date <- as.Date(NA)
max_date <- as.Date(NA)
for (k in 1:length(un_sites))
{
  #k <- 2
  temp <- filter(PW_data, site == un_sites[k])
  un_yrs <- unique(temp$season_year)
  
  yrs <- c(yrs, un_yrs)
  for (j in 1:length(un_yrs))
  {
    #j <- 2
    temp2 <- filter(temp, season_year == un_yrs[j])
    
    temp_dates <- as.Date(temp2$datetime, format = "%Y:%m:%d %H:%M:%S")
    t_min_date <- min(temp_dates)
    t_max_date <- max(temp_dates)
    max_date <- c(max_date, t_max_date)
    min_date <- c(min_date, t_min_date)
  }
}


#first date of season to use across all years (add 1 to start on full day)
#f_min_date <- format((min_date+1), '%m-%d')
#first_date <- max(f_min_date, na.rm = TRUE)
#start Jan 15
first_date <- format(as.Date('01-15', format = '%m-%d'), '%m-%d')

#use Feb 15 as last date of season to use across all years
last_date <- format(as.Date('02-15', format = '%m-%d'), '%m-%d')




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


#number of time steps (rows) in response data (+1 to account for start FIRST through LAST)
n_ts <- as.numeric(as.Date(paste0('2017', '-', last_date, format = '%Y-%m-%d')) - 
                     as.Date(paste0('2017', '-', first_date, format = '%Y-%m-%d')) + 1) * 24


#number of nests (columns) in response data - i
n_nests <- length(tog2)

#number of years (3rd dim) in response data - j
d_yrs <- sort(yrs[!duplicated(yrs)])
n_yrs <- length(d_yrs)

#number of sites (4th dim) in response data - k
n_sites <- length(un_sites)

#create blank array
nests_array <- array(NA, dim = c(n_ts, n_nests, n_yrs, n_sites))

#SITE 2, YEAR 2
#fill response data array
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
      
      #first and last days in season
      FIRST <- as.Date(paste0((d_yrs[j]), '-', first_date), format = "%Y-%m-%d")
      LAST <- as.Date(paste0(d_yrs[j], '-', last_date), format = "%Y-%m-%d")
      
      valid_dates <- which(temp_dates >= FIRST & temp_dates <= LAST)
      sel_dates <- temp_dates[which(temp_dates >= FIRST & temp_dates <= LAST)]
      
      if (min(sel_dates) > FIRST)
      {
        #add NA vals to front
        num_first <- length(which(sel_dates == min(sel_dates)))
        
        lna <- min(sel_dates) - FIRST
        na_first <- matrix(NA, ncol = length(tog2), nrow = (lna * 24) + (24 - num_first))
      } else {
        na_first <- NULL
      }
      
      if (max(sel_dates) < LAST)
      {
        #add NA vals to end
        num_last <- length(which(sel_dates == max(sel_dates)))
        
        lna <- LAST - max(sel_dates)
        na_last <- matrix(NA, ncol = length(tog2), nrow = (lna * 24) + (24 - num_last))
      } else {
        na_last <- NULL
      }
      
      #add buffers to front and back (if needed)
      vals <- as.matrix(temp2[valid_dates, tog2])
      n_vals <- rbind(na_first, vals, na_last)
      
      #appropriate date range and appropriate columns for nests
      nests_array[,,j,k] <- n_vals
    }
  }
}




# create real_nests matrix ------------------------------------------------

#create matrix that has number of nests at each site/year

#rows are years, columns are sites
real_nests <- matrix(NA, nrow = n_yrs, ncol = n_sites)
for (k in 1:dim(nests_array)[4])
{
  #k <- 1
  for (j in 1:dim(nests_array)[3])
  {
    #j <- 2
    #13 rows to make sure NA are due to night 
    temp <- nests_array[1:13,,j,k]
    
    tsum <- c()
    for (m in 1:13)
    {
      #m <- 7
      tsum <- c(tsum, sum(!is.na(temp[m,])))
    }
    #max is going to be number of nests in image
    real_nests[j,k] <- max(tsum)
  }
}





# create w_array ----------------------------------------------------------

#create w_array (day and night) - and before and after actual dates we have data for

w_array <- nests_array
#assign 1 to values with observations (either 0/1/2 for nests_array)
w_array[which(!is.na(w_array), arr.ind = TRUE)] <- 1

#only fill in 0 if this was a night (there is data on either side and not just missing nests)
for (k in 1:dim(nests_array)[4])
{
  #k <- 4
  for (j in 1:dim(nests_array)[3])
  {
    #j <- 1
    #if there are more than 0 nests in that year (there is data for that site/year)
    if (real_nests[j,k] > 0)
    {
      sub_i <- real_nests[j,k]
      w_array[cbind(which(is.na(w_array[,1:sub_i,j,k]), arr.ind = TRUE), j, k)] <- 0
    }
  }
}

#checks:
#nests_array[1:15, 2, 1, 4]
#w_array[1:15, 1:10, 1, 4]
#z_array[1:30, 1:10, 1, 4]


# observations with > 2 chicks --------------------------------------------

# 1% of nest observations have more than 2 chicks in them
# which(nests_array > 2, arr.ind = TRUE)
# num_o2 <- length(nests_array[which(nests_array > 2, arr.ind = TRUE)])
# total <- length(nests_array[which(!is.na(nests_array), arr.ind = TRUE)])
# num_o2/total

#determine which observation have more than two chicks observed and change them to 2
ind.g2 <- which(nests_array > 2, arr.ind = TRUE)
nests_array[ind.g2] <- 2




# create z_array ----------------------------------------------------------

#z_array - array of known true values

z_array <- nests_array

for (k in 1:dim(nests_array)[4])
{
  #k <- 1
  for (j in 1:dim(nests_array)[3])
  {
    #j <- 2
    for (i in 1:real_nests[j,k])
    {
      n1 <- 1
      if (sum(z_array[,i,j,k] == 2, na.rm = TRUE) > 0)
      {
        #last sight with two chicks
        n2 <- max(which(z_array[,i,j,k] == 2))
        #fill 2 for all between first val and last sight of 2
        z_array[n1:n2,i,j,k] <- 2
      }
      
      #NA at first state because model designates 2 chicks at time step 1
      z_array[n1,i,j,k] <- NA
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
  #k <- 1
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
#500km radius for June - Sep
#COULD ALSO ADD MAX OVER LAST 5 YEARS


setwd(dir[3])


sea_ice <- read.csv('SIC_500_W.csv')

i_SIC <- matrix(nrow = length(d_yrs), ncol = length(un_sites_ncc))
for (k in 1:length(un_sites_ncc))
{
  #k <- 1
  temp <- filter(sea_ice, SITE == un_sites_ncc[k])
  
  for (j in 1:length(d_yrs))
  {
    #j <- 1
    temp2 <- filter(temp, YEAR == d_yrs[j])
    i_SIC[j,k] <- temp2$WMN
  }
}


#standardize krill catch
t1 <- as.vector(i_SIC)
t2 <- scale(t1)[,1]
SIC <- matrix(t2, nrow = NROW(i_SIC))



# Create Data for JAGS ---------------------------------------------------------

#data object for JAGS model

#nests_array:
#dim1 (rows) [t] = time steps
#dim2 (cols) [i] = nests
#dim3 [j] = years (d_yrs)
#dim4 [k] = sites (un_sites)


DATA <- list(
  y = nests_array, #reponse
  NK = dim(nests_array)[4], #number of sites
  NJ = dim(nests_array)[3], #number of years covered for all sites
  NI = real_nests, #matrix with number of nests for each site/year NI[j,k]
  NT = dim(nests_array)[1], #number of time steps
  z = z_array, #known points of bird being alive
  w = w_array, #binary day (1)/night (0)
  x = scale(as.numeric(1:dim(nests_array)[1]), scale = FALSE)[,1], #time steps for increase in surv/detection over time 
  KRILL = KRILL, #standardized krill catch data (total krill caught over the previous winter and current breeding season)
  SIC = SIC) #standardized SIC for previous winter


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
      #both chicks alive at time step 1
      z[1,i,j,k] <- 2
      
      #time step
      for (t in 2:NT)
      {
      #state model
      z[t,i,j,k] ~ dbinom(p_alive[t,i,j,k], z[t-1,i,j,k])
      p_alive[t,i,j,k] <- ifelse(z[t-1,i,j,k] < 2, 
      phi[t,i,j,k] * z[t-1,i,j,k],
      phi[t,i,j,k])
      
      #observation model
      y[t,i,j,k] ~ dbinom(p_sight[t,i,j,k] * w[t,i,j,k], z[t,i,j,k]) #w binary day/night
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
      #beta_phi = slope for increasing surv over time (older chicks have higher surv)
      #eps_phi = residuals
      #pi_phi = effect of SIC on survival
      #rho_phi = effect of KRILL on survival
      logit(phi[t,i,j,k]) <- mu_phi + gamma_phi[j] + eta_phi[k] #+ beta_phi*x[t] + pi_phi * SIC[j,k] + rho_phi * KRILL[j,k]


      #p = detection prob
      #mu_p = grand mean for all sites/years
      #beta_phi = slope for increasing detection over time (older chicks have higher detection p)
      #nu_p = effect of nest
      logit(p[t,i,j,k]) <- mu_p + beta_p*x[t] #+ nu_p[i,j,k]

      } #t
      } #i
      } #j
      } #k
      
      
      
      #priors - phi
      mu_phi ~ dnorm(0, 0.01)

      #beta_phi ~ dnorm(0, 1000) T(0,0.03)


      #covariates
      #pi_phi ~ dnorm(0, 0.386)
      #rho_phi ~ dnorm(0, 0.386)
      
      for (k in 1:NK)
      {
      #eta_phi[k] ~ dnorm(0, tau_eta_phi)
      eta_phi[k] ~ dnorm(0, 0.386)    

      #for (j in 1:NJ)
      #{
      #for (t in 1:NT)
      #{
      #eps_phi[t,j,k] ~ dnorm(0, tau_eps_phi) #T(-10,10)
      #}
      #}
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
      
      #tau_eps_phi <- pow(sigma_eps_phi, -2)
      #sigma_eps_phi ~ dunif(0.25, 8)
      
      
      
      #priors - p
      mu_p ~ dnorm(0, 0.5)
      
      beta_p ~ dnorm(0, 1000) T(0,0.03)
      
      # for (k in 1:NK)
      # {
      # for (j in 1:NJ)
      # {
      # for (i in 1:NI[j,k])
      # {
      # #nu_p[i,j,k] ~ dnorm(0, tau_nu_p) #T(-10,10)
      # nu_p[i,j,k] ~ dnorm(0, 0.386) #T(-10,10)
      # }
      # }
      # }
      
      #tau_nu_p <- pow(sigma_nu_p, -2)
      #sigma_nu_p ~ dunif(0, 5)
      
      }",fill = TRUE)

  sink()
}



# Starting values ---------------------------------------------------------

#create and fill array with initial values for nu_p
narray <- array(NA, dim = c(max(real_nests), DATA$NJ, DATA$NK))
for (k in 1:DATA$NK)
{
  #k <- 1
  for (j in 1:DATA$NJ)
  {
    #j <- 1
    for (i in 1:DATA$NI[j,k])
    {
      #i <- 1
      if (DATA$NI[j,k] > 0)
      {
        narray[i,j,k] <- 0
      }
    }
  }
}


Inits_1 <- list(mu_phi = 7,
                #beta_phi = 0.01,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                #pi_phi = 0,
                #rho_phi = 0,
                mu_p = -2,
                beta_p = 0.002,
                #nu_p = narray,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                #sigma_nu_p = 1.06,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mu_phi = 7,
                #beta_phi = 0.01,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                #pi_phi = 0,
                #rho_phi = 0,
                mu_p = -2,
                beta_p = 0.002,
                #nu_p = narray,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                #sigma_nu_p = 1.06,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mu_phi = 7,
                #beta_phi = 0.01,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                #pi_phi = 0,
                #rho_phi = 0,
                mu_p = -2,
                beta_p = 0.002,
                #nu_p = narray,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                #sigma_nu_p = 1.06,
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

Inits_4 <- list(mu_phi = 7,
                #beta_phi = 0.01,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                #pi_phi = 0,
                #rho_phi = 0,
                mu_p = -2,
                beta_p = 0.002,
                #nu_p = narray,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                #sigma_nu_p = 1.06,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 4)

Inits_5 <- list(mu_phi = 7,
                #beta_phi = 0.01,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                #pi_phi = 0,
                #rho_phi = 0,
                mu_p = -2,
                beta_p = 0.002,
                #nu_p = narray,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                #sigma_nu_p = 1.06,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 5)

Inits_6 <- list(mu_phi = 7,
                #beta_phi = 0.01,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                #pi_phi = 0,
                #rho_phi = 0,
                mu_p = -2,
                beta_p = 0.002,
                #nu_p = narray,
                #sigma_eta_phi = 0.78,
                #sigma_gamma_phi = 0.84,
                #sigma_nu_p = 1.06,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 6)

F_Inits <- list(Inits_1, Inits_2, Inits_3, Inits_4, Inits_5, Inits_6)



# Parameters to track -----------------------------------------------------

Pars <- c('mu_phi',
          'eta_phi',
          'gamma_phi',
          #'beta_phi',
          #'pi_phi',
          #'rho_phi',
          #'sigma_eta_phi',
          #'sigma_gamma_phi',
          'mu_p',
          'beta_p'#,
          #'phi',
          #'p',
          #'nu_p'#,
          #'sigma_nu_p'
          )


# Run model ---------------------------------------------------------------

# #make sure model compiles
# jagsRun(jagsData = DATA,
#         jagsModel = 'pwatch_surv.jags',
#         jagsInits = F_Inits,
#         DEBUG = TRUE)


jagsRun(jagsData = DATA, 
               jagsModel = 'pwatch_surv.jags',
               jagsInits = F_Inits,
               params = Pars,
               jagsID = 'May_12_2018_gamma_eta',
               jagsDsc = 'Four sites. logit(phi) <- mu_phi + gamma_phi + eta_phi; logit(p) <- mu_p + beta_p*x; No partial pooling. Long queue.',
               db_hash = 'Markrecap_data_15.05.18.csv',
               n_chain = 6,
               n_adapt = 5000,
               n_burn = 80000,
               n_draw = 60000,
               n_thin = 20,
               DEBUG = FALSE,
               EXTRA = FALSE,
               Rhat_max = 1.1,
               n_max = 100000,
               save_data = TRUE)

