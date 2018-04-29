######################
#Mark-recapture for penguin watch 
#
#script to run model - model output saved to Results dir
#
#Authors: Casey Youngflesh
######################



# Clear environment -------------------------------------------------------

rm(list = ls())



# Load packages -----------------------------------------------------------

#devtools::install_github('caseyyoungflesh/jagsRun')

library(abind)
library(dplyr)
library(jagsRun)



# determine PW dates to use -----------------------------------------------

#ensures that row dimension (time steps within season) will have the same dimension

#setwd('~/Google_Drive/R/penguin_watch_model/Data/PW_data/RAW_Fiona_Apr_15_2018/')
setwd('../Data')

PW_data_p <- read.csv('Markrecap_data_15.05.18.csv', stringsAsFactors = FALSE)

#merge site and cam for independent cameras (treating each cam as a different site)
PW_data <- mutate(PW_data_p, site_p_cam = paste0(site, camera_letter))

#remove HALF bc it's not a full season
un_sites_p <- unique(PW_data$site_p_cam)
un_sites <- un_sites_p[which(un_sites_p != 'HALFb')]

#first date of each season
yrs <- c()
min_date <- as.Date(NA)
for (k in 1:length(un_sites))
{
  #k <- 1
  temp <- filter(PW_data, site_p_cam == un_sites[k])
  un_yrs <- unique(temp$season_year)
  
  yrs <- c(yrs, un_yrs)
  for (j in 1:length(un_yrs))
  {
    #j <- 1
    temp2 <- filter(temp, season_year == un_yrs[j])
    
    temp_dates <- as.Date(temp2$datetime, format = "%Y:%m:%d %H:%M:%S")
    t_min_date <- min(temp_dates)
    min_date <- c(min_date, t_min_date)
  }
}


#first date of season to use across all years (add 1 to start on full day)
f_min_date <- format((min_date+1), '%m-%d')
first_date <- max(f_min_date, na.rm = TRUE)

#use Feb 1 as last date of season to use across all years
last_date <- format(as.Date('02-01', format = '%m-%d'), '%m-%d')




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


#number of time steps (rows) in response data (-1 to account for greater than FIRST and less than LAST) - t
n_ts <- as.numeric(as.Date(paste0('2018', '-', last_date, format = '%Y-%m-%d')) - 
                     as.Date(paste0('2017', '-', first_date, format = '%Y-%m-%d')) - 1) * 24

#number of nests (columns) in response data - i
n_nests <- length(tog2)

#number of years (3rd dim) in response data - j
d_yrs <- sort(yrs[!duplicated(yrs)])
n_yrs <- length(d_yrs)

#number of sites (4th dim) in response data - k
n_sites <- length(un_sites)

#create blank array
nests_array <- array(NA, dim = c(n_ts, n_nests, n_yrs, n_sites))


#fill response data array
for (k in 1:n_sites)
{
  #k <- 1
  temp <- filter(PW_data, site_p_cam == un_sites[k])
  
  for (j in 1:n_yrs)
  {
    #j <- 2
    temp2 <- filter(temp, season_year == d_yrs[j])
    if (NROW(temp2) > 0)
    {
      temp_dates <- as.Date(temp2$datetime, format = "%Y:%m:%d %H:%M:%S")
      
      #first and last days in season
      FIRST <- as.Date(paste0((d_yrs[j]-1), '-', first_date), format = "%Y-%m-%d")
      LAST <- as.Date(paste0(d_yrs[j], '-', last_date), format = "%Y-%m-%d")
      valid_dates <- which(temp_dates > FIRST & temp_dates < LAST)
      
      #appropriate date range and appropriate columns for nests
      nests_array[,,j,k] <- as.matrix(temp2[valid_dates, tog2])
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

#create w_array (day and night)

w_array <- nests_array
#assign 1 to values with observations (either 0/1/2 for nests_array)
w_array[which(!is.na(w_array), arr.ind = TRUE)] <- 1

#only fille in 0 if this was a night (there is data on either side and not just missing nests)
for (k in 1:dim(nests_array)[4])
{
  #k <- 1
  for (j in 1:dim(nests_array)[3])
  {
    #j <- 2
    #if there are more than 0 nests in that year (there is data for that site/year)
    if (real_nests[j,k] > 0)
    {
      sub_i <- real_nests[j,k]
      w_array[cbind(which(is.na(w_array[,1:sub_i,j,k]), arr.ind = TRUE), j, k)] <- 0
    }
  }
}

#checks:
#nests_array[1:15, 1:10, 2, 1]
#w_array[1:15, 1:10, 2, 1]



# observations with > 2 chicks --------------------------------------------

# 1% of nest observations have more than 2 chicks in them
# which(DATA$y > 2, arr.ind = TRUE)
# num_o2 <- length(DATA$y[which(DATA$y > 2, arr.ind = TRUE)])
# total <- length(DATA$y[which(!is.na(DATA$y), arr.ind = TRUE)])
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
  x = as.numeric(1:dim(nests_array)[1])) #time steps for increase in surv/detection over time 

#ADD TO DATA
#sea ice for that year (matrix) - site/year
#krill catch in buffer for that year (matrix) - site/year
#ADD 'CAMERA' number in here somewhere (detection should not be the same for each camera at a site)


#setwd('~/Google_Drive/R/penguin_watch_model/Results/')
setwd('../Results')

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
      logit(phi[t,i,j,k]) <- mu_phi + eta_phi[k] + gamma_phi[j] + beta_phi*x[t]# + eps_phi[t,j,k] #+ pi_phi * SIC[j,k] + rho_phi * KRILL[j,k]


      #p = detection prob
      #mu_p = grand mean for all sites/years
      #beta_phi = slope for increasing detection over time (older chicks have higher detection p)
      #nu_p = effect of nest
      logit(p[t,i,j,k]) <- mu_p + beta_p*x[t] + nu_p[i,j,k]

      } #t
      } #i
      } #j
      } #k
      
      
      
      #priors - phi
      #Lunn prior - flat in probability space
      mu_phi ~ dnorm(0, 0.386)   

      beta_phi ~ dnorm(0, 1000) T(0,1) #[slope only pos] maybe variance 0.01 (precision 100) - plot histogram to get a look (will depend on time step length [i.e., one hour or one day])

      #covariates
      #pi_phi ~ dnorm(0, 0.01)
      #rho_phi ~ dnorm(0, 0.01)
      
      for (k in 1:NK)
      {
      eta_phi[k] ~ dnorm(0, tau_eta_phi)
      
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
      gamma_phi[j] ~ dnorm(0, tau_gamma_phi)
      }
      
      tau_eta_phi <- pow(sigma_eta_phi, -2)
      sigma_eta_phi ~ dunif(0.25, 3)
      
      tau_gamma_phi <- pow(sigma_gamma_phi, -2)
      sigma_gamma_phi ~ dunif(0.25, 3)
      
      #tau_eps_phi <- pow(sigma_eps_phi, -2)
      #sigma_eps_phi ~ dunif(0.25, 3)
      
      
      
      #priors - p
      #Lunn prior - flat in probability space
      mu_p ~ dnorm(0, 0.386)
      
      beta_p ~ dnorm(0, 10) T(0,1) #[slope only pos] maybe variance 0.1 (precision 10) - plot histogram to get a look (will depend on time step length [i.e., one hour or one day])
      
      for (k in 1:NK)
      {
      for (j in 1:NJ)
      {
      for (i in 1:NI[j,k])
      {
      nu_p[i,j,k] ~ dnorm(0, tau_nu_p) #T(-10,10)
      }
      }
      }
      
      tau_nu_p <- pow(sigma_nu_p, -2)
      sigma_nu_p ~ dunif(0.25, 3)
      
      }",fill = TRUE)

  sink()
}



# Starting values ---------------------------------------------------------


Inits_1 <- list(mu_phi = 0,
                beta_phi = 0.1,
                sigma_eta_phi = 0.26,
                sigma_gamma_phi = 0.26,
                #sigma_eps_phi = 0.26,
                mu_p = 0,
                beta_p = 0,
                sigma_nu_p = 0.26,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mu_phi = 0,
                beta_phi = 0.1,
                sigma_eta_phi = 0.27,
                sigma_gamma_phi = 0.27,
                #sigma_eps_phi = 0.27,
                mu_p = 0,
                beta_p = 0,
                sigma_nu_p = 0.27,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mu_phi = 0,
                beta_phi = 0.1,
                sigma_eta_phi = 0.28,
                sigma_gamma_phi = 0.28,
                #sigma_eps_phi = 0.28,
                mu_p = 0,
                beta_p = 0,
                sigma_nu_p = 0.28,
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

Inits_4 <- list(mu_phi = 0,
                beta_phi = 0.1,
                sigma_eta_phi = 0.29,
                sigma_gamma_phi = 0.29,
                #sigma_eps_phi = 0.28,
                mu_p = 0,
                beta_p = 0,
                sigma_nu_p = 0.29,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 4)

Inits_5 <- list(mu_phi = 0,
                beta_phi = 0.1,
                sigma_eta_phi = 0.30,
                sigma_gamma_phi = 0.30,
                #sigma_eps_phi = 0.28,
                mu_p = 0,
                beta_p = 0,
                sigma_nu_p = 0.30,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 5)

F_Inits <- list(Inits_1, Inits_2, Inits_3, Inits_4, Inits_5)



# Parameters to track -----------------------------------------------------

Pars <- c('mu_phi',
          'eta_phi',
          'gamma_phi',
          'beta_phi',
          #'eps_phi',
          'sigma_eta_phi',
          'sigma_gamma_phi',
          #'sigma_eps_phi',
          'mu_p',
          'beta_p',
          'nu_p',
          'sigma_nu_p')


# Run model ---------------------------------------------------------------


jagsRun(jagsData = DATA, 
               jagsModel = 'pwatch_surv.jags',
               jagsInits = F_Inits,
               params = Pars,
               jagsID = 'April_29_2018',
               jagsDsc = 'Remove eps_phi and eta_p',
               db_hash = 'Markrecap_data_15.05.18.csv',
               n_chain = 5,
               n_adapt = 8000,
               n_burn = 20000,
               n_draw = 20000,
               n_thin = 10,
               DEBUG = FALSE,
               EXTRA = FALSE,
               Rhat_max = 1.1,
               n_max = 100000)


