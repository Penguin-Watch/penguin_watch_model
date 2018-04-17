######################
#Mark-recapture for penguin watch 
#
#script to run model - model object saved as .rds
#
#Authors: Casey Youngflesh
######################


#QUESTIONS about model
#want to know how survival is changing. Add env covariates onto this
#no replication over period of closure (how to treat night - all images over one day as a single observation? should each hour be a time step with NAs over night?) - are the parameters identifiable? p219 Kerry and Schaub 2012
#should the random effect (time) for phi be in there? Is it necessary?


#how to model more than one site? maybe code each nest/site as a different nest (one site)


#nonidentifiability checks
#plot p[1,1] against phi[1,1] - CHECK
#check correlation p and phi - CHECK
#compare estimates p and phi against true p and phi - 
#plot posteriors for:
#1) betas - dnorm(0, 1) T(-1,1) - CHECK
#2) mean_phi and mean_p - dunif(0,0.1) - CHECK
#3) eps_phi and eps_p - dnorm(0, tau_p) T(-20,20) - p okay, phi maybe weakly identifiable
#and compare to priors


# Clear environment -------------------------------------------------------

rm(list = ls())



# Load packages -----------------------------------------------------------

# if('pacman' %in% rownames(installed.packages()) == FALSE)
# {
#   install.packages('pacman')
# }
# 
# pacman::p_load(rjags, MCMCvis, parallel)
#devtools::install_github('caseyyoungflesh/jagsRun')

require(jagsRun)
require(abind)




# determine PW dates to use -----------------------------------------------



#remove HALF bc it's not a full season
un_sites_p <- unique(PW_data$site)
un_sites <- un_sites_p[which(un_sites_p != 'HALF')]

min_date <- as.Date(NA)
for (k in 1:length(un_sites))
{
  #k <- 1
  temp <- filter(PW_data, site == un_sites[k])
  un_yrs <- unique(temp$season_year)
  
  for (j in 1:length(un_yrs))
  {
    #j <- 1
    temp2 <- filter(temp, season_year == un_yrs[j])
    
    temp_dates <- as.Date(temp2$datetime, format = "%Y:%m:%d %H:%M:%S")
    t_min_date <- min(temp_dates)
    min_date <- c(min_date, t_min_date)
  }
}


#first date of season to use (add 1 to start on full day)
f_min_date <- format((min_date+1), '%m-%d')
first_date <- max(f_min_date, na.rm = TRUE)

#use Feb 1 as last date of season to use
last_date <- format(as.Date('02-01', format = '%m-%d'), '%m-%d')




# Create PW array -------------------------------------------

#just nest time series columns 
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

PW_nests <- PW_data[,tog2]


dims <- c()
for (k in 1:length(un_sites))
{
  #k <- 5
  temp <- filter(PW_data, site == un_sites[k])
  un_yrs <- unique(temp$season_year)
  
  for (j in 1:length(un_yrs))
  {
    #j <- 2
    temp2 <- filter(temp, season_year == un_yrs[j])
    temp_dates <- as.Date(temp2$datetime, format = "%Y:%m:%d %H:%M:%S")
    
    #first and last days in season
    FIRST <- as.Date(paste0((un_yrs[j]-1), '-', first_date), format = "%Y-%m-%d")
    LAST <- as.Date(paste0(un_yrs[j], '-', last_date), format = "%Y-%m-%d")
    
    valid_dates <- which(temp_dates > FIRST & temp_dates < LAST)
    
    #appropriate date range and appropriate columns for nests
    temp_nests <- temp2[valid_dates, tog2]
    dims <- rbind(dims, dim(temp_nests))
  }
}

dims






max_rows <- max(dims[,1])

# Create PW array ---------------------------------------------------------


nests_array <- c()
un_sites <- unique(PW_data$site)
for (k in 1:length(un_sites))
{
  #k <- 2
  temp <- filter(PW_data, site == un_sites[k])
  un_yrs <- unique(temp$season_year)
  
  years_array <- c()
  for (j in 1:length(un_yrs))
  {
    #j <- 1
    temp2 <- filter(temp, season_year == un_yrs[j])
    temp_nests <- temp2[,tog2]
    years_array <- abind(years_array, temp_nests, along = 3)
  }
  
  nests_array <- abind(nests_array, years_array, along = 4)
}






#structure into 4 dimensional array ([X,,,] = i (nest_id), [,X,,] = t (time_step), 
#[,,X,] = j (year), [,,,X] = k (site))




DATA <- list(
  y = sim_data, #reponse
  N = NROW(sim_data), #number of nests
  L = NCOL(sim_data), #number of time points
  z = z_vals,
  x = 1:NCOL(sim_data),
  w = w_mat) #binary day (1)/night (0)



# Model -------------------------------------------------------------------

#L = length of time series (400)
#N = number of nests (5)

{
sink("pwatch_surv.jags")

cat("
    model {
  
    for (i in 1:N)
    {
    #both chicks alive at time step 1
    z[i,1] <- 2
    y.new[i,1] ~ dbinom(p[i,1], 2)
    
    for (t in 2:L)
    {
    
    #state model
    z[i,t] ~ dbinom(p_alive[i,t], z[i,t-1])
    p_alive[i,t] <- ifelse(z[i,t-1] < 2, 
                        phi[i,t] * z[i,t-1],
                        phi[i,t])
    
    #observation model
    y[i,t] ~ dbinom(p_sight[i,t] * w[i,t], z[i,t]) #binary day/night
    p_sight[i,t] <- ifelse(z[i,t] < 2,
                        p[i,t] * z[i,t],
                        p[i,t])
    
    
    #PPC
    y.new[i,t] ~ dbinom(p_sight[i,t] * w[i,t], z[i,t]) #binary day/night
    }
    }
    
    #PPC
    #mean
    mn.y <- mean(y)
    mn.y.new <- mean(y.new)
    pv.mn <- step(mn.y.new - mn.y)
    
    #sd
    sd.y <- sd(y)
    sd.y.new <- sd(y.new)
    pv.sd <- step(sd.y.new - sd.y)
    
    
    #transforms
    for (i in 1:N)
    {
    for (t in 1:L)
    {
    logit(phi[i,t]) <- mu_phi + beta_phi*x[t] + eps_phi[t]       #phi = survival prob
    logit(p[i,t]) <- mu_p + beta_p*x[t] + eps_p[i]            #p = detection prob
    }
    }
    
    
    #priors
    for (t in 1:L)
    {
    eps_phi[t] ~ dnorm(0, tau_phi) T(-10,10)
    }
    
    mean_phi ~ dbeta(1.5,1.5)                 #Mean survival
    mu_phi <- log(mean_phi / (1 - mean_phi))
    tau_phi <- pow(sigma_phi, -2)
    sigma_phi ~ dunif(0.25, 3)
    
    
    for (i in 1:N)
    {
    eps_p[i] ~ dnorm(0, tau_p) T(-10,10)
    }
    
    
    mean_p ~ dbeta(1.5,1.5)                    #Mean detection - could use alternative below
    mu_p <- log(mean_p / (1 - mean_p))         #Logit transform - could use alternative below
    #mean_p <- 1 / (1+exp(-mu_p))              #Mean detection - Inv-logit transform
    #mu_p ~ dnorm(0, 0.001)                    #Prior for logit of mean survival
    tau_p <- pow(sigma_p, -2)
    sigma_p ~ dunif(0.25, 3)
    
    beta_phi ~ dnorm(0, 1000) T(0,1) #[slope only pos] maybe variance 0.01 (precision 100) - plot histogram to get a look (will depend on time step length [i.e., one hour or one day])
    beta_p ~ dnorm(0, 100) T(0,1) #[slope only pos] maybe variance 0.1 (precision 10) - plot histogram to get a look (will depend on time step length [i.e., one hour or one day])
    
    
    }",fill = TRUE)

sink()
}



# Starting values ---------------------------------------------------------


Inits_1 <- list(mean_phi = 0.5,
                mean_p = 0.5,
                sigma_phi = 0.26,
                sigma_p = 0.26,
                beta_phi = 0.1,
                beta_p = 0,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mean_phi = 0.6,
                mean_p = 0.3,
                sigma_phi = 1,
                sigma_p = 1,
                beta_phi = 0.1,
                beta_p = 0,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mean_phi = 0.7,
                mean_p = 0.4,
                sigma_phi = 1.5,
                sigma_p = 1.5,
                beta_phi = 0.1,
                beta_p = 0,
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)



# Parameters to track -----------------------------------------------------

Pars <- c('mean_phi',
          'mean_p',
          'sigma_p',
          'sigma_phi',
          'beta_p',
          'beta_phi',
          'mu_phi',
          'mu_p',
          'eps_phi',
          'eps_p',
          'p',
          'phi',
          'pv.mn',
          'pv.sd')


# Run model ---------------------------------------------------------------

jagsRun(jagsData = DATA,
        jagsModel = 'pwatch_surv.jags',
        jagsInits = F_Inits,
        params = Pars,
        jagsID = 'April_17_2018',
        jagsDsc = 'First go with real data',
        n_chain = 3,
        n_adapt = 8000,
        n_burn = 50000,
        n_draw = 50000,
        n_thin = 10,
        DEBUG = TRUE)


