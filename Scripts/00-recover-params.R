#################
# Penguin Watch Model - 0 - Recovering gen values using one detection param vs many
#
# 0-detect-params.R | recovering generating values using one detection param vs many detection many
# 00-recover-params.R | recovering generating values using realistic data/model
# 1-process-krill-data.R | process krill data
# 2-process-SIC-data.R | process SIC data
# 3-model.R | penguin model
# 3-run-model.pbs | pbs script to run penguin model on HPC resources
# 4-analyze-output.R | analyze model output
#
# Author: Casey Youngflesh
#################

#simulate with realistic values for mu_p, mu_phi, beta_p


# Clear environment -------------------------------------------------------

rm(list = ls())


# DIR ---------------------------------------------------------------------


#laptop
setwd('~/Google_Drive/R/penguin_watch_model/Scripts/')


#HPC
dir <- '../Results'



# Load packages -----------------------------------------------------------

#devtools::install_github('caseyyoungflesh/jagsRun')

library(dplyr)
library(jagsRun)
library(MCMCvis)
library(boot)



# specify param values -----------------------------------------------

#blank data.frame
OUT <- data.frame()

#generate 5 sets of data to attempt to recover parameters
for (j in 1:10)
{
  #j <- 1
  set.seed(j)
  
  #specify values
  n_ts <- 40 #number of time steps
  x <- 1:n_ts
  nests <- 15 #number of nests
  
  # Values from May_25_2018 run - first day aggregation
  #                 mean     sd    2.5%     50%   97.5% Rhat n.eff
  # mu_phi        4.9776 2.0232  1.0291  4.9873  8.9211 1.00  4130
  # mu_p          1.6106 0.5973  0.4612  1.6045  2.7366 1.01  2142
  # beta_p        0.2979 0.0259  0.2493  0.2973  0.3511 1.00 14975
  
  #use values from model fit
  mu_phi <- 4.9873
  mu_p <- 1.6045
  beta_p <- 0.2973
  
  x_sc <- scale(x, scale = FALSE)[,1] - 1
  p_sim <- inv.logit(x_sc*beta_p + mu_p)
  phi_sim <- inv.logit(mu_phi)
  
  
  # survival ----------------------------------------------------------------
  
  #survival probabilities
  PHI <- matrix(rep(phi_sim, (n_ts - 1) * nests),
                ncol = nests)
  
  # detection ---------------------------------------------------------------
  
  #setection probabilities
  P <- matrix(rep(p_sim, times = nests),
              ncol = nests)
  
  
  # simulate time series ----------------------------------------------------
  
  #function to simulate time series
  sim_data_fun <- function(PHI_MAT, P_MAT, N_NESTS, TYPE)
  {
    #PHI_MAT <- PHI
    #P_MAT <- P
    #N_NESTS <- nests
    
    TS_LEN <- NROW(PHI_MAT) + 1
    CH <- matrix(0, 
                 ncol = N_NESTS, 
                 nrow = TS_LEN)
    
    for (i in 1:N_NESTS)
    {
      #i <- 1
      
      if (TYPE == 'TWO')
      {
        CH[i,1] <- 2 
        t_SP <- c(2, rep(NA, TS_LEN-1))
      }
      if (TYPE == 'ONE')
      {
        CH[i,1] <- 1
        t_SP <- c(1, rep(NA, TS_LEN-1))
      }
      
      for (t in 2:TS_LEN)
      {
        #t <- 2
        #TRUE STATE
        t_SP[t] <- rbinom(1, size = t_SP[t-1], prob = PHI_MAT[t-1,i])
        
        #OBSERVED STATE
        t_DP <- rbinom(1, size = t_SP[t], prob = P_MAT[t-1, i])
        CH[t,i] <- t_DP
      }
    }
    return(CH)
  }
  
  sim_data <- sim_data_fun(PHI, P, nests, TYPE = 'TWO')
  
  
  #get indexes for NA values - split half/half, NA/actual obs
  # ct <- 1
  # idx <- c()
  # for (i in 1:length(x))
  # {
  #   if (ct <= 12)
  #   {
  #     idx <- c(idx, i)
  #   }
  #   ct <- ct + 1
  #   
  #   if (ct == 25)
  #   {
  #     ct <- 1
  #   }
  # }
  
  #insert NA values for first 10 days
  
  #insert NA values
  sim_data[1:10,] <- NA
  
  
  # Create z-state matrix ---------------------------------------------------
  
  #fun modified from Kerry and Schaub 2012
  #assume that chicks are alive since time step 1
  known.state.fun <- function(INPUT, TYPE)
  {
    #INPUT <- sim_data
    state <- INPUT
  
    for (i in 1:NCOL(INPUT))
    {
      #i <- 4
      #two chicks alive at time step one
      state[1, i] <- 2
      
      if (TYPE == 'TWO')
      {
        if (sum(state[,i] == 2, na.rm = TRUE) > 1)
        {
          n2 <- max(which(INPUT[,i] == 2))
          state[1:n2, i] <- 2
        }
      }
      if (TYPE == 'ONE')
      {
        if (sum(state[,i] == 1, na.rm = TRUE) > 1)
        {
          n2 <- max(which(INPUT[,i] == 1))
          state[1:n2, i] <- 1
        }
      }
    }
    if (TYPE == 'TWO')
    {
      state[state == 0] <- NA
      state[state == 1] <- NA
    }
    if (TYPE == 'ONE')
    {
      state[state == 0] <- NA
    }
    return(state)
  }
  
  z_matrix <- known.state.fun(sim_data, TYPE = 'TWO')
  
  
  # Create Data for JAGS ---------------------------------------------------------
  
  #nests_array:
  #dim1 (rows) [t] = time steps
  #dim2 (cols) [i] = nests
  
  DATA <- list(
    y = sim_data, #response
    NI = dim(sim_data)[2], #number of nests
    NT = dim(sim_data)[1], #number of time steps
    z = z_matrix, #known points of bird being alive
    x = scale(as.numeric(1:dim(sim_data)[1]), scale = FALSE)[,1]) #time steps for increase in surv/detection over time   
  
  
  setwd(dir)
  
  # Model -------------------------------------------------------------------
  
  {
    sink("pwatch_sim.jags")
    
    cat("
        
        model {
        
        #nests
        for (i in 1:NI)
        {
        #both chicks alive at time step 1 (z[1,i] = 2)
        
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
        #phi = survival prob
        #mu_phi = grand mean for all sites/years
        
        logit(phi[t,i]) <- mu_phi
        
        #p = detection prob
        #mu_p = grand mean for all sites/years
        #beta_phi = slope for increasing detection over time (older chicks have higher detection p)
        
        logit(p[t,i]) <- mu_p + beta_p * x[t]
        
        } #t
        } #i
        
        #priors - phi
        mu_phi ~ dnorm(3, 0.1)
        
        #priors - p
        mu_p ~ dnorm(2, 0.1)
        beta_p ~ dnorm(0.1, 10) T(0, 0.5)
        
        }",fill = TRUE)
  
    sink()
  }
  
  
  
  # Starting values ---------------------------------------------------------
  
  
  Inits_1 <- list(mu_phi = 4,
                  mu_p = 1.5,
                  beta_p = 0.1,
                  .RNG.name = "base::Mersenne-Twister",
                  .RNG.seed = 1)
  
  Inits_2 <- list(mu_phi = 4,
                  mu_p = 1.5,
                  beta_p = 0.1,
                  .RNG.name = "base::Wichmann-Hill",
                  .RNG.seed = 2)
  
  Inits_3 <- list(mu_phi = 4,
                  mu_p = 1.5,
                  beta_p = 0.1,
                  .RNG.name = "base::Marsaglia-Multicarry",
                  .RNG.seed = 3)
  
  F_Inits <- list(Inits_1, Inits_2, Inits_3)
  
  
  
  # Parameters to track -----------------------------------------------------
  
  Pars <- c('mu_phi',
            'mu_p',
            'beta_p')
  
  
  # Run model ---------------------------------------------------------------
  
  #make sure model compiles
  # jagsRun(jagsData = DATA,
  #         jagsModel = 'pwatch_sim.jags',
  #         jagsInits = F_Inits,
  #         DEBUG = TRUE)
  
  out <- jagsRun(jagsData = DATA, 
          jagsModel = 'pwatch_sim.jags',
          jagsInits = F_Inits,
          params = Pars,
          n_chain = 3,
          n_adapt = 5000,
          n_burn = 10000,
          n_draw = 10000,
          n_thin = 10,
          report = FALSE)
  
  
  out_mu_phi <- MCMCsummary(out, params = 'mu_phi')
  out_mu_p <- MCMCsummary(out, params = 'mu_p')
  out_beta_p <- MCMCsummary(out, params = 'beta_p')
  m_rhat <- max(MCMCsummary(out)[,6])
  
  
  
  if (mu_phi > inv.logit(MCMCsummary(out, params = 'mu_phi')[3]) &
      mu_phi < inv.logit(MCMCsummary(out, params = 'mu_phi')[5]))
  {
    ret_phi <- TRUE
  } else {
    ret_phi <- FALSE
  }
  
  if (mu_p > inv.logit(MCMCsummary(out, params = 'mu_p')[3]) &
      mu_p < inv.logit(MCMCsummary(out, params = 'mu_p')[5]))
  {
    ret_p <- TRUE
  } else {
    ret_p <- FALSE
  }
  
  if (beta_p > inv.logit(MCMCsummary(out, params = 'beta_p')[3]) &
      beta_p < inv.logit(MCMCsummary(out, params = 'beta_p')[5]))
  {
    ret_beta_p <- TRUE
  } else {
    ret_beta_p <- FALSE
  }
  
  temp <- data.frame(RUN = j,
                     gen_mu_phi = mu_phi,
                     RET_phi = ret_phi,
                     model_mu_phi_LCI = out_mu_phi[3],
                     model_mu_phi_MED = out_mu_phi[4],
                     model_mu_phi_UCI = out_mu_phi[5],
                     gen_mu_p = mu_p,
                     RET_p = ret_p,
                     model_mu_p_LCI = out_mu_p[3],
                     model_mu_p_MED = out_mu_p[4],
                     model_mu_p_UCI = out_mu_p[5],
                     gen_beta_p = beta_p,
                     RET_beta_p = ret_beta_p,
                     model_beta_p_LCI = out_beta_p[3],
                     model_beta_p_MED = out_beta_p[4],
                     model_beta_p_UCI = out_beta_p[5],
                     max_Rhat = m_rhat
                     )
  
  OUT <- rbind(OUT, temp)
}


sink(paste0('simulation_results.txt'))
print(OUT)
sink()

