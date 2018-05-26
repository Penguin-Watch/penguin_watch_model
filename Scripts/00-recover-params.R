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

## Simulation_results_XXX.txt
# *number of time steps doesn't impact recovery of phi param (or CI on proportion of chicks alive at end of season)
#   -it does impact recovery of p and beta_p params

## p_beta_simulation_results.txt
# *gen values for mu_p and beta_p don't impact recovery of phi param (or CI on proportion of chicks alive at end of season)
#   -it does impact recovery of p param

## nests_simulation_results.txt
# *number of nests DOES NOT impact recovery of phi param (always recovers) but more nests do shrink CI on proportion of chicks alive at end of season
#   -more nests leads to worse recovery of mu_p for some reason (does not substantialy change the shape of the detection curve over the course of the season though - lines 126-128)



# Clear environment -------------------------------------------------------

rm(list = ls())


# DIR ---------------------------------------------------------------------


#laptop
#setwd('~/Google_Drive/R/penguin_watch_model/Scripts/')

#desktop
#setwd('~/gdrive/R/penguin_watch_model/Scripts/')

#HPC
dir <- '~/gdrive/R/penguin_watch_model/Results/'



# Load packages -----------------------------------------------------------

#devtools::install_github('caseyyoungflesh/jagsRun')

library(dplyr)
library(jagsRun)
library(MCMCvis)
library(boot)



# specify param values -----------------------------------------------

#values to cycle through
#detection
dp <- c(0.80, 0.90, 0.99)
#beta_p
bv <- c(0.125, 0.08, 0.05)
#nests
nn <- c(5, 15, 25, 35)

#blank data.frame
OUT <- data.frame()


#progress bar for det and beta_p
#pb <- txtProgressBar(min = 0, max = length(dp)*length(bv)*5, style = 3)
#progress bar for nests
pb <- txtProgressBar(min = 0, max = length(nn)*5, style = 3)
pbi <- 1

# for (m in 1:length(dp))
# {
#   #m <- 2
#   for (l in 1:length(bv))
#   {
    #l <- 3
  for (q in 1:length(nn))
  {
    out_ret_phi <- c()
    out_ret_p <- c()
    out_ret_beta_p <- c()
    out_surv_diff <- c()
    out_rhat <- c()
    
    #generate 5 sets of data to attempt to recover parameters
    for (j in 1:5)
    {
      #j <- 1
      set.seed(j)
    
      #number of time steps
      n_ts <- 160 #values$NTS[m] #number of time steps (two on, two off; first 40 are NA)
      x <- 1:n_ts
      period <- ifelse(n_ts == 40, 0, n_ts/80) #number of 'dark' time steps between observable periods
      ep <- ifelse(period == 0, 10, period * 2 * 10) #number of NA before data starts (egg period)
      
      #colony nest and survival params
      nests <- nn[q] #number of nests
      prop_surv <- 0.85 #proportion of chicks alive at creche
      
      #specify detection params
      detect <- 0.8 #dp[m] #actual detection prob
      beta_p <- 0.08 #bv[l] #values$BETA_P[m] #detection slope #40 = 0.4, 160 = 0.1, 320 = 0.05, 640 = 0.025
      
      mu_phi <- logit(prop_surv^(1/n_ts))
      mu_p <- logit(detect)
      
      x_sc <- scale(x, scale = FALSE)[,1] - 1
      p_sim <- inv.logit(x_sc*beta_p + mu_p) #mu_p shifts curve L/R; beta_p flattens curve
      ms <- max(diff(p_sim)) #max slope
      ms_idx <- which.max(diff(p_sim))
      # plot(p_sim, type = 'l')
      # #abline(v = ms_idx, col = 'red')
      # rect(0, 0, ep, 1, col = rgb(0,0,1,0.5))

      # Values from May_25_2018 run - first day aggregation
      #                 mean     sd    2.5%     50%   97.5% Rhat n.eff
      # mu_phi        4.9776 2.0232  1.0291  4.9873  8.9211 1.00  4130
      # mu_p          1.6106 0.5973  0.4612  1.6045  2.7366 1.01  2142
      # beta_p        0.2979 0.0259  0.2493  0.2973  0.3511 1.00 14975
      
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
      
      
      #insert NA values for 'egg period'
      
      #insert NA values
      sim_data[1:ep,] <- NA
      
      
      #get indexes for NA values (periods where camera was off)
    
      if (period > 0)
      {
        ct <- 1
        idx <- c()
        for (i in 1:length(x))
        {
          if (ct <= period)
          {
            idx <- c(idx, i)
          }
          ct <- ct + 1
      
          if (ct == ((period *2)+1))
          {
            ct <- 1
          }
        }
        sim_data[idx,] <- NA
      }
      
      
      
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
            mu_phi ~ dnorm(5, 0.1)
            
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
      

      
      if (mu_phi > MCMCsummary(out, params = 'mu_phi')[3] &
          mu_phi < MCMCsummary(out, params = 'mu_phi')[5])
      {
        ret_phi <- TRUE
      } else {
        ret_phi <- FALSE
      }
      
      if (mu_p > MCMCsummary(out, params = 'mu_p')[3] &
          mu_p < MCMCsummary(out, params = 'mu_p')[5])
      {
        ret_p <- TRUE
      } else {
        ret_p <- FALSE
      }
      
      if (beta_p > MCMCsummary(out, params = 'beta_p')[3] &
          beta_p < MCMCsummary(out, params = 'beta_p')[5])
      {
        ret_beta_p <- TRUE
      } else {
        ret_beta_p <- FALSE
      }
      
      
      surv_diff <- inv.logit(out_mu_phi[5])^n_ts - inv.logit(out_mu_phi[3])^n_ts
      
      out_ret_phi <- c(out_ret_phi, ret_phi)
      out_ret_p <- c(out_ret_p, ret_p)
      out_ret_beta_p <- c(out_ret_beta_p, ret_beta_p)
      out_surv_diff <- c(out_surv_diff, surv_diff)
      out_rhat <- c(out_rhat, m_rhat)
      
      #advance progress bar
      setTxtProgressBar(pb, pbi)
      pbi <- pbi + 1
    }#end repetition loop
    
    s_ret_phi <- sum(out_ret_phi)/5
    s_ret_p <- sum(out_ret_p)/5
    s_ret_beta_p <- sum(out_ret_beta_p)/5
    m_surv_diff <- mean(out_surv_diff)
    max_rhat <- max(out_rhat)
    
    temp <- data.frame(time_steps = n_ts,
                       period = period,
                       nests = nn[q],
                       GEN_mu_phi = mu_phi,
                       PRC_phi = s_ret_phi,
                       GEN_pca = prop_surv,
                       RNG_pca = m_surv_diff,
                       GEN_mu_p = mu_p,
                       GEN_detect = detect,
                       PRC_p = s_ret_p,
                       GEN_beta_p = beta_p,
                       PRC_beta_p = s_ret_beta_p,
                       max_Rhat = max_rhat)
    
    OUT <- rbind(OUT, temp)
  }
 # }
#} 

close(pb)

#PRC = proportion retrieved correctly
#GEN = generating value
#pca = proportion chicks alive (at end of season)
#RNG = mean range between 97.5 and 0.25 posterior quantiles

sink(paste0('nests_simulation_results.txt'))
print(OUT)
sink()




MCMCsummary(out)
ch_mu_phi <- MCMCchains(out, 
                        params = 'mu_phi')

pal <- apply(ch_mu_phi, 2, function(x) inv.logit(x)^160)
hist(pal)
quantile(pal, probs = c(0.025, 0.975))

#mu_p
PR <- rnorm(10000, 2, 1/sqrt(0.1))
MCMCtrace(out, 
          params = 'mu_p',
          ind = TRUE,
          priors = PR,
          post_zm = FALSE,
          pdf = FALSE)

#beta_p
PR <- rnorm(10000, 0.1, 1/sqrt(10))
MCMCtrace(out, 
          params = 'beta_p',
          ind = TRUE,
          priors = PR,
          post_zm = FALSE,
          pdf = FALSE)

#mu_phi
PR <- rnorm(10000, 5, 1/sqrt(0.1))
MCMCtrace(out, 
          params = 'mu_phi',
          ind = TRUE,
          priors = PR,
          post_zm = FALSE,
          pdf = FALSE)
