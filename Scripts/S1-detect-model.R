######################
#Mark-recapture for penguin watch 
#
#Model detection as a group and for individual nests using simulated data
#
#Authors: Casey Youngflesh
######################

#model one survival
#model one detection/model individual nest detection


# Clear environment -------------------------------------------------------

rm(list = ls())


# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman')
}

pacman::p_load(rjags, MCMCvis, parallel, jagsRun, boot)




# colony parameters -------------------------------------------------------

#simulate new data - script modified from Kerry and Schaub 2012
n_ts <- 150 #number of time steps
x <- 1:n_ts
nests <- 15 #number of nests



#run through J different parameter sets, generating J different datasets per parameter set
J <- 3
I <- 5

#progress bar
total <- I * J
pb <- txtProgressBar(min = 0, max = total, style = 3)
PB <- 0

#loop through
OUT <- data.frame()
for (j in 1:J)
{
  #j <- 1

  set.seed(j)
  #generating survival probabilities
  dphi <- runif(1, 0.99, 0.9999)
  #generating detection probabilities
  dp <- runif(nests, 0.3, 0.8)
  
  for (i in 1:I)
  {
    # survival ----------------------------------------------------------------
    
    #survival probabilities
    PHI <- matrix(rep(dphi, (n_ts - 1) * nests),
                  ncol = nests)
    
    
    # detection ---------------------------------------------------------------
    
    #setection probabilities
    P <- matrix(rep(dp, each = n_ts - 1),
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
    
    sim_data <- sim_data_fun(PHI, P, nests, TYPE = 'ONE')
    
    
    # Create z-state matrix ---------------------------------------------------
    
    #fun modified from Kerry and Schaub 2012
    #assume that chicks are alive since time step 1
    known.state.fun <- function(INPUT, TYPE)
    {
      #INPUT <- DATA$y
      state <- INPUT
      for (i in 1:NCOL(INPUT))
      {
        n1 <- 1
        
        if (TYPE == 'TWO')
        {
          if (sum(state[,i] == 2) > 0)
          {
            n2 <- max(which(INPUT[,i] == 2))
            state[n1:n2, i] <- 2
          }
          
          state[n1, i] <- NA
        }
        if (TYPE == 'ONE')
        {
          if (sum(state[,i] == 1) > 0)
          {
            n2 <- max(which(INPUT[,i] == 1))
            state[n1:n2, i] <- 1
          }
          
          state[n1, i] <- NA
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
    
    z_matrix <- known.state.fun(sim_data, TYPE = 'ONE')
    
    
    # construct data ----------------------------------------------------------
    
    DATA <- list(
      y = sim_data, #reponse
      NI = dim(sim_data)[2], #number of nests
      NT = dim(sim_data)[1], #number of time steps
      z = z_matrix) #known points of bird being alive
    
    
    # Model -------------------------------------------------------------------
    
    {
      sink("pwatch_sim.jags")
      
      cat("
          
          model {
          
          #nests
          for (i in 1:NI)
          {
          #chick alive at time step 1
          z[1,i] <- 1
          
          #time step
          for (t in 2:NT)
          {
          #state model
          z[t,i] ~ dbinom(p_alive[t,i], z[t-1,i])
          p_alive[t,i] <- phi[t,i] * z[t-1,i]
          
          #observation model
          y[t,i] ~ dbinom(p_sight[t,i], z[t,i])
          p_sight[t,i] <- p[t,i] * z[t,i]
          
          }
          }
          
          
          #transforms
          #nests
          for (i in 1:NI)
          {
          #time
          for (t in 1:NT)
          {
          logit(phi[t,i]) <- mu_phi
          logit(p[t,i]) <- mu_p + nu_p[i]
          } #t
          } #i
          
          #priors - phi
          mu_phi ~ dnorm(0, 0.01)
          
          #priors - p
          mu_p ~ dnorm(0, 0.5)
          
          for (i in 1:NI)
          {
          nu_p[i] ~ dnorm(0, tau_nu_p)
          }
          
          tau_nu_p ~ dunif(0, 10) 
          
          il_mu_phi <- ilogit(mu_phi)
  
  
          }",fill = TRUE)
  
      sink()
    }
    
    
    # Starting values ---------------------------------------------------------
    
    
    Inits_1 <- list(mu_phi = 5,
                    mu_p = 0,
                    tau_nu_p = 1,
                    .RNG.name = "base::Mersenne-Twister",
                    .RNG.seed = 1)
    
    Inits_2 <- list(mu_phi = 5,
                    mu_p = 0,
                    tau_nu_p = 1,
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = 2)
    
    Inits_3 <- list(mu_phi = 5,
                    mu_p = 0,
                    tau_nu_p = 1,
                    .RNG.name = "base::Marsaglia-Multicarry",
                    .RNG.seed = 3)
    
    F_Inits <- list(Inits_1, Inits_2, Inits_3)
    
    
    
    # Parameters to track -----------------------------------------------------
    
    Pars <- c('mu_phi',
              'il_mu_phi',
              'mu_p',
              'nu_p',
              'tau_nu_p')
    
    
    # Run model ---------------------------------------------------------------
    
    
    tt <- proc.time()[3]
    OUT_ind <- jagsRun(jagsData = DATA, 
            jagsModel = 'pwatch_sim.jags',
            jagsInits = F_Inits,
            params = Pars,
            n_chain = 3,
            n_adapt = 5000,
            n_burn = 30000,
            n_draw = 20000,
            report = FALSE)
    time_ind <- round((proc.time()[3] - tt)/60, digits = 1)
    
    
    
    # ONE DETECTION -------------------------------------------------------------------
    
    {
      sink("pwatch_sim_one.jags")
      
      cat("
          
          model {
          
          #nests
          for (i in 1:NI)
          {
          #chick alive at time step 1
          z[1,i] <- 1
          
          #time step
          for (t in 2:NT)
          {
          #state model
          z[t,i] ~ dbinom(p_alive[t,i], z[t-1,i])
          p_alive[t,i] <- phi[t,i] * z[t-1,i]
          
          #observation model
          y[t,i] ~ dbinom(p_sight[t,i], z[t,i])
          p_sight[t,i] <- p[t,i] * z[t,i]
          
          }
          }
          
          #transforms
          #nests
          for (i in 1:NI)
          {
          #time
          for (t in 1:NT)
          {
          logit(phi[t,i]) <- mu_phi
          logit(p[t,i]) <- mu_p
          } #t
          } #i
          
          #priors - phi
          mu_phi ~ dnorm(0, 0.01)
          
          #priors - p
          mu_p ~ dnorm(0, 0.5)

          il_mu_phi <- ilogit(mu_phi)
          
          }",fill = TRUE)
  
      sink()
    }
    
    
    
    # Starting values ---------------------------------------------------------
    
    
    Inits_1 <- list(mu_phi = 5,
                    mu_p = 0,
                    .RNG.name = "base::Mersenne-Twister",
                    .RNG.seed = 1)
    
    Inits_2 <- list(mu_phi = 5,
                    mu_p = 0,
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = 2)
    
    Inits_3 <- list(mu_phi = 5,
                    mu_p = 0,
                    .RNG.name = "base::Marsaglia-Multicarry",
                    .RNG.seed = 3)
    
    F_Inits <- list(Inits_1, Inits_2, Inits_3)
    
    
    # Parameters to track -----------------------------------------------------
    
    Pars <- c('mu_phi',
              'il_mu_phi',
              'mu_p')
    
    
    # Run model ---------------------------------------------------------------
    
    OUT_one <- jagsRun(jagsData = DATA, 
                       jagsModel = 'pwatch_sim.jags',
                       jagsInits = F_Inits,
                       params = Pars,
                       n_chain = 3,
                       n_adapt = 5000,
                       n_burn = 30000,
                       n_draw = 10000,
                       report = FALSE)
    
    
    # Analyze -----------------------------------------------------------------
    
    mn_phi <- dphi
    mn_p <- mean(dp)
    
    #individual detection
    ind_diff <- mn_phi - inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[4])
    #one detection
    one_diff <- mn_phi - inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[4])
    
    if (mn_phi > inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[3]) &
        mn_phi < inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[5]))
    {
      IND_log <- TRUE
    } else {
      IND_log <- FALSE
    }
    
    if (mn_phi > inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[3]) &
        mn_phi < inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[5]))
    {
      ONE_log <- TRUE
    } else {
      ONE_log <- FALSE
    }
    
    temp <- data.frame(GEN = j,
                       RUN = i,
                       MN_PHI = mn_phi, 
                       MN_P = mn_p,
                       I_phi_LCI = inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[3]),
                       I_phi_MED = inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[4]),
                       I_phi_UCI = inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[5]),
                       I_diff = ind_diff,
                       I_PASS = IND_log,
                       O_phi_LCI = inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[3]),
                       O_phi_MED = inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[4]),
                       O_phi_UCI = inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[5]),
                       O_diff = one_diff, 
                       O_PASS = ONE_log,
                       I_il_lci = MCMCsummary(OUT_ind, params = 'il_mu_phi')[3],
                       I_il_med = MCMCsummary(OUT_ind, params = 'il_mu_phi')[4],
                       I_il_uci = MCMCsummary(OUT_ind, params = 'il_mu_phi')[5],
                       O_il_lci = MCMCsummary(OUT_one, params = 'il_mu_phi')[3],
                       O_il_med = MCMCsummary(OUT_one, params = 'il_mu_phi')[4],
                       O_il_uci = MCMCsummary(OUT_one, params = 'il_mu_phi')[5],
                       I_p_LCI = inv.logit(MCMCsummary(OUT_ind, params = 'mu_p')[3]),
                       I_p_MED = inv.logit(MCMCsummary(OUT_ind, params = 'mu_p')[4]),
                       I_p_UCI = inv.logit(MCMCsummary(OUT_ind, params = 'mu_p')[5]),
                       O_p_LCI = inv.logit(MCMCsummary(OUT_one, params = 'mu_p')[3]),
                       O_p_MED = inv.logit(MCMCsummary(OUT_ind, params = 'mu_p')[4]),
                       O_p_UCI = inv.logit(MCMCsummary(OUT_one, params = 'mu_p')[5]),
                       TIME_ind = time_ind)
    
    OUT <- rbind(OUT, temp)
    PB <- PB + 1
    setTxtProgressBar(pb, PB)
  }
  row.names(OUT) <- NULL
}
close(pb)

