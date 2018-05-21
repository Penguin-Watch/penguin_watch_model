######################
#Mark-recapture for penguin watch 
#
#Model detection as a group and for individual nests using simulated data
#
#Authors: Casey Youngflesh
######################

#each nest own survival
#each nest own detection

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
n_ts <- 300 #number of time steps
x <- 1:n_ts
nests <- 15 #number of nests


OUT <- c()
for (i in 1:10)
{
  # survival ----------------------------------------------------------------
  
  #survival probability - change over time
  # sim_p_fun <- function(START, RATE = 0.008, TOP = 1)
  # {
  #   B <- rep(NA, n_ts-1)
  #   
  #   B[1] <- START
  #   K <- TOP
  #   r <- RATE
  #   for(t in 1:(n_ts-2))
  #   {
  #     B[t+1] <- B[t] + r*(1-B[t]/K)
  #   }
  #   return(B)
  # }
  # 
  # phi_data <- sim_p_fun(0.985)
  
  
  
  #i <- 1
  
  set.seed(i)
  dphi <- runif(nests, 0.990, 0.999)
  
  PHI <- matrix(rep(dphi, each = n_ts -1 ),
                ncol = nests)
  
  
  
  
  # detection ---------------------------------------------------------------
  
  #starting probs for each nest
  set.seed(i + 1)
  dp <- runif(nests, 0.3, 0.8)
  
  P <- matrix(rep(dp, each = n_ts - 1),
              ncol = nests)
  
  # for (i in 1:length(dp))
  # {
  #   P[i,] <- sim_p_fun(dp[i], RATE = 0.005, TOP = 0.95)
  #   #plot(p_data, type = 'l', ylim = c(0,1), main = paste0(i))
  # }
  
  #check trend
  #lm
  # cov_p <- 1:length(P[1,])
  # fit_p <- lm(P[1,] ~ cov_p)
  # summary(fit_p)
  
  
  
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
  
  set.seed(i + 2)
  sim_data <- sim_data_fun(PHI, P, nests, TYPE = 'ONE')
  
  
  
  
  # Create z-state matrix ---------------------------------------------------
  
  #fun modified from Kerry and Schaub 2012
  #assume that chicks are alive since time step 1 (even though they could have died as eggs)
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
        #phi = survival prob
        #mu_phi = grand mean for all sites/years
        
        logit(phi[t,i]) <- mu_phi
        
        
        #p = detection prob
        #mu_p = grand mean for all sites/years
        #beta_phi = slope for increasing detection over time (older chicks have higher detection p)
        #nu_p = effect of nest
        
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
        
        }",fill = TRUE)

    sink()
  }
  
  
  
  # Starting values ---------------------------------------------------------
  
  
  Inits_1 <- list(mu_phi = 7,
                  mu_p = -2,
                  tau_nu_p = 1,
                  .RNG.name = "base::Mersenne-Twister",
                  .RNG.seed = 1)
  
  Inits_2 <- list(mu_phi = 7,
                  mu_p = -2,
                  tau_nu_p = 1,
                  .RNG.name = "base::Wichmann-Hill",
                  .RNG.seed = 2)
  
  Inits_3 <- list(mu_phi = 7,
                  mu_p = -2,
                  tau_nu_p = 1,
                  .RNG.name = "base::Marsaglia-Multicarry",
                  .RNG.seed = 3)
  
  F_Inits <- list(Inits_1, Inits_2, Inits_3)
  
  
  
  # Parameters to track -----------------------------------------------------
  
  Pars <- c('mu_phi',
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
          n_burn = 300000,
          n_draw = 50000,
          report = FALSE)
  time_ind <- round(proc.time()[3] - tt, digits = 1)
  
  
  
  
  
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
        #phi = survival prob
        #mu_phi = grand mean for all sites/years
        
        logit(phi[t,i]) <- mu_phi
        
        
        #p = detection prob
        #mu_p = grand mean for all sites/years
        #beta_phi = slope for increasing detection over time (older chicks have higher detection p)
        #nu_p = effect of nest
        
        logit(p[t,i]) <- mu_p
        
        } #t
        } #i
        
        #priors - phi
        mu_phi ~ dnorm(0, 0.01)
        
        #priors - p
        mu_p ~ dnorm(0, 0.5)
        
        }",fill = TRUE)

    sink()
  }
  
  
  
  # Starting values ---------------------------------------------------------
  
  
  Inits_1 <- list(mu_phi = 7,
                  mu_p = -2,
                  .RNG.name = "base::Mersenne-Twister",
                  .RNG.seed = 1)
  
  Inits_2 <- list(mu_phi = 7,
                  mu_p = -2,
                  .RNG.name = "base::Wichmann-Hill",
                  .RNG.seed = 2)
  
  Inits_3 <- list(mu_phi = 7,
                  mu_p = -2,
                  .RNG.name = "base::Marsaglia-Multicarry",
                  .RNG.seed = 3)
  
  F_Inits <- list(Inits_1, Inits_2, Inits_3)
  
  
  
  # Parameters to track -----------------------------------------------------
  
  Pars <- c('mu_phi',
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
  
  mn_phi <- mean(dphi)
  
  #individual detection
  ind_diff <- mn_phi - inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[4])
  #one detection
  one_diff <- mn_phi - inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[4])
  
  temp <- data.frame(IND = ind_diff, ONE = one_diff, TIME_ind = time_ind)
  
  OUT <- rbind(OUT, temp)
}
