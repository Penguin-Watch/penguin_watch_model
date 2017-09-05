
#FIXED EFFECT FOR EPS_PHI - though should use random effect if we want to estimate temporal variabilty in EPS_PHI
#RANDOM EFFECT FOR EPS_P

#Probably want to remove EPS_PHI as a time step effect once we have covariates. For one site, something like logit(phi[i,j]) <- mu_phi + beta1*TIME[t] + beta2*SIC[t] + beta3*KRILL[t] + error[t]
#beta1 contrained time effect
#beta2 effect of SIC
#beta3 effect of KRILL
#error is variation not accounted for by these other covariates

#should also probably model (FOR PHI PARAMS) those beta parameters with site as a random effect (using group index) - mu should probably have a random effect as well - there may be simply be differences at sites that we aren't accounting for in our model

#probably don't need a random effect on any of the P PARAMS. Everything to do with detection probability should probably be kept the same. All nests can vary, so they are essentially treated as coming from one big colony for estimation of P.

#should just simulate two covariates that impact survival on a fine time scale (i.e., every time step) and remove EPS_PHI as time effect - do need to include EPS_PHI[t] as an error term though (or if modeling colony as a random effect, being sure to give each colony it's own EPS_PHI[t])


#PHI depends on time, but P depends on individual identity, so I think it's fine
#Look at correaltion between all PHI and P params, not just among PHI params and among P params
#Look at overlap betwen posterior and prior - put uniform (or broad) priors on betas and mus to look at this. Don't think this can be used for the error term, as this is modeled hierarchically - those parameters are themselves drawn from a distirbution.

#starting values of beta_p must be 0 when many time steps bc inv.logit(mu + beta_p*[large index for cov]) = 1 and if y = 0 at that index, JAGS doesn't like it

rm(list = ls())
require(rjags)

#CHANGE SIZE OF RESPONSE DATA HERE
n_ts <- 100 #number of time steps
x <- 1:n_ts
nests <- 30 #number of nests

#simulation function
sim_p_fun <- function(START, RATE = 0.01, TOP = 1)
{
  B <- rep(NA, n_ts-1)
  
  B[1] <- START
  K <- TOP
  r <- RATE
  for(t in 1:(n_ts-2))
  {
    B[t+1] <- B[t] + r*(1-B[t]/K)
  }
  return(B)
}

phi_data <- sim_p_fun(START = 0.99)

#simulate phi
PHI <- matrix(rep(phi_data, nests),
              nrow = nests,
              ncol = n_ts-1,
              byrow = TRUE)

#simulate P
P <- matrix(rep(NA, nests*(n_ts-1)),
            nrow = nests,
            ncol = n_ts-1)

set.seed(1)
dp <- runif(30, 0.3, 0.8)

for (i in 1:length(dp))
{
  P[i,] <- sim_p_fun(dp[i], RATE = 0.005, TOP = 0.95)
}

sim_data_fun <- function(PHI_MAT, P_MAT, N_NESTS)
{  
  TS_LEN <- NCOL(PHI_MAT) + 1
  CH <- matrix(0,
               ncol = TS_LEN,
               nrow = N_NESTS)
  
  for (i in 1:N_NESTS)
  {
    #i <- 1
    #both alive at start
    CH[i,1] <- 2
    
    t_SP <- c(2, rep(NA, TS_LEN-1))
    for (t in 2:TS_LEN)
    {
      #TRUE STATE
      t_SP[t] <- rbinom(1, size = t_SP[t-1], prob = PHI_MAT[i,t-1])
      
      #OBSERVED STATE
      t_DP <- rbinom(1, size = t_SP[t], prob = P_MAT[i,t-1])
      CH[i,t] <- t_DP
    }
  }
  return(CH)
}

sim_data <- sim_data_fun(PHI, P, nests)

#known info regarding z-state
known.state.fun <- function(INPUT)
{
  state <- INPUT
  for (i in 1:NROW(INPUT))
  {
    n1 <- 1
    
    if (sum(state[i,] == 2) > 0)
    {
      n2 <- max(which(INPUT[i,] == 2))
      state[i,n1:n2] <- 2
    }
    
    #NA at first state because model designates 2 alive at time step 1
    state[i,n1] <- NA
  }
  state[state == 0] <- NA
  state[state == 1] <- NA
  return(state)
}

z_vals <- known.state.fun(sim_data)



# Data for model ----------------------------------------------------------

DATA <- list(
  y = sim_data, #reponse
  N = NROW(sim_data), #number of nests
  L = NCOL(sim_data), #number of time points
  z = z_vals,
  x = 1:NCOL(sim_data))

# Model -------------------------------------------------------------------

#L = length of time series
#N = number of nests

{
  sink("mark_recapture.jags")
  
  cat("
      model {
      
      for (i in 1:N)
      {
      #both alive at time step 1
      z[i,1] <- 2
      
      for (t in 2:L)
      {
      
      #state model
      z[i,t] ~ dbinom(p_alive[i,t], z[i,t-1])
      p_alive[i,t] <- ifelse(z[i,t-1] < 2,
      phi[i,t] * z[i,t-1],
      phi[i,t])
      
      #observation model
      y[i,t] ~ dbinom(p_sight[i,t], z[i,t])
      p_sight[i,t] <- ifelse(z[i,t] < 2,
      p[i,t] * z[i,t],
      p[i,t])
      }
      }
      
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
      eps_phi[t] ~ dnorm(0, tau_phi)
      }
      
      mean_phi ~ dbeta(1.5,1.5)                 #Mean survival
      mu_phi <- log(mean_phi / (1 - mean_phi))
      tau_phi <- pow(sigma_phi, -2)
      sigma_phi ~ dunif(0, 10)
      
      for (i in 1:N)
      {
      eps_p[i] ~ dnorm(0, tau_p)
      }
      
      mean_p ~ dbeta(1.5,1.5)                    #Mean detection
      mu_p <- log(mean_p / (1 - mean_p))
      tau_p <- pow(sigma_p, -2)
      sigma_p ~ dunif(0, 10)
      
      beta_phi ~ dnorm(0, 100)
      beta_p ~ dnorm(0, 10)
      
      }",fill = TRUE)

  sink()
}

# Starting values ---------------------------------------------------------

Inits <- list(mean_phi = 0.5,
              mean_p = 0.5,
              sigma_phi = 0.5,
              sigma_p = 0.5,
              beta_phi = 0.1,
              beta_p = 0,
              .RNG.name = "base::Mersenne-Twister",
              .RNG.seed = 1)

# Parameters to track -----------------------------------------------------

Pars <- c('mean_phi',
          'mean_p',
          'sigma_p',
          'sigma_phi')

# Inputs for MCMC ---------------------------------------------------------

JAGS_FILE <- 'mark_recapture.jags'
n_adapt <- 8  # number for initial adapt
n_burn <- 10  # number burnin
n_draw <- 20  # number of final draws to make
n_thin <- 2   # thinning rate
n_chain <- 3  # number of chains

jm = jags.model(data = DATA,
                file = paste0(JAGS_FILE),
                inits = Inits,
                n.chains = 3,
                n.adapt = n_adapt)