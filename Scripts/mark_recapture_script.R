######################
#Mark-recapture for penguin watch 
#
#Assumption that chicks are alive at start of analysis (though could have failed as eggs)
#
#Authors: Casey Youngflesh
######################

#TODO
#add posterior predictive check - see Schrimpf script
#how to model more than one site?
#might want to know when the chicks die - how do we do that? should just fit surv param as random effect and look at that change across time rather than linear model?
#do the trends in both detection and survival make for nonidentifiability?

# Clear environment -------------------------------------------------------

rm(list = ls())


# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman')
}

pacman::p_load(rjags, MCMCvis)



# Load data ---------------------------------------------------------------


#setwd('Data')
#data <- read.csv('XXXX.csv', header=TRUE)

#ISSUES
#if don't see two chicks, can't say there were ever two chicks
#detection probability (and survival probability) will change over time - higher as they get older [maybe same survival probability for each nest, but they vary over time (linear function of age)? want survival prob to vary over time with env covariates maybe]
#also have possibly more than two chicks per cell in some cases when older
#what to do about no observation over night - can't be ignored right? Treating each hour as one time step here
#might not have to assume both chicks alive at start - first sighting of 2 chicks can be start - p182 Kerry and Schaub 2012
#all years should be run hierarchically for a site, and all sites hierarchically for each species?


#simulate new data - script modified from Kerry and Schaub 2012
n_ts <- 100 #number of time steps
x <- 1:n_ts
nests <- 5 #number of nests

surv_prob <- rep(0.985, n_ts-1)


#survival probability
a_s <- 0.985
b_s <- 0.0001

PHI <- rbind(a_s + b_s*x,
             a_s + b_s*x,
             a_s + b_s*x,
             a_s + b_s*x,
             a_s + b_s*x)


#detection probability
n1 <- 0.3
n2 <- 0.4
n3 <- 0.5
n4 <- 0.6
n5 <- 0.7
b <- 0.002

P <- rbind(n1 + b*x,
           n2 + b*x,
           n3 + b*x,
           n4 + b*x,
           n5 + b*x)


#function to simulate time series
sim_data_fun <- function(PHI_MAT, P_MAT, N_NESTS)
{
  TS_LEN <- NCOL(PHI_MAT) + 1
  CH <- matrix(0, 
               ncol = TS_LEN, 
               nrow = N_NESTS)
  
  for (i in 1:N_NESTS)
  {
    #both chicks alive at start
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

#fun modified from Kerry and Schaub 2012
#assume that chicks are alive since time step 1 (even though they could have died as eggs)
known.state.fun <- function(INPUT)
{
  #INPUT <- DATA$y
  state <- INPUT
  for (i in 1:NROW(INPUT))
  {
    n1 <- 1
    
    if (sum(state[i,] == 2) > 0)
    {
      n2 <- max(which(INPUT[i,] == 2))
      state[i,n1:n2] <- 2
    }
    
    if (sum(state[i,] == 1) > 0)
    {
      n3 <- max(which(state[i,] == 1))
      state[i,(n2+1):n3] <- 1
    }
    
    state[i,n1] <- NA
  }
  state[state == 0] <- NA
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

#L = length of time series (400)
#N = number of nests (5)

{
sink("mark_recapture.jags")

cat("
    model {
    
    for (i in 1:N)
    {
      #both chicks alive at time step 1
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
        logit(phi[i,t]) <- mu_phi + beta_phi*x[t]
        logit(p[i,t]) <- mu_p + beta_p*x[t] + eps_p[i]
      }
    }


    #priors
    for (i in 1:N)
    {
      eps_p[i] ~ dnorm(0, tau_p)
    }

    #phi = survival prob
    #p = detection prob
    
    mean_phi ~ dunif(0,1)
    mu_phi <- log(mean_phi / (1 - mean_phi))    
    
    mean_p ~ dunif(0,1)                 #Mean survival - could use alternative below
    mu_p <- log(mean_p / (1 - mean_p))  #Logit transform - could use alternative below
    tau_p <- pow(sigma, -2)
    sigma ~ dunif(0, 10)
    sigma2 <- pow(sigma, 2)
    
    beta_phi ~ dnorm(0, 0.1)
    beta_p ~ dnorm(0, 0.1)

    #mu_p ~ dnorm(0, 0.001)           #Prior for logit of mean survival
    #mean_p <- 1 / (1+exp(-mu_p))     #Inv-logit transform


    }",fill = TRUE)

sink()
}



# Starting values ---------------------------------------------------------


Inits_1 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma = runif(1, 0, 10),
                beta_p = 0,
                beta_phi = 0,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma = runif(1, 0, 10),
                beta_p = 0,
                beta_phi = 0,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma = runif(1, 0, 10),
                beta_p = 0,
                beta_phi = 0,
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)



# Parameters to track -----------------------------------------------------

Pars <- c('mean_phi',
          'mean_p',
          'sigma2',
          'beta_p',
          'beta_phi')


# Inputs for MCMC ---------------------------------------------------------

n_adapt <- 5000  # number for initial adapt
n_burn <- 6000 # number burnin
n_draw <- 20000  # number of final draws to make
n_thin <- 2    # thinning rate
n_chain <- 3  # number of chains

Rhat_max <- 1.02 # max allowable Rhat (close to 1 = convergence)
n_max <- 1e6 # max allowable iterations


# Run model (non-parallel) ---------------------------------------------------------------


#rjags
jm = jags.model(data = DATA, 
                file = "mark_recapture.jags", 
                inits = F_Inits, 
                n.chains = 3, 
                n.adapt = n_adapt)

update(jm, n.iter = n_burn)

out <- coda.samples(jm, 
                   n.iter = n_draw, 
                   variable.names = Pars, 
                   thin = n_thin)



#extra draws if didn't converge
n_total <- n_burn + n_draw
n_extra <- 0
while(max(MCMCsummary(out)[,5], na.rm = TRUE) > Rhat_max &
      n_total < n_max)
{
  
  out <- coda.samples(jm, 
                      n.iter = n_draw, 
                      variable.names = Pars,
                      n.thin = n_thin)
  
  n_extra <- n_extra + n_draw
  n_total <- n_total + n_draw
}

n_final <- floor((n_draw + n_extra)/n_thin)

#Inferences were derived from $`r n_final`$ samples drawn following an adaptation period of $`r n_adapt`$ draws, and a burn-in period of $`r (n_total - n_draw)`$ draws using $`r n_chain`$ chains and a thinning rate of $`r n_thin`$.




# Analyze posterior -------------------------------------------------------


#phi = survival prob
#p = detection prob

#summary
MCMCsummary(out, digits = 4)

MCMCtrace(out, pdf = TRUE)

p_ch <- MCMCchains(out, 'p')
p_v <- apply(p_ch, 2, median)[-c(1:3)]



#trace plots
MCMCtrace(out, params = 'beta', ind = TRUE)

#plots of beta parameters
MCMCplot(out, params = 'beta', rank = FALSE, labels = NULL,
         horiz = FALSE, ref_ovl = FALSE)

