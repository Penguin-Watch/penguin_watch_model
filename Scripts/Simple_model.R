##############################
# SIMPLE MODEL
#
# Penguin Watch Data simulated (no trend across time) - detection probability different for each nest (no day/night)
# Model p constant across time - random effect of nest
# Model phi constant across time and nest
#
# Justification for priors for parameters that contirbute to logit(p) and logit(phi)
#
# PPC (sd) results indicate mark-recapture with nest as random effect doesn't predict sd well - probably not an appropriate statistic to use in larger model
# PPC (mean) results suggest model is appropriate
# Kery and Schaub 2012 p 223 states that modeling using this format, you can't conduct the typical GOF stat, because we have individual data (they are doing something similar to a chi-square test for each iteration on both real and simualted data)
##############################

# Clear environment -------------------------------------------------------

rm(list = ls())


# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman', repos = "http://cran.case.edu")
}
pacman::p_load(rjags, parallel, MCMCvis)


#JAGS module
load.module("glm")


# Load data ---------------------------------------------------------------


#CHANGE SIZE OF RESPONSE DATA HERE
n_ts <- 100 #number of time steps
x <- 1:n_ts
nests <- 30 #number of nests


#simulate phi
PHI <- matrix(rep(0.995, nests*(n_ts-1)),
              nrow = nests,
              ncol = n_ts-1,
              byrow = TRUE)

#simulate P
set.seed(1)
dp <- runif(30, 0.3, 0.8)

P <- matrix(rep(NA, nests*(n_ts-1)),
            nrow = nests,
            ncol = n_ts-1)

for (i in 1:length(dp))
{
  P[i,] <- rep(dp[i], NCOL(P))
}



#simulate data
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
      y.new[i,1] ~ dbinom(p[i,1], 2)
      
      
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
      
      #PPC - simulate data
      y.new[i,t] ~ dbinom(p_sight[i,t], z[i,t])
      
      #calculate likelihood for each data point
      t.act.like[i,t] <- dbin(y[i,t], p_sight[i,t], z[i,t])
      t.sim.like[i,t] <- dbin(y.new[i,t], p_sight[i,t], z[i,t])
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
      
      #likelihood - King and Brook 2002 - p. 804 (not sure if this is acutally the full likelihood though)
      for (i in 1:N)
      {
      t2.act.like[i] <- prod(t.act.like)
      t2.sim.like[i] <- prod(t.sim.like)
      }

      #act.like <- prod(t2.act.like)
      #sim.like <- prod(t2.sim.like)

      
      #transforms
      for (i in 1:N)
      {
      for (t in 1:L)
      {
      logit(phi[i,t]) <- mu_phi                   #phi = survival prob
      logit(p[i,t]) <- mu_p + eps_p[i]            #p = detection prob
      }
      }
      
      mean_phi ~ dbeta(1,1)                 #Mean survival
      mu_phi <- logit(mean_phi)             #log(mean_phi / (1 - mean_phi))
      
      for (i in 1:N)
      {
      #want to shoot for sd of 1.7 for vague prior on p (i.e., constrain sigma_p)     
      eps_p[i] ~ dnorm(0, tau_p) T(-10, 10) 

      #n_p is detection probability for each nest
      n_p[i] <- ilogit(mu_p + eps_p[i])
      }
      
      mean_p ~ dbeta(1.5,1.5)                    #Mean detection
      mu_p <- logit(mean_p)                  #(mean_p / (1 - mean_p))
      tau_p <- pow(sigma_p, -2)
      sigma_p ~ dunif(0.25, 3)

      
      }",fill = TRUE)

  sink()
}




# Starting values ---------------------------------------------------------

Inits_1 <- list(mean_phi = 0.9,
                mean_p = 0.5,
                #sigma_phi = 0.1,
                sigma_p = 1.1,
                #beta_phi = 0.1,
                #beta_p = 0,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mean_phi = 0.9,
                mean_p = 0.6,
                #sigma_phi = 0.11,
                sigma_p = 1.5,
                #beta_phi = 0.1,
                #beta_p = 0,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mean_phi = 0.9,
                mean_p = 0.4,
                #sigma_phi = 0.09,
                sigma_p = 1.75,
                #beta_phi = 0.1,
                #beta_p = 0,
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)




# Parameters to track -----------------------------------------------------

Pars <- c('mean_phi',
          'mean_p',
          'sigma_p',
          'n_p',
          'sd.y',
          'sd.y.new',
          'mn.y',
          'mn.y.new',
          'act.like',
          'sim.like')




# Inputs for MCMC ---------------------------------------------------------

NAME <- 'out_Sep_05_2017_R1_simplified'

JAGS_FILE <- 'mark_recapture.jags'
n_adapt <- 5000  # number for initial adapt
n_burn <- 5000 # number burnin
n_draw <- 4000  # number of final draws to make
n_thin <- 2    # thinning rate
n_chain <- 3  # number of chains

EXTRA <- FALSE
Rhat_max <- 1.02 # max allowable Rhat (close to 1 = convergence)
n_max <- 100000 # max allowable iterations


# DEBUG -------------------------------------------------------------------

# jm = jags.model(data = DATA,
#                 file = paste0(JAGS_FILE),
#                 inits = F_Inits,
#                 n.chains = 3,
#                 n.adapt = n_adapt)
# 
# update(jm,n.iter = n_burn)
# 
# samples = coda.samples(jm,
#                        n.iter = n_draw,
#                        variable.names = Pars,
#                        thin = n_thin)
# 
# MCMCsummary(samples, ISB = FALSE)



# Run model (parallel) ---------------------------------------------------------------

#number of chains
cl <- parallel::makeCluster(n_chain)

pid <- NA
for(i in 1:n_chain)
{
  pidNum <- capture.output(cl[[i]])
  start <- regexpr("pid", pidNum)[[1]]
  end <- nchar(pidNum)
  pid[i] <- substr(pidNum, (start + 4), end)
}

parallel::clusterExport(cl,
                        c('DATA',
                          'n_adapt',
                          'n_burn',
                          'n_draw',
                          'n_thin',
                          'Pars',
                          'pid',
                          'F_Inits',
                          'JAGS_FILE'
                        ))


ptm <- proc.time()
out.1 <- parallel::clusterEvalQ(cl,
                                {
                                  require(rjags)
                                  processNum <- which(pid==Sys.getpid())
                                  m.inits <- F_Inits[[processNum]]
                                  
                                  jm = jags.model(data = DATA,
                                                  file = paste0(JAGS_FILE),
                                                  inits = m.inits,
                                                  n.chains = 1,
                                                  n.adapt = n_adapt)
                                  
                                  update(jm,
                                         n.iter = n_burn)
                                  
                                  samples = coda.samples(jm,
                                                         n.iter = n_draw,
                                                         variable.names = Pars,
                                                         thin = n_thin)
                                  return(samples)
                                })


out <- coda::mcmc.list(out.1[[1]][[1]],
                       out.1[[2]][[1]],
                       out.1[[3]][[1]])
stopCluster(cl)
tt <- (proc.time() - ptm)[3]/60 #minutes





#justification of priors on tau_p
require(boot)
#lower bounds
a <- rnorm(1000, 0, 1)
b <- inv.logit(a)
hist(b)
#upper bounds
a <- rnorm(1000, 2, 2)
b <- inv.logit(a)
hist(b)

#justification of priors on beta_p and beta_phi
#a <- rnorm(1000, 0, )




#summarize
MCMCsummary(out, digits = 4)
MCMCtrace(out)


#PPC
sd.y.ch <- MCMCchains(out, params = 'sd.y')
sd.y.new.ch <- MCMCchains(out, params = 'sd.y.new')

plot(sd.y.ch, sd.y.new.ch, pch = '.', asp = 1)
abline(0,1, col = 'red')
sum(sd.y.new.ch > sd.y.ch)/length(sd.y.ch)


mn.y.ch <- MCMCchains(out, params = 'mn.y')
mn.y.new.ch <- MCMCchains(out, params = 'mn.y.new')

plot(mn.y.ch, mn.y.new.ch, pch = '.', asp = 1)
abline(0,1, col = 'red')
sum(mn.y.ch > mn.y.new.ch)/length(mn.y.ch)
