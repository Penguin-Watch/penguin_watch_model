######################
#Mark-recapture for penguin watch 
#
#Assumption that chicks are alive at start of analysis (could have died as eggs)
#
#Authors: Casey Youngflesh
######################

#TODO
#simulate so params are actually known - see Kerry and Schaub 2012
#detection vary with nest
#detection vary with age (time)
#survival vary with age
#add posterior predictive check - see Schrimpf script


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


#simulate data
#model each site separately to start
#5 nests - 400 time steps - 2 possible chicks per nest

#ISSUES
#going to start with 0 at all nests (can't really see eggs)
#if don't see two chicks, can't say there were ever two chicks
#detection probability (and survival probability) will change over time - higher as they get older [maybe same survival probability for each nest, but they vary over time (linear function of age)? want survival prob to vary over time with env covariates maybe]
#also have possibly more than two chicks per cell in some cases when older



#simulate new data - script modified from Kerry and Schaub 2012
n_ts <- 100 #number of time steps
nests <- 5 #number of nests
surv_prob <- rep(0.99, n_ts-1)
detect_prob <- rep(0.5, n_ts-1)

SURV_PROB <- matrix(surv_prob, 
                    ncol = n_ts-1,
                    nrow = nests)

DETECT_PROB <- matrix(detect_prob,
                      ncol = n_ts-1,
                      nrow = nests)

#function to simulate time series
sim_data_fun <- function(SURV_MAT, DETECT_MAT, N_NESTS)
{
  TS_LEN <- NCOL(SURV_MAT) + 1
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
      t_SP[t] <- rbinom(1, size = t_SP[t-1], prob = SURV_MAT[i,t-1])
      
      #OBSERVED STATE
      t_DP <- rbinom(1, size = t_SP[t], prob = DETECT_MAT[i,t-1])
      CH[i,t] <- t_DP
    }
  }
  return(CH)
}

sim_data <- sim_data_fun(SURV_PROB, DETECT_PROB, nests)


# Data for model ----------------------------------------------------------


DATA <- list(
  y = sim_data, #reponse
  N = NROW(sim_data), #number of nests
  L = NCOL(sim_data)) #number of time points



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
                              surv_p[i,t] * z[i,t-1],
                              surv_p[i,t])

        #observation model
        y[i,t] ~ dbinom(p_sight[i,t], z[i,t])
        p_sight[i,t] <- ifelse(z[i,t] < 2,
                              detect_p[i,t] * z[i,t],
                              detect_p[i,t])

      }
    }


    #priors
    for (i in 1:N)
    {
      for (t in 1:L)
      {
        detect_p[i,t] <- mn_detect_p
        surv_p[i,t] <- mn_surv_p
      }
    }

    mn_detect_p ~ dunif(0,1)
    mn_surv_p ~ dunif(0,1)


    }",fill = TRUE)

sink()
}



# Starting values ---------------------------------------------------------

#produce inits for z-state - fill 2s and 1s
#fun modified from Kerry and Schaub 2012
#assume that chicks are alive since time step 1 (even though they could have died as eggs)
known.state.fun <- function(INPUT)
{
  #INPUT <- DATA$y
  state <- INPUT
  for (i in 1:NROW(INPUT))
  {
    n1 <- 1
    n2 <- max(which(INPUT[i,] == 2))
    
    state[i,n1:n2] <- 2
    
    n3 <- max(which(state[i,] == 1))

    state[i,(n2+1):n3] <- 1
    state[i,n1] <- NA
  }
  state[state == 0] <- NA
  return(state)
}

z_vals <- known.state.fun(DATA$y)

Inits_1 <- list(mn_surv_p = runif(1, 0, 1),
                mn_detect_p = runif(1, 0, 1),
                z = z_vals,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mn_surv_p = runif(1, 0, 1),
                mn_detect_p = runif(1, 0, 1),
                z = z_vals,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mn_surv_p = runif(1, 0, 1),
                mn_detect_p = runif(1, 0, 1),
                z = z_vals,
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)



# Parameters to track -----------------------------------------------------

Pars <- c('mn_surv_p',
          'mn_detect_p')


# Inputs for MCMC ---------------------------------------------------------

n_adapt <- 5000  # number for initial adapt
n_burn <- 2000 # number burnin
n_draw <- 2000  # number of final draws to make
n_thin <- 2    # thinning rate
n_chain <- 3  # number of chains

Rhat_max <- 1.1 # max allowable Rhat (close to 1 = convergence)
n_max <- 1e6 # max allowable iterations


# Run model (non-parallel) ---------------------------------------------------------------


#rjags
jm = jags.model(data = DATA, 
                file = "mark_recapture.jags", 
                inits = F_Inits, 
                n.chains = 3, 
                n.adapt = n_adapt)

update(jm, n.iter = n_burn)

out = coda.samples(jm, 
                   n.iter = n_draw, 
                   variable.names = Pars, 
                   thin = n_thin)



#extra draws if didn't converge
n_total <- n_burn + n_draw
n_extra <- 0
while(max(MCMCsummary(out)[,5]) > Rhat_max &
      n_total < n_max)
{
  
  out <- update(out,
                n.iter = n_draw,
                n.chains = n_chain,
                n.thin = n_thin)
  
  n_extra <- n_extra + n_draw
  n_total <- n_total + n_draw
}

n_final <- n_draw/n_thin

#Inferences were derived from $`r n_final`$ samples drawn following an adaptation period of $`r n_adapt`$ draws, and a burn-in period of $`r (n_total - n_draw)`$ draws using $`r n_chain`$ chains and a thinning rate of $`r n_thin`$.




# Analyze posterior -------------------------------------------------------


#summary
MCMCsummary(out)

#trace plots
MCMCtrace(out, params = 'beta', ind = TRUE)

#plots of beta parameters
MCMCplot(out, params = 'beta', rank = FALSE, labels = NULL,
         horiz = FALSE, ref_ovl = FALSE)


#for PPC - how does simulated error compare to actual error
s.sq.act <- MCMCchains(out, params = 'sq.act')
s.sq.sim <- MCMCchains(out, params = 'sq.sim')

plot(s.sq.act, s.sq.sim, pch = '.')
abline(a=0, b=1, lty=2, lwd = 5)
mean(s.sq.act > s.sq.sim) #should be about 0.5


#for r^2 and lamba (degree of pooling) - see Gelman and Pardoe 2006
e.y_out <- MCMCchains(out, params = 'e.y')
#id.y_out <- MCMCchains(out, params = 'id.y') - to calculate separate r^2 values

#r^2
r2.y <- 1 - mean(apply(e.y_out, 1, var))/var(data_in$y, na.rm=TRUE)

#lambda
lambda.y <- 1 - var(apply(e.y_out, 2, mean))/mean(apply(e.y_out, 1, var))
