######################
#Mark-recapture for penguin watch 
#
#
#Authors: Casey Youngflesh
######################



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


#1a - survive entire time
#1b - survive entire time
#2a - die day 100
#2b - die day 250
#3a - survive entire time
#3b - die day 300
#4a - die day 320
#4b - die day 150
#5a - survive entire time
#5b - die day 350


#true state
true_state <- data.frame(nest1 = rep(2, 400), 
                         nest2 = c(rep(2, 100), rep(1, 150), rep(0, 150)),
                         nest3 = c(rep(2, 300), rep(1, 100)),
                         nest4 = c(rep(2, 150), rep(1, 170), rep(0, 80)),
                         nest5 = c(rep(2, 350), rep(1, 50)))


#observed state - simulate imperfect detection

obs_fun <- function(input)
{
  #observation probability
  obs_prob <- 0.5
  o_state <- c()
  for (i in 1:length(input))
  {
    o_state[i] <- rbinom(1, size = input[i], prob = obs_prob)
  }
  return(o_state)
}

obs_state <- apply(true_state, 2, obs_fun)



DATA <- list(
  y = obs_state, #reponse
  N = NCOL(obs_state), #number of nests
  L = NROW(obs_state)) #number of time points



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
        #observation model
        
        yt1[t,i] ~ dbern(p_sight[t,i])
        yt2[t,i] ~ dbern(p_sight[t,i])
        
        y[t,i] <- ifelse(z[t,i]<2,
                        yt1[t,i],
                        yt1[t,i] + yt2[t,i])

        p_sight[t,i] <- ifelse(z[t,i]<2,
                        z[t,i] * detect_p[i],
                        detect_p[i])



        #state model

        zt1[t,i] ~ dbern(p_alive[t,i])
        zt2[t,i] ~ dbern(p_alive[t,i])        

        z[t,i] <- ifelse(z[t-1,i]<2,
                        zt1[t,i],
                        zt1[t,i] + zt2[t,i])

        p_alive[t,i] <- ifelse(z[t-1,i]<2,
                        z[t-1,i] * surv_p,
                        surv_p)
      }
    }


    #priors
    surv_p ~ dunif(0,1)

    for (i in 1:N)
    {
      detect_p[i] ~ dunif(0,1)
    }


#NOT SURE ABOUT ANY OF THIS IN THIS CASE
    #For posterior predictive check
    #simulate data and calc sq residual
    #y.rep[t,i] ~ dbern(p_sight.rep[t,i])
    #p_sight.rep[t,i] <- z.rep[t,i] * detect_p[i]
    #z.rep[t,i] ~ dbern(p_alive[t,i])

    #sim compared to fit      
    #sq.sim[i] <- pow(y.rep[t,i] - mu.y[i], 2)
    #actual compared to fit      
    #sq.act[i] <- pow((y[t,i] - mu.y[i]), 2)
    
    #error for R^2 calculation
    #e.y[i] <- y[i] - mu.y[i]
    #ids for each data point
    #id.y[i] <- id[i]
    
    }",fill = TRUE)

sink()
}



# Starting values ---------------------------------------------------------


Inits_1 <- list(surv_p = runif(1, min = 0, max = 1),
                detect_p = runif(DATA$N, min = 0, max = 1),
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(surv_p = runif(1, min = 0, max = 1),
                detect_p = runif(DATA$N, min = 0, max = 1),
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(surv_p = runif(1, min = 0, max = 1),
                detect_p = runif(DATA$N, min = 0, max = 1),
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)



# Parameters to track -----------------------------------------------------

Pars <- c('surv_p',
          'detect_p',
          'p_alive',
          'p_sight',
          'z')


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
                inits = Inits, 
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
