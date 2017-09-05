
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
#Look at overlap betwen posterior and prior - put uniform (or broad) priors on betas and mus to look at this. Don't think this can be used for the error term, as this is modeled hierarchically - those parameters are themselves drawn from a distribution.

#starting values of beta_p must be 0 when many time steps bc inv.logit(mu + beta_p*[large index for cov]) = 1 and if y = 0 at that index, JAGS doesn't like it

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

      y.new[i,t] ~ dbinom(p_sight[i,t], z[i,t])
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
      logit(phi[i,t]) <- mu_phi# + beta_phi*x[t] + eps_phi[t]       #phi = survival prob
      logit(p[i,t]) <- mu_p# + beta_p*x[t] + eps_p[i]            #p = detection prob
      }
      }
      

      #priors
      # for (t in 1:L)
      # {
      # eps_phi[t] ~ dnorm(0, tau_phi) T(-10, 10)
      # }
      
      mean_phi ~ dbeta(1,1)                 #Mean survival
      mu_phi <- log(mean_phi / (1 - mean_phi))
      # tau_phi <- pow(sigma_phi, -2)
      # sigma_phi ~ dunif(0, 10)
      
      # for (i in 1:N)
      # {
      # eps_p[i] ~ dnorm(0, tau_p) T(-10, 10)
      # }
      
      mean_p ~ dbeta(1,1)                    #Mean detection
      mu_p <- log(mean_p / (1 - mean_p))
      # tau_p <- pow(sigma_p, -2)
      # sigma_p ~ dunif(0, 10)
      
      # beta_phi ~ dnorm(0, 1000) T(0,1)
      # beta_p ~ dnorm(0, 100) T(0,1)
      
      }",fill = TRUE)

  sink()
}




# Starting values ---------------------------------------------------------

Inits_1 <- list(mean_phi = 0.9,
                mean_p = 0.5,
                #sigma_phi = 0.1,
                #sigma_p = 0.1,
                #beta_phi = 0.1,
                #beta_p = 0,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mean_phi = 0.9,
                mean_p = 0.6,
                #sigma_phi = 0.11,
                #sigma_p = 0.11,
                #beta_phi = 0.1,
                #beta_p = 0,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mean_phi = 0.9,
                mean_p = 0.4,
                #sigma_phi = 0.09,
                #sigma_p = 0.09,
                #beta_phi = 0.1,
                #beta_p = 0,
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)




# Parameters to track -----------------------------------------------------

Pars <- c('mean_phi',
          'mean_p',
          'sd.y',
          'sd.y.new',
          'mn.y',
          'mn.y.new')




# Inputs for MCMC ---------------------------------------------------------

NAME <- 'out_Sep_05_2017_R1_simplified'

JAGS_FILE <- 'mark_recapture.jags'
n_adapt <- 5000  # number for initial adapt
n_burn <- 1000 # number burnin
n_draw <- 2000  # number of final draws to make
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






#calculate rhats
var_names <- vapply(strsplit(colnames(out[[1]]), 
                             split = "[", fixed = TRUE), `[`, 1, FUN.VALUE=character(1))

params = c('mean_phi',
           'mean_p',
           'sigma_p',
           'sigma_phi',
           'beta_p',
           'beta_phi',
           'mu_phi',
           'mu_p')

grouped <- c()
for (i in 1:length(params))
{
  get.rows <- which(var_names %in% params[i])
  grouped <- c(grouped, get.rows)
}

#only params of interest - put back into mcmc.list object
nlist <- do.call(coda::mcmc.list, out[,grouped])

#calculate rhat
rhats <- round(gelman.diag(nlist, multivariate = FALSE)$psrf[,1], digits = 4)
rh_df <- data.frame(rhat = rhats)


#PPC
params = c('pv.mn', 'pv.sd')
grouped2 <- c()
for (i in 1:length(params))
{
  get.rows <- which(var_names %in% params[i])
  grouped2 <- c(grouped2, get.rows)
}

means <- apply(out[[1]][,grouped2], 2, mean)




#create directory
system(paste0('mkdir ', NAME))
setwd(paste0(NAME))

#set max number of rows to print to 5k
options(max.print = 5000)

#write results to text file
sink(paste0('results_', NAME,'.txt'))
print(paste0(NAME))
print(paste0('Number of nests: ', nests))
print(paste0('Number of time steps: ', n_ts))
print(paste0('Total minutes: ', round(tt, digits = 2)))
print(paste0('Total iterations: ', n_final))
print(paste0('n_adapt: ', n_adapt))
print(paste0('n_burn: ', n_burn))
print(paste0('n_draw: ', n_draw))
print(paste0('n_thin: ', n_thin))
print(paste0('n_chain: ', n_chain))
print(paste0('Extra: ', EXTRA))
print(paste0('Rhat_max: ', Rhat_max))
print(paste0('n_max: ', n_max))
print(rh_df)
print(paste0('Posterior Predictive Check:'))
print(means)
sink()




MCMCsummary(out, digits = 4)
sd.y.ch <- MCMCchains(out, params = 'sd.y')
sd.y.new.ch <- MCMCchains(out, params = 'sd.y.new')

plot(sd.y.ch, sd.y.new.ch, pch = '.',
     ylim = c(0.6, 0.8), xlim = c(0.6, 0.8))
abline(0,1, col = 'red')


mn.y.ch <- MCMCchains(out, params = 'mn.y')
mn.y.new.ch <- MCMCchains(out, params = 'mn.y.new')

plot(mn.y.ch, mn.y.new.ch, pch = '.',
     ylim = c(0.9, 1.1), xlim = c(0.9, 1.1))
abline(0,1, col = 'red')

