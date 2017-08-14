################
#Run pwatch mark-recapture script
################



# Clear environment -------------------------------------------------------

rm(list = ls())


# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman', repos = "http://cran.case.edu")
}
pacman::p_load(MCMCvis, rjags, parallel)


# Load data ---------------------------------------------------------------


#setwd('Data')
#data <- read.csv('XXXX.csv', header=TRUE)


#simulate new data - script modified from Kerry and Schaub 2012
n_ts <- 200 #number of time steps
x <- 1:n_ts
nests <- 30 #number of nests


#survival probability

sim_p_fun <- function(START, RATE = 0.008, TOP = 1)
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

phi_data <- sim_p_fun(0.985)

#check trend
#cov <- c(1:length(phi_data))
#plot(cov, phi_data, ylim = c(0.98,1), type ='l')

#nls
#fit_nls <- nls(phi_data ~ SSlogis(cov, Asym, xmid, scal))
#summary(fit_nls)
#lines(cov, predict(fit_nls), col = 'red')

#lm
#fit_lm <- lm(phi_data ~  cov)
#summary(fit_lm)
#abline(fit_lm, col = 'red')


PHI <- matrix(rep(phi_data, nests),
              nrow = nests,
              ncol = n_ts-1)


#detection probability
#starting probs for each nest
dp <- runif(30, 0.3, 0.8)

P <- matrix(rep(NA, nests*(n_ts-1)),
            nrow = nests,
            ncol = n_ts-1)

for (i in 1:length(dp))
{
  P[i,] <- sim_p_fun(dp[i], RATE = 0.005, TOP = 0.95)
  #plot(p_data, type = 'l', ylim = c(0,1), main = paste0(i))
}


#check trend
#lm
# cov_p <- 1:length(P[1,])
# fit_p <- lm(P[1,] ~ cov_p)
# summary(fit_p)






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
      logit(phi[i,t]) <- mu_phi + beta_phi*x[t] + eps_phi[t]       #phi = survival prob
      logit(p[i,t]) <- mu_p + beta_p*x[t] + eps_p[i]               #p = detection prob
      }
      }
      
      
      #priors
      for (t in 1:L)
      {
      eps_phi[t] ~ dnorm(0, tau_phi) T(-10,10)
      }
      
      mean_phi ~ dunif(0,1)
      mu_phi <- log(mean_phi / (1 - mean_phi))    
      tau_phi <- pow(sigma_phi, -2)
      sigma_phi ~ dunif(0, 10)
      sigma_phi2 <- pow(sigma_phi, 2)
      
      for (i in 1:N)
      {
      eps_p[i] ~ dnorm(0, tau_p) T(-10,10)
      }
      
      mean_p ~ dunif(0,1)                        #Mean survival - could use alternative below
      mu_p <- log(mean_p / (1 - mean_p))         #Logit transform - could use alternative below
      #mean_p <- 1 / (1+exp(-mu_p))              #Mean survival - Inv-logit transform    
      #mu_p ~ dnorm(0, 0.001)                    #Prior for logit of mean survival
      tau_p <- pow(sigma_p, -2)
      sigma_p ~ dunif(0, 10)
      sigma_p2 <- pow(sigma_p, 2)
      
      beta_phi ~ dnorm(0, 1) T(-1,1)
      beta_p ~ dnorm(0, 1) T(-1,1)
      
      
      }",fill = TRUE)

  sink()
}



# Starting values ---------------------------------------------------------


Inits_1 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma_phi = runif(1, 0, 10),
                sigma_p = runif(1, 0, 10),
                beta_phi = 0,
                beta_p = 0,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma_phi = runif(1, 0, 10),
                sigma_p = runif(1, 0, 10),
                beta_phi = 0,
                beta_p = 0,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma_phi = runif(1, 0, 10),
                sigma_p = runif(1, 0, 10),
                beta_phi = 0,
                beta_p = 0,
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)



# Parameters to track -----------------------------------------------------

Pars <- c('mean_phi',
          'mean_p',
          'sigma_p',
          'sigma_phi',
          'beta_p',
          'beta_phi',
          'mu_phi',
          'mu_p',
          'eps_phi',
          'eps_p')


# Inputs for MCMC ---------------------------------------------------------

JAGS_FILE <- 'mark_recapture.jags'
n_adapt <- 10000  # number for initial adapt
n_burn <- 100000 # number burnin
n_draw <- 20000  # number of final draws to make
n_thin <- 2    # thinning rate
n_chain <- 3  # number of chains

Rhat_max <- 1.02 # max allowable Rhat (close to 1 = convergence)
n_max <- 200000 # max allowable iterations


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




#more iterations if not converged
n_total <- n_burn + n_draw
n_extra <- 0
while(max(MCMCsummary(out)[,5], na.rm = TRUE) > Rhat_max &
      n_total < n_max)
{
  
  out.2 <- clusterEvalQ(cl, 
                        {
                          require(rjags)
                          processNum <- which(pid==Sys.getpid())
                          m.inits <- F_Inits[[processNum]]
                          
                          samples = coda.samples(jm, 
                                                 n.iter = n_draw, 
                                                 variable.names = Pars, 
                                                 thin = n_thin)
                          return(samples)
                        })
  
  out <- coda::mcmc.list(out.2[[1]][[1]], 
                         out.2[[2]][[1]], 
                         out.2[[3]][[1]])
  
  n_extra <- n_extra + n_draw
  n_total <- n_total + n_draw
}

stopCluster(cl)

n_final <- floor((n_draw + n_extra)/n_thin)
print(paste0('Total iterations: ', n_final))
(proc.time() - ptm)[3]/60 #minutes


#Inferences were derived from $`r n_final`$ samples drawn following an adaptation period of $`r n_adapt`$ draws, and a burn-in period of $`r (n_total - n_draw)`$ draws using $`r n_chain`$ chains and a thinning rate of $`r n_thin`$.




# Analyze posterior -------------------------------------------------------

#phi = survival prob
#p = detection prob

#adapt_burn_draw_time_rhatthreshold
saveRDS(out, 'out_10a_100b_20d_200t_102.rds')
#out <- readRDS('model_l_out.rds')
