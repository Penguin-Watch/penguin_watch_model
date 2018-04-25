##############################
# Kerry and Schaub 2012 - 7.5.3. Individual random effects
#
# Data simulated - no trend across time - survival different for each individual
# Mark-recapture for single individuals (1/0 rather than 2/1/0)
# Model p constant across time and individuals
# Model phi constant across time - random effect for individual
#
# Code modified to simulate data where all individuals are marked at time step 1
#
# Similar to No_trend_model.R, PPC (sd) indicates model doesn't replicate data sd well. sd not a good PPC metric in mark-recapture
##############################


rm(list = ls())

require(rjags)
require(MCMCvis)
require(parallel)

n.occasions <- 20                 # Number of capture occasions
marked <- 570
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
v.ind <- 0.5

# Draw annual survival probabilities
logit.phi <- rnorm(sum(marked), qlogis(mean.phi), v.ind^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = marked, byrow = FALSE)
P <- matrix(p, ncol = n.occasions-1, nrow = marked)


# Define function to simulate a capture-history (CH) matrix
set.seed(1)
sim_data_fun <- function(PHI_MAT, P_MAT, N_NESTS = 570)
{
  # PHI_MAT <- PHI
  # P_MAT <- P
  # N_NESTS <- 570
  
  TS_LEN <- NCOL(PHI_MAT) + 1
  CH <- matrix(0,
               ncol = TS_LEN,
               nrow = N_NESTS)
  
  for (i in 1:N_NESTS)
  {
    #i <- 1
    #both chicks alive at start
    t_SP <- c(1, rep(NA, TS_LEN-1))
    for (t in 2:TS_LEN)
    {
      #t <- 2
      #TRUE STATE
      t_SP[t] <- rbinom(1, size = t_SP[t-1], prob = PHI_MAT[i,t-1])
      
      #OBSERVED STATE
      t_DP <- rbinom(1, size = t_SP[t], prob = P_MAT[i,t-1])
      CH[i,t] <- t_DP
    }
    CH[i,1] <- 1
  }
  return(CH)
}


known.state.cjs <- function(ch)
{
  #ch <- CH
  state <- ch
  for (i in 1:dim(ch)[1])
  {
    #i <- 1
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}


# Simulate capture-histories
CH <- sim_data_fun(PHI, P)

# Bundle data
DATA <- list(y = CH, n.occasions = dim(CH)[2], z = known.state.cjs(CH))


# Specify model in BUGS language
sink("cjs-ind-raneff.jags")
cat("
    model {
    
    # Likelihood 
    for (i in 1:570)
    {
    # Define latent state at first capture
    z[i,1] <- 1
    y.new[i,1] ~ dbern(p[i,1])
    
    for (t in (2:20))
    {
    
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- phi[i,t-1] * z[i,t-1]
    
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <- p[i,t-1] * z[i,t]
    
    y.new[i,t] ~ dbern(mu2[i,t])
    
    } #t
    } #i
    
    #PPC
    #mean
    mn.y <- mean(y)
    mn.y.new <- mean(y.new)
    pv.mn <- step(mn.y.new - mn.y)
    
    #sd
    sd.y <- sd(y)
    sd.y.new <- sd(y.new)
    pv.sd <- step(sd.y.new - sd.y)   
    
    
    # Priors and constraints
    for (i in 1:570)
    {
    for (t in 1:(n.occasions-1))
    {
    logit(phi[i,t]) <- mu + epsilon[i]
    p[i,t] <- mean.p
    } #t
    } #i
    
    for (i in 1:570)
    {
    epsilon[i] ~ dnorm(0, tau)
    }
    
    mean.phi ~ dunif(0, 1)                   # Prior for mean survival 
    mu <- log(mean.phi / (1-mean.phi))       # Logit transformation
    sigma ~ dunif(0, 5)                      # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)
    mean.p ~ dunif(0, 1)                     # Prior for mean recapture 
    
    
    }
    ",fill = TRUE)
sink()


# Inits -------------------------------------------------------------------


# Initial values 

Inits_1 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma = runif(1, 0, 2),
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma = runif(1, 0, 2),
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma = runif(1, 0, 2),
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)


# params ------------------------------------------------------------------


# Parameters monitored
Pars <- c("mean.phi", 
                "mean.p", 
                "sigma",
                "pv.mn",
                "pv.sd",
                "sd.y",
                "sd.y.new",
                "mn.y",
                "mn.y.new")





# Inputs for MCMC ---------------------------------------------------------

JAGS_FILE <- 'cjs-ind-raneff.jags'
n_adapt <- 5000  # number for initial adapt
n_burn <- 20000 # number burnin
n_draw <- 10000  # number of final draws to make
n_thin <- 2    # thinning rate
n_chain <- 3  # number of chains


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





MCMCsummary(out, digits = 4)


sd.y.ch <- MCMCchains(out, params = 'sd.y')
sd.y.new.ch <- MCMCchains(out, params = 'sd.y.new')

plot(sd.y.ch, sd.y.new.ch, pch = '.', asp = 1)
abline(0,1, col = 'red')


mn.y.ch <- MCMCchains(out, params = 'mn.y')
mn.y.new.ch <- MCMCchains(out, params = 'mn.y.new')

plot(mn.y.ch, mn.y.new.ch, pch = '.', asp = 1)
abline(0,1, col = 'red')



