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
#do the trends in both detection and survival make for nonidentifiability? non-id of surv and detect in a particular year - plot chains against one another (for params think might not be identifiable) p219 Kerry and Schaub 2012
#maybe make function assymtotic rather than linear for varying surv and detection?
#change sim data to linear



#plot p[1,1] against phi[1,1] - CHECK
#check correlation p and phi - CHECK
#compare estimates p and phi against true p and phi - 
#plot posteriors for:
#1) betas - dnorm(0, 1) T(-1,1) - CHECK
#2) mean_phi and mean_p - dunif(0,0.1) - CHECK
#3) eps_phi and eps_p - dnorm(0, tau_p) T(-20,20) - p okay, phi maybe weakly identifiable
#and compare to priors




# Clear environment -------------------------------------------------------

rm(list = ls())


# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman')
}

pacman::p_load(rjags, MCMCvis, parallel)



# Load data ---------------------------------------------------------------


#setwd('Data')
#data <- read.csv('XXXX.csv', header=TRUE)

#ISSUES
#also have possibly more than two chicks per cell in some cases when older
#what to do about no observation over night - can't be ignored right? Treating each hour as one time step here
#might not have to assume both chicks alive at start - first sighting of 2 chicks can be start - p182 Kerry and Schaub 2012
#all years should be run hierarchically for a site, and all sites hierarchically for each species?


#simulate new data - script modified from Kerry and Schaub 2012
n_ts <- 400 #number of time steps
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
          'sigma_p2',
          'sigma_phi2',
          'beta_p',
          'beta_phi',
          'mu_phi',
          'mu_p',
          'eps_phi',
          'eps_p')


# Inputs for MCMC ---------------------------------------------------------

n_adapt <- 10000  # number for initial adapt
n_burn <- 40000 # number burnin
n_draw <- 20000  # number of final draws to make
n_thin <- 2    # thinning rate
n_chain <- 3  # number of chains

Rhat_max <- 1.02 # max allowable Rhat (close to 1 = convergence)
n_max <- 1e5 # max allowable iterations


# Run model (parallel) ---------------------------------------------------------------


#run_jags <- function(JAGS_FILE)

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
                          'F_Inits'
                        ))


ptm <- proc.time()
out.1 <- clusterEvalQ(cl, 
                      {
                        require(rjags)
                        processNum <- which(pid==Sys.getpid())
                        m.inits <- F_Inits[[processNum]]
                        
                        jm = jags.model(data = DATA, 
                                        file = "mark_recapture.jags", 
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


out <- mcmc.list(out.1[[1]][[1]], 
                 out.1[[2]][[1]], 
                 out.1[[3]][[1]])
#(proc.time() - ptm)[3]/60 #minutes



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
  
  out <- mcmc.list(out.2[[1]][[1]], 
                   out.2[[2]][[1]], 
                   out.2[[3]][[1]])
  
  n_extra <- n_extra + n_draw
  n_total <- n_total + n_draw
}
stopCluster(cl)
(proc.time() - ptm)[3]/60 #minutes

n_final <- floor((n_draw + n_extra)/n_thin)


# Run model - non-parallel ------------------------------------------------


# #rjags
# ptm <- proc.time()
# jm = jags.model(data = DATA,
#                 file = "mark_recapture.jags",
#                 inits = F_Inits,
#                 n.chains = 3,
#                 n.adapt = n_adapt)
# 
# update(jm, n.iter = n_burn)
# 
# out <- coda.samples(jm,
#                    n.iter = n_draw,
#                    variable.names = Pars,
#                    thin = n_thin)
# 
# 
# #extra draws if didn't converge
# n_total <- n_burn + n_draw
# n_extra <- 0
# while(max(MCMCsummary(out)[,5], na.rm = TRUE) > Rhat_max &
#       n_total < n_max)
# {
# 
#   out <- coda.samples(jm,
#                       n.iter = n_draw,
#                       variable.names = Pars,
#                       n.thin = n_thin)
# 
#   n_extra <- n_extra + n_draw
#   n_total <- n_total + n_draw
# }
# #(proc.time() - ptm)[3]/60 #minutes
# n_final <- floor((n_draw + n_extra)/n_thin)

#Inferences were derived from $`r n_final`$ samples drawn following an adaptation period of $`r n_adapt`$ draws, and a burn-in period of $`r (n_total - n_draw)`$ draws using $`r n_chain`$ chains and a thinning rate of $`r n_thin`$.




# Analyze posterior -------------------------------------------------------

#phi = survival prob
#p = detection prob

saveRDS(out, 'model_l_out.rds')
#out <- readRDS('model_out.rds')

#summary
MCMCtrace(out, ind = TRUE, pdf = TRUE)

MCMCsummary(out, digits = 4)

#cor of posteriors of p with posteriors of phi
pb <- txtProgressBar(min = 0, max = 30, style = 3)
CO <- matrix(nrow = 30, ncol = 100)
for (i in 1:NROW(DATA$y))
{
  for (j in 1:NCOL(DATA$y))
  {
    t_p <- MCMCchains(out, paste0('p[',i,',',j,']'))
    t_phi <- MCMCchains(out, paste0('phi[',i,',',j,']'))
    CO[i,j] <- cor(t_p, t_phi)
  }
  setTxtProgressBar(pb, i)
}
close(pb)

#saveRDS(CO, 'p_phi_cor.rds')




#correlation of beta_p and beta_phi
beta_p_ch <- MCMCchains(out, 'beta_p', excl = 'beta_phi')
beta_phi_ch <- MCMCchains(out, 'beta_phi')
plot(beta_p_ch, beta_phi_ch, pch = '.')
cor(beta_p_ch, beta_phi_ch)

#correlation of mean_p and mean_phi
mean_p_ch <- MCMCchains(out, 'mean_p', excl = 'mean_phi')
mean_phi_ch <- MCMCchains(out, 'mean_phi')
plot(mean_p_ch, mean_phi_ch, pch = '.')
cor(mean_p_ch, mean_phi_ch)





#posterior estimates for p
pb <- txtProgressBar(min = 0, max = 30, style = 3)
est_p <- matrix(nrow = 30, ncol = 100)
for (i in 1:NROW(DATA$y))
{
  for (j in 1:NCOL(DATA$y))
  {
    est_p[i,j] <- median(MCMCchains(out, paste0('p[',i,',',j,']')))
  }
  setTxtProgressBar(pb, i)
}
close(pb)


#posterior estimates for phi
pb <- txtProgressBar(min = 0, max = 30, style = 3)
est_phi <- matrix(nrow = 30, ncol = 100)
for (i in 1:NROW(DATA$y))
{
  for (j in 1:NCOL(DATA$y))
  {
    est_phi[i,j] <- median(MCMCchains(out, paste0('phi[',i,',',j,']')))
  }
  setTxtProgressBar(pb, i)
}
close(pb)



#compare posteriors to actual quantities






#plot posterior
#1) betas - dnorm(0, 1) T(-1,1) - check
#2) mean_phi and mean_p - dunif(0,0.1) - check
#3) eps_phi and eps_p - dnorm(0, tau_p) T(-20,20) - p okay, phi not converged?

beta_p_ch <- MCMCchains(out, 'beta_p', excl = 'beta_phi')
plot(density(beta_p_ch))
#prior
a <- rnorm(100000, 0, 1)
to.rm <- which(a > 1 | a < -1)
b <- a[-to.rm]
c <- b[1:7500] #first 5k
lines(density(c), col = 'red')

beta_phi_ch <- MCMCchains(out, 'beta_phi')
plot(density(beta_phi_ch))
#prior
a <- rnorm(100000, 0, 1)
to.rm <- which(a > 1 | a < -1)
b <- a[-to.rm]
c <- b[1:7500] #first 5k
lines(density(c), col = 'red')

mean_p_ch <- MCMCchains(out, 'mean_p', excl = 'mean_phi')
plot(density(mean_p_ch))
#prior
a2 <- runif(5000, 0, 1)
lines(density(a2), col = 'red')

mean_phi_ch <- MCMCchains(out, 'mean_phi')
plot(density(mean_phi_ch))
#prior
a2 <- runif(5000, 0, 1)
lines(density(a2), col = 'red')


sigma_p_ch <- sqrt(MCMCchains(out, 'sigma_p2'))
eps_p_1_ch <- MCMCchains(out, 'eps_p[1]')
plot(density(eps_p_1_ch))
eps_p_2_ch <- MCMCchains(out, 'eps_p[2]')
plot(density(eps_p_2_ch))
eps_p_3_ch <- MCMCchains(out, 'eps_p[3]')
plot(density(eps_p_3_ch))
eps_p_4_ch <- MCMCchains(out, 'eps_p[4]')
plot(density(eps_p_4_ch))
#prior
a3 <- rnorm(10000, 0, median(sigma_p_ch))
to.rm3 <- which(a3 > 20 | a3 < -20)
b3 <- a3#[-to.rm3]
c3 <- b3[1:7500]
lines(density(c3), col = 'red')

sigma_phi_ch <- sqrt(MCMCchains(out, 'sigma_phi2'))
eps_phi_1_ch <- MCMCchains(out, 'eps_phi[1]')
plot(density(eps_phi_1_ch))
eps_phi_2_ch <- MCMCchains(out, 'eps_phi[2]')
plot(density(eps_phi_2_ch))
eps_phi_3_ch <- MCMCchains(out, 'eps_phi[3]')
plot(density(eps_phi_3_ch))
eps_phi_4_ch <- MCMCchains(out, 'eps_phi[4]')
plot(density(eps_phi_4_ch))
eps_phi_5_ch <- MCMCchains(out, 'eps_phi[5]')
plot(density(eps_phi_5_ch))
#prior
a3 <- rnorm(10000, 0, median(sigma_phi_ch))
to.rm3 <- which(a3 > 20 | a3 < -20)
b3 <- a3#[-to.rm3]
c3 <- b3[1:7500]
lines(density(c3), col = 'red')




MCMCsummary(out, params = 'eps_phi', digits = 4)



