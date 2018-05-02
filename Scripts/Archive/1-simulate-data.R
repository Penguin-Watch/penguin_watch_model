######################
#Mark-recapture for penguin watch 
#
#script to create data
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

pacman::p_load(rjags, MCMCvis, parallel)



# Load data ---------------------------------------------------------------


#setwd('Data')
#data <- read.csv('XXXX.csv', header=TRUE)


#simulate new data - script modified from Kerry and Schaub 2012
n_ts <- 100 #number of time steps
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

#saveRDS(PHI, 'PHI.rds')




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


#saveRDS(P, 'P.rds')



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


#saveRDS(sim_data, 'sim_data.rds')





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


#saveRDS(z_vals, 'z_vals.rds')
