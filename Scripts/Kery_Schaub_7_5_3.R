
# 7.5.3. Individual random effects
# modified to simulate data where all individuals are marked at time step 1

rm(list = ls())

require(R2jags)


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


# # Function to create a matrix of initial values for latent state z
# cjs.init.z <- function(ch,f){
#   for (i in 1:dim(ch)[1]){
#     if (sum(ch[i,])==1) next
#     n2 <- max(which(ch[i,]==1))
#     ch[i,f[i]:n2] <- NA
#   }
#   for (i in 1:dim(ch)[1]){
#     ch[i,1:f[i]] <- NA
#   }
#   return(ch)
# }


# Simulate capture-histories
CH <- sim_data_fun(PHI, P)

# Bundle data
jags.data <- list(y = CH, n.occasions = dim(CH)[2], z = known.state.cjs(CH))


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


# Initial values 
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), sigma = runif(1, 0, 2))}  

# Parameters monitored
parameters <- c("mean.phi", 
                "mean.p", 
                "sigma2",
                "pv.mn",
                "pv.sd")

# MCMC settings
ni <- 5000
nt <- 6
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 73 min)
cjs.ind <- jags(jags.data, inits, parameters, "cjs-ind-raneff.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.ind, digits = 3)

# Produce graph
par(mfrow = c(1, 2), las = 1)
hist(cjs.ind$BUGSoutput$sims.list$mean.phi, nclass = 25, col = "gray", main = "", xlab = expression(bar(phi)), ylab = "Frequency")
abline(v = mean.phi, col = "red", lwd = 2)
hist(cjs.ind$BUGSoutput$sims.list$sigma2, nclass = 15, col = "gray", main = "", xlab = expression(sigma^2), ylab = "Frequency", xlim = c(0, 3))
abline(v = v.ind, col = "red", lwd = 2)
