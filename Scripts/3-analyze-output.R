######################
#Mark-recapture for penguin watch 
#
#script to analyze model output
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




# Analyze posterior -------------------------------------------------------

#phi = survival prob
#p = detection prob

setwd('HPC/Archive')
NAME <- 'out_8a_50b_20d_200t_PPC_diag'


out <- readRDS(paste0(NAME, '.rds'))



#summary
MCMCtrace(out, ind = TRUE, pdf = TRUE, iter = 10000)
MCMCsummary(out, digits = 4, params = c('mean_phi', 'mean_p'))
MCMCsummary(out, digits = 4, 
            params = c('mean_phi',
                       'mean_p',
                       'sigma_p',
                       'sigma_phi',
                       'beta_p',
                       'beta_phi',
                       'mu_phi',
                       'mu_p',
                       'eps_phi',
                       'eps_p',
                       'pv.mn',
                       'pv.sd'))


pv.mn_ch <- MCMCchains(out, params = 'pv.mn')
pv.sd_ch <- MCMCchains(out, params = 'pv.sd')
mean(pv.mn_ch)
mean(pv.sd_ch)







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








# Compare posterior p and phi estimates to generating params --------------

#simulate new data - to match with data input to model
n_ts <- 200 #number of time steps
x <- 1:n_ts
nests <- 30 #number of nests

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





#POSTERIOR p
est_p_matrix <- readRDS(paste0('est_p_', NAME, '.rds'))


#TRUE p - simualted in same way as model
P <- matrix(rep(NA, nests*(n_ts-1)),
            nrow = nests,
            ncol = n_ts-1)

set.seed(1) #to match with simulated data in HPC model script
dp <- runif(30, 0.3, 0.8)
for (i in 1:length(dp))
{
  P[i,] <- sim_p_fun(dp[i], RATE = 0.005, TOP = 0.95)
  #plot(p_data, type = 'l', ylim = c(0,1), main = paste0(i))
}
sim_P <- cbind(rep(NA, nrow(P)), P)

#compare two matrices - don't compare first time step, as those values are meaningless (do not affect simulated data bc those values are fixed in the fixed data - true state 2 chicks, 2 chick observed)

#estimated in BLACK - simulated model input in RED
plot(est_p_matrix[1,], type = 'l', ylim = c(0.3, 1))
for (i in 1:NROW(est_p_matrix))
{
  lines(est_p_matrix[i,])
}
for (i in 1:NROW(sim_P))
{
  lines(sim_P[i,], col = 'red')
}



#POSTERIOR PHI
est_phi_matrix <- readRDS(paste0('est_phi_', NAME, '.rds'))

#TRUE phi - simulated in same way as model input
phi_data <- sim_p_fun(START = 0.985)
PHI <- matrix(rep(phi_data, nests),
              nrow = nests,
              ncol = n_ts-1,
              byrow = TRUE)
sim_PHI <- cbind(rep(NA, nrow(PHI)), PHI)


#compare two matrices
plot(est_phi_matrix[1,], type = 'l', ylim = c(0.98, 1))
lines(sim_PHI[1,], col = 'red')



#correlation of posterior of p with posteriors of phi

p_phi_cor_matrix <- readRDS(paste0('p_phi_cor_', NAME, '.rds'))

hist(p_phi_cor_matrix) # no correlation










#plot posteriors on top of priors
# beta_p
# beta_phi
# mean_p
# mean_phi
# sigma_p
# sigma_phi
# eps_p
# eps_phi



#plot posterior
#1) betas - dnorm(0, 1) T(-1,1) - check
#2) mean_phi and mean_p - dunif(0,0.1) - check
#3) eps_phi and eps_p - dnorm(0, tau_p) T(-20,20) - p okay, phi not converged?

#https://stats.stackexchange.com/questions/97596/how-to-calculate-overlap-between-empirical-probability-densities
#define limits of where overlap is to be calculated
#prior (chain of prior)
#post (chain of posterior)
lower <- min(c(a, b)) - 1 
upper <- max(c(a, b)) + 1

# generate kernel densities
d_prior <- density(prior, from = lower, to = upper)
d_post <- density(post, from = lower, to = upper)
d <- data.frame(x = d_prior$x, prior = d_prior$y, post = d_post$y)

# calculate intersection densities
d$int <- pmin(d$prior, d$post)

#plot
plot(d$x, d$prior, xlim = c(lower, upper), ylim = c(0, 0.5), type = 'l') #prior
lines(d$x, d$post, col = 'red') #posterior
lines(d$x, d$int, col = 'green') #intersection

# integrate areas under curves
library(sfsmisc)
total <- integrate.xy(d$x, d$prior) + integrate.xy(d$x, d$posterior) #total area
#area of just to intersecting bit
intersection <- integrate.xy(d$x, d$int)

# compute overlap coefficient
(overlap <- 2 * intersection / total)






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



sigma_p_ch <- MCMCchains(out, 'sigma_p')
#prior
a3 <- rnorm(10000, 0, median(sigma_p_ch))
to.rm3 <- which(a3 > 20 | a3 < -20)
b3 <- a3#[-to.rm3]
c3 <- b3[1:7500]
eps_p_1_ch <- MCMCchains(out, 'eps_p[1]')
plot(density(eps_p_1_ch))
lines(density(c3), col = 'red')
eps_p_2_ch <- MCMCchains(out, 'eps_p[2]')
plot(density(eps_p_2_ch))
lines(density(c3), col = 'red')
eps_p_3_ch <- MCMCchains(out, 'eps_p[3]')
plot(density(eps_p_3_ch))
lines(density(c3), col = 'red')
eps_p_4_ch <- MCMCchains(out, 'eps_p[4]')
plot(density(eps_p_4_ch))
lines(density(c3), col = 'red')



#lots of overlap with prior - to be expected since variance was modeled hierarchically?
sigma_phi_ch <- MCMCchains(out, 'sigma_phi')
#prior
a3 <- rnorm(10000, 0, median(sigma_phi_ch))
to.rm3 <- which(a3 > 20 | a3 < -20)
b3 <- a3#[-to.rm3]
c3 <- b3[1:7500]
eps_phi_1_ch <- MCMCchains(out, 'eps_phi[1]')
plot(density(eps_phi_1_ch))
lines(density(c3), col = 'red')
eps_phi_2_ch <- MCMCchains(out, 'eps_phi[2]')
plot(density(eps_phi_2_ch))
lines(density(c3), col = 'red')
eps_phi_3_ch <- MCMCchains(out, 'eps_phi[3]')
plot(density(eps_phi_3_ch))
lines(density(c3), col = 'red')
eps_phi_4_ch <- MCMCchains(out, 'eps_phi[4]')
plot(density(eps_phi_4_ch))
lines(density(c3), col = 'red')
eps_phi_5_ch <- MCMCchains(out, 'eps_phi[5]')
plot(density(eps_phi_5_ch))
lines(density(c3), col = 'red')

