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
out <- readRDS('out_8a_20b_10d_200t_trackpphi_PPC.rds')


#summary
MCMCtrace(out, ind = TRUE, pdf = TRUE, iter = 10000)

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
pv.mn_sd <- MCMCchains(out, params = 'pv.sd')
mean(pv.mn_ch)



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





#need to track p and phi params
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

