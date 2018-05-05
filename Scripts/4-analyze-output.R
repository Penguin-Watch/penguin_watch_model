#################
# Penguin Watch Model - 4 - Analyze model output
#
# 1-process-krill-data.R | process krill data
# 2-process-SIC-data.R | process SIC data
# 3-model.R | penguin model
# 3-run-model.pbs | pbs script to run penguin model on HPC resources
# 4-analyze-output.R | analyze model output
#
# Author: Casey Youngflesh
#################



# Clear environment -------------------------------------------------------

rm(list = ls())



# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman')
}

pacman::p_load(MCMCvis)




# Load data -------------------------------------------------------

#phi = survival prob
#p = detection prob

NAME <- 'May_3_2018_covariates'

setwd(paste0('Results/', NAME))

out <- readRDS(paste0(NAME, '.rds'))



# Summarize ---------------------------------------------------------------

MCMCsummary(out)

MCMCsummary(out, params = 'mu_phi')
MCMCsummary(out, params = 'beta_phi')
MCMCsummary(out, params = 'sigma_eta_phi')
MCMCsummary(out, params = 'sigma_gamma_phi')
MCMCsummary(out, params = 'mu_p')
MCMCsummary(out, params = 'beta_p')
MCMCsummary(out, params = 'sigma_nu_p')
MCMCsummary(out, params = 'pi_phi')
MCMCsummary(out, params = 'rho_phi')





# PPO ---------------------------------------------------------------------

#' mu_phi ~ dnorm(0, 0.386)   
PR <- rnorm(15000, 0, 1/sqrt(0.386))
MCMCtrace(out, 
          params = 'mu_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' beta_phi ~ dnorm(0, 10) T(0,1)
tt <- rnorm(15000, 0, 1/sqrt(10))
PR <- tt[tt >= 0 & tt <= 1]
MCMCtrace(out, 
          params = 'beta_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' sigma_eta_phi ~ dunif(0, 100)
PR <- runif(15000, 0, 100)
MCMCtrace(out, 
          params = 'sigma_eta_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' sigma_gamma_phi ~ dunif(0, 100)
PR <- runif(15000, 0, 100)
MCMCtrace(out, 
          params = 'sigma_gamma_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)


#' mu_p ~ dnorm(0, 0.386)
PR <- rnorm(15000, 0, 1/sqrt(0.386))
MCMCtrace(out, 
          params = 'mu_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' beta_p ~ dnorm(0, 1) T(0,1)
tt <- rnorm(15000, 0, 1/sqrt(1))
PR <- tt[tt >= 0 & tt <= 1]
MCMCtrace(out, 
          params = 'beta_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' sigma_nu_p ~ dunif(0.25, 3)
PR <- runif(15000, 0, 10)
MCMCtrace(out, 
          params = 'sigma_nu_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' pi_phi ~ dnorm(0, 0.386)
PR <- rnorm(15000, 0, 1/sqrt(0.386))
MCMCtrace(out, 
          params = 'pi_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' rho_phi ~ dnorm(0, 0.386)
PR <- rnorm(15000, 0, 1/sqrt(0.386))
MCMCtrace(out, 
          params = 'rho_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)



# track detection and surv -----------------------------------------------------

#grand intercept - survival
mu_phi <- MCMCpstr(out, params = 'mu_phi', func = median, digits = 3)[[1]]
#site effect
eta_phi <- MCMCpstr(out, params = 'eta_phi', func = median, digits = 3)[[1]]
#year effect
gamma_phi <- MCMCpstr(out, params = 'gamma_phi', func = median, digits = 3)[[1]]
#slope (time) survival
beta_phi <- MCMCpstr(out, params = 'beta_phi', func = median, digits = 3)[[1]]*1:840

#grand intercept - detection
mu_p <- MCMCpstr(out, params = 'mu_p', func = median, digits = 3)[[1]]
#slope (time) detection
beta_p <- MCMCpstr(out, params = 'beta_p', func = median, digits = 3)[[1]]*1:840
#nest effect
nu_p <- MCMCpstr(out, params = 'nu_p', func = median, digits = 3)[[1]]


#surv <- data.frame(SITE = , YEAR = , NEST = )

#survival
#for (k in 1:NK)
#{
  k <- 2
#  for (j in 1:NJ)
#  {
    j <- 1
    surv <- inv.logit(mu_phi + eta_phi[k] + gamma_phi[j])
#  }
#}


#detect <- data.frame(SITE = , YEAR = , NEST = )

#detection
#for (k in 1:NK)
#{
   k <- 2
#  for (j in 1:NJ)
#  {
    j <- 1
    for (i in 1:NI[j,k])
    {
      i <- 1
      detect <- inv.logit(mu_p + nu_p[i,j,k])
    }
#  }
#}


    
    
    
# No correlations between parameters -----------------------------------------------

#correlation of beta_p and beta_phi
beta_p_ch <- MCMCchains(out, 'beta_p')
beta_phi_ch <- MCMCchains(out, 'beta_phi')
plot(beta_p_ch, beta_phi_ch, pch = '.')
cor(beta_p_ch, beta_phi_ch)

#correlation of mu_p and mu_phi
mean_p_ch <- MCMCchains(out, 'mu_p')
mean_phi_ch <- MCMCchains(out, 'mu_phi')
plot(mean_p_ch, mean_phi_ch, pch = '.')
cor(mean_p_ch, mean_phi_ch)




