#################
# Penguin Watch Model - 4 - Analyze model output
#
# 0-detect-params.R | recovering generating values using one detection param vs many detection params
# 00-recover-params.R | recovering generating values using realistic data/model
# 1-process-krill-data.R | process krill data
# 2-process-SIC-data.R | process SIC data
# 3-process-pw-data.R | process PW Pro data
# 4-model.R | penguin model
# 4-run-model.pbs | pbs script to run penguin model on HPC resources
# 5-analyze-output.R | analyze model output
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

pacman::p_load(MCMCvis, boot, dplyr)




# Load data -------------------------------------------------------

#phi = survival prob
#p = detection prob

NAME <- 'May_19_2018_one_detect'
NAME <- 'May_27_2018_one_site'

setwd(paste0('~/Google_Drive/R/penguin_watch_model/Results/', NAME))

out <- readRDS(paste0(NAME, '.rds'))
data <- readRDS('jagsData.rds')


  # Summarize ---------------------------------------------------------------

MCMCsummary(out, round = 4, n.eff = TRUE)
MCMCtrace(out,
          ind = TRUE)

inv.logit(MCMCsummary(out, params = 'mu_phi')[4])



MCMCsummary(out, excl = 'nu_p', round = 2, n.eff = TRUE)

#grand intercept - surv
MCMCsummary(out, params = 'mu_phi', round = 2)
#slope within a season
MCMCsummary(out, params = 'beta_phi', digits = 5)
#site effect
MCMCsummary(out, params = 'eta_phi')
MCMCsummary(out, params = 'sigma_eta_phi')
#year effect
MCMCsummary(out, params = 'gamma_phi')
MCMCsummary(out, params = 'sigma_gamma_phi')
#covariates
MCMCsummary(out, params = 'pi_phi')
MCMCsummary(out, params = 'rho_phi')


#grand intercept - detect
MCMCsummary(out, params = 'mu_p')
#slope within a season
MCMCsummary(out, params = 'beta_p')
#nest detection
MCMCsummary(out, params = 'nu_p')
MCMCsummary(out, params = 'sigma_nu_p')


#plots
MCMCplot(out, excl = 'nu_p')



# diagnose model errors -------------------------------------------------



Inits_1 <- list(mu_phi = 0,
                beta_phi = 0.01,
                eta_phi = rep(0, DATA$NK),
                gamma_phi = rep(0, DATA$NJ),
                pi_phi = 0,
                rho_phi = 0,
                mu_p = 0,
                beta_p = 0.01,
                nu_p = narray)

xvals <- scale(1:768, scale = FALSE)[,1]

# Error in node y[21,2,2,1] - 1
#DATA$y[1:21,2,2,1]
#DATA$w[1:21,2,2,1]
#DATA$z[1:21,2,2,1]
#phi
plot(inv.logit(5 + 0 + 0 + 0.005*xvals + 0 + 0))
#p
plot(inv.logit(-2 + 0.01*xvals + 0))





# PPO ---------------------------------------------------------------------

tf <- function(PR)
{
  hist(inv.logit(PR))
}


#' mu_phi ~ dnorm(0, 0.01)   
PR <- rnorm(15000, 3, 1/sqrt(0.01))
tf(PR)
MCMCtrace(out, 
          params = 'mu_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' eta_phi[k] ~ dnorm(0, 0.386)
PR <- rnorm(15000, 0, 1/sqrt(0.386))
tf(PR)
MCMCtrace(out, 
          params = 'eta_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' gamma_phi[k] ~ dnorm(0, 0.386)
PR <- rnorm(15000, 0, 1/sqrt(0.386))
tf(PR)
MCMCtrace(out, 
          params = 'gamma_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)


#' #' sigma_eta_phi ~ dunif(0, 100)
#' PR <- runif(15000, 0, 100)
#' MCMCtrace(out, 
#'           params = 'sigma_eta_phi',
#'           ind = TRUE, 
#'           priors = PR,
#'           pdf = FALSE,
#'           post_zm = FALSE)

#' #' sigma_gamma_phi ~ dunif(0, 100)
#' PR <- runif(15000, 0, 100)
#' MCMCtrace(out, 
#'           params = 'sigma_gamma_phi',
#'           ind = TRUE, 
#'           priors = PR,
#'           pdf = FALSE,
#'           post_zm = FALSE)


#' mu_p ~ dnorm(0, 0.1)
PR <- rnorm(15000, 2, 1/sqrt(0.1))
tf(PR)
MCMCtrace(out, 
          params = 'mu_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' beta_p ~ dnorm(0, 1000) T(0,0.03)
tt <- rnorm(15000, 0.1, 1/sqrt(10))
PR <- tt[tt >= 0 & tt <= 0.5]
tf(PR)
MCMCtrace(out, 
          params = 'beta_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

#' sigma_nu_p ~ dunif(0.25, 3)
# PR <- runif(15000, 0, 10)
# MCMCtrace(out, 
#           params = 'sigma_nu_p',
#           ind = TRUE, 
#           priors = PR,
#           pdf = FALSE,
#           post_zm = FALSE)

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

#' tau_nu_p ~ dunif(0, 25)
PR <- runif(15000, 0, 25)
MCMCtrace(out, 
          params = 'tau_nu_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)




# track detection and surv -----------------------------------------------------

#grand intercept - survival
mu_phi <- MCMCpstr(out, params = 'mu_phi', func = median)[[1]]
#site effect
eta_phi <- MCMCpstr(out, params = 'eta_phi', func = median)[[1]]
#year effect
gamma_phi <- MCMCpstr(out, params = 'gamma_phi', func = median)[[1]]


#grand intercept - detection
mu_p <- MCMCpstr(out, params = 'mu_p', func = median)[[1]]
#slope (time) detection
beta_p <- MCMCpstr(out, params = 'beta_p', func = median)[[1]] 
#beta_p <- 0.1
beta_p_vals <- beta_p * scale(as.numeric(1:63), scale = FALSE)[,1]
#site/year effect
nu_p <- MCMCpstr(out, params = 'nu_p', func = median)[[1]]


#surv <- data.frame(SITE = , YEAR = , NEST = )

#survival
#for (k in 1:NK)
#{
#k <- 2
#  for (j in 1:NJ)
#  {
#j <- 1
#surv <- inv.logit(mu_phi + eta_phi[k] + gamma_phi[j])
#  }
#}



#detection
#for (k in 1:NK)
#{
k <- 2
#  for (j in 1:NJ)
#  {
j <- 1
  #detect <- inv.logit(mu_p)
  #detect <- inv.logit(mu_p + nu_p[j,k])
  #detect <- inv.logit(mu_p + beta_p_vals)
  detect <- inv.logit(mu_p + beta_p_vals + nu_p[j,k])
  detect_nil <- mu_p + beta_p_vals + nu_p[j,k]
#  }
#}

plot(detect, type = 'l')
plot(detect_nil, type = 'l')
PR <- rnorm(15000, 0, 1/sqrt(0.386))
PR <- rnorm(15000, 0, 1/sqrt(5))
hist(PR)

data$y[1:30,,1,2]



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




# Relationship between krill and SIC --------------------------------------

#no clear relationship between SIC and krill

#50km radius for Dec - Feb
CCAMLR_krill <- read.csv('~/Google_Drive/R/penguin_watch_model/Data/Krill_data/CCAMLR/Processed_CCAMLR/CCAMLR_krill_breeding_season.csv')

#500km radius for June - Sep
SIC <- read.csv('~/Google_Drive/R/penguin_watch_model/Data/SIC_data/Processed/SIC_500_W.csv')

#500km radius MAX over last 5 years
#SIC_max <- read.csv('~/Google_Drive/R/penguin_watch_model/Data/SIC_data/Processed/SIC_500_MAX.csv')

j_kr_SIC <- left_join(CCAMLR_krill, SIC, by = c('YEAR', 'SITE'))

summary(lm(j_kr_SIC$T_KRILL ~ j_kr_SIC$WMN))




# last values from previous model -----------------------------------------

tt <- tail(MCMCchains(out), n = 1)
CHAIN <- 2
Inits2 <- list(beta_p = tail(MCMCchains(out, 
                                        params = 'beta_p', 
                                        chain_num = CHAIN), n = 1),
               eta_phi = c(tail(MCMCchains(out, 
                                           params = 'eta_phi\\[1\\]', 
                                           ISB = FALSE,
                                           chain_num = CHAIN), n = 1),
                           tail(MCMCchains(out, 
                                           params = 'eta_phi\\[2\\]', 
                                           ISB = FALSE,
                                           chain_num = CHAIN), n = 1),
                           tail(MCMCchains(out, 
                                           params = 'eta_phi\\[3\\]', 
                                           ISB = FALSE,
                                           chain_num = CHAIN), n = 1),
                           tail(MCMCchains(out, 
                                           params = 'eta_phi\\[4\\]', 
                                           ISB = FALSE,
                                           chain_num = CHAIN), n = 1)),
               gamma_phi = c(tail(MCMCchains(out, 
                                             params = 'gamma_phi\\[1\\]', 
                                             ISB = FALSE,
                                             chain_num = CHAIN), n = 1),
                             tail(MCMCchains(out, 
                                             params = 'gamma_phi\\[2\\]', 
                                             ISB = FALSE,
                                             chain_num = CHAIN), n = 1)),
               mu_p = tail(MCMCchains(out, 
                                      params = 'mu_p',
                                      chain_num = CHAIN), n = 1),
               mu_phi = tail(MCMCchains(out, 
                                        params = 'mu_phi',
                                        chain_num = CHAIN), n = 1),
               pi_phi = tail(MCMCchains(out, 
                                        params = 'pi_phi',
                                        chain_num = CHAIN), n = 1),
               rho_phi = tail(MCMCchains(out, 
                                         params = 'rho_phi',
                                         chain_num = CHAIN), n = 1))


Inits1$rho_phi - Inits2$rho_phi


MCMCsummary(out)

