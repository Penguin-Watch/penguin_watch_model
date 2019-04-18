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

NAME <- 'PW_60k_2019-04-12_FULL_beta[j,k]' #nu_p[i,j,k] + beta_p[j,k]
#NAME <- 'PW_60k_2019-04-14_FULL_nu_p[j,k]_beta[i,j,k]' #nu_p[j,k] + beta_p[i,j,k]

setwd(paste0('~/Google_Drive/R/penguin_watch_model/Results/', NAME))

# fit1 <- readRDS('PW_20k_2019-04-04_LOCK.rds')
# fit2 <- readRDS('PW_50k_2019-04-04_LOCK2.rds')
# fit3 <- readRDS('PW_60k_2019-04-04_LOCKNEKO.rds')
# fit4 <- readRDS('PW_60k_2019-04-05_LOCKNEKO_cov.rds')
# fit5 <- readRDS('PW_60k_2019-04-05_LOCKNEKOGEOR_cov.rds')
# fit5 <- readRDS('PW_60k_2019-04-05_full_no_missing_ALL.rds')
# fit5 <- readRDS('PW_60k_2019-04-09_FULL_z_out_p_out.rds')

fit5 <- readRDS(paste0(NAME, '.rds'))
data5 <- readRDS('jagsData.rds')


# MCMCvis::MCMCtrace(fit5, params = 'z_out', Rhat = TRUE, n.eff = TRUE,
#                    filename = 'z_out_trace.pdf')
# MCMCvis::MCMCtrace(fit5, params = 'nu_p', Rhat = TRUE, n.eff = TRUE,
#                    filename = 'nu_p_trace.pdf')
# MCMCvis::MCMCtrace(fit5, params = 'beta_p', Rhat = TRUE, n.eff = TRUE,
#                    filename = 'beta_p_trace.pdf')
# MCMCvis::MCMCtrace(fit5, params = 'mu_phi', Rhat = TRUE, n.eff = TRUE,
#                    filename = 'mu_phi_trace.pdf')


# plot BS ------------------------------------------

# #mu_phi
# mu_phi <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi')[[1]]
# mu_phi_LCI <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi', 
#                                 func = function(x) quantile(x, probs = c(0.025)))[[1]]
# mu_phi_UCI <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi', 
#                                 func = function(x) quantile(x, probs = c(0.975)))[[1]]
# 
# #mu_phi_bs
# mu_phi_bs <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi_bs')[[1]]
# mu_phi_bs_LCI <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi_bs', 
#                                 func = function(x) quantile(x, probs = c(0.025)))[[1]]
# mu_phi_bs_UCI <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi_bs', 
#                                 func = function(x) quantile(x, probs = c(0.975)))[[1]]
# 
# 
# setwd('~/Google_Drive/R/penguin_watch_model/Data/')
# SLL <- read.csv('site_ll.csv')
# 
# setwd('Krill_data/CCAMLR/Processed_CCAMLR/')
# krill <- read.csv('CCAMLR_krill_entire_season.csv')
# 
# setwd('../../../SIC_data/Processed')
# sea_ice <- read.csv('SIC_150_W.csv')
# sea_ice2 <- sea_ice[,c(1,2,7)]
# 
# mrg <- data.frame(SITE = rep(data5$unsites, each = 4),
#                   YEAR = as.vector(data5$yrs_array[-5,]), 
#                   mn_mu_phi = as.vector(mu_phi_bs),
#                   LCI_mu_phi = as.vector(mu_phi_bs_LCI),
#                   UCI_mu_phi = as.vector(mu_phi_bs_UCI))
# 
# 
# mrg2 <- dplyr::left_join(mrg, krill, by = c('SITE', 'YEAR'))
# mrg3 <- dplyr::left_join(mrg2, sea_ice2, by = c('SITE', 'YEAR'))
# mrg4 <- dplyr::left_join(mrg3, SLL, by = c('SITE' = 'site_id'))
# 
# to.rm <- which(is.na(mrg4$YEAR))
# mrg5 <- mrg4[-to.rm,]
# 
# 
# library(ggplot2)
# ggplot(mrg5, aes(YEAR, mn_mu_phi, color = SITE)) +
#   geom_point(size = 3, position = position_dodge(width = 0.7)) +
#   geom_errorbar(data = mrg5, 
#                 aes(ymin = LCI_mu_phi, ymax = UCI_mu_phi, 
#                     color = SITE), width = 1.1, 
#                 position = position_dodge(width = 0.7)) +
#   theme_bw() +
#   ylim(c(0, 2)) +
#   ylab('Breeding Success') +
#   xlab('Year')





# plot regression fit SIC/KRILL -----------------------------------------------------

# plot(mrg5$T_KRILL, mrg5$mn_mu_phi)
# plot(mrg5$W_MN, mrg5$mn_mu_phi)
# plot(mrg5$latitude, mrg5$mn_mu_phi)
# 
# 
# alpha_ch <- MCMCvis::MCMCchains(fit5, params = 'alpha_theta')[,1]
# pi_ch <- MCMCvis::MCMCchains(fit5, params = 'pi_theta')[,1]
# rho_ch <- MCMCvis::MCMCchains(fit5, params = 'rho_theta')[,1]
# 
# sim_KRILL <- seq(min(data5$KRILL, na.rm = TRUE)-1, 
#                  max(data5$KRILL, na.rm = TRUE)+1, length = 100)
# sim_SIC <- seq(min(data5$SIC, na.rm = TRUE)-1, 
#                max(data5$SIC, na.rm = TRUE)+1, length = 100)
# 
# mf_KRILL <- matrix(nrow = length(alpha_ch), ncol = 100)
# mf_SIC <- matrix(nrow = length(alpha_ch), ncol = 100)
# for (i in 1:length(sim_KRILL))
# {
#   mf_KRILL[,i] <- alpha_ch + rho_ch * sim_KRILL[i] #+ pi_ch * mean(sim_SIC) #mean of cov is 0
#   mf_SIC[,i] <- alpha_ch + pi_ch * sim_SIC[i] #+ rho_ch * mean(sim_KRILL) #mean of cov is 0
# }
# 
# med_mf_KRILL <- apply(mf_KRILL, 2, median)
# LCI_mf_KRILL <- apply(mf_KRILL, 2, function(x) quantile(x, probs = 0.025))
# UCI_mf_KRILL <- apply(mf_KRILL, 2, function(x) quantile(x, probs = 0.975))
# 
# FIT_PLOT <- data.frame(MN_FIT = med_mf_KRILL, 
#                        MN_X = sim_KRILL,
#                        LCI_FIT = LCI_mf_KRILL,
#                        UCI_FIT = UCI_mf_KRILL)
# 
# DATA_PLOT <- data.frame(MN_phi = as.vector(mu_phi), 
#                        LCI_phi = as.vector(mu_phi_LCI),
#                        UCI_phi = as.vector(mu_phi_UCI),
#                        KRILL = as.vector(data5$KRILL),
#                        SIC = as.vector(data5$SIC))
# 
# ggplot(data = DATA_PLOT, aes(KRILL, MN_phi), color = 'black', alpha = 0.6) +
#   geom_ribbon(data = FIT_PLOT, 
#               aes(x = MN_X, ymin = LCI_FIT, ymax = UCI_FIT),
#               fill = 'grey', alpha = 0.7,
#               inherit.aes = FALSE) +
#   geom_line(data = FIT_PLOT, aes(MN_X, MN_FIT), color = 'red',
#             alpha = 0.9,
#             inherit.aes = FALSE,
#             size = 1.4) +
#   geom_errorbar(data = DATA_PLOT, 
#                 aes(ymin = LCI_phi, ymax = UCI_phi), width = 0.3,
#                 color = 'black', alpha = 0.2) +
#   geom_point(data = DATA_PLOT, aes(KRILL, MN_phi), color = 'black',
#              inherit.aes = FALSE, size = 3, alpha = 0.7) +
#   theme_bw() +
#   xlab('Krill catch') +
#   ylab('survival rate') +
#   theme(
#     axis.text = element_text(size = 16),
#     axis.title = element_text(size = 18),
#     axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
#     axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
#     axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
# 
# 
# med_mf_SIC <- apply(mf_SIC, 2, median)
# LCI_mf_SIC <- apply(mf_SIC, 2, function(x) quantile(x, probs = 0.025))
# UCI_mf_SIC <- apply(mf_SIC, 2, function(x) quantile(x, probs = 0.975))
# 
# FIT_PLOT <- data.frame(MN_FIT = med_mf_SIC, 
#                        MN_X = sim_SIC,
#                        LCI_FIT = LCI_mf_SIC,
#                        UCI_FIT = UCI_mf_SIC)
# 
# ggplot(data = DATA_PLOT, aes(SIC, MN_phi), color = 'black', alpha = 0.6) +
#   geom_ribbon(data = FIT_PLOT, 
#               aes(x = MN_X, ymin = LCI_FIT, ymax = UCI_FIT),
#               fill = 'grey', alpha = 0.7,
#               inherit.aes = FALSE) +
#   geom_line(data = FIT_PLOT, aes(MN_X, MN_FIT), color = 'red',
#             alpha = 0.9,
#             inherit.aes = FALSE,
#             size = 1.4) +
#   geom_errorbar(data = DATA_PLOT, 
#                 aes(ymin = LCI_phi, ymax = UCI_phi), width = 0.3,
#                 color = 'black', alpha = 0.2) +
#   geom_point(data = DATA_PLOT, aes(SIC, MN_phi), color = 'black',
#              inherit.aes = FALSE, size = 3, alpha = 0.7) +
#   theme_bw() +
#   xlab('SIC') +
#   ylab('survival rate') +
#   theme(
#     axis.text = element_text(size = 16),
#     axis.title = element_text(size = 18),
#     axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
#     axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
#     axis.ticks.length= unit(0.2, 'cm')) #length of axis tick



# plotly 3D plot ----------------------------------------------------------

# #from here: https://stackoverflow.com/questions/36049595/mixing-surface-and-scatterplot-in-a-single-3d-plot
# 
# mf_KR_SIC <- matrix(nrow =  length(sim_KRILL), ncol = length(sim_SIC))
# d_pts <- data.frame(MP = as.vector(mu_phi),
#                     KR = as.vector(data5$KRILL),
#                     S = as.vector(data5$SIC))
# for (i in 1:length(sim_KRILL))
# {
#   for (j in 1:length(sim_SIC))
#   {
#     mf_KR_SIC[i,j] <- mean(alpha_ch) + 
#       mean(rho_ch) * sim_KRILL[i] + 
#       mean(pi_ch) * sim_SIC[j]
#   }
# }
# 
# 
# library(plotly)
# plot_ly(x = sim_KRILL,
#         y = sim_SIC,
#         z = mf_KR_SIC, type = "surface") %>% 
#   add_trace(data = d_pts, x = d_pts$KR, y = d_pts$S, z = d_pts$MP, 
#             mode = "markers", type = "scatter3d", 
#             marker = list(size = 5, color = "red", symbol = 104))




# time varying plots -----------------------------------------------------

z_out_mn <- MCMCvis::MCMCpstr(fit5, params = 'z_out', func = mean)[[1]]
z_out_sd <- MCMCvis::MCMCpstr(fit5, params = 'z_out', func = sd)[[1]]
z_out_LCI <- z_out_mn - z_out_sd
z_out_UCI <- z_out_mn + z_out_sd

p_out_mn <- MCMCvis::MCMCpstr(fit5, params = 'p_out', func = mean)[[1]]
p_out_sd <- MCMCvis::MCMCpstr(fit5, params = 'p_out', func = sd)[[1]]
p_out_LCI <- p_out_mn - p_out_sd
p_out_UCI <- p_out_mn + p_out_sd

for (i in 1:12)
{
  #i <- 1
  #remove years that don't have data (wasn't done in model script)
  t_date <- data5$date_array[,,i]
  to.rm.dt <- which(is.na(t_date[1,]))
  if (length(to.rm.dt) > 0)
  {
    t_date2 <- t_date[,-to.rm.dt]
  }
  for (j in 1:4)
  {
    #j <- 1
    
    if (sum(!is.na(z_out_mn[,j,i])) > 0)
    {
      #if t_date2 doesn't have dims (all NA cols were removed)
      if (is.null(dim(t_date2)))
      {
        dates <- as.Date(t_date2, origin = '1970-01-01')
      } else {
        dates <- as.Date(t_date2[,j], origin = '1970-01-01')
      }
      
      #keep every other date value
      n_dates <- dates[seq(1, 60, by = 5)]
      #breaks for x axis on plot
      n_breaks <- seq(1, 60, by = 5)
      #min number of chicks (from counts in images and known counts later)
      counts <- data5$c_array[,j,i] 
      #remove 0 vals
      to.na <- which(counts == 0)
      if (length(to.na) > 0)
      {
        counts[to.na] <- NA
      }
      
      #dotted line from start of season to first chick
      dt_seq <- seq(z_out_mn[1,j,i], z_out_mn[29,j,i], length = 30)
      
      SITE <- data5$unsites[i]
      YEAR <- data5$yrs_array[j,i]
      min_z_out <- min(z_out_LCI[,j,i])
      rng_z_out <- max(z_out_UCI[,j,i]) - min_z_out
      PLT_DF <- data.frame(time = 1:60,  
                           z_out_mn = c(rep(NA, 29), z_out_mn[30:60,j,i]),
                           z_out_LCI = c(rep(NA, 29), z_out_LCI[30:60,j,i]),
                           z_out_UCI = c(rep(NA, 29), z_out_UCI[30:60,j,i]),
                           z_out_dot = c(dt_seq, rep(NA, 30)),
                           p_out_mn = p_out_mn[,j,i], 
                           p_out_LCI = p_out_LCI[,j,i],
                           p_out_UCI = p_out_UCI[,j,i],
                           p_out_mn_sc = p_out_mn[,j,i] * rng_z_out + min_z_out,
                           p_out_LCI_sc = p_out_LCI[,j,i] * rng_z_out + min_z_out,
                           p_out_UCI_sc = p_out_UCI[,j,i] * rng_z_out + min_z_out,
                           count = counts)
      
      p <- ggplot(PLT_DF, aes(x = time)) + 
        geom_ribbon(aes(ymin = z_out_LCI, ymax = z_out_UCI),
                  fill = 'blue', alpha = 0.2) +
        geom_line(aes(y = z_out_mn), col = 'blue') +
        geom_line(aes(y = z_out_dot), col = 'blue', linetype = 2) +
        geom_ribbon(aes(ymin = p_out_LCI_sc, ymax = p_out_UCI_sc),
                    fill = 'red', alpha = 0.2) +
        geom_line(aes(y = p_out_mn_sc),
                  col = 'red') +
        #geom_line(aes(y = count),
        #          col = 'green') +
        theme_bw() +
        scale_y_continuous(sec.axis = sec_axis(~(.-min_z_out)/rng_z_out, 
                                               name = 'Detection probability')) +
        scale_x_continuous(labels = n_dates, breaks = n_breaks) +
        ylab('Number of chicks') +
        xlab('') +
        ggtitle(paste0(SITE, ' - ', YEAR)) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      
      #print(p)
      ggsave(p, filename = paste0(SITE, '-', YEAR, '-2019-04-12.pdf'))
    }
  }
}




# PPO ---------------------------------------------------------------------


tf <- function(PR)
{
  hist(inv.logit(PR))
}

# mu_p ~ dnorm(2, 0.1)
PR <- rnorm(15000, 0, 1/sqrt(0.1))
tf(PR)
MCMCtrace(fit5, 
          params = 'mu_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)


# mu_beta_p ~ dnorm(0.1, 10) T(0, 0.5)
PR_p <- rnorm(15000, 0.1, 1/sqrt(10))
PR <- PR_p[which(PR_p > 0 & PR_p < 0.5)]
tf(PR)
MCMCtrace(fit5, 
          params = 'mu_beta_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

# sigma_beta_p ~ dunif(0, 2)
PR <- runif(15000, 0, 2)
MCMCtrace(fit5, 
          params = 'sigma_beta_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

# theta_phi ~ dnorm(4, 0.25)
PR <- rnorm(15000, 4, 1/sqrt(0.25))
tf(PR)
MCMCtrace(fit5, 
          params = 'theta_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

# sigma_mu_phi ~ dunif(0, 3)
PR <- runif(15000, 0, 3)
MCMCtrace(fit5, 
          params = 'sigma_mu_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

# sigma_nu_p ~ dunif(0, 3)
PR <- runif(15000, 0, 3)
MCMCtrace(fit5, 
          params = 'sigma_nu_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

