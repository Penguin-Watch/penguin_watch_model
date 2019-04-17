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

#NAME <- 'PW_60k_2019-04-12_FULL_beta[j,k]' #nu_p[i,j,k] + beta_p[j,k]
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



# plot BS ------------------------------------------

#mu_phi
mu_phi <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi')[[1]]
mu_phi_LCI <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi', 
                                func = function(x) quantile(x, probs = c(0.025)))[[1]]
mu_phi_UCI <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi', 
                                func = function(x) quantile(x, probs = c(0.975)))[[1]]

#mu_phi_bs
mu_phi_bs <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi_bs')[[1]]
mu_phi_bs_LCI <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi_bs', 
                                func = function(x) quantile(x, probs = c(0.025)))[[1]]
mu_phi_bs_UCI <- MCMCvis::MCMCpstr(fit5, params = 'mu_phi_bs', 
                                func = function(x) quantile(x, probs = c(0.975)))[[1]]


setwd('~/Google_Drive/R/penguin_watch_model/Data/')
SLL <- read.csv('site_ll.csv')

setwd('Krill_data/CCAMLR/Processed_CCAMLR/')
krill <- read.csv('CCAMLR_krill_entire_season.csv')

setwd('../../../SIC_data/Processed')
sea_ice <- read.csv('SIC_150_W.csv')
sea_ice2 <- sea_ice[,c(1,2,7)]

mrg <- data.frame(SITE = rep(data5$unsites, each = 4),
                  YEAR = as.vector(data5$yrs_array[-5,]), 
                  mn_mu_phi = as.vector(mu_phi_bs),
                  LCI_mu_phi = as.vector(mu_phi_bs_LCI),
                  UCI_mu_phi = as.vector(mu_phi_bs_UCI))


mrg2 <- dplyr::left_join(mrg, krill, by = c('SITE', 'YEAR'))
mrg3 <- dplyr::left_join(mrg2, sea_ice2, by = c('SITE', 'YEAR'))
mrg4 <- dplyr::left_join(mrg3, SLL, by = c('SITE' = 'site_id'))

to.rm <- which(is.na(mrg4$YEAR))
mrg5 <- mrg4[-to.rm,]


library(ggplot2)
ggplot(mrg5, aes(YEAR, mn_mu_phi, color = SITE)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = mrg5, 
                aes(ymin = LCI_mu_phi, ymax = UCI_mu_phi, 
                    color = SITE), width = 1.1, 
                position = position_dodge(width = 0.7)) +
  theme_bw() +
  ylim(c(0, 2)) +
  ylab('Breeding Success') +
  xlab('Year')





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
      if (length(to.rm) > 0)
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
      ggsave(p, filename = paste0(SITE, '-', YEAR, '-2019-04-14.pdf'))
    }
  }
}




# posterior estimates -----------------------------------------------------


# setwd('~/Google_Drive/R/penguin_watch_model/Figures/')
# 
# pdf(file = 'eta_phi.pdf', width = 7, height = 6, useDingbats = FALSE)
# MCMCplot(out, params = 'eta_phi',
#          horiz = FALSE,
#          ylim = c(-4,4),
#          labels = NULL,
#          thick_sz = 8,
#          thin_sz = 4,
#          med_sz = 2.5,
#          main = 'ETA')
# dev.off()
# 
# pdf(file = 'gamma_phi.pdf', width = 7, height = 6, useDingbats = FALSE)
# MCMCplot(out, params = 'gamma_phi',
#          horiz = FALSE,
#          ylim = c(-4,4),
#          labels = NULL,
#          thick_sz = 8,
#          thin_sz = 4,
#          med_sz = 2.5,
#          main = 'GAMMA')
# dev.off()
# 
# pdf(file = 'pi_rho_phi.pdf', width = 7, height = 6, useDingbats = FALSE)
# MCMCplot(out, params = c('pi_phi', 'rho_phi'), 
#          horiz = FALSE, 
#          ylim = c(-4,4), 
#          labels = NULL,
#          thick_sz = 8,
#          thin_sz = 4,
#          med_sz = 2.5,
#          main = 'PI/RHO')
# dev.off()


pdf(file = 'param_estimates.pdf', width = 4, height = 6, useDingbats = FALSE)
MCMCplot(out, params = c('gamma_phi', 'eta_phi', 'pi_phi', 'rho_phi'),
         #horiz = FALSE,
         ylim = c(-4,4),
         labels = NULL,
         thick_sz = 8,
         thin_sz = 4,
         med_sz = 2.5,
         main = 'param estimates')
dev.off()


library(reshape2)
library(ggplot2)

tt <- matrix(c(0,1,0,0,1,0,1,1,1,1,1,0), nrow = 3)
ns_data <- melt(tt)
colnames(ns_data) <- c('year', 'site', 'data')

ggplot(data = ns_data, aes(x = year, y = site)) + 
  geom_tile(aes(fill = cut(data, breaks = 2, labels = 0:1)), 
            color = 'black', size = 0.5) +
  scale_fill_manual(values = colorRampPalette(c("grey90","light green"))(2), 
                    na.value = "#EEEEEE", name = "Data") +
  theme_void()


t4 <- rbinom(4, 2, 0.6)
CI <- 1.96*sqrt(2*(mean(t4)/2)*(1 - (mean(t4)/2))/4)
c((mean(t4)/2) - CI, (mean(t4)/2) + CI)

t10 <- rbinom(10, 2, 0.6)
(CI <- 1.96*sqrt(2*(mean(t10)/2)*(1 - (mean(t10)/2))/10))
c((mean(t10)/2) - CI, (mean(t10)/2) + CI)

t25 <- rbinom(25, 2, 0.6)
(CI <- 1.96*sqrt(2*(mean(t25)/2)*(1 - (mean(t25)/2))/25))
c((mean(t25)/2) - CI, (mean(t25)/2) + CI)

t50 <- rbinom(50, 2, 0.6)
(CI <- 1.96*sqrt(2*(mean(t50)/2)*(1 - (mean(t50)/2))/50))
c((mean(t50)/2) - CI, (mean(t50)/2) + CI)

t75 <- rbinom(75, 2, 0.6)
(CI <- 1.96*sqrt(2*(mean(t75)/2)*(1 - (mean(t75)/2))/75))
c((mean(t75)/2) - CI, (mean(t75)/2) + CI)

t100 <- rbinom(100, 2, 0.6)
(CI <- 1.96*sqrt(2*(mean(t100)/2)*(1 - (mean(t100)/2))/100))
c((mean(t100)/2) - CI, (mean(t100)/2) + CI)

t500 <- rbinom(500, 2, 0.6)
(CI <- 1.96*sqrt(2*(mean(t500)/2)*(1 - (mean(t500)/2))/500))
c((mean(t500)/2) - CI, (mean(t500)/2) + CI)

t1000 <- rbinom(1000, 2, 0.6)
(CI <- 1.96*sqrt(2*(mean(t1000)/2)*(1 - (mean(t1000)/2))/1000))
c((mean(t1000)/2) - CI, (mean(t1000)/2) + CI)


data$y
data$NI


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

xvals <- scale(1:60, scale = FALSE)[,1]

# Error in node y[21,2,2,1] - 1
#DATA$y[1:21,2,2,1]
#DATA$w[1:21,2,2,1]
#DATA$z[1:21,2,2,1]
#phi
plot(inv.logit(5 + 0 + 0 + 0.005*xvals + 0 + 0))
#p
plot(inv.logit(0 + 0.1*xvals))





# PPO ---------------------------------------------------------------------

tf <- function(PR)
{
  hist(inv.logit(PR))
}


#' mu_phi ~ dnorm(3, 0.1)   
PR <- rnorm(15000, 3, 1/sqrt(0.1))
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


#' mu_p ~ dnorm(2, 0.1)
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

#grand intercept - survival
ch_mu_phi <- MCMCchains(out, params = 'mu_phi')[[1]]
#site effect
ch_eta_phi <- MCMCpstr(out, params = 'eta_phi', type = 'chains')[[1]]
#year effect
ch_gamma_phi <- MCMCpstr(out, params = 'gamma_phi', type = 'chains')[[1]]


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
for (k in 1:data$NK)
{
#k <- 2
  for (j in 1:data$NJ)
  {
    #j <- 1
    surv <- inv.logit(mu_phi + eta_phi[k] + gamma_phi[j])
  }
}


#survival - uncertainty
for (k in 1:data$NK)
{
  #k <- 2
  for (j in 1:data$NJ)
  {
    #j <- 1
    surv <- inv.logit(ch_mu_phi + ch_eta_phi[k,] + ch_gamma_phi[j])
    hist(surv, main = paste0('k=', k, ' j=', j), xlim = c(0.9, 1))
  }
}


#covariates
pi <- MCMCchains(out, params = 'pi_phi')
hist(exp(pi))
rho <- MCMCchains(out, params = 'rho_phi')
hist(exp(rho))


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

