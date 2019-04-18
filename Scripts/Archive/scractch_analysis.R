#SCRATCH



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
