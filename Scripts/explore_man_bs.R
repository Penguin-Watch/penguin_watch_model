###############
#Analyze derived data
#
###############



dir <- setwd('~/Google_Drive/R/')

setwd(paste0(dir, '/penguin_watch_model/Data'))


man_bs <- read.csv('man_bs.csv')



# merge data with covariates ----------------------------------------------


setwd('Krill_data/CCAMLR/Processed_CCAMLR/')
krill <- read.csv('CCAMLR_krill_entire_season.csv')

setwd(paste0(dir, '/penguin_watch_model/Data/SIC_data/Processed'))
sic1_p <- read.csv('SIC_150_W.csv')
sic1 <- dplyr::filter(sic1_p, YEAR > 1995)

sic2_p <- read.csv('SIC_500_MAX.csv')
sic2 <- dplyr::filter(sic2_p, YEAR > 1995)

sic3_p <- read.csv('SIC_500_W.csv')
sic3 <- dplyr::filter(sic3_p, YEAR > 1995)
colnames(sic3)[7] <- 'W_500'

sic <- data.frame(sic1[,c('YEAR', 'SITE', 'W_MN')], 
                  MAX_SIC = sic2[,'MAX_SIC'], 
                  W_500 = sic3[,'W_500'])

md <- dplyr::left_join(man_bs, krill, by = c('SITE', 'YEAR'))
md2 <- dplyr::left_join(md, sic, by = c('SITE', 'YEAR'))




# plot data ---------------------------------------------------------------

plot(man_bs$LAT, man_bs$BS)
plot(md2$T_KRILL, md2$BS)
plot(md2$W_MN, md2$BS)
plot(md2$MAX_SIC, md2$BS)
plot(md2$W_500, md2$BS)

summary(lm(BS ~ T_KRILL, data = md2))
summary(lm(BS ~ W_MN, data = md2))
summary(lm(BS ~ MAX_SIC, data = md2))
summary(lm(BS ~ W_500, data = md2))

summary(lm(BS ~ T_KRILL + W_MN, data = md2))




# fit GAMM ----------------------------------------------------------------

library(rstanarm)
ITER <- 2000
CHAINS <- 3
DELTA <- 0.98

to.rm <- which(is.na(md2$T_KRILL))
md3 <- md2[-to.rm,]

fit <- rstanarm::stan_gamm4(BS ~ s(T_KRILL) + s(W_MN), 
                       data = md3,
                       iter = ITER,
                       chains = CHAINS,
                       cores = CHAINS,
                       adapt_delta = DELTA)
summary(fit)




# SEM ---------------------------------------------------------------------

# 14 LOCK 2013 1.6250000 -64.82     8
# 15 LOCK 2014 1.4000000 -64.82    15
# 16 LOCK 2015 1.6153846 -64.82    13


# t_phi[1,1]     5.5412 0.4513  4.7795  5.4977  6.5322 1.00  9000
# t_phi[2,1]     5.1822 0.3088  4.5949  5.1729  5.8174 1.00  8346
# t_phi[3,1]     5.5178 0.3805  4.8416  5.4928  6.3416 1.00  9000
# t_phi[1,2]     6.1074 0.4610  5.3276  6.0666  7.1265 1.00  8842
# t_phi[2,2]     5.1099 0.2543  4.6274  5.1029  5.6308 1.00  9234
# t_phi[3,2]     4.7428 0.2390  4.2974  4.7361  5.2285 1.00  8520

ch1 <- c(2,2,0,2,1,2,2,2)
ch2 <- c(1,1,1,2,2,2,2,1,1,1,1,1,2,1,2,1)
mean(ch1)
sd(ch1)/sqrt(8)
mean(ch2)
sd(ch1)/sqrt(15)

#derive breeding success = inverse logit of mu_phi (or equivalent), raise to t power, multiply by number of initial chicks
boot::inv.logit(4.78)^60*2
boot::inv.logit(5.50)^60*2
boot::inv.logit(6.53)^60*2

boot::inv.logit(5.33)^60*2
boot::inv.logit(6.0)^60*2
boot::inv.logit(7.1)^60*2
