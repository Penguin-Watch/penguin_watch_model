#################
# Analyze model output
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
NAME <- 'PW_100k_2019-07-10_nu_p[i,j,k]_beta[j,k]'
OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2019-07-10'

setwd(paste0('~/Google_Drive/R/penguin_watch_model/Results/', NAME))

fit <- readRDS(paste0(NAME, '.rds'))
data <- readRDS('jagsData.rds')

setwd('~/Google_Drive/R/penguin_watch_model/Data/PW_data/')
PW_data <- read.csv('PW_data_2019-07-10.csv', stringsAsFactors = FALSE)



# summaries ---------------------------------------------------------------

sm <- MCMCvis::MCMCsummary(fit, excl = c('z_out', 'p_out'))
max(sm[,'Rhat'])
hist(sm[,'Rhat'])
min(sm[,'n.eff'])
hist(sm[,'n.eff'])



# trace plots -------------------------------------------------------------

#create dir for figs if doesn't exist
ifelse(!dir.exists(OUTPUT), 
       dir.create(OUTPUT), 
       FALSE)

setwd(OUTPUT)

# MCMCvis::MCMCtrace(fit, params = 'z_out', Rhat = TRUE, n.eff = TRUE,
#                    filename = 'z_out_trace.pdf')
MCMCvis::MCMCtrace(fit, params = 'nu_p', Rhat = TRUE, n.eff = TRUE,
                   ind = TRUE, filename = 'nu_p_trace.pdf')
MCMCvis::MCMCtrace(fit, params = 'beta_p', Rhat = TRUE, n.eff = TRUE,
                   ind = TRUE, filename = 'beta_p_trace.pdf')
MCMCvis::MCMCtrace(fit, params = 'mu_phi', Rhat = TRUE, n.eff = TRUE,
                   ind = TRUE, filename = 'mu_phi_trace.pdf')


# plot BS ------------------------------------------

#mu_phi
mu_phi <- MCMCvis::MCMCpstr(fit, params = 'mu_phi')[[1]]
mu_phi_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_phi',
                                func = function(x) quantile(x, probs = c(0.025)))[[1]]
mu_phi_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_phi',
                                func = function(x) quantile(x, probs = c(0.975)))[[1]]

#mu_phi_bs
mu_phi_bs <- MCMCvis::MCMCpstr(fit, params = 'mu_phi_bs')[[1]]
mu_phi_bs_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_phi_bs',
                                func = function(x) quantile(x, probs = c(0.025)))[[1]]
mu_phi_bs_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_phi_bs',
                                func = function(x) quantile(x, probs = c(0.975)))[[1]]



#USE d_mrg in data object
yrs_rng <- range(data$yrs_array, na.rm = TRUE)

site_vec <- c()
year_vec <- c()
for (i in 1:length(data$unsites))
{
  #i <- 1
  tv <- rep(data$unsites[i], sum(!is.na(data$yrs_array[,i])))
  yv <- data$yrs_array[which(!is.na(data$yrs_array[,i])),i]
  
  site_vec <- c(site_vec, tv)
  year_vec <- c(year_vec, yv)
}

mrg2 <- data.frame(SITE = site_vec,
                  YEAR = year_vec,
                  mn_mu_phi = as.vector(mu_phi_bs)[!is.na(as.vector(mu_phi_bs))],
                  LCI_mu_phi = as.vector(mu_phi_bs_LCI)[!is.na(as.vector(mu_phi_bs))],
                  UCI_mu_phi = as.vector(mu_phi_bs_UCI)[!is.na(as.vector(mu_phi_bs))])


library(ggplot2)
p <- ggplot(mrg2, aes(YEAR, mn_mu_phi, color = SITE)) + 
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbar(data = mrg2,
                aes(ymin = LCI_mu_phi, ymax = UCI_mu_phi,
                    color = SITE), width = 1.1,
                position = position_dodge(width = 0.7)) +
  theme_bw() +
  ylim(c(0, 2)) +
  ylab('Breeding Success') +
  xlab('Year')

ggsave('site_bs.pdf', p)



# plotly 3D plot ----------------------------------------------------------

# #from here: https://stackoverflow.com/questions/36049595/mixing-surface-and-scatterplot-in-a-single-3d-plot
# 
# mf_KR_SIC <- matrix(nrow =  length(sim_KRILL), ncol = length(sim_SIC))
# d_pts <- data.frame(MP = as.vector(mu_phi),
#                     KR = as.vector(data$KRILL),
#                     S = as.vector(data$SIC))
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


# precip events -----------------------------------------------------------

#aggregate precip events to day (sum rain and max snow)

setwd('~/Google_Drive/R/penguin_watch_model/Data/precip_data')

fls <- list.files()

precip_df <- data.frame()
for (i in 1:length(data$unsites))
{
  #i <- 1
  fls2 <- fls[grep(data$unsites[i], fls)]
  
  years <- as.numeric(substr(fls2, start = 6, stop = 9))
  for (j in 1:length(years))
  {
    #j <- 1
    
    #read in csv for that site/year
    tt <- read.csv(fls2[grep(years[j], fls2)])
    
    #remove time and transform to date
    dates <- as.Date(sapply(strsplit(as.character(tt$datetime), 
                                    split = ' '), '[', 1), 
                    format = '%Y:%m:%d')
    
    udates <- unique(dates)
    for (k in 1:length(udates))
    {
      #k <- 1
      t2 <- tt[which(dates == udates[k]),]
      
      #max snow score
      t_snow <- t2$score[grep('S', t2$score)]
      if (length(t_snow) > 0)
      {
        m_snow <- max(as.numeric(substr(t_snow, start = 2, stop = 2)))
      } else {
        m_snow <- 0
      }
      
      #rain
      t_rain <- t2$score[grep('R', t2$score)]
      if (length(t_rain) > 0)
      {
        s_rain <- sum(as.numeric(substr(t_rain, start = 2, stop = 2)))
      } else {
        s_rain <- 0
      }
      
      t_precip <- data.frame(site = data$unsites[i],
                             season_year = years[j],
                             date = udates[k],
                             m_snow,
                             s_rain)
      
      precip_df <- rbind(precip_df, t_precip)
    }
  }
}




# time varying plots -----------------------------------------------------

z_out_mn <- MCMCvis::MCMCpstr(fit, params = 'z_out', func = mean)[[1]]
z_out_sd <- MCMCvis::MCMCpstr(fit, params = 'z_out', func = sd)[[1]]
z_out_LCI <- z_out_mn - z_out_sd
z_out_UCI <- z_out_mn + z_out_sd

p_out_mn <- MCMCvis::MCMCpstr(fit, params = 'p_out', func = mean)[[1]]
p_out_sd <- MCMCvis::MCMCpstr(fit, params = 'p_out', func = sd)[[1]]
p_out_LCI <- p_out_mn - p_out_sd
p_out_UCI <- p_out_mn + p_out_sd

#plot snow events in purple (2 or greater) and rain events in orange (2 or greater)
#width of vertical lines indicate severity of precipitation event

#create dir for figs if doesn't exist
ifelse(!dir.exists(paste0(OUTPUT, '/time_plots')), 
       dir.create(paste0(OUTPUT, '/time_plots')), 
       FALSE)

setwd(paste0(OUTPUT, '/time_plots'))

for (i in 1:dim(data$date_array)[3])
{
  #i <- 1
  #remove years that don't have data (wasn't done in model script)
  t_date <- data$date_array[,,i]
  to.rm.dt <- which(is.na(t_date[1,]))
  if (length(to.rm.dt) > 0)
  {
    t_date2 <- t_date[,-to.rm.dt]
  }
  for (j in 1:dim(z_out_mn)[2])
  {
    #j <- 2
    
    if (sum(!is.na(z_out_mn[,j,i])) > 0)
    {
      #if t_date2 doesn't have dims (all NA cols were removed)
      if (is.null(dim(t_date2)))
      {
        dates <- as.Date(t_date2, origin = '1970-01-01')
      } else {
        dates <- as.Date(t_date2[,j], origin = '1970-01-01')
      }
      
      #filter for daily precip events for site/year
      t_precip <- dplyr::filter(precip_df, 
                                site == data$unsites[i],
                                season_year == data$yrs_array[j,i],
                                m_snow >= 2 | s_rain >= 2,
                                date < max(dates),
                                date > min(dates))
      t_precip$ts <- as.numeric(t_precip[,'date'] - min(dates))
      
      
      #keep every 5th date value
      n_dates <- dates[seq(1, 60, by = 5)]
      #breaks for x axis on plot
      n_breaks <- seq(1, 60, by = 5)
      #min number of chicks (from counts in images and known counts later)
      counts <- data$c_array[,j,i] 
      #remove 0 vals
      to.na <- which(counts == 0)
      if (length(to.na) > 0)
      {
        counts[to.na] <- NA
      }
      
      #NEEDS TO BE CHANGED AFTER NEXT MODEL RUN TO REFLECT DIFFERENT STARTING POINTS
      #dotted line from start of season to first chick
      dt_seq <- seq(z_out_mn[1,j,i], z_out_mn[31,j,i], length = 31)
      
      SITE <- data$unsites[i]
      YEAR <- data$yrs_array[j,i]
      min_z_out <- min(z_out_LCI[,j,i])
      rng_z_out <- max(z_out_UCI[,j,i]) - min_z_out
      rng_dt_seq <- (range(dt_seq)[2])
      PLT_DF <- data.frame(time = 1:60,  
                           z_out_mn = c(rep(NA, 30), z_out_mn[31:60,j,i]),
                           z_out_LCI = c(rep(NA, 30), z_out_LCI[31:60,j,i]),
                           z_out_UCI = c(rep(NA, 30), z_out_UCI[31:60,j,i]),
                           z_out_dot = c(dt_seq, rep(NA, 29)),
                           p_out_mn = p_out_mn[,j,i], 
                           p_out_LCI = p_out_LCI[,j,i],
                           p_out_UCI = p_out_UCI[,j,i],
                           # p_out_mn_sc = p_out_mn[,j,i] * rng_z_out + min_z_out,
                           # p_out_LCI_sc = p_out_LCI[,j,i] * rng_z_out + min_z_out,
                           # p_out_UCI_sc = p_out_UCI[,j,i] * rng_z_out + min_z_out,
                           p_out_mn_sc = p_out_mn[,j,i] * rng_dt_seq,
                           p_out_LCI_sc = p_out_LCI[,j,i] * rng_dt_seq,
                           p_out_UCI_sc = p_out_UCI[,j,i] * rng_dt_seq,
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
        #precipitation
        geom_vline(xintercept = t_precip$ts, 
                       color = 'purple', size = (t_precip$m_snow/3),
                   alpha = 0.5) +
        geom_vline(xintercept = t_precip$ts, 
                   color = 'orange', size = (t_precip$s_rain/3),
                   alpha = 0.5) +
        # geom_line(aes(y = count),
        #           col = 'green') +
        theme_bw() +
        scale_y_continuous(limits = c(-1, (max(dt_seq) + 1)),
                           sec.axis = sec_axis(~(.)/max(dt_seq),
                                               name = 'Detection probability')) +
        # scale_y_continuous(sec.axis = sec_axis(~(.-min_z_out)/rng_z_out,
        #                                        name = 'Detection probability')) +
        scale_x_continuous(labels = n_dates, breaks = n_breaks) +
        ylab('Number of chicks') +
        xlab('') +
        ggtitle(paste0(SITE, ' - ', YEAR)) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      
      #print(p)
      ggsave(p, filename = paste0(SITE, '-', YEAR, '-2019-07-10.pdf'))
    }
  }
}




# analyze -----------------------------------------------------------------

#to check where there might be problems - reference plots created above and look at time series
data$unsites
data$yrs_array

site <- 12
year <- 2
tdat <- data$y[-c(1:30), , year, site]
n_id <- which(apply(tdat, 2, function(x) sum(!is.na(x))) > 0)

tdat[, n_id]



# compare death rates at different period of season -----------------------

mort <- data.frame()
for (i in 1:length(data$unsites))
{
  #i <- 1
  for (j in 1:length(data$yrs_array[,i]))
  {
    #j <- 5
    if (!is.na(data$yrs_array[j,i]))
    {
      #lay to hatch
      lh <- round((z_out_mn[1,j,i] - z_out_mn[31,j,i]) / 31, 3)
      
      #young hatch to older hatch
      yo <- round((z_out_mn[31,j,i] - z_out_mn[45,j,i]) / 15, 3)
      
      #older hatch to creche
      oc <- round((z_out_mn[45,j,i] - z_out_mn[60,j,i]) / 15, 3)
      
      t_m <- data.frame(site = data$unsites[i], 
                        season_year = data$yrs_array[j,i],
                        lh, yo, oc)
      
      mort <- rbind(mort, t_m)
    }
  }
}

num <- cbind(rep(1, NROW(mort)), rep(2, NROW(mort)), rep(3, NROW(mort)))
mval <- mort[,c('lh', 'yo', 'oc')]




setwd(OUTPUT)

pdf('death_rates_season.pdf', width = 5, height = 5)
plot(density(mval[,1]), xlim = c(0, 0.6), ylim = c(0,5), 
     col = rgb(1,0,0,0.3), lwd = 3, 
     main = 'Chick death rates over nesting season', 
     xlab = 'Chick deaths per day')
lines(density(mval[,2]), col = rgb(0,1,0,0.3), lwd = 3)
lines(density(mval[,3]), col = rgb(0,0,1,0.3), lwd = 3)
legend('topright',
       legend = c('lay -> hatch', 'hatch -> young chick', 
                  'young chick -> creche'),
       col = c(rgb(1,0,0,0.3), rgb(0,1,0,0.3), rgb(0,0,1,0.3)),
       lwd = c(2,2,2,2), cex = 1)
dev.off()

#friedman test to check for differences in means (repeated measures - does not assume normality)
friedman.test(as.matrix(mval))
#post-hoc test
require(PMCMR)
posthoc.friedman.nemenyi.test(as.matrix(mval))




# merge data ----------------------------------------------------------

p2 <- data.frame()
for (i in 1:length(data$unsites))
{
  #i <- 1
  tsite <- dplyr::filter(precip_df, site == data$unsites[i])
  
  u_year <- unique(tsite$season_year)
  for (j in 1:length(u_year))
  {
    #j <- 3
    tyear <- dplyr::filter(tsite, season_year == u_year[j])
    tsnow <- length(which(tyear$m_snow >= 2))
    train <- length(which(tyear$s_rain >= 2))
    
    tt <- data.frame(SITE = data$unsites[i], 
                     YEAR = u_year[j],
                     tsnow,
                     train)
    p2 <- rbind(p2, tt)
  }
}

#merge lat/ln with precip, BS, and creche date
cll <- unique(PW_data[,c('site', 'col_lat', 'col_lon')])
mrg3 <- dplyr::left_join(mrg2, cll, by = c('SITE' = 'site'))
mrg4 <- dplyr::left_join(mrg3, p2, by = c('SITE', 'YEAR'))

mrg5 <- dplyr::left_join(mrg4, data$d_mrg, by = c('SITE' = 'site', 'YEAR' = 'season_year'))


#using lat/lon from Penguin Watch for PETE
#Merge with Hinke et al. 2017 (MEE) data
hinke_2017 <- data.frame(SITE = c('CAPE', 'CIER', 'COPA', 'GALI', 'LION', 'PETE'), 
                         YEAR = rep(2017, 6),
                         mn_mu_phi = c(1.63, 1.47, 1.53, 1.46, 1.26, 1.51),
                         LCI_mu_phi = rep(NA, 6),
                         UCI_mu_phi = rep(NA, 6),
                         col_lat = c(-62.46, -64.143, -62.175, 
                                     -65.244, -62.135, -65.17),
                         col_lon = c(-60.789, -60.984, -58.456, 
                                     -64.247, -58.126, -64.14),
                         tsnow = rep(NA, 6),
                         train = rep(NA, 6),
                         j = rep(NA, 6),
                         k = rep(NA, 6),
                         num_nests = c(8, 15, 58, 28, 19, 37),
                         chick_date = rep(NA, 6),
                         creche_date = c('2017-01-18', '2017-01-20', '2016-12-22',
                                         '2017-02-04', '2016-12-20', '2017-02-08'),
                         days = rep(NA, 6))

#merge with Lynch et al. 2009 (Polar Bio) data
lynch_2009 <- data.frame(SITE = 'PETE', 
                         YEAR = c('2004', '2005', '2006', '2007', '2008'),
                         mn_mu_phi = c((3260/2145), (2781/2265), (3453/2438),
                                       (3343/2293), (3348/2719)),
                         LCI_mu_phi = rep(NA, 5),
                         UCI_mu_phi = rep(NA, 5),
                         col_lat = rep(-65.17, 5),
                         col_lon = rep(-64.14, 5),
                         tsnow = rep(NA, 5),
                         train = rep(NA, 5),
                         j = rep(NA, 5),
                         k = rep(NA, 5),
                         num_nests = c(2145, 2265, 2438, 2293, 2719),
                         chick_date = rep(NA, 5),
                         creche_date = rep(NA, 5),
                         days = rep(NA, 5))

#From Lynch et al. 2009 (Polar Bio)
#LOCK - 1.24 - 1.39 - Cobley and Shears 1999
#BIRD ISLAND - 0 - 1.2 - Croxall and Prince 1979
#BIRD ISLAND - 0.9 - 1.02 - Williams 1990
#MACQUARIE - 0.93 +- 0.45 - Holmes et al. 2006
#MACQUARIE - 0.36 - 1.14 - Reilly and Kerle 1981
#MACQUARIE - 0 - 1.52 - Robertson 1986

mrg6 <- rbind(mrg5, hinke_2017, lynch_2009)

mrg6$SOURCE <- c(rep('PW', length(which(!is.na(mrg6$train)))), 
                 rep('Hinke', length(which(is.na(mrg6$train) & !is.na(mrg6$creche_date)))),
                 rep('Lynch', length(which(is.na(mrg6$creche_date)))))


#save  object
setwd(OUTPUT)
saveRDS(mrg6, 'mrg6.rds')





# PPO ---------------------------------------------------------------------

tf <- function(PR)
{
  hist(inv.logit(PR))
}

# mu_p ~ dnorm(2, 0.1)
PR <- rnorm(15000, 0, 1/sqrt(0.1))
tf(PR)
MCMCtrace(fit, 
          params = 'mu_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)


# mu_beta_p ~ dnorm(0.1, 10) T(0, 0.5)
PR_p <- rnorm(15000, 0.1, 1/sqrt(10))
PR <- PR_p[which(PR_p > 0 & PR_p < 0.5)]
tf(PR)
MCMCtrace(fit, 
          params = 'mu_beta_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

# sigma_beta_p ~ dunif(0, 2)
PR <- runif(15000, 0, 2)
MCMCtrace(fit, 
          params = 'sigma_beta_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

# theta_phi ~ dnorm(4, 0.25)
PR <- rnorm(15000, 4, 1/sqrt(0.25))
tf(PR)
MCMCtrace(fit, 
          params = 'theta_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

# sigma_mu_phi ~ dunif(0, 3)
PR <- runif(15000, 0, 3)
MCMCtrace(fit, 
          params = 'sigma_mu_phi',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

# sigma_nu_p ~ dunif(0, 3)
PR <- runif(15000, 0, 3)
MCMCtrace(fit, 
          params = 'sigma_nu_p',
          ind = TRUE, 
          priors = PR,
          pdf = FALSE,
          post_zm = FALSE)

