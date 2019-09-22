#################
# Mortality timing plots
#
# Author: Casey Youngflesh
#################



# Clear environment -------------------------------------------------------

rm(list = ls())



# dir ---------------------------------------------------------------------

dir <- '~/Google_Drive/R/penguin_watch_model/'
OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2019-09-08'


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(dplyr)
library(ggplot2)



# load data ---------------------------------------------------------------

setwd(OUTPUT)

#read in latent true state (number of chicks) from model output
z_out <- readRDS('z_out.rds')

#read in detection prob from model output
p_out <- readRDS('p_out.rds')

#read in model input data
data <- readRDS('jagsData.rds')

#read in precipitation data
precip_df <- readRDS('precip_df.rds')



# process data ------------------------------------------------------------

#extract mean and 1 sd for latent state
z_out_mn <- MCMCvis::MCMCpstr(z_out, params = 'z_out', func = mean)[[1]]
z_out_sd <- MCMCvis::MCMCpstr(z_out, params = 'z_out', func = sd)[[1]]
z_out_LCI <- z_out_mn - z_out_sd
z_out_UCI <- z_out_mn + z_out_sd

p_out_mn <- MCMCvis::MCMCpstr(p_out, params = 'p_out', func = mean)[[1]]
p_out_sd <- MCMCvis::MCMCpstr(p_out, params = 'p_out', func = sd)[[1]]
p_out_LCI <- p_out_mn - p_out_sd
p_out_UCI <- p_out_mn + p_out_sd



# compare death rates at different period of season -----------------------

mort <- data.frame()
#site
for (i in 1:length(data$unsites))
{
  #i <- 1
  #year
  for (j in 1:length(data$yrs_array[,i]))
  {
    #j <- 1
    if (!is.na(data$yrs_array[j,i]))
    {
      #lay to hatch
      lh <- round((z_out_mn[1,j,i] - z_out_mn[30,j,i]) / 30, 3)
      
      #young hatch to older hatch
      yo <- round((z_out_mn[31,j,i] - z_out_mn[45,j,i]) / 15, 3)
      
      #older hatch to creche
      oc <- round((z_out_mn[46,j,i] - z_out_mn[60,j,i]) / 15, 3)
      
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

#palette from here: colormind.io
colors <- c('#225857', '#ECC151', '#E54E3D')
cols <- grDevices::col2rgb(colors)/255

#figure of mortality rates by egg/chick stage
pdf('mort_rates_season.pdf', width = 5, height = 5)
plot(density(mval[,1]), xlim = c(0, 0.7), ylim = c(0,6), 
     col = rgb(cols[1,3], cols[2,3], cols[3,3], 1), lwd = 3, 
     main = 'Chick death rates over nesting season', 
     xlab = 'Chick deaths per day')
lines(density(mval[,2]), 
      col = rgb(cols[1,2], cols[2,2], cols[3,2], 1), lwd = 3)
lines(density(mval[,3]), 
      col = rgb(cols[1,1], cols[2,1], cols[3,1], 1), lwd = 3)
legend('topright',
       legend = c('lay -> hatch', 'hatch -> young chick', 
                  'young chick -> creche'),
       col = c(rgb(cols[1,3], cols[2,3], cols[3,3], 1), 
               rgb(cols[1,2], cols[2,2], cols[3,2], 1), 
               rgb(cols[1,1], cols[2,1], cols[3,1], 1)),
       lwd = c(2,2,2,2), cex = 1)
dev.off()

#friedman test to check for differences in means (repeated measures - does not assume normality)
friedman.test(as.matrix(mval))




# time varying plots -----------------------------------------------------

#plot mean true latent state (total number of chicks in colony) as blue line (95% CI in blue ribbon)
#plot detection probability (prob of detecting chick at a given point in time)
#plot snow events in purple (2 or greater) and rain events in orange (2 or greater)
#width of vertical lines indicate severity of precipitation event

#create dir for figs if doesn't exist
ifelse(!dir.exists(paste0(OUTPUT, '/time_plots')), 
       dir.create(paste0(OUTPUT, '/time_plots')), 
       FALSE)

setwd(paste0(OUTPUT, '/time_plots'))

#sites
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
  #years
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
      
      #filter for daily precip events for site/year - only events 2+
      t_precip <- dplyr::filter(precip_df, 
                                site == data$unsites[i],
                                season_year == data$yrs_array[j,i],
                                m_snow >= 2 | s_rain >= 2,
                                date < max(dates),
                                date > min(dates))
      t_precip$ts <- as.numeric(t_precip[,'date'] - min(dates) + 1)
      
      
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
      
      #when do observations start being modeled at first nest for that site/year
      st_obs <- min(which(!is.na(data$y[,1,j,i])))
      
      #dotted line from start of season to first chick sighting
      dt_seq <- seq(z_out_mn[1,j,i], z_out_mn[st_obs,j,i], length = st_obs)
      
      
      SITE <- data$unsites[i]
      YEAR <- data$yrs_array[j,i]
      min_z_out <- min(z_out_LCI[,j,i])
      rng_z_out <- max(z_out_UCI[,j,i]) - min_z_out
      rng_dt_seq <- (range(dt_seq)[2])
      PLT_DF <- data.frame(time = 1:60,  
                           z_out_mn = c(rep(NA, (st_obs - 1)), 
                                        z_out_mn[st_obs:60,j,i]),
                           z_out_LCI = c(rep(NA, (st_obs - 1)), 
                                         z_out_LCI[st_obs:60,j,i]),
                           z_out_UCI = c(rep(NA, (st_obs - 1)), 
                                         z_out_UCI[st_obs:60,j,i]),
                           z_out_dot = c(dt_seq, rep(NA, 60 - st_obs)),
                           p_out_mn = p_out_mn[,j,i], 
                           p_out_LCI = p_out_LCI[,j,i],
                           p_out_UCI = p_out_UCI[,j,i],
                           p_out_mn_sc = p_out_mn[,j,i] * rng_z_out + min_z_out,
                           p_out_LCI_sc = p_out_LCI[,j,i] * rng_z_out + min_z_out,
                           p_out_UCI_sc = p_out_UCI[,j,i] * rng_z_out + min_z_out,
                           # p_out_mn_sc = p_out_mn[,j,i] * rng_dt_seq,
                           # p_out_LCI_sc = p_out_LCI[,j,i] * rng_dt_seq,
                           # p_out_UCI_sc = p_out_UCI[,j,i] * rng_dt_seq,
                           count = counts)
      
      p <- ggplot(PLT_DF, aes(x = time)) + 
        #THIS IS FOR THE LATENT TRUE STATE (number of chicks)
        geom_ribbon(aes(ymin = z_out_LCI, ymax = z_out_UCI),
                    fill = 'blue', alpha = 0.2) +
        geom_line(aes(y = z_out_mn), col = 'blue', size = 3) +
        geom_line(aes(y = z_out_dot), col = 'blue', linetype = 2, size = 3) +
        #THIS IS FOR THE DETECTION PROBABILITY
        geom_ribbon(aes(ymin = p_out_LCI_sc, ymax = p_out_UCI_sc),
                    fill = 'red', alpha = 0.2) +
        geom_line(aes(y = p_out_mn_sc),
                  col = 'red', size = 3) +
        # #Y-axis options
        # #second y-axis with number of chicks from 0 to max
        # scale_y_continuous(limits = c(-1, (max(dt_seq) + 1)),
        #                    sec.axis = sec_axis(~(.)/max(dt_seq),
        #                                        name = 'Detection probability')) +
        # #second y-axis with number of chicks from min to max
        scale_y_continuous(sec.axis = sec_axis(~(.-min_z_out)/rng_z_out,
                                               name = 'Detection probability')) +
        #one y-axis with number of chicks from 0 to max
        # scale_y_continuous(limits = c(0, (max(dt_seq) + 1))) +
        #comment out all y_continuous to plot min(chicks) : max(chicks)
        #THIS IS FOR PRECIPITATION
        geom_vline(xintercept = t_precip$ts,
                   color = 'purple', size = (t_precip$m_snow/3)*2,
                   alpha = 0.5) +
        geom_vline(xintercept = t_precip$ts,
                   color = 'orange', size = (t_precip$s_rain/3)*2,
                   alpha = 0.5) +
        #JUST 2018-01-13 BROW
        # geom_vline(xintercept = 48,
        #            color = 'red', size = 3,
        #            alpha = 0.9) +
        #THIS IS FOR THE NUMBER OF CONFIRMED CHICKS
        # geom_line(aes(y = count),
        #           col = 'green') +
        theme_bw() +
        scale_x_continuous(labels = n_dates, breaks = n_breaks) +
        ylab('Number of chicks') +
        xlab('') +
        ggtitle(paste0(SITE, ' - ', YEAR)) + 
        theme(
          plot.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title = element_text(size = 24),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          #ONLY FOR SECOND AXIS
          axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 15)),
          axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
          axis.ticks.length= unit(0.2, 'cm'))
      
      #print(p)
      # setwd('..')
      # ggsave(p, filename = 'Fig_5_b.pdf',
      #        width = 8,
      #        height = 8)
      
      ggsave(p, filename = paste0(SITE, '-', YEAR, '-2019-08-23.jpg'), 
             width = 8,
             height = 8)
    }
  }
}
