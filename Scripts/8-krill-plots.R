#################
# Analyze krill
#
# Author: Casey Youngflesh
#################


# Clear environment -------------------------------------------------------

rm(list = ls())


# dir ---------------------------------------------------------------------

dir <- '~/Google_Drive/R/penguin_watch_model/'
OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2020-03-16'


# Load packages -----------------------------------------------------------

library(ggplot2)


# Load data -------------------------------------------------------

setwd(paste0(dir, 'Data/Krill_data/CCAMLR/Processed_CCAMLR'))

#read in weighted krill catch
krill_weighted_150 <- readRDS('krill_weighted_150.rds')

#read in average (across years) krill catch over entire season
CCAMLR_kr_AY <- read.csv('CCAMLR_krill_average.csv')


# krill figs --------------------------------------------------------------

# which sites have the most fishing at them?

#CCAMLR
p <- ggplot(CCAMLR_kr_AY, aes(SITE, T_KRILL)) +
  geom_col() +
  ylab('Mean krill catch (tonnes)') +
  theme_bw() +
  ggtitle('CCAMLR - Mean krill catch March - Jan (2000-2018)')

setwd(OUTPUT)
ggsave(p, filename = 'krill_catch_by_site.jpg')


# When is krill fishing most intense in this region? ----------------------

#CCAMLR
krill_time <- data.frame(MONTH = as.integer(1:12), KRILL = rep(NA, 12))
for (i in 1:12)
{
  #i <- 1
  temp <- dplyr::filter(krill_weighted_150, MONTH == i)
  swk <- sum(temp$WEIGHTED_KRILL, na.rm = TRUE)
  krill_time[i,2] <- swk/1000
}

months <- c('Jan', 'Mar', 'May', 'Jul',
            'Sep', 'Nov')

p2 <- ggplot(krill_time, aes(x = MONTH, y = KRILL)) +
  geom_col() + 
  theme_bw() +
  xlab('Month') +
  ylab('Krill catch (thousands of tons)') +
  scale_x_continuous(breaks = seq(1,12, by = 2), labels = months) +
  theme(
    plot.title = element_text(size = 22),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
    axis.ticks.length= unit(0.2, 'cm')) + #length of axis tick
  ggtitle('')
  #ggtitle('CCAMLR - Krill catch by month')

setwd(OUTPUT)
ggsave(p2, filename = 'krill_catch_by_month.jpg')
