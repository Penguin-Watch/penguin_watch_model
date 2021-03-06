#################
# gifs of PW imagery and model results
#
# Author: Casey Youngflesh
#################


# Clear environment -------------------------------------------------------

rm(list = ls())


# dir ---------------------------------------------------------------------

dir <- '~/Google_Drive/R/penguin_watch_model/'
OUTPUT <- '~/Google_Drive/R/penguin_watch_model/Results/OUTPUT-2020-03-16'


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(dplyr)
library(ggplot2)
#requires installation of imagemagick - script written to work on macOS


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

#model input
setwd('~/Google_Drive/R/penguin_watch_model/Data/PW_data/')
PW_data <- read.csv('PW_data_2019-10-07.csv', stringsAsFactors = FALSE)


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


# gifs --------------------------------------------------------------------

#takes a few minutes per site - may need to be altered for sites other than BROW 2018 (because of different frequency of images, etc.)

#for (i in 1:12)
#{
i <- 1
#remove years that don't have data (wasn't done in model script)
t_date <- data$date_array[,,i]
to.rm.dt <- which(is.na(t_date[1,]))
if (length(to.rm.dt) > 0)
{
  t_date2 <- t_date[,-to.rm.dt]
}

#for (j in 1:4)
#{
j <- 2

# if (sum(!is.na(z_out_mn[,j,i])) > 0)
# {
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

#hourly time steps (only 8 hours)
times <- seq(from = st_obs, to = 61, by = 1/8)
ltimes <- length(times)
times2 <- times[-ltimes]

SITE <- data$unsites[i]
YEAR <- data$yrs_array[j,i]
min_z_out <- min(z_out_LCI[,j,i])
rng_z_out <- max(z_out_UCI[,j,i]) - min_z_out
rng_dt_seq <- (range(dt_seq)[2])
PLT_DF <- data.frame(time = 1:62,  
                     z_out_mn = c(rep(NA, (st_obs - 1)), 
                                  z_out_mn[st_obs:60,j,i], rep(NA, 2)),
                     z_out_LCI = c(rep(NA, (st_obs - 1)), 
                                   z_out_LCI[st_obs:60,j,i], rep(NA, 2)),
                     z_out_UCI = c(rep(NA, (st_obs - 1)), 
                                   z_out_UCI[st_obs:60,j,i], rep(NA, 2)),
                     z_out_dot = c(dt_seq, rep(NA, 60 - st_obs), rep(NA, 2)),
                     p_out_mn = c(p_out_mn[,j,i], rep(NA, 2)),
                     p_out_LCI = c(p_out_LCI[,j,i], rep(NA, 2)),
                     p_out_UCI = c(p_out_UCI[,j,i], rep(NA, 2)),
                     p_out_mn_sc = c(p_out_mn[,j,i] * rng_z_out + min_z_out, rep(NA, 2)),
                     p_out_LCI_sc = c(p_out_LCI[,j,i] * rng_z_out + min_z_out, rep(NA, 2)),
                     p_out_UCI_sc = c(p_out_UCI[,j,i] * rng_z_out + min_z_out, rep(NA, 2)),
                     count = c(counts, rep(NA, 2)))


#create gif dir if doesn't exist and change to that dir
giff_dir <- paste0(OUTPUT, '/gif/')

ifelse(!dir.exists(giff_dir), 
       dir.create(giff_dir), 
       FALSE)

#create dir within gif/
gdir <- paste0(giff_dir,
                   SITE, '-', YEAR, '-2019-10-07-plot')

#create dir for figs if doesn't exist and change to that dir
ifelse(!dir.exists(gdir), 
       dir.create(gdir), 
       FALSE)

setwd(gdir)

#navigate to directory with camera images
imgdir <- paste0('/Users/caseyyoungflesh/Google_Drive/Research/Projects/Penguin_watch/PW_surv_model_data/Full_res_images')
setwd(imgdir)

#names of all image directories
dirs <- list.files()

#grep for site and year
idx <- c(grep(SITE, dirs), grep(YEAR, dirs))

#directory for that specific site and year - navgate to that dir
sydir <- dirs[idx[duplicated(idx)]]
setwd(sydir)
files <- list.files()

#filter PW data by site and year
sydata <- dplyr::filter(PW_data, season_year == YEAR, site == SITE)

#convert to comprehensible dates
sub_dates <- substr(sydata$datetime, start = 1, stop = 10)
sub_times <- substr(sydata$datetime, start = 12, stop = 19)
conv_dates <- as.Date(sub_dates, format = '%Y:%m:%d')
conv_times <- as.numeric(substr(sub_times, start = 1, stop = 2))

#times for first day
fd_times <- as.numeric(substr(sub_times[which(conv_dates == conv_dates[1])], 
                              start = 1, stop = 2))

#starting image name
img_name_st <- sydata$imagename[1]

#just jpg files
files2 <- files[grep('.JPG', files)]

#find index first relevant image
img_idx <- grep(img_name_st, files2)

#jpg files that are relevant
files3 <- files2[img_idx:length(files2)]

#how many images remain for the first day
r_img <- 17 - min(fd_times)
#unique dates
ucd <- unique(conv_dates)

#vector of dates
full_dt <- c(rep(ucd[1], r_img), rep(ucd[-1], each = 8))

#trim off tail (images end before end of day)
tr_date <- full_dt[1:(length(files3) - 1)]
full_time <- c((conv_times[1] + 1):17, rep(10:17, times = 29))
tr_time <- full_time[1:(length(files3) - 1)]

#cbind filenames, dates, and times
img_td <- data.frame(filename = paste0(imgdir, '/', sydir, '/', files3), 
                     date = c(conv_dates[1], tr_date), 
                     time = c(conv_times[1], tr_time))


# ##############################
# #insert black images into day
# 
# #each day is 10:00 - 17:00
# hours <- 10:17
# night <- c(18:24, 1:9)
# #length day
# lh <- length(as.numeric(hours))
# #length night
# ln <- length(night)
# 
# #number of images before night on first day
# nidx <- lh - which(hours == min(fd_times)) + 1
# #files for day 1
# day_1 <- c(files3[1:nidx], rep('~/Google_Drive/R/penguin_watch_model/Data/black.JPG', ln))
# 
# #vector of images
# counter <- (nidx + 1)
# day_img <- day_1
# for (m in 2:30)
# {
#   day_img <- c(day_img, files3[counter:(counter + lh - 1)], 
#                rep('~/Google_Drive/R/penguin_watch_model/Data/black.JPG', ln))
#   counter <- counter + lh - 1
# }
# 
# #keep only images that are represented in plot
# final_imgs <- day_img[1:length(times)]
# ##############################



##############################
#no black images

#each day is 10:00 - 17:00 (8 of 24 hours)

#pad with black before starting and after ending
#number of black img before start
n_blk_img_st <- img_td[1,]$time - 10
tb <- data.frame(filename = 
                   rep('~/Google_Drive/R/penguin_watch_model/Data/black.JPG', 
                       n_blk_img_st), 
                 date = rep(img_td$date[1], n_blk_img_st),
                 time = 10:(10 + n_blk_img_st - 1))

img_td2 <- rbind(tb, img_td)

#number of black img after end
#total number of end images - current number of img - don't fill last day

n_blk_img_end <- (length(st_obs:60) * 8) - NROW(img_td2)

#number of hours left in day
eday <- 17 - tail(img_td2, n = 1)$time

#last day from previous df
ld <- tail(img_td2$date, n = 1)

tb2 <- data.frame(filename = 
                    rep('~/Google_Drive/R/penguin_watch_model/Data/black.JPG', 
                        n_blk_img_end), 
                  date = rep(ld, eday),
                  time = (17 - eday + 1):17)

img_td3 <- rbind(img_td2, tb2)

#remove last 4 images (second half of last day)
limg <- NROW(img_td3)
img_td4 <- img_td3[-c((limg - 4):limg),]
##############################

#8 hours per day
#create fig for each hourly time step
setwd(gdir)
for (d in 1:length(times2))
{
  p <- ggplot(PLT_DF, aes(x = time)) + 
    geom_ribbon(aes(ymin = z_out_LCI, ymax = z_out_UCI),
                fill = 'blue', alpha = 0.2) +
    geom_line(aes(y = z_out_mn), col = 'blue', size = 3) +
    geom_line(aes(y = z_out_dot), col = 'blue', linetype = 2, size = 3) +
    #THIS IS FOR THE DETECTION PROBABILITY
    # geom_ribbon(aes(ymin = p_out_LCI_sc, ymax = p_out_UCI_sc),
    #             fill = 'red', alpha = 0.2) +
    # geom_line(aes(y = p_out_mn_sc),
    #           col = 'red', size = 3) +
    #geom_line(aes(y = count),
    #          col = 'green') +
    geom_vline(xintercept = (times2[d] - 0.5), color = 'green', size = 2) +
    theme_bw() +
    # scale_y_continuous(sec.axis = sec_axis(~(.-min_z_out)/rng_z_out,
    #                                        name = 'Detection probability')) +
    scale_x_continuous(labels = n_dates, breaks = n_breaks) +
    ylab('Number of chicks') +
    xlab('') +
    #ggtitle(paste0(SITE, ' - ', YEAR)) + 
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 24),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title = element_text(size = 24),
      axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
      axis.ticks.length= unit(0.2, 'cm'))
  
  #convert d to string
  st_d <- toString(d)
  #number of characters for d
  nch_d <- nchar(st_d)
  #if less than three, pad with zeros
  if (nch_d == 2)
  {
    st_d <- paste0('0', st_d)
  }
  if (nch_d == 1)
  {
    st_d <- paste0('00', st_d)
  }
  
  #print(p)
  ggsave(p, filename = paste0('GIF-', st_d, '-', 
                              SITE, '-', YEAR, '-2019-10-07.jpg'),
         width = 8,
         height = 8)
}

#create dir for figs if doesn't exist and change to that dir
cdir <- paste0(OUTPUT, '/gif/', 
               SITE, '-', YEAR, '-2019-10-07-cam')

ifelse(!dir.exists(cdir), 
       dir.create(cdir), 
       FALSE)

#scale plots and cam images
#2048 x 1536 - 16%
W <- 2048 * 0.75
H <- 1536/2048 * W


#rescale plot images
system(paste0('for f in *.jpg; do convert $f -scale ', 
              W, 'X', H, ' ', cdir, '/$f-scale.jpg; done'))

setwd(cdir)
plt_img <- list.files()[grep('GIF', list.files())]

#rescale cam images
setwd(paste0(imgdir, '/', sydir))
system(paste0('for f in *.JPG; do convert $f -scale ', 
              W, 'X', H, ' ', cdir, '/$f-scale.jpg; done'))

#rescale black
system(paste0('convert ~/Google_Drive/R/penguin_watch_model/Data/black.JPG -scale ', W, 'X', H, ' ', cdir, '/black.JPG-scale.jpg'))

#cam img paths
cam_img_path <- paste0(as.character(img_td3$filename), '-scale.jpg')
#cam img names
cam_img <- sapply(strsplit(cam_img_path, split = '/'), tail, 1)

# #every other image to reduce size
# plt_img2 <- plt_img[c(TRUE, FALSE)]
# cam_img2 <- cam_img[c(TRUE, FALSE)]

#every image
plt_img2 <- plt_img
cam_img2 <- cam_img


#create copies of video caption
cap_dir <- '/Users/caseyyoungflesh/Google_Drive/R/penguin_watch_model/Appendices/'
#number of caption frames
ncap <- 175
setwd(cdir)

#need prefix of zeros
zn <- 1:ncap
d1 <- paste0('00', as.character(zn[which(nchar(zn) == 1)]))
d2 <- paste0('0', as.character(zn[which(nchar(zn) == 2)]))
d3 <- paste0(as.character(zn[which(nchar(zn) == 3)]))

dj <- c(d1, d2, d3)

for (i in 1:ncap)
{
  #i <- 1
  system(paste0('cp ', cap_dir, 'Video_S1_caption.jpg "g-', dj[i],'.jpg"')) 
}

#append images to one another
for (i in 1:length(plt_img2))
{
  #i <- 1
  i2 <- i + ncap
  
  if (nchar(i2) == 2)
  {
    n <- paste0('0', i2)
  }
  if (nchar(i2) == 3)
  {
    n <- paste0(i2)
  }
  
  system(paste0('convert ', cdir, '/', cam_img2[i], 
                ' ', cdir, '/', plt_img2[i],
                ' +append ', cdir, '/g-', i2, '.jpg'))
}


#create mp4
setwd(cdir)
system(paste0("ffmpeg -framerate 7 -pattern_type glob -i 'g-*.jpg' -c:v libx264 ", 
              sydir, "-cam-plot-combine.mp4"))

#create gif from mp4
setwd(cdir)
system(paste0("ffmpeg -i ", sydir, "-cam-plot-combine.mp4 -r 15 -vf scale=800:-1 -ss 00:00:25 -to 00:00:57 ", sydir, "-cam-plot-combine.gif"))

