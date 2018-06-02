######################
#Script to plot simulated data time series
#
#0,1,2 chicks for specified number of nests in tile plot
#
######################


# load packages -----------------------------------------------------------


library(ggplot2)
library(reshape2)


# function to sim data ----------------------------------------------------



#function to simulate time series
sim_data_fun <- function(PHI_MAT, P_MAT, N_NESTS, TYPE)
{
  #PHI_MAT <- PHI
  #P_MAT <- P
  #N_NESTS <- nests
  
  TS_LEN <- NROW(PHI_MAT) + 1
  CH <- matrix(0, 
               ncol = N_NESTS, 
               nrow = TS_LEN)
  
  for (i in 1:N_NESTS)
  {
    #i <- 1
    
    if (TYPE == 'TWO')
    {
      CH[i,1] <- 2 
      t_SP <- c(2, rep(NA, TS_LEN-1))
    }
    if (TYPE == 'ONE')
    {
      CH[i,1] <- 1
      t_SP <- c(1, rep(NA, TS_LEN-1))
    }
    
    for (t in 2:TS_LEN)
    {
      #t <- 2
      #TRUE STATE
      t_SP[t] <- rbinom(1, size = t_SP[t-1], prob = PHI_MAT[t-1,i])
      
      #OBSERVED STATE
      t_DP <- rbinom(1, size = t_SP[t], prob = P_MAT[t-1, i])
      CH[t,i] <- t_DP
    }
  }
  return(CH)
}


# specify params ----------------------------------------------------------



nests <- 5
n_ts <- 25


P <- matrix(rep(0.7, times = nests*n_ts),
            ncol = nests)


PHI <- matrix(rep(0.97, (n_ts - 1) * nests),
              ncol = nests)



# simulate data -----------------------------------------------------------


sim_data <- sim_data_fun(PHI, P, nests, TYPE = 'TWO')
ns_data <- melt(sim_data)
colnames(ns_data) <- c('time', 'nest', 'chicks')



# plot data ---------------------------------------------------------------

ggplot(data = ns_data, aes(x = time, y = nest)) + 
  geom_tile(aes(fill = cut(chicks, breaks = 3, labels = 0:2)), 
            color = 'black', size = 0.5) +
  scale_fill_manual(values = colorRampPalette(c("white","red"))(5), 
                    na.value = "#EEEEEE", name = "Chicks") +
  theme_void()
