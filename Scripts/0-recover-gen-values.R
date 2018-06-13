#################
# Penguin Watch Model - 0 - Recovering gen values using one detection param vs many
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

#determine whether one detection params for all nests or individual detection params recovers generating survival value better
#Run time < 2.5 hours

#RESULTS: Using 150 time steps and 'runif(1, 0.99, 0.9999)' for survival, one detection and ind detection params obth recover generating values



# Clear environment -------------------------------------------------------

rm(list = ls())


# Load packages -----------------------------------------------------------

if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman')
}

pacman::p_load(rjags, MCMCvis, parallel, jagsRun, boot)




# colony parameters -------------------------------------------------------

#simulate new data - script modified from Kerry and Schaub 2012
n_ts <- 150 #number of time steps
x <- 1:n_ts
nests <- 15 #number of nests



#run through J different parameter sets, generating J different datasets per parameter set
J <- 3
I <- 5

#progress bar
total <- I * J
pb <- txtProgressBar(min = 0, max = total, style = 3)
PB <- 0

#loop through
OUT <- data.frame()
for (j in 1:J)
{
  #j <- 1

  set.seed(j)
  #generating survival probabilities
  dphi <- runif(1, 0.99, 0.9999)
  #generating detection probabilities
  dp <- runif(nests, 0.3, 0.8)
  
  for (i in 1:I)
  {
    # survival ----------------------------------------------------------------
    
    #survival probabilities
    PHI <- matrix(rep(dphi, (n_ts - 1) * nests),
                  ncol = nests)
    
    
    # detection ---------------------------------------------------------------
    
    #setection probabilities
    P <- matrix(rep(dp, each = n_ts - 1),
                ncol = nests)
    
    
    # simulate time series ----------------------------------------------------
    
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
    
    sim_data <- sim_data_fun(PHI, P, nests, TYPE = 'ONE')
    
    
    # Create z-state matrix ---------------------------------------------------
    
    #fun modified from Kerry and Schaub 2012
    #assume that chicks are alive since time step 1
    known.state.fun <- function(INPUT, TYPE)
    {
      #INPUT <- DATA$y
      state <- INPUT
      for (i in 1:NCOL(INPUT))
      {
        n1 <- 1
        
        if (TYPE == 'TWO')
        {
          if (sum(state[,i] == 2) > 0)
          {
            n2 <- max(which(INPUT[,i] == 2))
            state[n1:n2, i] <- 2
          }
          
          state[n1, i] <- NA
        }
        if (TYPE == 'ONE')
        {
          if (sum(state[,i] == 1) > 0)
          {
            n2 <- max(which(INPUT[,i] == 1))
            state[n1:n2, i] <- 1
          }
          
          state[n1, i] <- NA
        }
      }
      if (TYPE == 'TWO')
      {
        state[state == 0] <- NA
        state[state == 1] <- NA
      }
      if (TYPE == 'ONE')
      {
        state[state == 0] <- NA
      }
      return(state)
    }
    
    z_matrix <- known.state.fun(sim_data, TYPE = 'ONE')
    
    
    # construct data ----------------------------------------------------------
    
    DATA <- list(
      y = sim_data, #reponse
      NI = dim(sim_data)[2], #number of nests
      NT = dim(sim_data)[1], #number of time steps
      z = z_matrix) #known points of bird being alive
    
    
    # Model -------------------------------------------------------------------
    
    {
      sink("pwatch_sim.jags")
      
      cat("
          
          model {
          
          #nests
          for (i in 1:NI)
          {
          #chick alive at time step 1
          z[1,i] <- 1
          
          #time step
          for (t in 2:NT)
          {
          #state model
          z[t,i] ~ dbinom(p_alive[t,i], z[t-1,i])
          p_alive[t,i] <- phi[t,i] * z[t-1,i]
          
          #observation model
          y[t,i] ~ dbinom(p_sight[t,i], z[t,i])
          p_sight[t,i] <- p[t,i] * z[t,i]
          
          }
          }
          
          
          #transforms
          #nests
          for (i in 1:NI)
          {
          #time
          for (t in 1:NT)
          {
          logit(phi[t,i]) <- mu_phi
          logit(p[t,i]) <- mu_p + nu_p[i]
          } #t
          } #i
          
          #priors - phi
          mu_phi ~ dnorm(0, 0.01)
          
          #priors - p
          mu_p ~ dnorm(0, 0.5)
          
          for (i in 1:NI)
          {
          nu_p[i] ~ dnorm(0, tau_nu_p)
          }
          
          tau_nu_p ~ dunif(0, 10) 
          
          il_mu_phi <- ilogit(mu_phi)
  
  
          }",fill = TRUE)
  
      sink()
    }
    
    
    # Starting values ---------------------------------------------------------
    
    
    Inits_1 <- list(mu_phi = 5,
                    mu_p = 0,
                    tau_nu_p = 1,
                    .RNG.name = "base::Mersenne-Twister",
                    .RNG.seed = 1)
    
    Inits_2 <- list(mu_phi = 5,
                    mu_p = 0,
                    tau_nu_p = 1,
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = 2)
    
    Inits_3 <- list(mu_phi = 5,
                    mu_p = 0,
                    tau_nu_p = 1,
                    .RNG.name = "base::Marsaglia-Multicarry",
                    .RNG.seed = 3)
    
    F_Inits <- list(Inits_1, Inits_2, Inits_3)
    
    
    
    # Parameters to track -----------------------------------------------------
    
    Pars <- c('mu_phi',
              'il_mu_phi',
              'mu_p',
              'nu_p',
              'tau_nu_p')
    
    
    # Run model ---------------------------------------------------------------
    
    
    tt <- proc.time()[3]
    OUT_ind <- jagsRun(jagsData = DATA, 
            jagsModel = 'pwatch_sim.jags',
            jagsInits = F_Inits,
            params = Pars,
            n_chain = 3,
            n_adapt = 5000,
            n_burn = 30000,
            n_draw = 20000,
            report = FALSE)
    time_ind <- round((proc.time()[3] - tt)/60, digits = 1)
    
    
    
    # ONE DETECTION -------------------------------------------------------------------
    
    {
      sink("pwatch_sim_one.jags")
      
      cat("
          
          model {
          
          #nests
          for (i in 1:NI)
          {
          #chick alive at time step 1
          z[1,i] <- 1
          
          #time step
          for (t in 2:NT)
          {
          #state model
          z[t,i] ~ dbinom(p_alive[t,i], z[t-1,i])
          p_alive[t,i] <- phi[t,i] * z[t-1,i]
          
          #observation model
          y[t,i] ~ dbinom(p_sight[t,i], z[t,i])
          p_sight[t,i] <- p[t,i] * z[t,i]
          
          }
          }
          
          #transforms
          #nests
          for (i in 1:NI)
          {
          #time
          for (t in 1:NT)
          {
          logit(phi[t,i]) <- mu_phi
          logit(p[t,i]) <- mu_p
          } #t
          } #i
          
          #priors - phi
          mu_phi ~ dnorm(0, 0.01)
          
          #priors - p
          mu_p ~ dnorm(0, 0.5)

          il_mu_phi <- ilogit(mu_phi)
          
          }",fill = TRUE)
  
      sink()
    }
    
    
    
    # Starting values ---------------------------------------------------------
    
    
    Inits_1 <- list(mu_phi = 5,
                    mu_p = 0,
                    .RNG.name = "base::Mersenne-Twister",
                    .RNG.seed = 1)
    
    Inits_2 <- list(mu_phi = 5,
                    mu_p = 0,
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = 2)
    
    Inits_3 <- list(mu_phi = 5,
                    mu_p = 0,
                    .RNG.name = "base::Marsaglia-Multicarry",
                    .RNG.seed = 3)
    
    F_Inits <- list(Inits_1, Inits_2, Inits_3)
    
    
    # Parameters to track -----------------------------------------------------
    
    Pars <- c('mu_phi',
              'il_mu_phi',
              'mu_p')
    
    
    # Run model ---------------------------------------------------------------
    
    OUT_one <- jagsRun(jagsData = DATA, 
                       jagsModel = 'pwatch_sim.jags',
                       jagsInits = F_Inits,
                       params = Pars,
                       n_chain = 3,
                       n_adapt = 5000,
                       n_burn = 30000,
                       n_draw = 10000,
                       report = FALSE)
    
    
    # Analyze -----------------------------------------------------------------
    
    mn_phi <- dphi
    mn_p <- mean(dp)
    
    #individual detection
    ind_diff <- mn_phi - inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[4])
    #one detection
    one_diff <- mn_phi - inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[4])
    
    if (mn_phi > inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[3]) &
        mn_phi < inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[5]))
    {
      IND_log <- TRUE
    } else {
      IND_log <- FALSE
    }
    
    if (mn_phi > inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[3]) &
        mn_phi < inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[5]))
    {
      ONE_log <- TRUE
    } else {
      ONE_log <- FALSE
    }
    
    temp <- data.frame(GEN = j,
                       RUN = i,
                       MN_PHI = mn_phi, 
                       MN_P = mn_p,
                       I_phi_LCI = inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[3]),
                       I_phi_MED = inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[4]),
                       I_phi_UCI = inv.logit(MCMCsummary(OUT_ind, params = 'mu_phi')[5]),
                       I_diff = ind_diff,
                       I_PASS = IND_log,
                       O_phi_LCI = inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[3]),
                       O_phi_MED = inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[4]),
                       O_phi_UCI = inv.logit(MCMCsummary(OUT_one, params = 'mu_phi')[5]),
                       O_diff = one_diff, 
                       O_PASS = ONE_log,
                       I_il_lci = MCMCsummary(OUT_ind, params = 'il_mu_phi')[3],
                       I_il_med = MCMCsummary(OUT_ind, params = 'il_mu_phi')[4],
                       I_il_uci = MCMCsummary(OUT_ind, params = 'il_mu_phi')[5],
                       O_il_lci = MCMCsummary(OUT_one, params = 'il_mu_phi')[3],
                       O_il_med = MCMCsummary(OUT_one, params = 'il_mu_phi')[4],
                       O_il_uci = MCMCsummary(OUT_one, params = 'il_mu_phi')[5],
                       I_p_LCI = inv.logit(MCMCsummary(OUT_ind, params = 'mu_p')[3]),
                       I_p_MED = inv.logit(MCMCsummary(OUT_ind, params = 'mu_p')[4]),
                       I_p_UCI = inv.logit(MCMCsummary(OUT_ind, params = 'mu_p')[5]),
                       O_p_LCI = inv.logit(MCMCsummary(OUT_one, params = 'mu_p')[3]),
                       O_p_MED = inv.logit(MCMCsummary(OUT_ind, params = 'mu_p')[4]),
                       O_p_UCI = inv.logit(MCMCsummary(OUT_one, params = 'mu_p')[5]),
                       TIME_ind = time_ind)
    
    OUT <- rbind(OUT, temp)
    PB <- PB + 1
    setTxtProgressBar(pb, PB)
  }
  row.names(OUT) <- NULL
}
close(pb)



# > OUT
#     GEN RUN    MN_PHI      MN_P I_phi_LCI I_phi_MED I_phi_UCI        I_diff I_PASS O_phi_LCI O_phi_MED O_phi_UCI
# 1    1   1 0.9926285 0.5656943 0.9856388 0.9916209 0.9956767  0.0010076544   TRUE 0.9856528 0.9916679 0.9956944
# 2    1   2 0.9926285 0.5656943 0.9880138 0.9931887 0.9965988 -0.0005601178   TRUE 0.9881232 0.9931956 0.9966262
# 3    1   3 0.9926285 0.5656943 0.9847386 0.9911398 0.9954506  0.0014886947   TRUE 0.9847858 0.9911530 0.9954089
# 4    1   4 0.9926285 0.5656943 0.9918196 0.9957759 0.9981866 -0.0031473413   TRUE 0.9918205 0.9957649 0.9982090
# 5    1   5 0.9926285 0.5656943 0.9784914 0.9869643 0.9928255  0.0056642240   TRUE 0.9786032 0.9870042 0.9928752
# 6    2   1 0.9918303 0.5767803 0.9825792 0.9896006 0.9944623  0.0022297657   TRUE 0.9826267 0.9895855 0.9944578
# 7    2   2 0.9918303 0.5767803 0.9835018 0.9901085 0.9947201  0.0017217955   TRUE 0.9835604 0.9901075 0.9946995
# 8    2   3 0.9918303 0.5767803 0.9847138 0.9913995 0.9957818  0.0004308565   TRUE 0.9844647 0.9911355 0.9955753
# 9    2   4 0.9918303 0.5767803 0.9834469 0.9903817 0.9950231  0.0014486602   TRUE 0.9842946 0.9922175 1.0000000
# 10   2   5 0.9918303 0.5767803 0.9868840 0.9924003 0.9960662 -0.0005699720   TRUE 0.9868980 0.9923818 0.9960787
# 11   3   1 0.9916636 0.5720155 0.9855776 0.9918097 0.9960037 -0.0001461351   TRUE 0.9854969 0.9918195 0.9959513
# 12   3   2 0.9916636 0.5720155 0.9855184 0.9922847 0.9992594 -0.0006211387   TRUE 0.9851407 0.9915149 0.9957534
# 13   3   3 0.9916636 0.5720155 0.9841175 0.9909750 0.9955266  0.0006885639   TRUE 0.9840805 0.9909796 0.9955695
# 14   3   4 0.9916636 0.5720155 0.9892045 0.9946397 0.9999989 -0.0029761364   TRUE 0.9887949 0.9937826 0.9970830
# 15   3   5 0.9916636 0.5720155 0.9876435 0.9929589 0.9964935 -0.0012952511   TRUE 0.9876267 0.9929390 0.9964936
#           O_diff O_PASS  I_il_lci  I_il_med  I_il_uci  O_il_lci  O_il_med  O_il_uci   I_p_LCI   I_p_MED   I_p_UCI
# 1   0.0009606724   TRUE 0.9856388 0.9916209 0.9956767 0.9856528 0.9916679 0.9956944 0.4851134 0.5819708 0.6678209
# 2  -0.0005671126   TRUE 0.9880138 0.9931887 0.9965988 0.9881232 0.9931956 0.9966262 0.4971405 0.5645756 0.6328864
# 3   0.0014755520   TRUE 0.9847386 0.9911398 0.9954506 0.9847858 0.9911530 0.9954089 0.4850582 0.5627962 0.6409059
# 4  -0.0031364041   TRUE 0.9918196 0.9957759 0.9981866 0.9918205 0.9957649 0.9982090 0.4959189 0.5777225 0.6543024
# 5   0.0056242943   TRUE 0.9784914 0.9869643 0.9928255 0.9786032 0.9870042 0.9928752 0.5277686 0.5996446 0.6697682
# 6   0.0022448398   TRUE 0.9825792 0.9896006 0.9944623 0.9826267 0.9895855 0.9944578 0.5174093 0.6072875 0.6866078
# 7   0.0017227949   TRUE 0.9835018 0.9901085 0.9947201 0.9835604 0.9901075 0.9946995 0.5290054 0.6130926 0.6890525
# 8   0.0006947868   TRUE 0.9847138 0.9913995 0.9957818 0.9844647 0.9911355 0.9955753 0.4407429 0.5888300 0.6776632
# 9  -0.0003871373   TRUE 0.9834469 0.9903817 0.9950231 0.9842946 0.9922175 1.0000000 0.4754054 0.5758282 0.6734665
# 10 -0.0005514293   TRUE 0.9868840 0.9924003 0.9960662 0.9868980 0.9923818 0.9960787 0.4767220 0.5850599 0.6850635
# 11 -0.0001559138   TRUE 0.9855776 0.9918097 0.9960037 0.9854969 0.9918195 0.9959513 0.5312644 0.6111212 0.6918195
# 12  0.0001486855   TRUE 0.9855184 0.9922847 0.9992594 0.9851407 0.9915149 0.9957534 0.2147440 0.5496394 0.6358817
# 13  0.0006840067   TRUE 0.9841175 0.9909750 0.9955266 0.9840805 0.9909796 0.9955695 0.5429464 0.6150007 0.6831790
# 14 -0.0021190193   TRUE 0.9892045 0.9946397 0.9999989 0.9887949 0.9937826 0.9970830 0.2797293 0.5809285 0.6861618
# 15 -0.0012754325   TRUE 0.9876435 0.9929589 0.9964935 0.9876267 0.9929390 0.9964936 0.5164488 0.5849106 0.6490374
#       O_p_LCI   O_p_MED   O_p_UCI TIME_ind
# 1  0.4898901 0.5819708 0.6628838      6.4
# 2  0.4943305 0.5645756 0.6366552      5.9
# 3  0.4884144 0.5627962 0.6420444      6.7
# 4  0.5057977 0.5777225 0.6560570      5.4
# 5  0.5232023 0.5996446 0.6683863      7.5
# 6  0.5157991 0.6072875 0.6873760      7.6
# 7  0.5333326 0.6130926 0.6885692      6.9
# 8  0.5148639 0.5888300 0.6728887      6.9
# 9  0.1655992 0.5758282 0.6634529      7.0
# 10 0.4779864 0.5850599 0.6780807      5.7
# 11 0.5294043 0.6111212 0.6888228      7.0
# 12 0.4823084 0.5496394 0.6372791      7.1
# 13 0.5422229 0.6150007 0.6805903      6.9
# 14 0.5144682 0.5809285 0.6948145      6.1
# 15 0.5156998 0.5849106 0.6511619      5.8




