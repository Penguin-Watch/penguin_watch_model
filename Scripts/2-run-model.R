######################
#Mark-recapture for penguin watch 
#
#script to run model - model object saved as .rds
#
#Authors: Casey Youngflesh
######################


#QUESTIONS about model
#want to know how survival is changing. Add env covariates onto this
#no replication over period of closure (how to treat night - all images over one day as a single observation? should each hour be a time step with NAs over night?) - are the parameters identifiable? p219 Kerry and Schaub 2012
#should the random effect (time) for phi be in there? Is it necessary?


#how to model more than one site? maybe code each nest/site as a different nest (one site)


#nonidentifiability checks
#plot p[1,1] against phi[1,1] - CHECK
#check correlation p and phi - CHECK
#compare estimates p and phi against true p and phi - 
#plot posteriors for:
#1) betas - dnorm(0, 1) T(-1,1) - CHECK
#2) mean_phi and mean_p - dunif(0,0.1) - CHECK
#3) eps_phi and eps_p - dnorm(0, tau_p) T(-20,20) - p okay, phi maybe weakly identifiable
#and compare to priors


# Clear environment -------------------------------------------------------

rm(list = ls())



# Load packages -----------------------------------------------------------

# if('pacman' %in% rownames(installed.packages()) == FALSE)
# {
#   install.packages('pacman')
# }
# 
# pacman::p_load(rjags, MCMCvis, parallel)
#devtools::install_github('caseyyoungflesh/jagsRun')

require(jagsRun)


# Data for model ----------------------------------------------------------


un_sites <- unique(PW_data$site)
for (i in 1:length(un_sites))
{
  #i <- 1
  temp <- filter(PW_data, site == un_sites[i])
  un_yrs <- unique(temp$season_year)
  
  for (j in 1:length(un_yrs))
  {
    #j <- 1
    temp2 <- filter(temp, season_year == un_yrs[j])
    
  }
  
}




sim_data <- readRDS('sim_data.rds')
z_vals <- readRDS('z_vals.rds')


DATA <- list(
  y = sim_data, #reponse
  N = NROW(sim_data), #number of nests
  L = NCOL(sim_data), #number of time points
  z = z_vals,
  x = 1:NCOL(sim_data)) 



# Model -------------------------------------------------------------------

#L = length of time series (400)
#N = number of nests (5)

{
sink("mark_recapture.jags")

cat("
    model {
  
      for (i in 1:N)
    {
    #both chicks alive at time step 1
    z[i,1] <- 2
    y.new[i,1] ~ dbinom(p[i,1], 2)
    
    for (t in 2:L)
    {
    
    #state model
    z[i,t] ~ dbinom(p_alive[i,t], z[i,t-1])
    p_alive[i,t] <- ifelse(z[i,t-1] < 2, 
                        phi[i,t] * z[i,t-1],
                        phi[i,t])
    
    #observation model
    y[i,t] ~ dbinom(p_sight[i,t] * w[i,t], z[i,t]) #binary day/night
    p_sight[i,t] <- ifelse(z[i,t] < 2,
                        p[i,t] * z[i,t],
                        p[i,t])
    
    
    #PPC
    y.new[i,t] ~ dbinom(p_sight[i,t] * w[i,t], z[i,t]) #binary day/night
    }
    }
    
    #PPC
    #mean
    mn.y <- mean(y)
    mn.y.new <- mean(y.new)
    pv.mn <- step(mn.y.new - mn.y)
    
    #sd
    sd.y <- sd(y)
    sd.y.new <- sd(y.new)
    pv.sd <- step(sd.y.new - sd.y)
    
    
    #transforms
    for (i in 1:N)
    {
    for (t in 1:L)
    {
    logit(phi[i,t]) <- mu_phi + beta_phi*x[t] + eps_phi[t]       #phi = survival prob
    logit(p[i,t]) <- mu_p + beta_p*x[t] + eps_p[i]            #p = detection prob
    }
    }
    
    
    #priors
    for (t in 1:L)
    {
    eps_phi[t] ~ dnorm(0, tau_phi) T(-10,10)
    }
    
    mean_phi ~ dbeta(1.5,1.5)                 #Mean survival
    mu_phi <- log(mean_phi / (1 - mean_phi))
    tau_phi <- pow(sigma_phi, -2)
    sigma_phi ~ dunif(0.25, 3)
    
    
    for (i in 1:N)
    {
    eps_p[i] ~ dnorm(0, tau_p) T(-10,10)
    }
    
    
    mean_p ~ dbeta(1.5,1.5)                    #Mean detection - could use alternative below
    mu_p <- log(mean_p / (1 - mean_p))         #Logit transform - could use alternative below
    #mean_p <- 1 / (1+exp(-mu_p))              #Mean detection - Inv-logit transform
    #mu_p ~ dnorm(0, 0.001)                    #Prior for logit of mean survival
    tau_p <- pow(sigma_p, -2)
    sigma_p ~ dunif(0.25, 3)
    
    beta_phi ~ dnorm(0, 1000) T(0,1) #[slope only pos] maybe variance 0.01 (precision 100) - plot histogram to get a look (will depend on time step length [i.e., one hour or one day])
    beta_p ~ dnorm(0, 100) T(0,1) #[slope only pos] maybe variance 0.1 (precision 10) - plot histogram to get a look (will depend on time step length [i.e., one hour or one day])
    
    
    }",fill = TRUE)

sink()
}



# Starting values ---------------------------------------------------------


Inits_1 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma_phi = runif(1, 0, 10),
                sigma_p = runif(1, 0, 10),
                beta_phi = 0,
                beta_p = 0,
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = 1)

Inits_2 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma_phi = runif(1, 0, 10),
                sigma_p = runif(1, 0, 10),
                beta_phi = 0,
                beta_p = 0,
                .RNG.name = "base::Wichmann-Hill",
                .RNG.seed = 2)

Inits_3 <- list(mean_phi = runif(1, 0, 1),
                mean_p = runif(1, 0, 1),
                sigma_phi = runif(1, 0, 10),
                sigma_p = runif(1, 0, 10),
                beta_phi = 0,
                beta_p = 0,
                .RNG.name = "base::Marsaglia-Multicarry",
                .RNG.seed = 3)

F_Inits <- list(Inits_1, Inits_2, Inits_3)



# Parameters to track -----------------------------------------------------

Pars <- c('mean_phi',
          'mean_p',
          'sigma_p',
          'sigma_phi',
          'beta_p',
          'beta_phi',
          'mu_phi',
          'mu_p',
          'eps_phi',
          'eps_p')


# Inputs for MCMC ---------------------------------------------------------

JAGS_FILE <- 'mark_recapture.jags'
n_adapt <- 10000  # number for initial adapt
n_burn <- 80000 # number burnin
n_draw <- 20000  # number of final draws to make
n_thin <- 2    # thinning rate
n_chain <- 3  # number of chains

Rhat_max <- 1.01 # max allowable Rhat (close to 1 = convergence)
n_max <- 70000 # max allowable iterations


# Run model (parallel) ---------------------------------------------------------------

#number of chains
cl <- parallel::makeCluster(n_chain)

pid <- NA
for(i in 1:n_chain)
{
  pidNum <- capture.output(cl[[i]])
  start <- regexpr("pid", pidNum)[[1]]
  end <- nchar(pidNum)
  pid[i] <- substr(pidNum, (start + 4), end)
}
 
parallel::clusterExport(cl, 
                        c('DATA', 
                          'n_adapt', 
                          'n_burn', 
                          'n_draw',
                          'n_thin',
                          'Pars', 
                          'pid',
                          'F_Inits',
                          'JAGS_FILE'
                        ))

ptm <- proc.time()
out.1 <- parallel::clusterEvalQ(cl, 
                      {
                        require(rjags)
                        processNum <- which(pid==Sys.getpid())
                        m.inits <- F_Inits[[processNum]]
                        
                        jm = jags.model(data = DATA, 
                                        file = paste0(JAGS_FILE), 
                                        inits = m.inits, 
                                        n.chains = 1, 
                                        n.adapt = n_adapt)
                        
                        update(jm, 
                               n.iter = n_burn)
                        
                        samples = coda.samples(jm, 
                                               n.iter = n_draw, 
                                               variable.names = Pars, 
                                               thin = n_thin)
                        return(samples)
                        })


out <- coda::mcmc.list(out.1[[1]][[1]], 
                       out.1[[2]][[1]], 
                       out.1[[3]][[1]])




#more iterations if not converged
n_total <- n_burn + n_draw
n_extra <- 0
while(max(MCMCsummary(out)[,5], na.rm = TRUE) > Rhat_max &
      n_total < n_max)
{
  
  out.2 <- clusterEvalQ(cl, 
                        {
                          require(rjags)
                          processNum <- which(pid==Sys.getpid())
                          m.inits <- F_Inits[[processNum]]
                          
                          samples = coda.samples(jm, 
                                                 n.iter = n_draw, 
                                                 variable.names = Pars, 
                                                 thin = n_thin)
                          return(samples)
                        })
  
  out <- coda::mcmc.list(out.2[[1]][[1]], 
                         out.2[[2]][[1]], 
                         out.2[[3]][[1]])
  
  n_extra <- n_extra + n_draw
  n_total <- n_total + n_draw
}

stopCluster(cl)

n_final <- floor((n_draw + n_extra)/n_thin)
print(paste0('Total iterations: ', n_final))
(proc.time() - ptm)[3]/60 #minutes


#saveRDS(out, 'out_10a_80b_20d_101.rds')


#Inferences were derived from $`r n_final`$ samples drawn following an adaptation period of $`r n_adapt`$ draws, and a burn-in period of $`r (n_total - n_draw)`$ draws using $`r n_chain`$ chains and a thinning rate of $`r n_thin`$.


