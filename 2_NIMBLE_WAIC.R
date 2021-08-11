




# STEP 1 : WAIC analysis


#########################################################################



library(nimble)
library(dplyr)
library(tidyr)

options(scipen = 999)

setwd("C:/Users/gkonstan/Desktop/COPD_temperature/")



# Load the simulated data


COPD_dat <- readRDS("dat")


# calculate the quantiles and scale
quantile(COPD_dat$temperature, prob = seq(from = .50, to = .95, by = .05))
(quantile(COPD_dat$temperature, prob = seq(from = .50, to = .95, by = .05)) - 
mean(COPD_dat$temperature))/sd(COPD_dat$temperature) -> tmp.quant


# Number in parallel

# the following lines are to run the model in parallel from terminal. If this is not intended one can fix thr 
# thr takes values from 1 to 10 (number of temperature thresholds to be assessed)
thr <- commandArgs()[6]
print(thr)
thr <- as.numeric(thr)



x.change <- tmp.quant[thr]


adj <- TRUE # TO RUN THE ADJUSTED MODELS
adj <- FALSE # TO RUN THE UNADJUSTED MODELS


#----------------------------------------- Model STEP 1



# code
model_LT <- nimbleCode(
  {
    for(i in 1:N){

      O[i] ~ dpois(mu[i])
      
      if(adj){
        mu[i] <- exp(beta_tmp[Q[i]]*temperature[i] + inprod(beta[1:K.b], X[i,1:K.b]) + u[IDCC[i]])
      }else{
        mu[i] <- exp(beta_tmp[Q[i]]*temperature[i] + u[IDCC[i]])
      }
      

      Q[i] <- 1 + step(temperature[i] - x.change)

    }

    for(j in 1:J){
      u[j] ~ dnorm(0, sd = 100) # the fixed effect to make the Poisson conditional logistic regression
    }

    for (k in 1:2) {
      beta_tmp[k] ~ dnorm(0, sd = 1)
    }
    
    beta_tmp_unscaled[1] <- beta_tmp[1]/x.sd
    beta_tmp_unscaled[2] <- beta_tmp[2]/x.sd
    
    # confounders
    if(adj){
      for(kb in 1:K.b){
        beta[kb] ~ dnorm(0, sd = 1)
      }
    }else{
      
    }

    # monitor some nodes of lp and re
    mu_keep[1:3] <- c(mu[1], mu[100000], mu[153])
  }
)






# data

N <- nrow(COPD_dat)
J <- max(COPD_dat$ID)


conf.mat <- cbind(

     COPD_dat$hol,
     as.vector(scale(COPD_dat$O3)),
     as.vector(scale(COPD_dat$PM25)),
     as.vector(scale(COPD_dat$RH))

)

K.b <-  ncol(conf.mat)


COPD_NIMBLE_data <- list(

  O = COPD_dat$O,
  temperature = as.vector(scale(COPD_dat$temperature)),
  X = conf.mat
  
)


COPD_NIMBLE_constants <- list(
  
  N = N,
  IDCC = COPD_dat$ID,
  J = J,
  K.b = K.b,
  x.change = x.change, 
  x.sd = sd(COPD_dat$temperature), 
  adj = adj
  
)




# ------------------------- NIMBLE RUN



# initials

if(adj == TRUE){
  
  initials =
    list(
      beta_tmp = rep(0, times = 2),
      beta = rep(0, K.b),
      u = rep(-1, times = J)
    )
  
  # parameters to monitor
  parameters = c("beta_tmp", "beta_tmp_unscaled", "mu_keep", "u", "beta")
  nam.store <- "adjusted"
  
}else{
  
  initials =
    list(
      beta_tmp = rep(0, times = 2),
      u = rep(-1, times = J)
    )
  
  # parameters to monitor
  parameters = c("beta_tmp", "beta_tmp_unscaled", "mu_keep", "u")
  nam.store <- "unadjusted"
  
}





# NIMBLE call

# define the model
t_0 <- Sys.time()
nimble_model <- nimbleModel(
  code = model_LT,
  data = COPD_NIMBLE_data,
  constants = COPD_NIMBLE_constants,
  inits = initials,
  name = "model_LT",
  calculate = FALSE
)

# compile model
Cmodel <- compileNimble(nimble_model)
Cmodel$calculate()


nimble_model$origInits


# configure and build the MCMC
conf <- configureMCMC(nimble_model, monitors = parameters, enableWAIC = TRUE)
MCMC <- buildMCMC(conf)

# compile the MCMC
cMCMC <- compileNimble(MCMC, project = Cmodel)
t_1 <- Sys.time()
t_1 - t_0



# MCMC setting

ni <- 200000   # nb iterations
nt <- 100      # thinning interval
nb <- 100000   # nb iterations as burn-in
nc <- 1       # nb chains


# run the MCMC
t_0 <- Sys.time()
mod_LT_res <- runMCMC(cMCMC, niter = ni, nburnin = nb, thin = nt, samples = TRUE, summary = FALSE, WAIC = TRUE)
t_1 <- Sys.time()
t_1 - t_0

# and store

readRDS(mod_LT_res, paste0("mod_LT_res_", nam.store, thr))


#########################################################################
#########################################################################
#########################################################################
#########################################################################

