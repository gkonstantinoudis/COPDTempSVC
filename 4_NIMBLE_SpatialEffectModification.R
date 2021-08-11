



# STEP 3 : Spatial effect modification


#########################################################################



library(nimble)
library(dplyr)
library(tidyr)
library(sf)
library(spdep)
library(fastDummies)
library(doBy)


options(scipen = 999)

setwd("C:/Users/gkonstan/Desktop/COPD_temperature/")



# Load the simulated data


COPD_dat <- readRDS("dat")



# calculate the quantiles and scale
quantile(COPD_dat$temperature, prob = seq(from = .50, to = .95, by = .05))
(quantile(COPD_dat$temperature, prob = seq(from = .50, to = .95, by = .05)) - 
    mean(COPD_dat$temperature))/sd(COPD_dat$temperature) -> tmp.quant




# Here we define thr to be the threshold of the model that minimizes the WAIC on step 1. 
thr <- 8
x.change <- tmp.quant[thr]




# NOTE!! The adj here refers to the spatial effect modification
adj_EM <- TRUE # TO RUN THE ADJUSTED MODELS
adj_EM <- FALSE # TO RUN THE UNADJUSTED MODELS


# Read the shp with the spatial effect modifiers.
shp <- readRDS("SpatialEffectModifiers")

shp.EM <- shp
shp.EM$geometry <- NULL
shp.EM <- shp.EM[,c("IMD_quintiles", "Broad.RUC11", "gr.space", "average.tmp")]
tmpIMD <- dummy_cols(shp.EM$IMD_quintiles)[,-c(1:2)]
tmpURBN <- dummy_cols(shp.EM$Broad.RUC11)[,-c(1:2)]

shp.EM <- cbind(shp.EM[,c("gr.space", "average.tmp")],
                tmpIMD,
                tmpURBN)

# lets scale the green space too
shp.EM$gr.space <- as.vector(scale(shp.EM$gr.space))
shp.EM$average.tmp <- as.vector(scale(shp.EM$average.tmp))




#----------------------------------------- Model: STEP 3


model_SVC <- nimbleCode(
  {
    for(i in 1:N){

      O[i] ~ dpois(mu[i])
      mu[i] <- exp(beta_tmp[i]*temperature[i] + inprod(beta[1:K.b], X[i,1:K.b]) + u[IDCC[i]])
      Q[i] <- step(temperature[i] - x.change)
      beta_tmp[i] <- (1 - Q[i])*b_tmp_low + Q[i]*b_tmp_SVC[ID.sp[i]]

    }

    for(j in 1:J){

      u[j] ~ dnorm(0, sd = 100)

    }

    for(m in 1:M){
      
      b[m] <- v[m] + w[m]
      
      if(adj_EM){
        b_tmp_SVC[m] <- b[m] + inprod(gamma[1:EM.b], Z[m,1:EM.b])
      }else{
        b_tmp_SVC[m] <- b[m]
      }
      v[m] ~ dnorm(b_tmp_high, sd = sd.v)

    }

    w[1:M] ~ dcar_normal(adj[1:LM], weights[1:LM], num[1:M], tau.w, zero_mean = 1)
    # set the priors of the fixed effects

    # confounders
    for(kb in 1:K.b){

    beta[kb] ~ dnorm(0, sd = 1)

    }


    # spatial effect modifiers
    if(adj_EM){
      for(emb in 1:EM.b){
        gamma[emb] ~ dnorm(0, sd = 1)
      }
    }else{
    }
    
    b_tmp_low ~ dnorm(0, sd = 1)
    b_tmp_low_unscaled <- b_tmp_low/x.sd
    b_tmp_high ~ dnorm(0.01570643, sd = 0.004612425) # the posterior of the WAIC analysis to help convergence
    b_tmp_high_unscaled <- b_tmp_high/x.sd

    # set the priors of the hyperparameter of the random effects
    sd.v ~ dgamma(shape = 1, rate = 2)
    sd.w ~ dgamma(shape = 1, rate = 2)
    tau.w <- 1/(sd.w*sd.w)

    # monitor some us and mus
    mu_keep[1:3] <- c(mu[1], mu[100], mu[153])
    u_keep[1:3] <- c(u[100], u[62], u[105])

  }
)



# data
N <- nrow(COPD_dat)
J <- max(COPD_dat$ID)
N.LTLA <- max(COPD_dat$ladid)







# Modify spatial data so there are no discontinuities:

shp_nb <- poly2nb(shp)
W.scale <- nb2mat(shp_nb, zero.policy = TRUE, style = "B")
which(apply(W.scale, 2, sum)==0)
# the two truths are two islands, namely in the isle of Scilly and isle of Wight. 

# I will connect the isle of Scilly with Cornwal and the isle of Wight with Southampton and Portsmouth
# For the isle of Scilly:
W.scale[52, 51] <- W.scale[51, 52] <- 1
W.scale[46, 45] <- W.scale[45, 46] <- W.scale[46, 44] <- W.scale[44, 46] <- 1

tmp_nb <- mat2listw(W.scale)
nbWB_B <- listw2WB(tmp_nb)




# confounder data

conf.mat <- cbind(

     COPD_dat$hol,
     as.vector(scale(COPD_dat$O3)),
     as.vector(scale(COPD_dat$PM25)),
     as.vector(scale(COPD_dat$RH))

)

K.b <-  ncol(conf.mat)
EM.b <- ncol(shp.EM)

COPD_NIMBLE_data <- list(

  O = COPD_dat$O,
  temperature = as.vector(scale(COPD_dat$temperature)),
  X = conf.mat, 
  Z = shp.EM

)




COPD_NIMBLE_constants <- list(

  N = N,
  J = J,
  K.b =  K.b,
  IDCC = COPD_dat$ID,
  ID.sp = COPD_dat$ladid,

  LM = length(nbWB_B$weights),
  M = nrow(shp),

  adj = nbWB_B$adj,
  num = nbWB_B$num,
  weights = nbWB_B$weights,

  x.sd = sd(COPD_dat$temperature),
  EM.b = EM.b,
  x.change = x.change, 
  adj_EM = adj_EM

)






# ------------------------- NIMBLE RUN



# initials

if(adj_EM == TRUE){
  
  initials =
    list(
      beta = rep(0, K.b),
      u = rep(0, times = J),
      v = rep(0, times = nrow(shp)),
      w = rep(0, times = nrow(shp)),
      gamma = rep(0, EM.b),
      sd.v = 0.1, 
      sd.w = 0.01,
      b_tmp_low = 0, 
      b_tmp_high = 0
    )
  
  
  # parameters to monitor
  parameters = c("beta", "mu_keep", "u_keep", 
                 "b_tmp_low", "b_tmp_high", 
                 "b_tmp_low_unscaled","b_tmp_high_unscaled",
                 "sd.v", "sd.w", "b", "gamma")
  
  nam.store <- "adjusted"
  
}else{
  
  initials =
    list(
      beta = rep(0, K.b),
      u = rep(0, times = J),
      v = rep(0, times = nrow(shp)),
      w = rep(0, times = nrow(shp)),
      sd.v = 0.1, 
      sd.w = 0.01,
      b_tmp_low = 0, 
      b_tmp_high = 0
    )
  
  
  # parameters to monitor
  parameters = c("beta", "mu_keep", "u_keep", 
                 "b_tmp_low", "b_tmp_high", 
                 "b_tmp_low_unscaled","b_tmp_high_unscaled",
                 "sd.v", "sd.w", "b")
  
  nam.store <- "unadjusted"
  
}







# NIMBLE call

# define the model
t_0 <- Sys.time()
nimble_model <- nimbleModel(
  code = model_SVC,
  data = COPD_NIMBLE_data,
  constants = COPD_NIMBLE_constants,
  inits = initials,
  name = "model_SVC",
  calculate = FALSE
)

nimble_model$initializeInfo()

# compile model
Cmodel <- compileNimble(nimble_model)
Cmodel$calculate()

# nimble_model$origInits


# configure and build the MCMC
conf <- configureMCMC(nimble_model, monitors = parameters)
MCMC <- buildMCMC(conf)

# compile the MCMC
cMCMC <- compileNimble(MCMC, project = Cmodel)
t_1 <- Sys.time()
t_1 - t_0



# MCMC setting
ni <- 250000  # nb iterations
nt <- 100     # thinning interval
nb <- 150000  # nb iterations as burn-in
nc <- 1       # nb chains


# run the MCMC
t_0 <- Sys.time()
modre_SVC <- runMCMC(cMCMC, niter = ni , nburnin = nb, thin = nt, samples = TRUE, summary = FALSE)
t_1 <- Sys.time()
t_1 - t_0

saveRDS(mod_SVC_res, file = paste0("mod_SVC_res_", nam.store))


#########################################################################
#########################################################################
#########################################################################
#########################################################################


