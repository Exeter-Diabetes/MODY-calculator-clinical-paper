#:--------------------------------------------------------
# 
# In this file, we fit a Bayesian shrinkage model to non-insulin treated patients
#
#:--------------------------------------------------------

# load libraries
library(nimble)
library(rms)
library(tidyverse)

# load functions needed for generating data
source("data/create_data.R")

# load datasets
## Load control-case dataset
mody <- create_data(dataset = "case-control t2d")

## Load population representative dataset
united <- create_data(dataset = "united t2d")

# Model fitting

### Define model structure
model_code <- nimbleCode({
  
  ## case-control likelihood component
  for (i in 1:nCC) {
    
    ## regression mean
    logit(pCC[i]) <- beta0 + inprod(beta[1:np], xCC[i, 1:np])
    
    ## likelihood
    MCC[i] ~ dbern(pCC[i])
  }
  
  ## likelihood component
  for (j in 1:nind) {
    
    ## prediction from UNITED
    Z[j] <- beta0 + inprod(beta[1:np], x[j, 1:np])
    
    ## regression mean - with shrinkage applied
    logit(p[j]) <- gamma0 + gamma1 * Z[j]
    
    ## likelihood
    M[j] ~ dbern(p[j])
    
  }
  
  ## regression priors
  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  for(j in 1:np) {
    beta[j] ~ dnorm(0, 0.01)
  }
  
})


### Set up data for model
#:--- Case control covariates
xCC <- as.matrix(select(mody, agedx, bmi, hba1c, pardm, agerec, insoroha, sex))
MCC <- as.numeric(as.character(mody$mody))

#:--- UNITED covariates
x <- as.matrix(select(united, agedx, bmi, hba1c, pardm, agerec, insoroha, sex))

# ## remove observations with missing covariates in UNITED
M <- united$M

### set up data list for NIMBLE
data <- list(M = M, MCC = MCC)

### set up other components of model
consts <- list(nind = length(M),             # number of individuals in UNITED
               nCC = length(MCC),            # number of individuals in Case control
               np = ncol(x),                 # number of covariates
               xCC = xCC,                    # covariate matrix for Case control
               x = x)                        # covariate matrix for UNITED\

### set up initial values for parameters in the model
initFn <- function(M, x, np) {
  
  ## simulate from priors
  gamma0 <- rnorm(1)
  gamma1 <- rnorm(1)
  beta0 <- rnorm(1)
  beta <- rnorm(np)
  
  ## simulate missing information
  Z <- t(beta %*% t(x))
  M1 <- gamma0 + gamma1 * Z
  M1 <- exp(M1) / (1 + exp(M1))
  M1 <- rbinom(length(M1), 1, M1)
  
  M1[!is.na(M)] <- NA
  
  ## return initial values
  list(gamma0 = gamma0,
       gamma1 = gamma1,
       beta0 = beta0,
       beta = beta,
       M = M1
  )
}


# adding this so that the logProb isn't -Inf
logProb <- "-Inf"

while (logProb == "-Inf") {
  
  ### set list of initial values for NIMBLE
  inits <- initFn(data$M, consts$x, consts$np)
  
  ### define the model, data, inits and constants
  model <- nimbleModel(code = model_code, constants = consts, data = data, inits = inits)
  
  logProb <- model$calculate()
  
}

### compile the model
cModel <- compileNimble(model)

### configure MCMC and monitor parameters needed for prediction
config <- configureMCMC(cModel)
config$printMonitors()


### build the model
built <- buildMCMC(config)
cBuilt <- compileNimble(built)

### run the model
run <- runMCMC(cBuilt,
               niter = 500000,
               nburnin = 300000,
               nchains = 4,
               progressBar = TRUE,
               summary = TRUE,
               samplesAsCodaMCMC = TRUE,
               thin = 1)

### save the model for future comparisons
saveRDS(run, "model_development/type_2_model_posteriors.rds")
