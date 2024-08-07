#########################################################################
#:-----------------------------------------------------------------------
#:----: Simulation study
#:----:
#:----: In this file we run a simulation study where we bootstrap the 
#:----:   population representative dataset and refit the model
#########################################################################


### Load libraries
library(tidyverse)
library(modelr)


### Number of simulations
number_simulations = 1000

### Number of iterations in the Bayesian model
niter = 20000
nburn = 10000

#:------------------------------------------------------------------------
#:------------------------------------------------------------------------
#:------------------------------------------------------------------------

### set up Parallel cores
library(foreach)

# collect number of cores
n.cores <- 10

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()



#:------------------------------------------------------------------------
#:------------------------------------------------------------------------
#:------------------------------------------------------------------------

### load original datasets

# load functions needed for generating data
source("data/create_data.R")

# load datasets
## Load control-case dataset
mody <- create_data(dataset = "case-control t2d")

## Load population representative dataset
united <- create_data(dataset = "united t2d", biomarkers = "reduced", commonmody = FALSE)


### set up seeds for the draws
seeds <- sample(1:1000000, number_simulations)

set.seed(123)
dataset_bootstrap <- modelr::bootstrap(united, number_simulations)



### start the loop - this happens in its own core (session due to parallel computing)
simulation <- foreach(iterations = 1:number_simulations) %dopar% {
  
  # load libraries (again as if new session)
  library(posterior)
  library(modelr)
  library(rms)
  library(nimble)
  library(tidyverse)
  
  # get the chosen dataset
  set.seed(seeds[iterations])
  bootstrap_data <- dataset_bootstrap$strap[[iterations]]
  test_cal <- bootstrap_data$data[bootstrap_data$idx,]
  
  
  ## load functions
  source("new_data_predictions/prediction_functions.R")
  
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
  x <- as.matrix(select(test_cal, agedx, bmi, hba1c, pardm, agerec, insoroha, sex))
  
  # ## remove observations with missing covariates in UNITED
  M <- test_cal$M
  
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
                 niter = niter,
                 nburnin = nburn,
                 nchains = 4,
                 inits = inits,
                 progressBar = TRUE,
                 summary = TRUE,
                 samplesAsCodaMCMC = TRUE,
                 thin = 1)
  
  # calculate Rhat values for this simulation
  rhat_run_recalibration_mixture <- run$samples$chain1 %>%
    as.data.frame() %>%
    rbind(
      run$samples$chain2 %>%
        as.data.frame(),
      run$samples$chain3 %>%
        as.data.frame(),
      run$samples$chain4 %>%
        as.data.frame()
    ) %>%
    gather() %>%
    group_by(key) %>%
    mutate(rhat = rhat(value)) %>%
    ungroup() %>%
    select(-value) %>%
    unique()
  
  
  # modify posteriors to use prediction functions
  post_recalibration_mixture <- list(post = run$samples)
  class(post_recalibration_mixture) <- "T2D"
  
  # make predictions
  newdata_x <- as_tibble(as.matrix(select(test_cal, pardm, agerec, hba1c, agedx, sex, bmi, insoroha)))
  
  preds_recalibration_mixture <- predict(post_recalibration_mixture, newdata_x, rcs_parms)
  predictions_recalibration_mixture <- apply(preds_recalibration_mixture, 2, mean)
  
  
  combined_dataset <- data.frame(
    id = bootstrap_data$idx,
    recalibration = predictions_recalibration_mixture
  )
  
  
  output <- list(
    # seed for iterations
    seed = seeds[iterations],
    # dataset same covariate distribution
    bootstrap_dataset = combined_dataset,
    # rhat for recalibration
    rhat_recalibration = rhat_run_recalibration_mixture
  )
  
  return(output)
  
}


simulation <- list(
  iterations = number_simulations,
  # number of iterations
  niter = niter,
  # number of burn in
  nburn = nburn,
  # simulation values
  simulations = simulation
)


### save output
saveRDS(simulation, "bootstrap_riley_simulation/output/simulation_t2d.rds")



### turn off Parallel cores
parallel::stopCluster(cl = my.cluster)

