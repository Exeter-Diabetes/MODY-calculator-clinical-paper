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
mody <- create_data(dataset = "case-control t1d")

## Load population representative dataset
united <- create_data(dataset = "united t1d", biomarkers = "reduced", commonmody = FALSE)


### set up seeds for the draws
seeds <- sample(1:1000000, number_simulations)

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
  
  ## Estimation of prior distribution for C- or A+ patients
  pM <- rbeta(100000, 2 + 1, 247 - 2 + 1)
  pT_M <- rbeta(100000, 6 + 1, 282 - 6 + 1)
  p <- pM * pT_M
  
  ## fit beta distribution to this
  fn <- function(pars, psamp) {
    alpha <- pars[1]
    beta <- pars[2]
    if(alpha < 0 | beta < 0) return(NA)
    sum(dbeta(psamp, alpha, beta, log = TRUE))
  }
  pbeta <- optim(c(1, 1), fn, control = list(fnscale = -1), psamp = p)
  
  ### Define model structure
  model_code <- nimbleCode({
    
    ## case-control likelihood component
    for (i in 1:nCC) {
      
      ## regression mean
      logit(pCC[i]) <- beta0 + inprod(beta[1:np], xCC[i, 1:np])
      
      ## likelihood
      MCC[i] ~ dbern(pCC[i])
    }
    
    ## real world likelihood
    for (j in 1:nind) { # iterate through patients from real world data
      
      ## prediction from UNITED
      Z[j] <- beta0 + inprod(beta[1:np], x[j, 1:np])
      
      ## regression of real world data - shrinkage applied
      logit(p_shrinkage[j]) <- gamma0 + gamma1 * Z[j]
      
      ## MODY modelling
      M[j] ~ dbern(p[j]) # probability of MODY
      p[j] <- T[j] * pMp_Cn_or_Ap + ( 1 - T[j] ) * p_shrinkage[j] # decomposition of MODY based on T
      T[j] ~ dbern(pT[j]) # explaining how T should be defined
      
      logit(pT[j]) <- beta_t0 + inprod(beta_t[1:np_T], xT[j, 1:np_T]) + inprod(beta_spline[1:np_spline], x_spline[j, 1:np_spline])
    }
    
    ## testing priors
    pMp_Cn_or_Ap ~ dbeta(pMp_Cn_or_Ap_alpha, pMp_Cn_or_Ap_beta)       # probability of MODY given C- or A+ is < 0.01
    
    beta0 ~ dnorm(0, 0.01)
    for(j in 1:np) {
      beta[j] ~ dnorm(0, 0.01)
    }
    
    beta_t0 ~ dnorm(0, 0.01)
    for(j in 1:np_T) {
      beta_t[j] ~ dnorm(0, 0.01)
    }
    for(j in 1:np_spline) {
      beta_spline[j] ~ dnorm(0, 0.01)
    }
    
    ## regression priors
    gamma0 ~ dnorm(0, sd = 10)
    gamma1 ~ dnorm(0, sd = 10)
    
  })
  
  ### Set up data for model
  #:--- Case control covariates
  xCC <- as.matrix(select(mody, pardm, agerec, hba1c, agedx, sex))
  MCC <- as.numeric(as.character(mody$mody))
  
  #:--- UNITED covariates
  x <- as.matrix(select(test_cal, pardm, agerec, hba1c, agedx, sex))
  xT <- as.matrix(select(test_cal, bmi, agedx, pardm, agerec))
  x_spline <- test_cal %>%
    select(bmi, agedx, agerec) %>%
    mutate(bmi_spline = as.numeric(rms::rcs(bmi, 3)[, 2]),
           agedx_spline = as.numeric(rms::rcs(agedx, 3)[, 2]),
           agerec_spline = as.numeric(rms::rcs(agerec, 3)[, 2])) %>%
    select(-bmi, -agedx, -agerec) %>%
    set_names(c("bmi_spline", "agedx_spline", "agerec_spline")) %>%
    as.matrix()
  
  rcs_parms <- data.frame(
    bmi = attr(rms::rcs(test_cal$bmi, 3)[, 2], "parms"),
    agedx = attr(rms::rcs(test_cal$agedx, 3)[, 2], "parms"),
    agerec = attr(rms::rcs(test_cal$agerec, 3)[, 2], "parms")
  )
  
  ## remove observations with missing covariates in UNITED
  M <- test_cal$M
  T <- test_cal$T
  
  ### set up data list for NIMBLE
  data <- list(M = M,
               MCC = MCC,
               T = T)
  
  ### set up other components of model
  consts <- list(nind = length(M),             # number of individuals in UNITED
                 nCC = length(MCC),            # number of individuals in Case control
                 np = ncol(x),                 # number of covariates
                 np_T = ncol(xT),                 # number of covariates T
                 np_spline = ncol(x_spline),   # number of covariates splines
                 xCC = xCC,                    # covariate matrix for Case control
                 x = x,                        # covariate matrix for UNITED
                 xT = xT,                        # covariate matrix for UNITED T
                 x_spline = x_spline,          # covariate matrix for UNITED splines
                 pMp_Cn_or_Ap_alpha = pbeta$par[1],
                 pMp_Cn_or_Ap_beta = pbeta$par[2])
  
  ### set up initial values for parameters in the model
  initFn <- function(M, x, T, np, np_T, np_spline, pbeta) {
    
    ## simulate from priors
    gamma0 <- rnorm(1)
    gamma1 <- rnorm(1)
    pMp_Cn_or_Ap <- rbeta(1, pbeta[1], pbeta[2])
    beta0 <- rnorm(1)
    beta <- rnorm(np)
    beta_t0 <- rnorm(1)
    beta_t <- rnorm(np_T)
    beta_spline <- rnorm(np_spline)
    
    ## simulate missing information for MODY
    Z <- t(beta %*% t(x))
    M1 <- gamma0 + gamma1 * Z
    M1 <- exp(M1) / (1 + exp(M1))
    M1 <- rbinom(length(M1), 1, M1)
    M1[!is.na(M)] <- NA
    
    # simulate latent missing variable
    T1 <- T
    T1[is.na(T)] <- 1
    T1[!is.na(T)] <- NA
    
    ## return initial values
    list(gamma0 = gamma0,
         gamma1 = gamma1,
         beta0 = beta0,
         beta = beta,
         beta_t0 = beta_t0,
         beta_t = beta_t,
         beta_spline = beta_spline,
         pMp_Cn_or_Ap = pMp_Cn_or_Ap,
         T = T1,
         M = M1
    )
  }
  
  
  # adding this so that the logProb isn't -Inf
  logProb = "-Inf"
  
  while (logProb == "-Inf") {
    
    ### set list of initial values for NIMBLE
    inits <- initFn(data$M, consts$x, data$T, consts$np, consts$np_T, consts$np_spline, pbeta$par)
    
    ### define the model, data, inits and constants
    model <- nimbleModel(code = model_code, constants = consts, data = data, inits = inits)
    
    logProb = model$calculate()
    
    # print(logProb)
    
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
  class(post_recalibration_mixture) <- "T1D"
  
  # make predictions
  newdata_x <- as_tibble(as.matrix(select(test_cal, pardm, agerec, hba1c, agedx, sex, bmi, C, A)))
  
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
saveRDS(simulation, "bootstrap_riley_simulation/output/simulation_t1d.rds")


### turn off Parallel cores
parallel::stopCluster(cl = my.cluster)

