rm(list = ls())
################ LBP data generation scenario 1 ###################
scenario = 1 # rho = 0.1

# Load libraries
library(fields)
library(mvnfast)
library(rstan)
library(coda)
library(foreach)

# just for demonstration, one core
cl <- parallel::makeCluster(1); doParallel::registerDoParallel(1)

# <<<<<<<<<<<<<<<<<<<<   uncomment below for actucal parallel run
#ncores = 34
#if(ncores >= parallel::detectCores()) stop("ncores is greater than the number of detected cores")
#cl <- parallel::makeCluster(ncores)
#doParallel::registerDoParallel(cl)


# Source files to run the simulation
# data generation
source("simul/simul_data_lbp.R")
# lbp
source("R/polya.R")
source("R/LBBer_Matern.R")
source("R/LBBer_Matern_unblocked.R")
source("R/LBBer_Matern_predict.R")
# copula
source("simul/copula_codes/copula_Matern_predict.R")
# running with prediction and summary
source("simul/run_copula.R")
source("simul/run_LBBer.R")

# Parameters of the beta marginals
a <- 1
b <- 2

# Burn-in and number of iterations
nburn <- 1000
nsamples <- 1000
nsims = 1 # <<<<<<<<<<<<<<<<<<<<  increase to 100 in actual run

# Set the parameters for the simulation
n <- 400
ntest <- 100


save = FALSE # <<<<<<<<<<<<<<<<< change to TRUE in actual run

results = foreach(isim = 1:nsims,
                  .export = ls(globalenv()),
                  .errorhandling = "pass",
                  .packages = c("coda","rstan","mvnfast","fields"),
                  .combine = "c") %dopar% {
                    # step 0: Create the directory to store the output
                    out_dir <- paste0("simul/output_lbpdatafit/scenario",scenario,"/simulation_",isim)
                    if(!dir.exists(out_dir)){
                      dir.create(out_dir)
                    }
                    # step 1: data generation and set seed
                    data <- simulate_data_lbp(n = n, ntest = ntest, seed = isim, scenario = scenario)
                    if(save) saveRDS(data, file = paste0(out_dir, "/data_scenario_", scenario, ".rds"))
                    
                    # Set the model parameters
                    model_pars <- list(a = 1,
                                       b = 2,
                                       nsamples = nsamples,
                                       nburn = nburn,
                                       range_lb = 0.01,
                                       range_ub = 0.5,
                                       ngrid = 50)
                    
                    # Input of the stan function
                    stan_input <- list(n = n,
                                       z = data$z,
                                       D = data$Dxx,
                                       a = 1, # this is just a placeholder, a=1 fixed anyways
                                       b = 2,
                                       range_lb = 0.01,
                                       range_ub = 0.5,
                                       m = nrow(data$Duu),
                                       D_star = data$Duu,
                                       D_site_star = data$Dxu)
                    pars_save <- c("theta", "rangepar")
                    
                    # Run the logistic beta, setting 1, block, adapt, MH
                    results_LBP <- run_LBBer(data = data, model_pars = model_pars,
                                             nparticle = 0, adapt = TRUE, block = TRUE)
                    if(save) saveRDS(results_LBP, file = paste0(out_dir, "/LBP_scenario",scenario,"_setting1.rds"))
                    
                    # Run the logistic beta, setting 2, block, nonadapt, MH
                    results_LBP <- run_LBBer(data = data, model_pars = model_pars,
                                             nparticle = 0, adapt = FALSE, block = TRUE)
                    if(save) saveRDS(results_LBP, file =paste0(out_dir, "/LBP_scenario",scenario,"_setting2.rds"))
                    
                    # Run the logistic beta, setting 3, NONblock, adapt, MH
                    results_LBP <- run_LBBer(data = data, model_pars = model_pars,
                                             nparticle = 0, adapt = TRUE, block = FALSE)
                    if(save) saveRDS(results_LBP, file = paste0(out_dir, "/LBP_scenario",scenario,"_setting3.rds"))
                    
                    # Run the logistic beta, setting 4, NONblock, nonadapt, MH
                    results_LBP <- run_LBBer(data = data, model_pars = model_pars,
                                             nparticle = 0, adapt = FALSE, block = FALSE)
                    if(save) saveRDS(results_LBP, file =paste0(out_dir, "/LBP_scenario",scenario,"_setting4.rds"))
                    
                    # Run the logistic beta, setting 5, block, adapt, PARTICLE
                    results_LBP <- run_LBBer(data = data, model_pars = model_pars,
                                             nparticle = 10, adapt = TRUE, block = TRUE)
                    if(save) saveRDS(results_LBP, file = paste0(out_dir, "/LBP_scenario",scenario,"_setting5.rds"))
                    
                    # Run the logistic beta, setting 6, block, nonadapt, PARTICLE
                    results_LBP <- run_LBBer(data = data, model_pars = model_pars,
                                             nparticle = 10, adapt = FALSE, block = TRUE)
                    if(save) saveRDS(results_LBP, file =paste0(out_dir, "/LBP_scenario",scenario,"_setting6.rds"))
                    
                    # Run the logistic beta, setting 7, NONblock, adapt, PARTICLE
                    results_LBP <- run_LBBer(data = data, model_pars = model_pars,
                                             nparticle = 10, adapt = TRUE, block = FALSE)
                    if(save) saveRDS(results_LBP, file = paste0(out_dir, "/LBP_scenario",scenario,"_setting7.rds"))
                    
                    # Run the logistic beta, setting 8, NONblock, nonadapt, PARTICLE
                    results_LBP <- run_LBBer(data = data, model_pars = model_pars,
                                             nparticle = 10, adapt = FALSE, block = FALSE)
                    if(save) saveRDS(results_LBP, file =paste0(out_dir, "/LBP_scenario",scenario,"_setting8.rds"))
                    
                    # Run the copula
                    results_copula <- run_copula(data = data, stan_input = stan_input, model_pars = model_pars, seed = isim)
                    if(save) saveRDS(results_copula, file =paste0(out_dir, "/copula_scenario",scenario,".rds"))
                    
                    
                    isim
                  }
results

# end with saving all objects in RDS. Make sure save = TRUE

