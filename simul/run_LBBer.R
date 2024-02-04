library(scoringRules)
library(coda)
library(mcmcse)

run_LBBer <- function(data, model_pars, nparticle, adapt, block){
  results <- list()
  results$nparticle = nparticle
  results$adapt = adapt
  results$block = block
  #------------------------------------------------- Step 1 - Fit the model
  if(block){
    fit <- LBBer_Matern(z = data$z,
                        x = data$coords,
                        a = model_pars$a,
                        b = model_pars$b,
                        range_grid = seq(model_pars$range_lb,
                                         model_pars$range_ub,
                                         length = model_pars$ngrid),
                        smoothness = 1.5,
                        Dxx = data$Dxx,
                        nburn = model_pars$nburn,
                        nsave = model_pars$nsamples,
                        nthin = 1,
                        nparticle = nparticle,
                        adapt = adapt)
  } else {
    fit <- LBBer_Matern_unblocked(z = data$z,
                                  x = data$coords,
                                  a = model_pars$a,
                                  b = model_pars$b,
                                  range_grid = seq(model_pars$range_lb,
                                                   model_pars$range_ub,
                                                   length = model_pars$ngrid),
                                  smoothness = 1.5,
                                  Dxx = data$Dxx,
                                  nburn = model_pars$nburn,
                                  nsave = model_pars$nsamples,
                                  nthin = 1,
                                  nparticle = nparticle,
                                  adapt = adapt)
  }
  # Save the result
  results$fit <- fit
  
  #------------------------------------------------- Step 2 - Predict in- and out-of-sample
  
  print("\n Predict")
  # In-sample
  theta_LBP <- colMeans(1/(1+exp(-fit$eta_save)))
  MAE_LBP <- round(mean(abs(theta_LBP - data$true_probs)), 6)
  RMSE_LBP <- round(sqrt(mean((theta_LBP - data$true_probs)^2)),6)
  crps_LB = scoringRules::crps_sample(data$true_probs, t(1/(1+exp(-fit$eta_save))))
  
  # Out-of-sample
  theta_pred_LBP <- LBBer_Matern_predict(fit = fit, xpred = data$coords_test)
  
  prob_test_LBP <- colMeans(theta_pred_LBP)
  MAE_test_LBP <- round(mean(abs(prob_test_LBP - data$true_probs_test)), 6)
  RMSE_test_LBP <- round(sqrt(mean((prob_test_LBP - data$true_probs_test)^2)),6)
  crps_test_LBP = scoringRules::crps_sample(data$true_probs_test, t(theta_pred_LBP))
  
  # Save the result
  results$predict <- theta_pred_LBP
  
  #------------------------------------------------- Step 3 - calculate performance metrics
  # Effective sample sizes
  avg_ess <- mean(coda::effectiveSize(1/(1+exp(-fit$eta_save))))
  mult_ess <- try(mcmcse::multiESS(1/(1+exp(-fit$eta_save))), silent = TRUE)
  
  # time (in seconds)
  time_sampling <- as.numeric(fit$t_mcmc)
  
  # Save the result
  results$prediction_metrics <- c("MAE" = MAE_LBP,
                                  "RMSE" = RMSE_LBP,
                                  "crps" = mean(crps_LB),
                                  "MAE_test" = MAE_test_LBP,
                                  "RMSE_test" = RMSE_test_LBP,
                                  "crps_test" = mean(crps_test_LBP),
                                  "Ess" = avg_ess,
                                  "mult_ess" = mult_ess,
                                  "time_sampling" = time_sampling)
  return(results)
}

