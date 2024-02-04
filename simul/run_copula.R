library(scoringRules)
library(rstan)
library(coda)
library(mcmcse)

run_copula <- function(data, stan_input, model_pars, seed){
  results <- list()
  fit <- stan(file = "simul/copula_codes/copula_Matern.stan",
              data = stan_input,
              warmup = model_pars$nburn, # to see convergence behavior
              iter = model_pars$nsamples + model_pars$nburn, chains = 1, par = pars_save, seed = seed)
  results$fit <- fit
  #------------------------------------------------- Step 2 - Predict in- and out-of-sample
  print("\n Predict")
  pars <- rstan::extract(fit)
  # In-sample
  theta_copula <- colMeans(pars$theta)
  MAE_copula <- round(mean(abs(theta_copula - data$true_probs)), 6)
  RMSE_copula <- round(sqrt(mean((theta_copula - data$true_probs)^2)),6)
  crps_copula = scoringRules::crps_sample(data$true_probs, t(pars$theta))
  
  # Out-of-sample
  theta_pred_copula <- copula_Matern_predict(n = length(data$z),
                                      ntest = length(data$z_test),
                                      D = fields::rdist(rbind(data$coords, data$coords_test)),
                                      theta_samples = pars$theta,
                                      rangepar_samples = pars$rangepar,
                                      b = stan_input$b,
                                      a = stan_input$a)
  
  prob_test_copula <- colMeans(theta_pred_copula)
  MAE_test_copula <- round(mean(abs(prob_test_copula - data$true_probs_test)), 6)
  RMSE_test_copula <- round(sqrt(mean((prob_test_copula - data$true_probs_test)^2)),6)
  crps_test_copula = scoringRules::crps_sample(data$true_probs_test, t(theta_pred_copula))
  
  # Save the result
  results$predict <- theta_pred_copula
  #------------------------------------------------- Step 3 - calculate performance metrics
  # Effective sample sizes
  avg_ess <- mean(coda::effectiveSize(pars$theta))
  mult_ess <- try(mcmcse::multiESS(pars$theta), silent = TRUE)
  
  # time (in seconds)
  time_sampling <- sum(get_elapsed_time(fit))
  
  # Save the result
  results$prediction_metrics <- c("MAE" = MAE_copula,
                                  "RMSE" = RMSE_copula,
                                  "crps" = mean(crps_copula),
                                  "MAE_test" = MAE_test_copula,
                                  "RMSE_test" = RMSE_test_copula,
                                  "crps_test" = mean(crps_test_copula),
                                  "Ess" = avg_ess,
                                  "mult_ess" = mult_ess,
                                  "time_sampling" = time_sampling)
  return(results)
}
