# postprocessing codes for LB-DDP mixture model
# fit is LB-DDP output
# X, Phi are at new evaluation points

LBDDP_density_old <- function(y_eval_density, X, Phi, fit){

  nX = nrow(X)
  nsave = fit$nsave
  H = fit$H
  a = 1
  b = fit$b

  pdf_mcmc = matrix(0, nrow = nsave, ncol = nX)

  for(i in 1:nX){
    for(isave in 1:nsave){
      sticks = 1/(1+exp(-0.5*fit$lambda_save[isave,-H]*(a-b)-sqrt(fit$lambda_save[isave,-H])*Phi[i,]%*%fit$alpha_save[isave,,-H])) # V_h(x_i), h=1,...H for each column
      mixture_weights = as.numeric(sticks_to_probs(sticks))
#      browser()
      pdf_mcmc[isave, i] = sum(mixture_weights*dnorm(y_eval_density, as.numeric(X[i,]%*%fit$beta_save[isave,,]) ,sqrt(1/fit$tau_save[isave,])))
    }
  }
  return(pdf_mcmc)

}
# for multiple y_eval
LBDDP_density <- function(y_eval, X_new, Phi_new, fit){
  
  n_eval = length(y_eval)
  X_new = as.numeric(X_new)
  Phi_new = as.numeric(Phi_new)
  nsave = fit$nsave
  H = fit$H
  a = 1
  b = fit$b
  
  Y_eval = matrix(y_eval, nrow = n_eval, ncol = H)
  pdf_mcmc = matrix(0, nrow = nsave, ncol = n_eval)
  condmean_mcmc = matrix(0, nrow = nsave, ncol = 1)
  
  
  for(isave in 1:nsave){
    sticks = 1/(1+exp(-0.5*fit$lambda_save[isave,-H]*(a-b)-sqrt(fit$lambda_save[isave,-H])*Phi_new%*%fit$alpha_save[isave,,-H])) # V_h(x_i), h=1,...H for each column
    mixture_weights = matrix(sticks_to_probs(sticks), nrow = n_eval, ncol = H, byrow = TRUE)
    pdf_mcmc[isave, ] = rowSums(mixture_weights*
                                  t(dnorm(t(Y_eval), t(X_new%*%fit$beta_save[isave,,]) ,sqrt(1/fit$tau_save[isave,]))))
    condmean_mcmc[isave] = sum(mixture_weights[1,]*as.numeric(X_new%*%fit$beta_save[isave,,]))
  }
  out = list(pdf_mcmc = pdf_mcmc, condmean_mcmc = condmean_mcmc)
  return(out)
}








LBDDP_cdf <- function(threshold, X, Phi, fit){

  nX = nrow(X)
  nsave = fit$nsave
  H = fit$H
  a = 1
  b = fit$b

  cdf_mcmc = matrix(0, nrow = nsave, ncol = nX)

  for(i in 1:nX){
    for(isave in 1:nsave){
      sticks = 1/(1+exp(-0.5*fit$lambda_save[isave,-H]*(a-b)-sqrt(fit$lambda_save[isave,-H])*Phi[i,]%*%fit$alpha_save[isave,,-H])) # V_h(x_i), h=1,...H for each column
      mixture_weights = as.numeric(sticks_to_probs(sticks))
      #      browser()
      cdf_mcmc[isave, i] = sum(mixture_weights*pnorm(threshold, as.numeric(X[i,]%*%fit$beta_save[isave,,]) ,sqrt(1/fit$tau_save[isave,])))
    }
  }
  return(cdf_mcmc)

}





#
#
# lbddp_density_old <- function(y_eval_density, Xnew, Phinew, beta, alpha, lambda, tau, M){
#
#   nX = nrow(Xnew)
#   H = ncol(Phinew)
#   a = 1
#   b = M
#
#   pdf = numeric(nX)
#   for(i in 1:nX){ # for each x
#     for(h in 1:(H-1)){
#       eta = 0.5*lambda[h]*(a-b) + sqrt(lambda[h])*Phinew[i,]%*%alpha[,h]
#     }
#     mixture_weights = sticks_to_probs(1/(1+exp(-eta))) # length h prob
#     means = Xnew[i,]%*%beta
#     sds = sqrt(1/tau)
#     browser()
#     pdf[i] = sum(mixture_weights*dnorm(y_eval_density, means, sds))
#   }
# }
#
#
#
