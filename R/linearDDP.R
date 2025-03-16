library(BayesLogit)
library(fields)
library(mvnfast)

#' Logistic-beta DDP mixture model with linearly dependent atoms
#' based on normalized feature map kernel for LBP
#' Using blocked Gibbs, adapted polya proposal, independent M-H 
#'
#' @param y response
#' @param X design matrix for atoms
#' @param b DP concentration parameter
#' @param H max number of components, upper bound on the number of clusters
#' @param prior list of prior parameters, prior$mu_beta, prior$Sigma_beta, (normal priors on coefficient in atom processes) prior$a_tau, prior$b_tau, (IG priors on variacne in atom processes)
#' @param nburn number of burn-in iterations 
#' @param nsave number of saved posterior samples. Total number of MCMC iteration is \code{nburn + nsave * nthin}.
#' @param nthin thin-in rate
#'
#' @return
#' @export
#'
#' @examples
linearDDP = function(y, X, b = 1, H = 20, prior = NULL, nburn = 100, nsave = 1000, nthin = 1)
{
  
  n = length(y)
  p = ncol(X)
#  q = ncol(Phi)
  #if(any(abs(rowSums(Phi^2)-1)>1e-10)) stop("Phi must have row normalized")
  
  # logistic-beta parameter
  a = 1
  b = b # DP concentration
  
  # priors
  if(is.null(prior$mu_beta)) mu_beta = rep(0, p) else mu_beta = prior$mu_beta
  if(is.null(prior$Sigma_beta)) Sigma_beta = diag(100, nrow = p, ncol = p) else Sigma_beta = prior$Sigma_beta
  Sigma_beta_inv = solve(Sigma_beta)
  
  # gamma prior of tau
  if(is.null(prior$a_tau)) a_tau = 1 else a_tau = prior$a_tau
  if(is.null(prior$b_tau)) b_tau = 1 else b_tau = prior$b_tau
  
  # initialize 
  r = rep(1,n)
  v = numeric(H-1)
  beta = matrix(0, nrow = p, ncol = H)
  tau = rep(1,H)
  #alpha = matrix(0, nrow = q, ncol = H)
  lambda = rep(1,H)
  lambda_acc = numeric(H)
  running_avg_lambda = rep(0, H)
  
  # saving objects 
  beta_save = array(0, dim = c(nsave, p, H))
  tau_save = array(0, dim = c(nsave, H))
  #alpha_save = array(rnorm(q*H), dim = c(nsave, q, H))
  weights_save = array(0, dim = c(nsave, H))
  #lambda_acc_save = array(0, dim = c(nsave, H))
  ncluster_save = array(0, dim = c(nsave))
  
  # Run MCMC       
  pb <- txtProgressBar(style=3)
  isave = 1
  ##### Start of MCMC #####
  t_start = Sys.time()
  for(imcmc in 1:(nburn + nsave*nthin)){
    
    setTxtProgressBar(pb, imcmc/(nburn + nsave*nthin))
    
    # (Step 1 is pushed to the back)
    ##### Step 2 #####
    ns <- sapply(1:H, function(x, v) sum(v == x), v = r)
    
    
    v[1:(H-1)] <- rbeta(H-1, 1 + ns[1:(H-1)], b + rev(cumsum(rev(ns[-1]))))
    sticks_all = matrix(v, nrow = n, ncol = H-1, byrow = TRUE)
    # 
     for (h in 1:H) {
    #   if (h < H) { # except last component H
    #     # Subsetting observations
    #     index <- (r > h - 1)
    #     nhbar <- sum(index)
    #     zh <- (r[index] == h) # act as binary response
    #     Phih <- matrix(Phi[index,,drop =F], nhbar, q)
    #     
    #     ###### Step 2-1: update Polya-gamma ######
    #     
    #     if(nhbar > 0){
    #       linh <- as.numeric(0.5*lambda[h]*(a-b) + sqrt(lambda[h])*Phih %*% alpha[,h])
    #       omegah <- BayesLogit::rpg.devroye(num = nhbar, h = 1, z = linh) # length nh
    #     }
    #     
    #     ###### Step 2-2: update Polya mixing variable ######
    #     
    #     # collapsed Gibbs: update lambda_h when alpha_h is marginalized out
    #     if(nhbar > 0){
    #       # indep MH, adapted proposal
    #       # recursive calculation of average
    #       running_avg_lambda[h] = running_avg_lambda[h] + (lambda[h] - running_avg_lambda[h])/imcmc
    #       abprime = find_polya_proposal(a, b, lambda_curr = running_avg_lambda[h])
    #       lambda_new = rpolya(1, abprime[1], abprime[2])
    #       logprob_new =  dmvn_lowrankstr(y = (zh - 0.5)/omegah,
    #                                      mu = rep(0.5*lambda_new*(a-b), length(omegah)),
    #                                      X = Phih,
    #                                      Dinv = diag(1/lambda_new, ncol(Phih)),
    #                                      Ediag = 1/omegah, log = T)
    #       
    #       logprob_old =  dmvn_lowrankstr(y = (zh - 0.5)/omegah,
    #                                      mu = rep(0.5*lambda[h]*(a-b), length(omegah)),
    #                                      X = Phih,
    #                                      Dinv = diag(1/lambda[h], ncol(Phih)),
    #                                      Ediag = 1/omegah, log = T)
    #       
    #       logratio = logprob_new - logprob_old + lambda[h]*(a*b-prod(abprime))/2 - lambda_new*(a*b-prod(abprime))/2
    #       # accept or reject in metropolis step
    #       if(log(runif(1))<logratio){
    #         lambda[h] = lambda_new
    #         lambda_acc[h] = T
    #       }else{
    #         lambda_acc[h] = F
    #       }
    #     }else{
    #       lambda[h] = rpolya(1, a, b) # draw from prior
    #       lambda_acc[h] = NA
    #     }
    #     
    #     ###### Step 2-3: update logistic-beta ######
    #     
    #     if(nhbar > 0){
    #       # eta_h | lambda = 0.5lambda(a-b) + sqrt(lambda) Phih alpha_h
    #       y_surrogate = (zh - 0.5)/omegah - 0.5*lambda[h]*(a-b)
    #       X_surrogate = sqrt(lambda[h])*Phih
    #       # updating from as if y_surrogate = X_surrogate*alpha + ePhilon, epslion ~ N(0, Omega^{-1})
    #       
    #       # prior on alpha_h: normal with mean zero, diagonal covariance 11
    #       Qtemp = crossprod(X_surrogate*sqrt(omegah)) + diag(1, nrow = q, ncol = q)
    #       btemp = crossprod(X_surrogate*omegah, y_surrogate)
    #       alpha[,h] = spam::rmvnorm.canonical(1, btemp, Qtemp)
    #     }else{
    #       alpha[,h] = rnorm(q) # draw from prior
    #     }
    #   }
      
      ##### Step 3: update component-specific parameters (atom process) #####
  
      yh = y[r==h]
      nh = length(yh)
      Xh = X[r==h,, drop = F]
      
      Qtemp = Sigma_beta_inv + crossprod(Xh)*tau[h] # when nh = 0, this is Sigma0_inv
      btemp = Sigma_beta_inv%*%mu_beta + crossprod(Xh,yh)*tau[h] # when nh = 0, this is prior
      beta[,h] = spam::rmvnorm.canonical(1, btemp, Qtemp)
      
      residual <- as.numeric(yh - Xh%*%beta[,h]) # when nh = 0, this is numeric(0)
      tau[h] <- rgamma(1, a_tau + nh / 2, b_tau + sum(residual^2) / 2)
      
    }
    #browser()
    ##### Step 1: update component allocation #####
    
    # sticks_all = matrix(0, nrow = n, ncol = H-1)
    # for(h in 1:(H-1)){
    #   etah = 0.5*lambda[h]*(a-b) + sqrt(lambda[h])*Phi%*%alpha[,h]
    #   sticks_all[,h] = 1/(1+exp(-etah))
    # }
    #sticks_all = 1/(1+exp(-t(0.5*lambda[-H]*(a-b)-sqrt(lambda[-H])*t(Phi%*%alpha[,-H])))) # V_h(x_i), h=1,...H for each column
    mixture_weights_all = sticks_to_probs(sticks_all)
    
    sdmat = matrix(1/sqrt(tau), nrow = n, ncol = H, byrow = T)
    mumat = X%*%beta
    ymat = matrix(y, nrow =n, ncol = H, byrow = F)
    
    select_probs_all = mixture_weights_all*matrix(dnorm(ymat, mumat, sdmat), nrow=n, ncol = H)
    # https://stackoverflow.com/questions/20508658/sampling-repeatedly-with-different-probability
    r = sample.rowwise(select_probs_all)
    #r = apply(select_probs_all, 1, function(x) sample.int(H, size = 1, prob = x))
    #browser()
    # save
    if((imcmc > nburn)&&((imcmc-nburn)%%nthin==0)){
      #browser()
      beta_save[isave,,] = beta
      tau_save[isave,] = tau
      #alpha_save[isave,,] = alpha
      weights_save[isave,] = mixture_weights_all[1,]
      #lambda_acc_save[isave,] = lambda_acc
       
      ncluster_save[isave] = length(unique(r))
      isave = isave + 1
    }
  }
  t_mcmc = Sys.time() - t_start
  out = list()
  
  # input
  out$b = b
  out$H = H
  
  out$beta_save = beta_save
  out$tau_save = tau_save
#  out$alpha_save = alpha_save
  out$weights_save = weights_save
 # out$lambda_acc_save = lambda_acc_save
  out$ncluster_save = ncluster_save
  out$nsave = nsave
  
  out$t_mcmc = t_mcmc
  return(out)
}


#' Transform probability vector to stick weights
#'
#' @param p probability vector, or matrix with each row being probability vector
#' @param drop_lastcol default T, drop last column of 1's?
#'
#' @return matrix of sticks weights (with one less column by default), each row corresponds to stick weights
#' @export
#' @importFrom matrixStats rowCumsums
#'
#' @examples
#'
#' probs_to_sticks(c(0.5, 0.3, 0.2))
#' probmat = matrix(c(.5, .3, .2, .6, .2, .2), ncol = 3, byrow = T)
#' probs_to_sticks(probmat)
#'
probs_to_sticks <- function(p, drop_lastcol = T){
  if(is.vector(p)) p = matrix(p, nrow = 1, ncol = length(p))
  sticks = p
  sticks[,-1] = p[,-1]/(1-matrixStats::rowCumsums(p[,-ncol(p),drop = F]))
  sticks[which(is.nan(sticks))] <- 0 # when 0/0
  if(drop_lastcol){
    return(sticks[,-ncol(sticks), drop = F])
  }else{
    return(sticks)
  }
}

#' Transform stick weights to probability vector
#'
#' Assuming last element (of the column) does not have 1, which is dummy variable for stick weights
#'
#' @param sticks stick vector, or matrix with each row being stick weights.
#'
#' @return matrix of probabilities, each row corresponds to probabilities.
#' @importFrom matrixStats rowCumprods
#' @export
#'
#' @examples
#'
#' sticks_to_probs(c(0.1, 0.2, 0.3)) # output: length 4 probability vector
#' stickmat = matrix(c(.1, .2, .3, .4, .5, .6), ncol = 3, byrow = T)
#' sticks_to_probs(stickmat)
#'
sticks_to_probs <- function(sticks){
  if(is.vector(sticks)) sticks = matrix(sticks, nrow = 1, ncol = length(sticks))
  #if(any(abs(sticks[,ncol(sticks)] - 1) < 1e-12)) stop("last column has dummy 1, please remove")
  p = sticks
  p[,-1] = sticks[,-1]*matrixStats::rowCumprods(1-sticks[,-ncol(p), drop = F])
  p = cbind(p, 1-rowSums(p))
  return(p)
}


# prob_to_sticks <- function(p){
#   if(is.matrix(p)) warning("p should be numeric, not a vector")
#   sticks = p
#   sticks[-1] = p[-1]/(1-cumsum(p)[-length(p)])
#   sticks[which(is.nan(sticks))] <- 0 # when 0/0
#   round(sticks, digits = 16)
#   sticks
# }


#
# ptilde_to_p <- function(ptilde){
#   p = ptilde
#   p[-1] = ptilde[-1]*cumprod(1-ptilde[-length(p)])
#   p
# }

#' Vectorized sampling (by row by row) with probability stored in rows of matrix
#' source: https://stackoverflow.com/questions/20508658/sampling-repeatedly-with-different-probability
#'
#' @param Weight n by K matrix, each row corresponds to unnormalized probabilities of sampling weights
#'
#' @return integers with length n, each element in {1,...,K}
#' @export
#'
#' @examples
sample.rowwise <- function(Weight) {
  x <- runif(nrow(Weight))
  cumul.w <- Weight %*% upper.tri(diag(ncol(Weight)), diag = TRUE)/rowSums(Weight)
  i <- rowSums(x > cumul.w) + 1L
  i
}
# #
# temp1 = dpois(1:100, lambda = 30)*10
# temp2 = dpois(1:100, lambda = 40)
# temp12 = rbind(temp1,temp2)
#
# nsave = 1000
# tempsave = matrix(0,nrow = nsave,ncol = 2)
# for(isave in 1:nsave) tempsave[isave,] =  sample.rowwise(temp12)
#
# hist(tempsave[,1]) # mean 30
# hist(tempsave[,2]) # mean 40
#
# var(tempsave[,1])
# var(tempsave[,2])
# #
# #




# for multiple y_eval
linearDDP_density <- function(y_eval, X_new, fit){
  
  n_eval = length(y_eval)
  X_new = as.numeric(X_new)
  nsave = fit$nsave
  H = fit$H
  
  Y_eval = matrix(y_eval, nrow = n_eval, ncol = H)
  pdf_mcmc = matrix(0, nrow = nsave, ncol = n_eval)
  
  condmean_mcmc = matrix(0, nrow = nsave, ncol = 1)
  for(isave in 1:nsave){
    mixture_weights = matrix(fit$weights_save[isave,], nrow = n_eval, ncol = H, byrow = TRUE)
    #browser()
    pdf_mcmc[isave, ] = rowSums(mixture_weights*
                                  t(dnorm(t(Y_eval), t(X_new%*%fit$beta_save[isave,,]) ,sqrt(1/fit$tau_save[isave,]))))
    condmean_mcmc[isave] = sum(mixture_weights[1,]*as.numeric(X_new%*%fit$beta_save[isave,,]))
  }
  
  out = list(pdf_mcmc = pdf_mcmc, condmean_mcmc = condmean_mcmc)
  return(out)
}

