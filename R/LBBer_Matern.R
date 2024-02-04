library(BayesLogit)
library(fields)
library(mvnfast)

#' Latent LBP model for binary data with Matern correlation
#' Implements Algorithm 1, blocked Gibbs sampler
#'
#' @param z n by 1 binary data
#' @param x n by d matrix, spatial coordinates, d = 1 or 2
#' @param a LBP shape parameter 1, fixed
#' @param b LBP shape parameter 2, fixed
#' @param range_grid Matern range parameter, grid for discrete uniform prior
#' @param smoothness Matern smoothness parameter, fixed
#' @param Dxx optional, n by n matrix, distance between data (x)
#' @param nburn number of burn-in iterations 
#' @param nsave number of saved posterior samples. Total number of MCMC iteration is \code{nburn + nsave * nthin}.
#' @param nthin thin-in rate
#' @param nparticle number of particles. Default 0, which corresponds to independent M-H
#' @param adapt whether to use adaptive proposal. Default TRUE
#' @param verbose show progress. Default TRUE
#'
#' @export
#'
#' @examples
LBBer_Matern <- function(z, x, a, b, range_grid, smoothness = 1.5, Dxx =NULL,
                             nburn = 100, nsave = 1000, nthin = 1, nparticle = 0, adapt = T, verbose=TRUE){
  if(nparticle==0) indepMH = TRUE else indepMH = FALSE
  t_start = Sys.time()
  # sanity check for a and b which should be positive reals
  if(a <= 0 || b <= 0) stop("a and b must be positive reals")
  # sanity check for range and smoothness which should be positive reals
  if(smoothness <= 0) stop("smoothness must be positive reals")
  # sanity check for z which should be a vector of binary
  if(!is.vector(z) || !all(z %in% c(0,1))) stop("z must be a vector of binary")
  
  
  x = as.matrix(x)
  n = nrow(x)
  d = ncol(x)
  if(is.null(Dxx)) Dxx = fields::rdist(x)
  
  running_avg_lambda = 0
  ngrid = length(range_grid)
  Rxx_grid = list()
  Rxxinv_grid = list()
  Rxxinvchol_grid = list()
  Rxxlogdet_grid = numeric(ngrid)
  
  igrid = 0
  print("Pre-MCMC calculation... ")
  for(range in range_grid){
    igrid = igrid + 1
    Rxx_grid[[igrid]] = fields::Matern(Dxx, range = range, smoothness = smoothness) # save
    Rxxinv_grid[[igrid]] = solve(Rxx_grid[[igrid]]) # save
    Rxxinvchol_grid[[igrid]] = chol(Rxxinv_grid[[igrid]])
    Rxxlogdet_grid[igrid] =as.numeric(determinant(Rxx_grid[[igrid]], log = T)$modulus)
  }
  # ngrid*n times n block matrix
  # each n by n block is cholesky of Rxxinv_grid[[igrid]]
  # for fast update of range parameter
  sqrt_quadform = matrix(0,ngrid*n, n)
  for(igrid in 1:ngrid){
    sqrt_quadform[((igrid-1)*n+1):(igrid*n),] = Rxxinvchol_grid[[igrid]]
  }
  
  # Initialize
  range  = range_grid[ceiling(ngrid/2)]
  Rxx = fields::Matern(Dxx, range = range, smoothness = smoothness)
  Rxxinv = Rxxinv_grid[[ceiling(ngrid/2)]]
  
  lambda = epolya(a,b)
  eta = as.numeric(spam::rmvnorm(1, mu = rep(0.5*(a-b), n), Sigma = lambda*Rxx))
  
  # Saving objects
  range_save = array(0, dim = c(nsave))
  eta_save = array(0, dim = c(nsave, n))
  eta_pred_save = array(0, dim = c(nsave, ngrid))
  lambda_save = array(0, dim = c(nsave))
  lambda_save_all = array(NA, dim = c(nburn + nsave*nthin)) # for adaptation
  lambda_acc_save = array(0, dim = c(nsave))
  
  t_end = Sys.time()
  t_premcmc = difftime(t_end, t_start, units = "secs")
  
  # Run MCMC
  if(verbose){
    pb <- txtProgressBar(style=3)
  }
  isave = 1
  ##### Start of MCMC #####
  t_start = Sys.time()
  for(imcmc in 1:(nburn + nsave*nthin)){
    if(verbose){
      setTxtProgressBar(pb, imcmc/(nburn + nsave*nthin))
    }
    
    ##### Step 1: update omega, Polya-Gamma auxiliary variable #####
    
    omega = BayesLogit::rpg.devroye(n, h = 1, z = eta)
    
    ##### Step 2: update lambda (Polya mixing parameter) #####
    
    # Adaptive polya proposal?
    if(adapt){
      # running average of lambda adaptive MCMC
      running_avg_lambda = running_avg_lambda + (lambda - running_avg_lambda)/imcmc
      abprime = find_polya_proposal(a, b, lambda_curr = running_avg_lambda)
    }else{
      abprime = c(a,b)
    }
    
    # Independent M-H or particle Gibbs
    if(indepMH){
      # proposal
      lambda_new = rpolya(1, abprime[1], abprime[2])
      
      logprob_new = mvnfast::dmvn((z-0.5)/omega,
                                  mu = rep(0.5*lambda_new*(a-b), n),
                                  sigma = diag(1/omega)+lambda_new*Rxx, log = T)
      logprob_old = mvnfast::dmvn((z-0.5)/omega,
                                  mu = rep(0.5*lambda*(a-b), n),
                                  sigma = diag(1/omega)+lambda*Rxx, log = T)
      logratio = logprob_new - logprob_old + (lambda - lambda_new)*(a*b - prod(abprime))/2
      # accept or reject in metropolis step
      if(log(runif(1)) < logratio){
        lambda = lambda_new; lambda_acc = T;
      }else{
        lambda_acc = F
      }
    }else{ # particle Gibbs
      # proposal
      lambda_particles = c(rpolya(nparticle, abprime[1], abprime[2]),lambda)
      
      # collapsed likelihood evaluations
      logweights = numeric(nparticle+1)
      for(iparticle in 1:(nparticle+1)){
        logweights[iparticle] = mvnfast::dmvn((z-0.5)/omega,
                                              mu = rep(0.5*lambda_particles[iparticle]*(a-b), n),
                                              sigma = diag(1/omega)+lambda_particles[iparticle]*Rxx,
                                              log = T)
      }
      
      logweights = logweights - lambda_particles*(a*b-prod(abprime))/2
      
      logweights = exp(logweights - matrixStats::logSumExp(logweights)) # normalize
      idx = sample.int(nparticle + 1, size = 1, prob = logweights)
      # update lambda
      lambda = lambda_particles[idx]
      if(idx <= nparticle)  lambda_acc = T else lambda_acc =F
    }
    
    ##### Step 3: update eta (logistic-beta)   #####
    
    Q_eta = diag(omega) + 1/lambda*Rxxinv
    b_eta = (z - 1/2) + Rxxinv%*%rep(0.5*(a-b), n)
    eta = as.numeric(spam::rmvnorm.canonical(1, b_eta, Q_eta))
    
    ##### Step 4: update range parameter  #####
    
    # eta | lambda ~ N(0.5 lambda(a-b) , lambda*Rxx)
    # normal density evaluation
    loglik = numeric(ngrid)
    quadform = colSums(matrix(sqrt_quadform%*%(eta - 0.5*lambda*(a-b)), nrow = n, ncol = ngrid)^2)
    loglik = -0.5*Rxxlogdet_grid - 0.5*(1/lambda)*quadform
    # sample
    grididx = sample.int(ngrid, size = 1, prob = exp(loglik - max(loglik)))
    range = range_grid[grididx]
    
    # update
    Rxx = Rxx_grid[[grididx]]
    Rxxinv = Rxxinv_grid[[grididx]]
    Rxxlogdet = Rxxlogdet_grid[grididx]
    sqrt_quadform[((grididx-1)*n+1):(grididx*n),] = Rxxinvchol_grid[[grididx]]
    
    # save
    lambda_save_all[imcmc] = lambda
    if((imcmc > nburn)&&((imcmc-nburn)%%nthin==0)){
      eta_save[isave,] = eta
      range_save[isave] = range
      lambda_save[isave] = lambda
      lambda_acc_save[isave] = lambda_acc
      isave = isave + 1
    }
  }
  t_end = Sys.time()
  t_mcmc = difftime(t_end, t_start, units = "secs")
  
  out = list()
  
  out$x = x
  
  out$eta_save = eta_save
  out$lambda_save = lambda_save
  out$lambda_acc_save = lambda_acc_save
  out$nsave = nsave
  out$range_save = range_save
  out$smoothness = smoothness
  out$t_mcmc = t_mcmc
  out$t_premcmc = t_premcmc
  out$indepMH = indepMH
  
  return(out)
}
