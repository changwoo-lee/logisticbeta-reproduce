library(fields)
library(dqrng)

rmlb <- function(n, a, b, R, kmax = 200){
  p = ncol(R)

  lambda = rpolya(n, a, b, kmax = kmax)
  theta = (a - b)/2 # scalar

  #x = mvnfast::rmvn(n, mu = rep(0,p), sigma = R)
  x = matrix(rnorm(n*p), n, p)%*%chol(R)
  #rmvn(n,lambda*theta, sqrt(lambda)*t(Sigma*sqrt(lambda)))
  return(lambda*theta + sqrt(lambda)*x)
}

rpolya = function(n, a, b, kmax = 200){
  k = 0:kmax # start from 0
  denom = (a+k)*(b+k)/2 # length kmax + 1
  #temp = matrix(rexp(n*(kmax + 1))/rep(denom, n), nrow = kmax + 1, ncol = n, byrow = F)
  temp = matrix(dqrng::dqrexp(n*(kmax + 1))/rep(denom, n), nrow = kmax + 1, ncol = n, byrow = F)
  colSums(temp)
}

simulate_data_lbp <- function(n, ntest, seed, scenario){

  set.seed(seed)
  # Sample the coordinates at random
  coords <- cbind(runif(n), runif(n)) # spatial locations
  coords_test <- cbind(runif(ntest), runif(ntest))
  # Pre-calculate the distances for the Matern covariances
  D <- fields::rdist(rbind(coords, coords_test)) # (n + ntest) x (n + ntest)
  Dxx <- D[1:n, 1:n] # n x n

  if (scenario == 1){
    Sig = fields::Matern(D, smoothness = 1.5, range = 0.1)
  } else if (scenario == 2){
    Sig = fields::Matern(D, smoothness = 1.5, range = 0.2)
  } else if (scenario == 3){
    Sig = fields::Matern(D, smoothness = 1.5, range = 0.4)
  }

  # generate true_probs from lbp with matern
  # here a=1, b=2
  eta = rmlb(1, a=1, b=2, Sig)
  true_probs_all = as.numeric(1/(1+exp(-eta)))# a = 1, b = 2
  #hist(true_probs_all, xlim = c(0,1))
  true_probs <- true_probs_all[1:n]
  true_probs_test <- true_probs_all[(n+1):(n+ntest)]

  z <- rbinom(n, size = 1, true_probs)
  z_test <- rbinom(ntest, size = 1, true_probs_test)

  # Save the data
  data <- list(z = z,
               z_test = z_test,
               true_probs = true_probs,
               true_probs_test = true_probs_test,
               coords = coords,
               coords_test = coords_test,
               Dxx = Dxx)
  return(data)
}

