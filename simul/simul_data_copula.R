library(fields)

simulate_data_copula <- function(n, ntest, seed, scenario){

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

  # generate true_probs from gaussian copula with matern
  zeta = as.numeric(mvnfast::rmvn(1, mu = rep(0, n + ntest),sigma = Sig))
  true_probs_all = qbeta(pnorm(zeta), 1, 2) # a = 1, b = 2

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

