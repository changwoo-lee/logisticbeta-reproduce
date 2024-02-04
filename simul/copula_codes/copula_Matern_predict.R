# This function samples from the posterior predictive distribution of the
# gaussian copula with beta marginals
copula_Matern_predict <- function(n, ntest, D, theta_samples, rangepar_samples, a = 1, b){
  nsamples <- nrow(theta_samples)
  pb <- txtProgressBar(style=3)
  out <- matrix(NA, nsamples, ntest)
  for(i in 1:nsamples){
    # Check progress
    setTxtProgressBar(pb, i/(nsamples))
    # Parameter
    theta <- theta_samples[i, ]
    rangepar <- rangepar_samples[i]
    # Retrieve w from Stan
    w <- qnorm(1 - (1 - theta)^b)
    # Matern covariance
    K <- exp(-D / rangepar) * (1 + D / rangepar)
    # Calculate the posterior of the kernel
    Kxx_inv <- solve(K[1:n, 1:n])
    Kxstar_x <- K[(n + 1):(n + ntest), 1:n]
    Kxstar_xstar <- K[(n + 1):(n + ntest), (n + 1):(n + ntest)]
    # Predict
    M <- Kxstar_x %*% Kxx_inv
    m_pred <- M %*% w
    cov_pred <- Kxstar_xstar -  tcrossprod(M, Kxstar_x)
    # Adjust for symmetry
    cov_pred[lower.tri(cov_pred)] <- t(cov_pred)[lower.tri(cov_pred)]
    # Sample from the multivariate normal
    w_test <-as.numeric(mvnfast::rmvn(1, mu = m_pred, sigma = cov_pred))
    out[i, ] <- 1 - (1 - pnorm(w_test))^(1/b)
  }
  return(out)
}

