# fit: output from LBBer_Matern
# xpred: npred by d matrix (d = 1 or 2) of new locations 
LBBer_Matern_predict <- function(fit, xpred) {
  x <- fit$x
  y <- as.matrix(xpred)
  nsave <- fit$nsave
  eta_save <- fit$eta_save
  lambda_save <- fit$lambda_save
  eta_pred_save <- matrix(0, nrow = nsave, ncol = nrow(y))
  dist_x <- fields::rdist(x)
  dist_x_y <- fields::rdist(x, y)
  dist_y <- fields::rdist(y)
  # Progress bar
  pb <- txtProgressBar(style = 3)
  for (isave in 1:nsave) {
    # Check progress
    setTxtProgressBar(pb, isave / (nsave))
    # Calculate kernels
    Rxx <- fields::Matern(dist_x, range = fit$range[isave], smoothness = fit$smoothness)
    Rxy <- fields::Matern(dist_x_y, range = fit$range[isave], smoothness = fit$smoothness)
    Ryy <- fields::Matern(dist_y, range = fit$range[isave], smoothness = fit$smoothness)
    Rxxinv_Rxy <- solve(Rxx, Rxy)
    # Prediction
    eta <- eta_save[isave, ]
    lambda <- lambda_save[isave]
    temp_mean <- t(Rxxinv_Rxy) %*% eta
    temp_var <- lambda * (Ryy - t(Rxxinv_Rxy) %*% Rxy)
    eta_pred_save[isave, ] <- mvnfast::rmvn(1, mu = temp_mean, sigma = temp_var)
  }
  theta <- 1 / (1 + exp(-eta_pred_save))
  return(theta)
}