# DDE exposure analysis with smoking group

rm(list = ls())

library(ggplot2)
library(patchwork)
library(reshape2)
library(coda)    # For MCMC analysis
library(splines) # For computing the natural B-splines basis
library(logitnorm) # prior property comparison


source('R/polya.R')
source('R/dmvn_lowrankstr.R')
# LBDDP
source("R/LBDDP_basis.R")
source("R/LSBP_basis.R")

######## 1. Prepare data ###############

#library(BNPmix)
#data("CPP") # from BNPmix package
# or alternatively from CRAN mirror repository
load(url("https://github.com/cran/BNPmix/raw/master/data/CPP.RData"))

y_original = CPP$gest
x_original = CPP$dde
# nonsmoking group
y_original = y_original[CPP$smoke==1]
x_original = x_original[CPP$smoke==1]
# scale data
y = scale(y_original)
x = scale(x_original)
xsd = sd(x_original)
xmean = mean(x_original)


n = length(y)
summary(x_original)
summary(y_original)




p = 2   # number of column of design matrix of linearly dependent atoms
nbasis = 6   # Number of basis functions for dependent weight

# normalized basis 
Basis = ns(x,nbasis, intercept = T)
Basis_normalized = Basis/sqrt(rowSums(Basis^2))

nsave = 30000
nburn = 5000
run = FALSE
if(run){

  ########## prior setting 1 ############################################
priorLBDDP       <- list(mu_beta = rep(0,p),
                         Sigma_beta = diag(1,p),
                         a_tau = 1, b_tau = 1) # b = 0.2 (direct input of function)

priorLSBP_set1       <- list(mu_beta = rep(0,p),
                        Sigma_beta = diag(1,p),
                        a_tau = 1, b_tau = 1, sigma2_alpha = 100) # sigma2_alpha = 100 <<<<

# comparison of prior co-clustering probability
# LB-DDP, b = 0.2
1/(0.2+1) # 0.8333
# LSBP, sigma2_alpha = 100
mu1 = logitnorm::momentsLogitnorm(0, sigma = 10)[1]
mu2 = logitnorm::momentsLogitnorm(0, sigma = 10)[2] + mu1^2 
mu2/(2*mu1 - mu2) # 0.854

set.seed(1)
fit_LBDDP_set1 = LBDDP_basis(y,
                             cbind(1, x),
                             Basis_normalized, b = 0.2, H = 20,
                             prior = priorLBDDP,
                             nburn = nburn, nsave = nsave, nthin = 1)

set.seed(11)
fit_LSBP_set1 = LSBP_basis(y,
                           cbind(1, x),
                           Basis_normalized, H = 20,
                           prior = priorLSBP_set1,
                           nburn = nburn, nsave = nsave, nthin = 1)

########## prior setting 2 ############################################
priorLBDDP       <- list(mu_beta = rep(0,p),
                         Sigma_beta = diag(1,p),
                         a_tau = 1, b_tau = 1) # b = 1 (direct input of function)

priorLSBP_set2       <- list(mu_beta = rep(0,p),
                             Sigma_beta = diag(1,p),
                             a_tau = 1, b_tau = 1, sigma2_alpha = pi^2/3) # sigma2_alpha = 100 <<<<
# comparison of prior co-clustering probability
# LB-DDP, b = 1
1/(1+1) # 0.5
# LSBP, sigma2_alpha = 100
mu1 = logitnorm::momentsLogitnorm(0, sigma = pi/sqrt(3))[1]
mu2 = logitnorm::momentsLogitnorm(0, sigma = pi/sqrt(3))[2] + mu1^2 
mu2/(2*mu1 - mu2) # 0.515

set.seed(2)
fit_LBDDP_set2 = LBDDP_basis(y,
                             cbind(1, x),
                             Basis_normalized, b = 1, H = 20,
                             prior = priorLBDDP,
                             nburn = nburn, nsave = nsave, nthin = 1)

set.seed(22)
fit_LSBP_set2 = LSBP_basis(y,
                           cbind(1, x),
                           Basis_normalized, H = 20,
                           prior = priorLSBP_set2,
                           nburn = nburn, nsave = nsave, nthin = 1)


########## prior setting 3 ############################################
priorLBDDP       <- list(mu_beta = rep(0,p),
                         Sigma_beta = diag(1,p),
                         a_tau = 1, b_tau = 1) # b = 2 (direct input of function)

priorLSBP_set3       <- list(mu_beta = rep(0,p),
                             Sigma_beta = diag(1,p),
                             a_tau = 1, b_tau = 1, sigma2_alpha = 0.2^2)
# comparison of prior co-clustering probability
# LB-DDP, b = 2
1/(2+1) # 0.333
# LSBP, sigma2_alpha = 100
mu1 = logitnorm::momentsLogitnorm(0, sigma = 0.2)[1]
mu2 = logitnorm::momentsLogitnorm(0, sigma = 0.2)[2] + mu1^2 
mu2/(2*mu1 - mu2) # 0.338

set.seed(3)
fit_LBDDP_set3 = LBDDP_basis(y,
                             cbind(1, x),
                             Basis_normalized, b = 2, H = 20,
                             prior = priorLBDDP,
                             nburn = nburn, nsave = nsave, nthin = 1)

set.seed(33)
fit_LSBP_set3 = LSBP_basis(y,
                           cbind(1, x),
                           Basis_normalized, H = 20,
                           prior = priorLSBP_set3,
                           nburn = nburn, nsave = nsave, nthin = 1)

save.image("DDEexposure_nonsmoking.RData")

} else {
  load("DDEexposure_nonsmoking.RData")
}
# time
fit_LBDDP_set2$t_mcmc
fit_LSBP_set2$t_mcmc

source("R/LBDDP_postprocess.R")
source("R/LSBP_postprocess.R")

######## 2. Posterior predictive samples ###############

######## 2.1 Posterior predictive samples of preterm birth ##############

xgrid = seq(-1.2, max(x), length = 100)
xgrid_original = xgrid*xsd + xmean
xgrid_original

X1_grid = cbind(1, xgrid)
X2_grid     <- ns(xgrid,          # Design matrix for the stick-breaking weights
                  knots=attr(Basis,"knots"),
                  Boundary.knots=attr(Basis,"Boundary.knots"), intercept = T)
X2_grid_normalized = X2_grid/sqrt(rowSums(X2_grid^2))


cdf_LSBP_set1 <- LSBP_cdf((37 - mean(y_original))/sd(y_original), X1_grid, X2_grid_normalized, fit_LSBP_set1)
cdf_LSBP_set2 <- LSBP_cdf((37 - mean(y_original))/sd(y_original), X1_grid, X2_grid_normalized, fit_LSBP_set2)
cdf_LSBP_set3 <- LSBP_cdf((37 - mean(y_original))/sd(y_original), X1_grid, X2_grid_normalized, fit_LSBP_set3)

cdf_LBDDP_set1 <- LBDDP_cdf((37 - mean(y_original))/sd(y_original), X1_grid, X2_grid_normalized, fit_LBDDP_set1)
cdf_LBDDP_set2 <- LBDDP_cdf((37 - mean(y_original))/sd(y_original), X1_grid, X2_grid_normalized, fit_LBDDP_set2)
cdf_LBDDP_set3 <- LBDDP_cdf((37 - mean(y_original))/sd(y_original), X1_grid, X2_grid_normalized, fit_LBDDP_set3)




############# 3. plot ##############################



cdf_LSBP_set1_df <- data.frame(xgrid_original,
                               y = apply(cdf_LSBP_set1, 2, mean),
                               ylo = apply(cdf_LSBP_set1, 2, function(x) quantile(x, prob  = 0.025)),
                               yup = apply(cdf_LSBP_set1, 2, function(x) quantile(x, prob  = 0.975)))
cdf_LBDDP_set1_df <- data.frame(xgrid_original,
                                y = apply(cdf_LBDDP_set1, 2, mean),
                                ylo = apply(cdf_LBDDP_set1, 2, function(x) quantile(x, prob  = 0.025)),
                                yup = apply(cdf_LBDDP_set1, 2, function(x) quantile(x, prob  = 0.975)))


cdf_LSBP_set2_df <- data.frame(xgrid_original,
                               y = apply(cdf_LSBP_set2, 2, mean),
                               ylo = apply(cdf_LSBP_set2, 2, function(x) quantile(x, prob  = 0.025)),
                               yup = apply(cdf_LSBP_set2, 2, function(x) quantile(x, prob  = 0.975)))
cdf_LBDDP_set2_df <- data.frame(xgrid_original,
                                y = apply(cdf_LBDDP_set2, 2, mean),
                                ylo = apply(cdf_LBDDP_set2, 2, function(x) quantile(x, prob  = 0.025)),
                                yup = apply(cdf_LBDDP_set2, 2, function(x) quantile(x, prob  = 0.975)))

cdf_LSBP_set3_df <- data.frame(xgrid_original,
                               y = apply(cdf_LSBP_set3, 2, mean),
                               ylo = apply(cdf_LSBP_set3, 2, function(x) quantile(x, prob  = 0.025)),
                               yup = apply(cdf_LSBP_set3, 2, function(x) quantile(x, prob  = 0.975)))

cdf_LBDDP_set3_df <- data.frame(xgrid_original,
                                y = apply(cdf_LBDDP_set3, 2, mean),
                                ylo = apply(cdf_LBDDP_set3, 2, function(x) quantile(x, prob  = 0.025)),
                                yup = apply(cdf_LBDDP_set3, 2, function(x) quantile(x, prob  = 0.975)))




cdf_LBDDP_set1_df$setting = "Logistic-beta, setting 1"
cdf_LBDDP_set2_df$setting = "Logistic-beta, setting 2"
cdf_LBDDP_set3_df$setting = "Logistic-beta, setting 3"
cdf_LBDDP_set1_df$model = "LBDDP"
cdf_LBDDP_set2_df$model = "LBDDP"
cdf_LBDDP_set3_df$model = "LBDDP"


cdf_LSBP_set1_df$setting = "Logit stick-breaking, setting 1"
cdf_LSBP_set2_df$setting = "Logit stick-breaking, setting 2"
cdf_LSBP_set3_df$setting = "Logit stick-breaking, setting 3"
cdf_LSBP_set1_df$model = "LSBP"
cdf_LSBP_set2_df$model = "LSBP"
cdf_LSBP_set3_df$model = "LSBP"

cdf_LBDDP_df = rbind(cdf_LBDDP_set1_df, cdf_LBDDP_set2_df, cdf_LBDDP_set3_df)

cdf_LSBP_df = rbind(cdf_LSBP_set1_df, cdf_LSBP_set2_df, cdf_LSBP_set3_df)
cdf_df = rbind(cdf_LBDDP_df, cdf_LSBP_df)


cdf_plot_1=
  ggplot(cdf_LBDDP_df, aes(x = xgrid_original, y = y)) +
  geom_line(color = "#00BFC4") +
  geom_point(aes(x, y), shape = 108, alpha = 1, size = 1, data = data.frame(x = x_original, y= rep(0, length(x_original))) ) +
  geom_ribbon(aes(ymin = ylo, ymax = yup), alpha = 0.2, fill = "#00BFC4") + ylim(0,1) +
  xlab("DDE (µg/L)") + ylab("Probability of preterm birth") + theme_light() +
  facet_wrap(~setting)
cdf_plot_2 =
  ggplot(cdf_LSBP_df, aes(x = xgrid_original, y = y)) +
  geom_line(color = "#F8766D") +
  geom_point(aes(x, y), shape = 108, alpha = 1, size = 1, data = data.frame(x = x_original, y= rep(0, length(x_original))) ) +
  geom_ribbon(aes(ymin = ylo, ymax = yup), alpha = 0.2, fill = "#F8766D") + ylim(0,1) +
  xlab("DDE (µg/L)") + ylab("Probability of preterm birth") + theme_light() +
  facet_wrap(~setting)

cdf_plot_12 = cdf_plot_1 / cdf_plot_2
cdf_plot_12
#ggsave("pretermbirth_nonsmoking.pdf", cdf_plot_12, width = 9, height = 5)















