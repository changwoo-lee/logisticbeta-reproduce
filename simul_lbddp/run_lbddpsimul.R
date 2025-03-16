# simulation code for n = 500, datagen = A and b = 1 settings

rm(list = ls())
slurm <- TRUE
if(slurm){
  slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
} else{
  slurm_id <- 0
}

library(DDPstar)
library(splines)
source('R/polya.R')
source('R/dmvn_lowrankstr.R')
source("R/LBDDP_basis.R")
source("R/LSBP_basis.R")
source("R/linearDDP.R")

isim = slurm_id
n = 500
set.seed(isim)
x <- runif(n)
y1 <- x + rnorm(n, 0, 0.1)
y2 <- x^4 + rnorm(n, 0, 0.2)
u <- runif(n)
prob <- exp(-2*x)
y <- ifelse(u < prob, y1, y2)

df <- data.frame(x = x, y = y)


############### b = 1 ###############
# 1. linear ddp
priorLBDDP       <- list(mu_beta = rep(0,2),
                         Sigma_beta = diag(100,2),
                         a_tau = 1, b_tau = 1)
# measure time
start_time = Sys.time()
fit_linearddp = linearDDP(y,
                       cbind(1,x),
                       b = 1, H = 20,
                       prior = priorLBDDP,
                       nburn = 20000, nsave = 5000, nthin = 1)
end_time = Sys.time()
time_linearddp = end_time - start_time
# 2. 
nbasis = 6
Basis = ns(x, nbasis, intercept = T)
Basis_normalized = Basis/sqrt(rowSums(Basis^2))

start_time = Sys.time()
fit_lbddp = LBDDP_basis(y,
                        cbind(1,x),
                        Basis_normalized, b = 1, H = 20,
                        prior = priorLBDDP,
                        nburn = 20000, nsave = 5000, nthin = 1)
end_time = Sys.time()
time_lbddp = end_time - start_time
fit_lbddp$Basis = Basis
fit_lbddp$Basis_normalized = Basis_normalized
 
priorLSBP       <- list(mu_beta = rep(0,2),
                             Sigma_beta = diag(100,2),
                             a_tau = 1, b_tau = 1, sigma2_alpha = pi^2/3) # corresponding to b = 1
start_time = Sys.time()
fit_lsbp = LSBP_basis(y,
                           cbind(1, x),
                           Basis_normalized, H = 20,
                           prior = priorLSBP,
                           nburn = 20000, nsave = 5000, nthin = 1)
end_time = Sys.time()
time_lsbp = end_time - start_time
#######################################################################################

n_new = 100
x_new= seq(0,1,length = n_new) 

# True regression function
m_true_new = exp(-2*x_new)*x_new + (1-exp(-2*x_new))*x_new^4	

# True density
y_grid=seq(-1,2,length = 500)
y_delta = y_grid[2]-y_grid[1]
m2=length(y_grid)
f_true_new = matrix(0,m2,n_new)
for(i in 1:n_new){
  f_true_new[,i] = exp(-2*x_new[i])*dnorm(y_grid, mean = x_new[i], sd = 0.1) + (1-exp(-2*x_new[i]))*dnorm(y_grid, mean = x_new[i]^4, sd = 0.2)
}

str(f_true_new) # row corresponds to ygrid, col corresponds to x_new



###########################################################################
source("R/LSBP_postprocess.R")
source("R/LBDDP_postprocess.R")


newbasis           <- ns(x_new,          # Design matrix for the stick-breaking weights
                         knots=attr(Basis,"knots"),
                         Boundary.knots=attr(Basis,"Boundary.knots"), intercept = T)
newbasis_normalized = newbasis/sqrt(rowSums(newbasis^2))

fhat_linearddp = matrix(0,length(y_grid),n_new)
fhat_lbddp = matrix(0,length(y_grid),n_new)
fhat_lsbp = matrix(0,length(y_grid),n_new)

ehat_linearddp = numeric(n_new)
ehat_lbddp = numeric(n_new)
ehat_lsbp = numeric(n_new)


# take some time......
for(i_new in 1:n_new){
  temp_linearddp = linearDDP_density(y_grid, cbind(1,x_new[i_new]), fit_linearddp)
  fhat_linearddp[,i_new] = colMeans(temp_linearddp$pdf_mcmc)
  ehat_linearddp[i_new] = mean(temp_linearddp$condmean_mcmc)
 
  ###################################################################################################
  temp_lbddp = LBDDP_density(y_grid, cbind(1,x_new[i_new]), newbasis_normalized[i_new,], fit_lbddp)
  fhat_lbddp[,i_new] = colMeans(temp_lbddp$pdf_mcmc)
  ehat_lbddp[i_new] = mean(temp_lbddp$condmean_mcmc)
  
  ###################################################################################################
  temp_lsbp = LSBP_density(y_grid, cbind(1,x_new[i_new]), newbasis_normalized[i_new,], fit_lsbp)
  fhat_lsbp[,i_new] = colMeans(temp_lsbp$pdf_mcmc)
  ehat_lsbp[i_new] = mean(temp_lsbp$condmean_mcmc)
}


regerr_linearddp = sqrt(mean((ehat_linearddp - m_true_new)^2))
regerr_lbddp = sqrt(mean((ehat_lbddp - m_true_new)^2))
regerr_lsbp = sqrt(mean((ehat_lsbp - m_true_new)^2))

denserr_linearddp = sum(abs(fhat_linearddp - f_true_new)*y_delta)/n_new
denserr_lbddp = sum(abs(fhat_lbddp - f_true_new)*y_delta)/n_new
denserr_lsbp = sum(abs(fhat_lsbp - f_true_new)*y_delta)/n_new











