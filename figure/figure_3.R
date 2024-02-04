# reproduce figure 3
rm(list = ls())
source("R/polya.R")
source("R/lb.R")
library(ggplot2) # plotting

set.seed(1)
x = c(seq(0,0.99,length = 295), seq(1,2,length = 10), seq(2.01,3,length = 295))
n = length(x)

ngrid= 200 # for prediction
xgrid = seq(0.01, 2.99, length = ngrid) # for prediction

p_true = 1/(1+exp(-cos(x*pi)))
z = rbinom(n, 1, prob = p_true)
pgrid_true = 1/(1+exp(-cos(xgrid*pi)))



set.seed(1234)
dqrng::dqset.seed(1234) # rpolya uses dqrng package
R22 = fields::Matern(fields::rdist(xgrid), range = 0.3, smoothness = 1.5)
eta = rmlb(8, a= 2, b = 4, R22) # eight LBP samples

dfprior = data.frame(xgrid = xgrid, pi_draws = t(1/(1+exp(-eta))))
g1 = ggplot(dfprior) +
  geom_line(aes(x=xgrid, y = pi_draws.1), col = 1) +
  geom_line(aes(x=xgrid, y = pi_draws.2), col = 2) +
  geom_line(aes(x=xgrid, y = pi_draws.3), col = 3) +
  geom_line(aes(x=xgrid, y = pi_draws.4), col = 4) +
  geom_line(aes(x=xgrid, y = pi_draws.5), col = 5) +
  geom_line(aes(x=xgrid, y = pi_draws.6), col = 6) +
  geom_line(aes(x=xgrid, y = pi_draws.7), col = 7) +
  geom_line(aes(x=xgrid, y = pi_draws.8), col = 8) +
  ylim(c(0,1)) + ylab("") + xlab("") + theme_bw()
g1


# take some time, without low-rank structure
source("R/LBBer_Matern.R")
source("R/LBBer_Matern_predict.R")
set.seed(1)
fit = LBBer_Matern(z, x, a = 2, b = 4, 
                   range_grid = 0.3, # discrete uniform prior on range_grid. Fixed if only one value provided 
                   smoothness = 1.5, 
                   nburn = 1000, nsave = 1000, nthin = 1)
# prediction of success probability
fit_pred = LBBer_Matern_predict(fit, xgrid)



df = data.frame(xgrid = xgrid,
                pgrid_true = pgrid_true,
                pi_est = apply(fit_pred, 2, mean),
                pi_lb = apply(fit_pred, 2, function(x) quantile(x, 0.025)),
                pi_ub = apply(fit_pred, 2, function(x) quantile(x, 0.975)),
                x = x,
                z = z)

g2 = ggplot(df) +
  geom_point(aes(x = x, y = z), alpha = 1, shape= 108) +
  geom_ribbon(aes(x = xgrid, ymin = pi_lb, ymax = pi_ub), fill = "grey85") +
  geom_line(aes(x=xgrid, y = pi_est)) +
  geom_line(aes(x=xgrid, y = pgrid_true), color = "blue", linetype = "dashed")+
  ylab("") + xlab("") +
  theme_bw()

g2
library(patchwork)
g12 = g1 + g2
#ggsave("lbp_1d_matern_range0.3_nu1.5.pdf", g12, height = 3, width = 10)
