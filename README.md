# logisticbeta-reproduce

Code to reproduce figures, simulation results, and real data analysis results from the paper 

> Lee, C. J., Zito, A., Sang, H., & Dunson, D. B. (2025). Logistic-beta processes for dependent random probabilities with beta marginals. Bayesian Analysis (in press), https://doi.org/10.1214/25-BA1541


## Codes for reproducing figures

* Figure 1 (density of univariate logistic-beta and Polya):
  - `figure/figure_1.R`
* Figure 2 (density of bivariate logistic-beta):
  - `figure/figure_2.R`
* Figure 3 (LBP prior example and latent LBP):
  - `figure/figure_3.R`
* Figure 4 (flexibility ranges):
  - `figure/figure_4_and_scatter.R`
* Figure 5 (preterm birth probability plot for smoking group): 
  - `real_analysis/DDEexposure_smoking.R`
* Figure S.1.1 (correlation lower bounds and bivariate
  scatterplots): 
  - `figure/figure_4_and_scatter.R`
* Figure S.3.2 (preterm birth probability plot for nonsmoking group):
  - `real_analysis/DDEexposure_nonsmoking.R`

## Codes for reproducing simulation results

### Section 5.1: nonparametric binary regression

Scenario 1, Scenario 2, Scenario 3 correspond to rho = 0.1, 0.2, 0.4, respectively, for data generation settings with different Matern range parameters.

- See `simul/execute_simul_copula*.R` for running all simulations with data generation based on Gaussian copula.
- See `simul/execute_simul_lbp*.R` for running all simulations data generation based on latent LBP.
- `simul/simul_data_copula.R` and `simul/simul_data_lbp.R` are scripts to generate data under Gaussian copula and latent LBP models, respectively.
- `run_copula.R` and `run_lbp.R` are scripts to run a single simulation for Gaussian copula and latent LBP models, respectively.
- Stan code and associated prediction code for the Gaussian copula model can be found in the folder `simul/copula_codes.`

### Section 5.2: Bayesian density regression

See `simul_lbddp/run_lbddpsimul.R`.

## Codes for reproducing real data analysis results

- Smoking group: `real_analysis/DDEexposure_smoking.R`
- Nonsmoking group: `real_analysis/DDEexposure_nonsmoking.R`

## R codes

1.  `R/lb.R`: basic functions for logistic-beta

- Random variate generation: univariate `rulb()`, multivariate `rmlb()`.
- Density evaluation: univariate `dulb()`, multivariate `dmlb()` (warning: unstable!).
- Others: `mlb_cov()`, `pulb()`

2.  `R/polya.R`: basic functions for Polya distribution

- Random variate generation: `rpolya()` (requires `dqrng` package for fast random number generation)
- Density evaluation: `dpolya()` (warning: unstable!)
- Adaptive Polya proposal: `find_polya_proposal()` to find a’, b’ with moment matching criterion
- Others: mean `epolya()`, variance `vpolya()`.

3. Functions for latent LBP model for binary data

- `R/LBBer_Matern.R` implements Algorithm 1 (blocked Gibbs sampler) with Matern correlation kernel.
   * Can handle both 1 and 2-dimensional covariates
   * Fixed smoothness parameter. Discrete uniform prior is assumed on the Matern range parameter.
   * If `adapt = TRUE`, it utilizes the adaptive Polya proposal (see `find_polya_proposal()` in `R/polya.R`)
   * If `nparticle = 0`, it runs independent M-H for sampling the Polya mixing parameter. Otherwise, it runs particle Gibbs. 
   * For the unblocked version, see `R/LBBer_Matern_unblocked.R`
- `R/LBBer_Matern_predict_LB.R` is a function to get posterior predictive samples of success probabilities.

4. `R/LBDDP_basis.R` and `R/LBDDP_postprocess.R`: functions for LB-DDP mixture model

- `R/LBDDP_basis.R` implements Algorithm 2 (blocked Gibbs sampler) with normalized basis.
- `R/LBDDP_postprocess.R` is a function to get posterior predictive samples of conditional density and cumulative probability.

5.  `R/LSBP_basis.R` and `R/LSBP_postprocess.R`: functions for logit stick-breaking mixture model, implemented in R.

6. `R/linearDDP.R`: functions for linear DDP, implemented in R.
    
7. `R/dmvn_lowrankstr.R`: functions for MVN density evaluation with low-rank structure

## Example 

```R
# reproduce figure 3, latent LBP model with binary data
rm(list = ls())
source("R/polya.R")
source("R/lb.R")
library(ggplot2) # plotting
library(dqrng)
library(fields)
library(mvnfast)
library(spam)

set.seed(1)
x = c(seq(0,0.99,length = 295), seq(1,2,length = 10), seq(2.01,3,length = 295))
n = length(x)

ngrid= 200 # for prediction
xgrid = seq(0.01, 2.99, length = ngrid) # for prediction

p_true = 1/(1+exp(-cos(x*pi)))
z = rbinom(n, 1, prob = p_true) # binary data
pgrid_true = 1/(1+exp(-cos(xgrid*pi)))


# prior draw with Matern kernel
set.seed(1234)
dqrng::dqset.seed(1234) # rpolya uses dqrng package
R22 = fields::Matern(fields::rdist(xgrid), range = 0.3, smoothness = 1.5)
eta = rmlb(8, a= 2, b = 4, R22) # eight LBP samples

dfprior = data.frame(xgrid = xgrid, pi_draws = t(1/(1+exp(-eta)))) # transform to [0,1]
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

# posterior
# take some time, full-rank Matern correlation kernel
source("R/LBBer_Matern.R")
source("R/LBBer_Matern_predict.R")
set.seed(1)
fit = LBBer_Matern(z, x, a = 2, b = 4, 
                   range_grid = 0.3, # discrete uniform prior on range_grid. Fixed if only one value is provided 
                   smoothness = 1.5, 
                   nburn = 1000, nsave = 1000, nthin = 1)
# prediction of success probabilities
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
```


