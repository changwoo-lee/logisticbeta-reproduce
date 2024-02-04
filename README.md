# logisticbeta-reproduce

Code to reproduce figures, simulation results, and real data analysis results from the paper 

“Logistic-beta processes for modeling dependent random probabilities
with beta marginals” 

## Codes for reproducing figures

* Figure 1 (density of univariate logistic-beta and Polya):
  - `figure/figure_1.R`
* Figure 2 (density of bivariate logistic-beta):
  - `figure/figure_2.R`
* Figure 3 (LBP prior example and latent LBP):
  - `figure/figure_3.R`
* Figure 4 (correlation lower bounds of bivariate betas and LB-DDP):
  - `figure/figure_4_and_A1.R`
* Figure 5 (preterm birth probability plot for smoking group): 
  - `real_analysis/DDEexposure_smoking.R`
* Figure 6 (preterm birth probability plot for nonsmoking group):
  - `real_analysis/DDEexposure_nonsmoking.R`
* Figure A.1 (correlation lower bounds of LB-DDP and bivariate
  scatterplots): 
  - `figure/figure_4_and_A1.R`
* Figure C.1 (effect of feature map normalization):
  - `figure/figure_C1.R`

## Codes for reproducing simulation results

Scenario 1, Scenario 2, Scenario 3 corresponds to rho = 0.1, 0.2, 0.4, respectively for data generation.

- See `simul/execute_simul_copula*.R` for running all simulations with data generation based on Gaussian copula.
- See `simul/execute_simul_lbp*.R` for running all simulations data generation based on latent LBP.
- `simul/simul_data_copula.R` and `simul/simul_data_lbp.R` are scripts to generate data under Gaussian copula and latent LBP models, respectively.
- `run_copula.R` and `run_lbp.R` are scripts to run a single simulation for Gaussian copula and latent LBP models, respectively.
- Stan code and associated prediction code for Gaussian copula model can be found at `simul/copula_codes.R`.

## Codes for reproducing real data analysis results

- Smoking group: `real_analysis/DDEexposure_smoking.R`
- Nonsmoking group: `real_analysis/DDEexposure_nonsmoking.R`

## R codes

1.  `R/lb.R`: basic functions for logistic-beta

- Random variate generation: univariate `rulb()`, multivariate `rmlb()`.
- Density evaluation: univariate `dulb()`, multivariate `dmlb()`
  (warning: unstable!).
- Others: `mlb_cov()`, `pulb()`

2.  `R/polya.R`: basic functions for Polya distribution

- Random variate generation: `rpolya()` (requires `dqrng` package for
  fast random number generation)
- Density evaluation: `dpolya()` (warning: unstable!)
- Adaptive Polya proposal: `find_polya_proposal()` to find a’, b’ with
  moment matching
- Others: mean `epolya()`, variance `vpolya()`.

3.  Functions for latent LBP model for binary data

- `R/LBBer_Matern.R` implements Algorithm 1 (blocked Gibbs sampler)
  with Matern correlation kernel. For unblocked version (only for simulation), see `R/LBBer_Matern_unblocked.R`
- `R/LBBer_Matern_predict_LB.R` is a function to get posterior predictive samples of
  success probabilities
- Can handle both 1 and 2-dimensional, based on Matern with fixed
  smoothness parameter
- Discrete uniform prior is assumed on Matern range parameter
- If `nparticle = 0`, it runs independnt M-H for sampling Polya mixing
  parameter, otherwise particle Gibbs.
- If `adapt = TRUE`, it utilizes adaptive Polya proposal (see
  `find_polya_proposal()` in `R/polya.R`)

4.  `R/LBDDP_basis.R` and `R/LBDDP_postprocess.R`: functions for LB-DDP mixture
    model

- `R/LBDDP_basis.R` implements Algorithm 2 (blocked Gibbs sampler) with
  normalized basis
- `R/LBDDP_postprocess.R` is a function to get posterior predictive
  samples of conditional density and cumulative distribution

5.  `R/LSBP_basis.R` and `R/LSBP_postprocess.R`: functions for logit
    stick-breaking mixture model, implemented in R.
    
6. `R/dmvn_lowrankstr.R`: functions for MVN density evaluation with low-rank structure


