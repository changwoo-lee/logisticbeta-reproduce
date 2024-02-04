// based on the code from https://github.com/mbjoseph/gpp-speed-test
// gaussian copula with beta(1, b) mariginals. For beta(a,1) marginal, modify below

data {
  real a; 
  real b;
  real<lower = 0> range_lb;
  real<lower = 0> range_ub;
  int<lower = 1> n; // number of data
  int<lower = 0, upper=1> z[n]; // binary data
  matrix[n, n] D;// distance matrix for sample locations
}

parameters {
//  real<lower = 0> sig_sq; // sd of latent gp : correlation!
  real<lower = 0> rangepar; // range parameter of exponential kernel
  vector[n] zzz; // stnormal
}

transformed parameters {
  vector[n] theta; // prob
  vector[n] w;
  cov_matrix[n] Sigma;

  for (i in 1:(n-1)) {
    for (j in (i + 1):n) {
      // there are multiple parametrizaton of Matern, such as scale range by a factor or not...
      Sigma[i, j] = exp(-D[i, j] / rangepar) * (1 + D[i, j] / rangepar);
      Sigma[j, i] = Sigma[i, j];
    }
  }
  for (k in 1:n) Sigma[k, k] = 1;

  w = cholesky_decompose(Sigma) * zzz;

  // # Transformaton from latent GP to Bernoulli probability theta
  // # Gaussian copula: theta = B^(-1)(Phi(w)) where
  // # Phi is cdf of standard normal and B^(-1) is inverse cdf of beta(a,b)
  // # a, b are parameters of beta distribution

  // # When a= 1, b = 1
  // theta = Phi(w);
  // # general a, b=1
  // theta = pow(Phi(w), 1/a);
  // # a=1, general b
  theta = 1 - pow(1-Phi(w), 1/b);
  // # for general a and b, need inverse cdf of beta, need stanheader version 2.30 or more....
  // # I got compatibility issues with rstan and dev version of StanHeaders package...
  // # inv_inc_beta does not support vector input...
  // for(i in 1:n){
  //   theta[i] = inv_inc_beta(a, b, Phi(w[i])); // inverse cdf of beta(a,b),
  // }
}

model {
  rangepar ~ uniform(range_lb, range_ub);
  zzz ~ normal(0, 1); // gp
  z ~ bernoulli(theta);
}
