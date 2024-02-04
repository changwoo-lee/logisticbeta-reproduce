library(dqrng) # faster RNG
library(matrixStats) # log-sum-exp

#' Generate Polya random variables with parameters a, b
#'
#' @param n integer, number of samples
#' @param a positive real, first parameter
#' @param b positive real, second parameter
#' @param kmax integer (default 1000), truncation point
#'
#' @return length n vector, random samples from Polya distribution
#' @export
#' @importFrom stats rexp dqrng dqrexp
#' @examples
#'
#' rpolya(10, a = 1, b = 2)
#'
rpolya = function(n, a, b, kmax = 200){
  k = 0:kmax # start from 0
  denom = (a+k)*(b+k)/2 # length kmax + 1
  #temp = matrix(rexp(n*(kmax + 1))/rep(denom, n), nrow = kmax + 1, ncol = n, byrow = F)
  temp = matrix(dqrng::dqrexp(n*(kmax + 1))/rep(denom, n), nrow = kmax + 1, ncol = n, byrow = F)
  colSums(temp)
}

#' Evaluate density of the Polya distribution
#'
#' @param x scalar or vector of positive real, evaluation point
#' @param a positive real, first parameter
#' @param b positive real, second parameter
#' @param log default F, return log density?
#' @param kmax default 100, number of summands
#'
#' @return (log) density evaulated at x
#' @export
#' @importFrom matrixStats colLogSumExps
#'
#' @examples
#'
#' a = 1
#' b = 2
#' xdraw = rpolya(1000, a, b)
#' hist(xdraw, breaks=  40, freq = F)
#' grid = seq(0.1, 10, length = 300)
#' lines(grid, dpolya(grid, a, b))
#'
dpolya <- function(x, a, b, log = F, kmax = 200){
  kseq = 0:kmax
  logwk = lchoose(-(a + b), 0:kmax) + log((a + b)/2 + kseq) - lbeta(a, b)
  exponent_mat = matrix(-0.5*(a + kseq)*(b + kseq)) %*% x # (kmax + 1) by n matrix
  logsummand_mat = logwk + exponent_mat # each column corresponds to x

  signs = (-1)^kseq
  logsums_positive = matrixStats::colLogSumExps(logsummand_mat[signs == 1,, drop = F])
  logsums_negative = matrixStats::colLogSumExps(logsummand_mat[signs == -1,, drop = F])

  if(any(logsums_positive < logsums_negative)){
    warning("numerical error, return 0 density value")
    logsums_positive = logsums_negative + pmax(logsums_positive-logsums_negative, 0)

    logdensity = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }else{
    logdensity = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }
  if(log){
    return(logdensity)
  }else{
    return(exp(logdensity))
  }
}

# find a',b' that satisfies a'+b' = a+b and has mean closest to lambda_curr
# WLOG, a' < b' and a' is in the interval (eps, (a+b)/2 - eps)
# the order does not matter, since Polya(a,b) is same as Polya(b,a)
find_polya_proposal <- function(a,b,lambda_curr, eps = 0.001){
  c = a + b
  aprime = optimize(function(x) (2*(digamma(x)-digamma(c-x))/(2*x-c) - lambda_curr)^2, interval = c(eps, c/2-eps))$minimum
  bprime = c - aprime
  return(c(aprime, bprime))
}


# mean
epolya <- function(a,b){
  if(a==b){
    return(2*trigamma(a))
  }else{
    return(2*(digamma(a)-digamma(b))/(a-b))
  }
}

# variance
vpolya <- function(a,b){
  if(a==b){
    return(2*psigamma(a, 3)/3)
  }else{
    return( 4/(a-b)^2*(trigamma(a)+trigamma(b)- (2*(digamma(a)-digamma(b))/(a-b))) )
  }
}
