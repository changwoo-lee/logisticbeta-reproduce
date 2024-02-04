#' Generate multivariate logistic-beta random variables
#'
#' @param n number of samples
#' @param a positive real, logistic beta parameter
#' @param b positive real, logistic beta parameter
#' @param R p by p matrix, correlation matrix
#'
#' @return n by p matrix, each row corresponds to random samples
#' @importFrom mvnfast rmvn
#' @export
#'
#' @examples
#'
#'
#' R = matrix(c(1,0.99,0.99,1), 2, 2)
#' draws = rmlb(1000, a = 1, b = 3, R)
#'
#' plot(draws)
#'
#' grid = seq(-10,10, length = 200)
#'
#' hist(draws[,1], freq = F)
#' lines(grid, dlb(grid, 1, 3))
#'
#' hist(draws[,2], freq = F)
#' lines(grid, dlb(grid, 1, 3))
#'
#'
rmlb <- function(n, a, b, R, kmax = 200){
  p = ncol(R)

  lambda = rpolya(n, a, b, kmax = kmax)
  theta = (a - b)/2 # scalar

  #x = mvnfast::rmvn(n, mu = rep(0,p), sigma = R)
  x = matrix(rnorm(n*p), n, p)%*%chol(R)
  #rmvn(n,lambda*theta, sqrt(lambda)*t(Sigma*sqrt(lambda)))
  return(lambda*theta + sqrt(lambda)*x)
}


#' Density of multivariate logistic beta distribution (with shared parameter)
#'
#' @param x d-dimensional vector
#' @param a
#' @param b
#' @param R
#' @param kmax
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dmlb <- function(x, a, b, R, kmax = 1000, log = F){
  x = as.numeric(x)
  kseq = 0:kmax
  d = ncol(R)

  logwk = lchoose(-(a + b), 0:kmax) + log((a + b)/2 + kseq) - lbeta(a, b)
  signs = (-1)^kseq
  Rinv = solve(R)
  Rinvx = Rinv%*%x # vector
  zetak = (kseq+a)*(kseq+b) + ((a - b)/2)^2*sum(Rinv)
  xRx = sum(x*Rinvx)
  logsummand = logwk + (d-2)/4*(log(zetak) - log(xRx)) + log(base::besselK(sqrt(xRx*zetak), nu = d/2 - 1))
  logsums_positive = matrixStats::logSumExp(logsummand[signs == 1])
  logsums_negative = matrixStats::logSumExp(logsummand[signs == -1])

  if(any(logsums_positive < logsums_negative)){
    stop("numerical error, negative density")
  }else{
    logsum = logsums_positive + log1p(-exp(logsums_negative-logsums_positive)) # log-minus-exp
  }
  logdensity = logsum + log(2) + (a - b)/2*sum(Rinvx) -(d/2)*log(2*pi) - 0.5*determinant(R, logarithm = T)$modulus
  logdensity = as.numeric(logdensity)
  if(log){
    return(logdensity)
  }else{
    return(exp(logdensity))
  }
}

mlb_cov <- function(a,b, r12){
  if(a==b){
    return(r12*(trigamma(a) + trigamma(b)))
  }else{
    return(trigamma(a) + trigamma(b) + (r12-1)*2*(digamma(a)-digamma(b))/(a-b))
  }
}



#' Univariate logistic beta distribution
#'
#' @param x real, evaluation point
#' @param n integer, number of samples for random variable generation.
#' @param a positive real, logistic beta parameter
#' @param b positive real, logistic beta parameter
#' @param log logical(default F), return log value?
#'
#' @name lb
#' @return density of logistic beta distribution evaluated at x
#' @export
#'
#' @examples
#'
#' grid = seq(-5, 5, length = 200)
#' plot(grid, dulb(grid, a = 1, beta = 2))
#'
dulb = function(x, a, b, log = F){
  logdensity = -lbeta(a, b) - b*x - (a + b)*log(1+exp(-x))
  if(log){
    return(logdensity)
  }else{
    return(exp(logdensity))
  }
}

#' @rdname lb
#' @importFrom stats rbeta
#' @return random samples from logistic beta distribution
#' @export
#' @examples
#'
#' rulb(100, a = 1, b = 2)
#'
rulb = function(n, a, b){
  p = rbeta(n, a, b)
  return(log(p/(1-p)))
}


#' @rdname lb
#' @importFrom stats rbeta
#' @return random samples from logistic beta distribution
#' @export
#' @examples
#'
#' pulb(1, a = 1, b = 2)
#'
pulb = function(x, a, b, lower.tail = T, log.p = F){
  pbeta(1/(1+exp(-x)),a, b, lower.tail = lower.tail, log.p = log.p)
}

