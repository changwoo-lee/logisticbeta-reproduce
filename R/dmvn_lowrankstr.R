#' calculate N_d(y | mu, XDX^t+ E), E = diag(evec)
#'
#' utilizing Woodbury matrix identity (XDX^t + E)^(-1) = E^(-1) - E^(-1)X(D^(-1)+ X^t E^(-1) X)^(-1) X^t E^(-1)
#' and matrix determinant lemma det(XDX^t + E) = det(D^(-1)+ X^t E^(-1) X) det(D) det(E)
#'
#' @param y n by d
#' @param mu length d
#' @param X d by k, k<<d
#' @param Dinv inverse of D
#' @param Dlogdet log determinant of D, optional
#' @param Ediag length d, diagonal of E
#' @param log logical
#'
#' @return
#' @export
#'
#' @examples
dmvn_lowrankstr = function(y, mu, X, Dinv, Dlogdet = NULL, Ediag, log = T){
  dimen = length(Ediag)
  if(is.vector(y)) y = matrix(y, ncol = dimen)
  if(ncol(y)!=dimen) stop("wrong dimension of y, ncol(y) should be dimension")
  resid = t(y)-mu
  Xt_Einv_y = crossprod(X/Ediag, resid)

  U = chol( Dinv+ crossprod(X/sqrt(Ediag)) )
  if(is.null(Dlogdet)) Dlogdet = -determinant(Dinv, log = T)$modulus
  logdetSig = Dlogdet + sum(log(Ediag)) + 2*sum(log(diag(U)))

  quad1 = colSums((resid/sqrt(Ediag))^2) # t(resid)%*%diag(1/Ediag)%*%resid
  #quad2 = colSums(solve(t(U), Xt_Einv_y)^2)
  quad2 = colSums(backsolve(U, Xt_Einv_y, transpose = T)^2)

  quadform = quad1 - quad2
  #browser()
  loglik = -(length(Ediag)/2)*log(2*pi) - 0.5*logdetSig - quadform/2
  if(log) return(loglik) else return(exp(loglik))
}
