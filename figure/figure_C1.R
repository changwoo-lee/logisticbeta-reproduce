
rm(list = ls())

# reproducing figure C.1

# misc functions
#################################################
#' Transform probability vector to stick weights
#'
#' @param p probability vector, or matrix with each row being probability vector
#' @param drop_lastcol default T, drop last column of 1's?
#'
#' @return matrix of sticks weights (with one less column by default), each row corresponds to stick weights
#' @export
#' @importFrom matrixStats rowCumsums
#'
#' @examples
#'
#' probs_to_sticks(c(0.5, 0.3, 0.2))
#' probmat = matrix(c(.5, .3, .2, .6, .2, .2), ncol = 3, byrow = T)
#' probs_to_sticks(probmat)
#'
probs_to_sticks <- function(p, drop_lastcol = T){
  if(is.vector(p)) p = matrix(p, nrow = 1, ncol = length(p))
  sticks = p
  sticks[,-1] = p[,-1]/(1-matrixStats::rowCumsums(p[,-ncol(p),drop = F]))
  sticks[which(is.nan(sticks))] <- 0 # when 0/0
  if(drop_lastcol){
    return(sticks[,-ncol(sticks), drop = F])
  }else{
    return(sticks)
  }
}

#' Transform stick weights to probability vector
#'
#' Assuming last element (of the column) does not have 1, which is dummy variable for stick weights
#'
#' @param sticks stick vector, or matrix with each row being stick weights.
#'
#' @return matrix of probabilities, each row corresponds to probabilities.
#' @importFrom matrixStats rowCumprods
#' @export
#'
#' @examples
#'
#' sticks_to_probs(c(0.1, 0.2, 0.3)) # output: length 4 probability vector
#' stickmat = matrix(c(.1, .2, .3, .4, .5, .6), ncol = 3, byrow = T)
#' sticks_to_probs(stickmat)
#'
sticks_to_probs <- function(sticks){
  if(is.vector(sticks)) sticks = matrix(sticks, nrow = 1, ncol = length(sticks))
  #if(any(abs(sticks[,ncol(sticks)] - 1) < 1e-12)) stop("last column has dummy 1, please remove")
  p = sticks
  p[,-1] = sticks[,-1]*matrixStats::rowCumprods(1-sticks[,-ncol(p), drop = F])
  p = cbind(p, 1-rowSums(p))
  return(p)
}


#' Vectorized sampling (by row by row) with probability stored in rows of matrix
#' source: https://stackoverflow.com/questions/20508658/sampling-repeatedly-with-different-probability
#'
#' @param Weight n by K matrix, each row corresponds to unnormalized probabilities of sampling weights
#'
#' @return integers with length n, each element in {1,...,K}
#' @export
#'
#' @examples
sample.rowwise <- function(Weight) {
  x <- runif(nrow(Weight))
  cumul.w <- Weight %*% upper.tri(diag(ncol(Weight)), diag = TRUE)/rowSums(Weight)
  i <- rowSums(x > cumul.w) + 1L
  i
}
# #
# temp1 = dpois(1:100, lambda = 30)*10
# temp2 = dpois(1:100, lambda = 40)
# temp12 = rbind(temp1,temp2)
#
# nsave = 1000
# tempsave = matrix(0,nrow = nsave,ncol = 2)
# for(isave in 1:nsave) tempsave[isave,] =  sample.rowwise(temp12)
#
# hist(tempsave[,1]) # mean 30
# hist(tempsave[,2]) # mean 40
#
# var(tempsave[,1])
# var(tempsave[,2])
# #
# #

##########################################################





library(BNPmix)
data("CPP")
y_original = CPP$gest
x_original = CPP$dde
# smoking group
y_original = y_original[CPP$smoke==2]
x_original = x_original[CPP$smoke==2]

y = scale(y_original)
x = scale(x_original)

xmean = attr(x, "scaled:center")
xsd = attr(x, "scaled:scale")

n = length(y)

nbasis = 6


# unnormalized basis, based on natural cubic splines
Phi = splines::ns(x, nbasis, intercept = T)
Phi_normalized = Phi/sqrt(rowSums(Phi^2))

# 1. Prior analysis
# difference in terms of co-clustering probability

xgrid = seq(min(x), max(x), length = 100)
Phi_grid = splines::ns(xgrid,
                       knots=attr(Phi,"knots"),
                       Boundary.knots=attr(Phi,"Boundary.knots"), intercept = T)

Phi_grid_normalized = Phi_grid/sqrt(rowSums(Phi_grid^2))

sigma2_gamma = pi^2/3

sigma2_grid = diag(Phi_grid %*%diag(sigma2_gamma, nrow = nbasis, ncol = nbasis) %*% t(Phi_grid))
sigma2_grid_normalized = diag(Phi_grid_normalized %*%diag(sigma2_gamma, nrow = nbasis, ncol = nbasis) %*% t(Phi_grid_normalized))


coclust_logitnormalstick <- function(sigma){
  out = numeric(length(sigma))
  for(i in 1:length(sigma)){
    temp = logitnorm::momentsLogitnorm(0, sigma[i])
    E1 = temp[1]
    E2 = as.numeric(temp[2] + temp[1]^2)
    out[i] = E2/(2*E1 - E2)
  }
  out
}




Kmax = 50
nsim = 100000
# simulate stick weights, takes some time....
set.seed(1)
ncluster = matrix(0, nsim, length(sigma2_grid))
for(i in 1:length(sigma2_grid)){
  temp = matrix(rnorm(nsim*(Kmax-1), mean = 0, sd = sqrt(sigma2_grid[i])), nsim, Kmax-1)
  sticks = 1/(1+exp(-temp))
  Weight = sticks_to_probs(sticks)
  for(isim in 1:nsim){
    ncluster[isim,i] = length(unique(sample.int(Kmax, size = n, replace = T, prob = pmax(Weight[isim,], 0))))
  }
}

ncluster_normalized = matrix(0, nsim, length(sigma2_grid_normalized))
for(i in 1:length(sigma2_grid)){
  temp = matrix(rnorm(nsim*(Kmax-1), mean = 0, sd = sqrt(sigma2_grid_normalized[i])), nsim, Kmax-1)
  sticks = 1/(1+exp(-temp))
  Weight = sticks_to_probs(sticks)
  for(isim in 1:nsim){
    ncluster_normalized[isim,i] = length(unique(sample.int(Kmax, size = n, replace = T, prob = pmax(Weight[isim,], 0))))
  }
}


M = 1

par(mfrow = c(1,2))

ddegrid = seq(min(x), max(x), length = 100)*xsd + xmean



plot(ddegrid, coclust_logitnormalstick(sqrt(sigma2_grid)), type="l", ylim = c(0.4, 0.55))
abline(h = 1/(M+1), col = "red")

plot(ddegrid, apply(ncluster, 2, mean), type="l", ylim= c(7,9))
abline(h= M*digamma(M+n) - M*digamma(M), col = 2)


library(ggplot2)
library(reshape2)
df = data.frame(ddegrid, coclust = coclust_logitnormalstick(sqrt(sigma2_grid)), ncluster = apply(ncluster, 2, mean),
                coclust_normalized = coclust_logitnormalstick(sqrt(sigma2_grid_normalized)), ncluster_normalized = apply(ncluster_normalized, 2, mean))
colnames(df) = c("dde", "coclust", "ncluster", "coclust_normalized", "ncluster_normalized")

coclust_normalized = coclust_logitnormalstick(sqrt(sigma2_grid_normalized))
ncluster_normalized = apply(ncluster_normalized, 2, mean)

# plot ddegrid versus coclust
g1 = ggplot(df, aes(x=dde, y=coclust)) +
  geom_line(linetype="dashed") +
  #geom_line(aes(y=coclust_normalized), color = "red") +
  geom_hline(yintercept = mean(coclust_normalized)) +
#  geom_hline(yintercept = 1/(M+1)) +
  ylim(c(0.4, 0.55)) + theme_bw()
g1 = g1 + ylab(expression("Pr("~s[i]~"="~s[j]~")")) + xlab("DDE (mg/L)")

g1

# plot ddegrid versus ncluster
g2 = ggplot(df, aes(x=dde, y=ncluster)) +
  geom_line(linetype="dashed") +
#  geom_line(aes(y=ncluster_normalized), color = "red") +
  geom_hline(yintercept = mean(ncluster_normalized)) +
#  geom_hline(yintercept = M*digamma(M+n) - M*digamma(M)) +
  ylim(c(7,9.2)) + theme_bw()
g2 = g2 + ylab(expression("E["~K[1023]~"]")) + xlab("DDE (mg/L)")
g2

library(patchwork)
g12 = (g1 / g2)
g12
#ggsave("figures/lbddp_priors_new.pdf", g12, width = 3, height = 3)


##

library(ggplot2)
library(reshape2)
df = data.frame(ddegrid, Phi_grid, Phi_grid_normalized)
colnames(df) = c("dde", paste0("Phi", 1:nbasis), paste0("Phi_normalized", 1:nbasis))
df = melt(df, id.vars = "dde")

df2 = df[df$variable %in% paste0("Phi_normalized", 1:nbasis),]
# linetype dashed
g3temp = ggplot(df2, aes(x=dde, y=value, col=variable)) + geom_line() + theme_bw()
# only draw Phi
df3 = df[df$variable %in% paste0("Phi", 1:nbasis),]
# add geom_line from df3 to g2
# draw dde versus Phi_normalized,
g3 = g3temp + geom_line(data = df3, aes(x=dde, y=value, col=variable), linetype = "dashed")

# specify geom_line colors
g3 = g3 + scale_color_manual(values =
                             rep(c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"), 2),
                             labels = expression(phi[1], phi[2], phi[3], phi[4], phi[5], phi[6],
                                                 varphi[1], varphi[2], varphi[3], varphi[4], varphi[5], varphi[6]))
# legend position at top
g3 = g3 + theme(legend.position = "top")
# place legend in 2 columns, not in 1 column
g3 = g3 + guides(colour = guide_legend(ncol = 6, byrow = T, override.aes = list(linetype = c(1,1,1,1,1,1,2,2,2,2,2,2))))
g3

# legend title change
g3 = g3 + labs(colour = "Feature maps \n(dashed: unnormalized)") + ylab("Feature map") + xlab("DDE (mg/L)")
g3
#ggsave("figures/lbddp_basis.pdf", g3, width = 6, height = 3)

