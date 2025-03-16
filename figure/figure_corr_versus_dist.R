### correlation as a function of distance
### reproducing figure S.1.2


rm(list=ls())

source("R/polya.R")
source("R/lb.R")

# calculate correlation of logisitc-transformed LB (thus having beta marginals) using monte carlo
# rhovec: numeric vector of correlation parameter between two points (R(x_1, x_2))
corr_logiLB<- function(rhogrid, a, b, nsim = 5000, nsize = 1000){
  # check all rhogrid is within -1 and 1
  if(any(abs(rhogrid) > 1)){
    stop("all rhogrid must be within -1 and 1")
  }
  ngrid = length(rhogrid)
  corrsave = array(0, dim = c(nsim, ngrid))
  for(i in 1:ngrid){
    rho = rhogrid[i]
    if(rho>-1 & rho < 1){
      R = matrix(c(1,rho, rho, 1), 2, 2)
      for(isim in 1:nsim){
        Zsamples = rmlb(nsize, a, b, R)
        Vsamples = 1/(1+exp(-Zsamples))
        corrsave[isim,i] = cor(Vsamples[,1],Vsamples[,2])
      }
    }else if(rho==-1){
      for(isim in 1:nsim){
        lambda = rpolya(nsize, a, b)
        if(any(lambda < 0)) browser()
        Zsamples1 = rnorm(nsize, sd = sqrt(lambda))
        Zsamples2 = -Zsamples1
        Zsamples = cbind(Zsamples1 + 0.5*lambda*(a-b), Zsamples2+ 0.5*lambda*(a-b))
        Vsamples = 1/(1+exp(-Zsamples))
        corrsave[isim,i] = cor(Vsamples[,1],Vsamples[,2])
      }
    }else{
      for(isim in 1:nsim){
        lambda = rpolya(nsize, a, b)
        Zsamples1 = rnorm(nsize, sd = sqrt(lambda))
        Zsamples2 = Zsamples1
        Zsamples = cbind(Zsamples1 + 0.5*lambda*(a-b), Zsamples2+ 0.5*lambda*(a-b))
        Vsamples = 1/(1+exp(-Zsamples))
        corrsave[isim,i] = cor(Vsamples[,1],Vsamples[,2])
      }
    }
  }
  out = colMeans(corrsave)
  # monte carlo standard error
  attr(out, "MonteCarloSE") = apply(corrsave, 2, sd)/sqrt(nsim)
  return(out)
}



dgrid = seq(0.0001, 1.5, length = 200)

# 1. Exponential kernel, range 1/3,  effective range = 1
fields::Matern(1, range = 1/3, smoothness = 0.5) # approx. 0.05
# Matern 1.5 kernel, range 1/4.75 = 4/19, effective range 1
fields::Matern(1, range = 1/4.75, smoothness = 1.5) # approx. 0.05
# Matern 2.5 kernel, range 1/5.9 = 10/59, effective range 1
fields::Matern(1, range = 1/5.9, smoothness = 2.5) # approx. 0.05


fields::Matern(dgrid, range = 1/3, smoothness = 0.5)

plot(dgrid,fields::Matern(dgrid, range = 1/3, smoothness = 0.5), type="l")
lines(dgrid,fields::Matern(dgrid, range = 1/4.75, smoothness = 1.5), type= "l")
lines(dgrid,fields::Matern(dgrid, range = 1/5.9, smoothness = 2.5), type="l")


rhogrid1 = fields::Matern(dgrid, range = 1/3, smoothness = 0.5)
rhogrid2 = fields::Matern(dgrid, range = 1/4.75, smoothness = 1.5)
rhogrid3 = fields::Matern(dgrid, range = 1/5.9, smoothness = 2.5)

corrV_102 = list()
corrV_102[[1]] = corr_logiLB(rhogrid1, a=1, b=0.2, nsim = 1000)
corrV_102[[2]] = corr_logiLB(rhogrid2, a=1, b=0.2, nsim = 1000)
corrV_102[[3]] = corr_logiLB(rhogrid3, a=1, b=0.2, nsim = 1000)

corrV_11 = list()
corrV_11[[1]] = corr_logiLB(rhogrid1, a=1, b=1, nsim = 1000)
corrV_11[[2]] = corr_logiLB(rhogrid2, a=1, b=1, nsim = 1000)
corrV_11[[3]] = corr_logiLB(rhogrid3, a=1, b=1, nsim = 1000)

corrV_12 = list()
corrV_12[[1]] = corr_logiLB(rhogrid1, a=1, b=2, nsim = 1000)
corrV_12[[2]] = corr_logiLB(rhogrid2, a=1, b=2, nsim = 1000)
corrV_12[[3]] = corr_logiLB(rhogrid3, a=1, b=2, nsim = 1000)

corrV_110 = list()
corrV_110[[1]] = corr_logiLB(rhogrid1, a=1, b=10, nsim = 1000)
corrV_110[[2]] = corr_logiLB(rhogrid2, a=1, b=10, nsim = 1000)
corrV_110[[3]] = corr_logiLB(rhogrid3, a=1, b=10, nsim = 1000)


mlb_cov <- function(a,b, r12){
  if(a==b){
    return(r12*(trigamma(a) + trigamma(b)))
  }else{
    return(trigamma(a) + trigamma(b) + (r12-1)*2*(digamma(a)-digamma(b))/(a-b))
  }
}

correta_102 = list()
correta_11 = list()
correta_12 = list()
correta_110 = list()

correta_102[[1]] = numeric(length(rhogrid1))
correta_11[[1]] = numeric(length(rhogrid1))
correta_12[[1]] = numeric(length(rhogrid1))
correta_110[[1]] = numeric(length(rhogrid1))
for(i in 1:length(rhogrid1)){
  correta_102[[1]][i] = mlb_cov(a=1,b=0.2, r12 = rhogrid1[i])/mlb_cov(a=1,b=0.2, r12 = 1)
  correta_11[[1]][i] = mlb_cov(a=1,b=1, r12 = rhogrid1[i])/mlb_cov(a=1,b=1, r12 = 1)
  correta_12[[1]][i] = mlb_cov(a=1,b=2, r12 = rhogrid1[i])/mlb_cov(a=1,b=2, r12 = 1)
  correta_110[[1]][i] = mlb_cov(a=1,b=10, r12 = rhogrid1[i])/mlb_cov(a=1,b=10, r12 = 1)
}
correta_102[[2]] = numeric(length(rhogrid2))
correta_11[[2]] = numeric(length(rhogrid2))
correta_12[[2]] = numeric(length(rhogrid2))
correta_110[[2]] = numeric(length(rhogrid2))
for(i in 1:length(rhogrid2)){
  correta_102[[2]][i] = mlb_cov(a=1,b=0.2, r12 = rhogrid2[i])/mlb_cov(a=1,b=0.2, r12 = 1)
  correta_11[[2]][i] = mlb_cov(a=1,b=1, r12 = rhogrid2[i])/mlb_cov(a=1,b=1, r12 = 1)
  correta_12[[2]][i] = mlb_cov(a=1,b=2, r12 = rhogrid2[i])/mlb_cov(a=1,b=2, r12 = 1)
  correta_110[[2]][i] = mlb_cov(a=1,b=10, r12 = rhogrid2[i])/mlb_cov(a=1,b=10, r12 = 1)
}
correta_102[[3]] = numeric(length(rhogrid3))
correta_11[[3]] = numeric(length(rhogrid3))
correta_12[[3]] = numeric(length(rhogrid3))
correta_110[[3]] = numeric(length(rhogrid3))
for(i in 1:length(rhogrid3)){
  correta_102[[3]][i] = mlb_cov(a=1,b=0.2, r12 = rhogrid3[i])/mlb_cov(a=1,b=0.2, r12 = 1)
  correta_11[[3]][i] = mlb_cov(a=1,b=1, r12 = rhogrid3[i])/mlb_cov(a=1,b=1, r12 = 1)
  correta_12[[3]][i] = mlb_cov(a=1,b=2, r12 = rhogrid3[i])/mlb_cov(a=1,b=2, r12 = 1)
  correta_110[[3]][i] = mlb_cov(a=1,b=10, r12 = rhogrid3[i])/mlb_cov(a=1,b=10, r12 = 1)
}

corrV_to_corrG <- function(mycorr, b){
  beta_mean = 1/(1+b)
  beta_var = b/((1+b)^2*(b+2))
  
  E12 = mycorr*beta_var + beta_mean^2
  return((1+b)^2/(2/E12 - 1-b ))
}

corrG_102 = lapply(corrV_102, function(x) corrV_to_corrG(x, 0.2))
corrG_11 = lapply(corrV_11, function(x) corrV_to_corrG(x, 1))
corrG_12 = lapply(corrV_12, function(x) corrV_to_corrG(x, 2))
corrG_110 = lapply(corrV_110, function(x) corrV_to_corrG(x, 10))


corrGlowerbound_102 = corrV_to_corrG(corr_logiLB(-1, a=1, b=0.2, nsim = 1000), 0.2)
corrGlowerbound_11 = corrV_to_corrG(corr_logiLB(-1, a=1, b=1, nsim = 1000), 1)
corrGlowerbound_12 = corrV_to_corrG(corr_logiLB(-1, a=1, b=2, nsim = 1000), 2)
corrGlowerbound_110 = corrV_to_corrG(corr_logiLB(-1, a=1, b=10, nsim = 1000), 10)

par(mfrow = c(4,3))
plot(dgrid, correta_102[[1]], ylim = c(0,1), type="l")
lines(dgrid, correta_102[[2]], col = 2, lty = 2)
lines(dgrid, correta_102[[3]], col = 3, lty = 3)

plot(dgrid, corrV_102[[1]], ylim = c(0,1), type="l")
lines(dgrid, corrV_102[[2]], col = 2, lty = 2)
lines(dgrid, corrV_102[[3]], col = 3, lty = 3)

plot(dgrid, corrG_102[[1]], ylim = c(0,1), type="l")
lines(dgrid, corrG_102[[2]], col = 2, lty = 2)
lines(dgrid, corrG_102[[3]], col = 3, lty = 3)
abline(h = corrGlowerbound_102, col = 4)

######################################################

plot(dgrid, correta_11[[1]], ylim = c(0,1), type="l")
lines(dgrid, correta_11[[2]], col = 2, lty = 2)
lines(dgrid, correta_11[[3]], col = 3, lty = 3)

plot(dgrid, corrV_11[[1]], ylim = c(0,1), type="l")
lines(dgrid, corrV_11[[2]], col = 2, lty = 2)
lines(dgrid, corrV_11[[3]], col = 3, lty = 3)

plot(dgrid, corrG_11[[1]], ylim = c(0,1), type="l")
lines(dgrid, corrG_11[[2]], col = 2, lty = 2)
lines(dgrid, corrG_11[[3]], col = 3, lty = 3)
abline(h = corrGlowerbound_11, col = 4)

######################################################
plot(dgrid, correta_12[[1]], ylim = c(0,1), type="l")
lines(dgrid, correta_12[[2]], col = 2, lty = 2)
lines(dgrid, correta_12[[3]], col = 3, lty = 3)


plot(dgrid, corrV_12[[1]], ylim = c(0,1), type="l")
lines(dgrid, corrV_12[[2]], col = 2, lty = 2)
lines(dgrid, corrV_12[[3]], col = 3, lty = 3)

plot(dgrid, corrG_12[[1]], ylim = c(0,1), type="l")
lines(dgrid, corrG_12[[2]], col = 2, lty = 2)
lines(dgrid, corrG_12[[3]], col = 3, lty = 3)
abline(h = corrGlowerbound_12, col = 4)

##################################################

plot(dgrid, correta_110[[1]], ylim = c(0,1), type="l")
lines(dgrid, correta_110[[2]], col = 2, lty = 2)
lines(dgrid, correta_110[[3]], col = 3, lty = 3)

plot(dgrid, corrV_110[[1]], ylim = c(0,1), type="l")
lines(dgrid, corrV_110[[2]], col = 2, lty = 2)
lines(dgrid, corrV_110[[3]], col = 3, lty = 3)

plot(dgrid, corrG_110[[1]], ylim = c(0,1), type="l")
lines(dgrid, corrG_110[[2]], col = 2, lty = 2)
lines(dgrid, corrG_110[[3]], col = 3, lty = 3)
abline(h = corrGlowerbound_110, col = 4)


# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

list_to_df <- function(curve_list, dgrid, type, group) {
  # curve_list: a list (e.g. correta_102) with 3 vectors (for the three curves)
  # dgrid: the common x-axis values
  # type: a string indicating the curve type ("correta", "corrV", or "corrG")
  # group: a string indicating the group ("102", "11", "12", "110")
  map_dfr(seq_along(curve_list), function(i) {
    tibble(
      dgrid  = dgrid,
      value  = curve_list[[i]],
      line   = factor(i),   # factor indicating the curve index (1, 2, or 3)
      type   = type,
      group  = group
    )
  })
}


df_all <- bind_rows(
  list_to_df(correta_102, dgrid, " correta", "a=1,  b=0.2"),
  list_to_df(corrV_102,    dgrid, " corrV",    "a=1,  b=0.2"),
  list_to_df(corrG_102,    dgrid, "corrG",    "a=1,  b=0.2"),
  
  list_to_df(correta_11,  dgrid, " correta", "a=1,  b=1"),
  list_to_df(corrV_11,     dgrid, " corrV",    "a=1,  b=1"),
  list_to_df(corrG_11,     dgrid, "corrG",    "a=1,  b=1"),
  
  list_to_df(correta_12,  dgrid, " correta", "a=1,  b=2"),
  list_to_df(corrV_12,     dgrid, " corrV",    "a=1,  b=2"),
  list_to_df(corrG_12,     dgrid, "corrG",    "a=1,  b=2"),
  
  list_to_df(correta_110, dgrid, " correta", "a=1, b=10"),
  list_to_df(corrV_110,    dgrid, " corrV",    "a=1, b=10"),
  list_to_df(corrG_110,    dgrid, "corrG",    "a=1, b=10")
)

df_hline <- tibble(
  group       = factor(c("a=1,  b=0.2", "a=1,  b=1", "a=1,  b=2", "a=1, b=10"), 
                       levels = c("a=1,  b=0.2", "a=1,  b=1", "a=1,  b=2", "a=1, b=10")),
  type        = factor("corrG", levels = c(" correta", " corrV", "corrG")),
  lower_bound = c(corrGlowerbound_102, corrGlowerbound_11, corrGlowerbound_12, corrGlowerbound_110)
)

annotation_df <- tibble(
  group = factor("a=1,  b=1", levels = c("a=1,  b=0.2", "a=1,  b=1", "a=1,  b=2", "a=1, b=10")),
  type  = factor(" correta", levels = c(" correta", " corrV", "corrG")),
  x     = 1.15,
  y     = 0.2,
  xend  = 1.05,
  yend  = 0.1
)

g = ggplot(df_all, aes(x = dgrid, y = value, group = line, color = line, linetype = line)) +
  geom_line() +
  geom_hline(data = df_hline, aes(yintercept = lower_bound), color = "darkgrey", linetype = 4) +
  facet_grid(group ~ type, 
             labeller = labeller(
               type = as_labeller(
                 c(" correta" = "corr(eta(x),eta(x*\"'\"))",
                   " corrV"   = "corr(V(x),V(x*\"'\"))",
                   "corrG"   = "corr(G[x](B),G[x*\"'\"](B))"),
                 label_parsed
               )
             )) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("black", "red", "blue")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  labs(x = "dist(x,x')", y = "Correlation") +
  scale_color_manual(
    name   = "Correlation as a function of dist(x,x')",
    values = c("black", "red", "blue"),
    labels = c("Matern_0.5", "Matern_1.5", "Matern_2.5")
  ) +
  scale_linetype_manual(
    name   = "Correlation as a function of dist(x,x')",
    values = c("solid", "dashed", "dotted"),
    labels = c("Matern_0.5", "Matern_1.5", "Matern_2.5")
  ) +
  # Add the annotation only in the desired panel
  # Add the arrow annotation in the specified panel
  # geom_segment(data = annotation_df,
  #              aes(x = x, y = y, xend = xend, yend = yend),
  #              arrow = arrow(length = unit(0.3, "cm")),
  #              inherit.aes = FALSE,
  #              color = "black", size = 1) +
  theme_bw() + theme(legend.position = "top")
g

#ggsave("corr_distplot_withoutarrow.pdf", g, width = 6.5, height = 7)
