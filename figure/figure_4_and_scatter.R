### correlation lower bound analysis
### reproducing figures 4 and A1

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


# 1. Logistic-beta

# calculate correlation of logisitc-transformed LB (thus having beta marginals) using monte carlo
# rhovec: numeric vector of correlation parameter between two points (R(x_1, x_2))
E12_logiLB<- function(rhogrid, a, b, nsim = 5000, nsize = 1000){
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
        corrsave[isim,i] = mean(Vsamples[,1]*Vsamples[,2])
      }
    }else if(rho==-1){
      for(isim in 1:nsim){
        lambda = rpolya(nsize, a, b)
        if(any(lambda < 0)) browser()
        Zsamples1 = rnorm(nsize, sd = sqrt(lambda))
        Zsamples2 = -Zsamples1
        Zsamples = cbind(Zsamples1 + 0.5*lambda*(a-b), Zsamples2+ 0.5*lambda*(a-b))
        Vsamples = 1/(1+exp(-Zsamples))
        corrsave[isim,i] = mean(Vsamples[,1]*Vsamples[,2])
      }
    }else{
      for(isim in 1:nsim){
        lambda = rpolya(nsize, a, b)
        Zsamples1 = rnorm(nsize, sd = sqrt(lambda))
        Zsamples2 = Zsamples1
        Zsamples = cbind(Zsamples1 + 0.5*lambda*(a-b), Zsamples2+ 0.5*lambda*(a-b))
        Vsamples = 1/(1+exp(-Zsamples))
        corrsave[isim,i] = mean(Vsamples[,1]*Vsamples[,2])
      }
    }
  }
  out = colMeans(corrsave)
  # monte carlo standard error
  attr(out, "MonteCarloSE") = apply(corrsave, 2, sd)/sqrt(nsim)
  return(out)
}

# 2. copula (frechet lower bound)


# calculate correlation of logisitc-transformed LB (thus having beta marginals) using monte carlo
# rhovec: numeric vector of correlation parameter between two points (R(x_1, x_2))
E12_copula<- function(rhogrid, a, b, nsim = 5000, nsize = 1000){
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
        Zsamples = mvnfast::rmvn(nsize, sigma = R)
        Vsamples = qbeta(pnorm(Zsamples), a, b)
        corrsave[isim,i] = mean(Vsamples[,1]*Vsamples[,2])
      }
    }else if(rho==-1){
      for(isim in 1:nsim){
        Zsamples1 = rnorm(nsize)
        Zsamples2 = -Zsamples1
        Zsamples = cbind(Zsamples1, Zsamples2)
        Vsamples = qbeta(pnorm(Zsamples), a, b)
        corrsave[isim,i] = mean(Vsamples[,1]*Vsamples[,2])
      }
    }else{
      for(isim in 1:nsim){
        Zsamples1 = rnorm(nsize)
        Zsamples2 = Zsamples1
        Zsamples = cbind(Zsamples1, Zsamples2)
        Vsamples = qbeta(pnorm(Zsamples), a, b)
        corrsave[isim,i] = mean(Vsamples[,1]*Vsamples[,2])
      }
    }
  }
  out = colMeans(corrsave)
  # monte carlo standard error
  attr(out, "MonteCarloSE") = apply(corrsave, 2, sd)/sqrt(nsim)
  return(out)
}

#
DKddp_mixmomentlowerbd <- function(b){
  b^(3/2)/((b+1)*sqrt(2+b)) + 2/(1+b) - 1
}


#### reproduce figure 4, lb-ddp range #####


# grid of b, same distance in log scale, length 100, must include b=1
bgrid_1 = 2^(seq(log2(1/16), log2(16), length = 101))

bgrid = bgrid_1[13:101]
bgrid # length 89
ngrid = length(bgrid)

corbeta1 = corbeta2 = corbeta3 = corbeta4 = numeric(ngrid)
hypertie1 = hypertie2 = hypertie3 = hypertie4 = numeric(ngrid)
corddp1 = corddp2 = corddp3 = corddp4 = numeric(ngrid)


# take some time...
for(i in 1:ngrid){
  temp1 = E12_logiLB(rho=-1, a=1, b=bgrid[i])
  corbeta1[i] = (temp1 - 1/(1+bgrid[i])^2) / (bgrid[i]/((1+bgrid[i])^2*(2+bgrid[i])))
  hypertie1[i] = (1+bgrid[i])/ ( 2/temp1 - 1- bgrid[i])
  corddp1[i] = (1+bgrid[i])^2/ ( 2/temp1 - 1- bgrid[i])
  
  temp2 = E12_copula(rho=-1, a=1, b=bgrid[i])
  corbeta2[i] = (temp2 - 1/(1+bgrid[i])^2) / (bgrid[i]/((1+bgrid[i])^2*(2+bgrid[i])))
  hypertie2[i] = (1+bgrid[i])/ ( 2/temp2 - 1- bgrid[i])
  corddp2[i] = (1+bgrid[i])^2/ ( 2/temp2 - 1- bgrid[i])
  
  corbeta3[i] = 0
  hypertie3[i] = 1/(1+2*bgrid[i])
  corddp3[i] = (1+bgrid[i])/(1+2*bgrid[i])
  
  temp4 = DKddp_mixmomentlowerbd(bgrid[i])
  corbeta4[i] = (temp4 - 1/(1+bgrid[i])^2) / (bgrid[i]/((1+bgrid[i])^2*(2+bgrid[i])))
  # corbeta4 = sqrt(bgrid)*(bgrid+1)*sqrt(bgrid+ 2) - bgrid*(bgrid+2)
  hypertie4[i] = (1+bgrid[i])/ ( 2/temp4 - 1- bgrid[i])
  corddp4[i] = (1+bgrid[i])^2/ ( 2/temp4 - 1- bgrid[i])
}

plot(bgrid,corbeta1, type= "l", log = "x", ylim = c(-1,1))
lines(bgrid,corbeta2, col = 2, type="b")
lines(bgrid,corbeta3, col = 3)
lines(bgrid,corbeta4, col = 4)

plot(bgrid, hypertie1, type= "l", log = "x", ylim = c(0,1))
lines(bgrid,hypertie2, col = 2, type="b")
lines(bgrid,hypertie3, col = 3)
lines(bgrid,hypertie4, col = 4)

plot(bgrid,corddp1, type= "l", log = "x", ylim = c(0,1))
lines(bgrid,corddp2, col = 2, type="b")
lines(bgrid,corddp3, col = 3)
lines(bgrid,corddp4, col = 4)

corbeta4_1 = sqrt(bgrid)*(bgrid+1)*sqrt(bgrid+ 2) - bgrid*(bgrid+2)
plot(bgrid, corbeta4_1)
lines(bgrid, corbeta4)

# make it ggplot
library(reshape)
library(ggplot2)
library(latex2exp)

df = data.frame(bgrid, temp1 = corbeta1, temp2 = corbeta4, temp3 = corbeta3, temp4 = corbeta2, temp5 = rep(2, length(corbeta2))) # change ordering between temp2 and temp4
df = melt(df, id.vars = "bgrid")

#ggplot(df, aes(x = bgrid, y = value, col = variable)) + geom_line() + scale_x_log10() + ylim(-1,1)
# log scale in x
# 4 linetypes: solid, dotted, twodash, longdash
# g = ggplot(df, aes(x = bgrid, y = value, col = variable, linetype=variable)) + geom_line() +
#   theme_bw() + theme(legend.position = "bottom") +
#   scale_color_manual(values = c("black", "red", "purple", "blue"), labels = c("(i) LB-DDP", "(ii) AR1sq-DDP","(iii) tsDDP (inf: indep. beta)","(iv) Copula-based DDP (inf: Fréchet lower bound)"), breaks = c("temp1", "temp2", "temp3", "temp4")) +
#   scale_linetype_manual(values = c(1,2,3,4), labels = c("(i) LB-DDP", "(ii) AR1sq-DDP","(iii) tsDDP (inf: indep. beta)","(iv) Copula-based DDP (inf: Fréchet lower bound)"), breaks = c("temp1", "temp2", "temp3", "temp4"))+
#   scale_x_continuous(trans = "log2", breaks =c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16), labels = c("1/8", "1/4","1/2","1","2","4","8","16")) + ylim(-1,1) + xlab("b") + ylab(TeX("inf Corr$(V_h(x_i), V_h(x_j))$"))+
#   labs(color = "", linetype="") + theme_bw() + theme(text = element_text(family = "serif"))#+ theme(legend.position = "top")
#g

mylabel = c("Infimum of M1 (LB-DDP)", "Infimum of M2", "Infimum of M3", "Infimum of M4", "Supremum (common for M1-M4)")

g = ggplot(df, aes(x = bgrid, y = value, col = variable, linetype=variable)) + geom_line() +
  theme_bw() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("black", "red", "purple", "blue", "grey60"), labels = mylabel, breaks = c("temp1", "temp2", "temp3", "temp4","temp5")) +
  scale_linetype_manual(values = c(1,2,3,4,5), labels = mylabel, breaks = c("temp1", "temp2", "temp3", "temp4","temp5"))+
  scale_x_continuous(trans = "log2", breaks =c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16), labels = c("1/8", "1/4","1/2","1","2","4","8","16")) + ylim(-1,1) + xlab("b") + ylab(TeX("Corr$(V_h(x), V_h(x'))$"))+
  labs(color = "", linetype="") + theme_bw() + theme(text = element_text(family = "serif"))#+ theme(legend.position = "top")
# add horizontal line at y=1
# add (bgrid, 1)
g = g + geom_line(aes(x = bgrid, y = 1), linetype = "aa", col = "grey60")
g

df1 = data.frame(bgrid, temp1 = hypertie1, temp2 = hypertie4, temp3 = hypertie3, temp4 = hypertie2) # change ordering between temp2 and temp4
df1 = melt(df1, id.vars = "bgrid")
ggplot(df1, aes(x = bgrid, y = value, col = variable)) + geom_line() + scale_x_log10() + ylim(0,1)
# log scale in x
# 4 linetypes: solid, dotted, twodash, longdash
g1 = ggplot(df1, aes(x = bgrid, y = value, col = variable, linetype=variable)) + geom_line() +
  theme_bw() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("black", "red", "purple", "blue"), labels = c("(i) LB-DDP", "(ii) AR1sq-DDP","(iii) tsDDP (inf: indep. beta)","(iv) Copula-based DDP (inf: Fréchet lower bound)"), breaks = c("temp1", "temp2", "temp3", "temp4")) +
  scale_linetype_manual(values = c(1,2,3,4), labels = c("(i) LB-DDP", "(ii) AR1sq-DDP","(iii) tsDDP (inf: indep. beta)","(iv) Copula-based DDP (inf: Fréchet lower bound)"), breaks = c("temp1", "temp2", "temp3", "temp4"))+
  scale_x_continuous(trans = "log2", breaks =c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16), labels = c("1/8", "1/4","1/2","1","2","4","8","16")) + ylim(0,1) + xlab("b") + ylab(TeX("tie probability"))+
  labs(color = "", linetype="") + theme_bw() + theme(text = element_text(family = "serif"))
#+ theme(legend.position = "top")
# 4 linetypes: solid, dotted, twodash, longdash
# g1 = ggplot(df1, aes(x = bgrid, y = value, col = variable, linetype=variable)) + geom_line() +
#   theme_bw() + theme(legend.position = "bottom") +
#   scale_color_manual(values = c("black", "red", "purple", "blue"), labels = mylabel, breaks = c("temp1", "temp2", "temp3", "temp4")) +
#   scale_linetype_manual(values = c(1,2,3,4), labels = mylabel, breaks = c("temp1", "temp2", "temp3", "temp4"))+
#   scale_x_continuous(trans = "log2", breaks =c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16), labels = c("1/8", "1/4","1/2","1","2","4","8","16")) + ylim(0,1) + xlab("b") + ylab("tie probability")+
#   labs(color = "", linetype="") + theme_bw() + theme(text = element_text(family = "serif"))
g1
#+ theme(legend.position = "top")
g1 = g1 + geom_line(aes(x = bgrid, y = 1/(bgrid+1)), linetype = "aa", col = "grey60")

g1

df2 = data.frame(bgrid, temp1 = corddp1, temp2 = corddp4, temp3 = corddp3, temp4 = corddp2) # change ordering between temp2 and temp4
df2 = melt(df2, id.vars = "bgrid")

ggplot(df2, aes(x = bgrid, y = value, col = variable)) + geom_line() + scale_x_log10() + ylim(-1,1)
# log scale in x
# 4 linetypes: solid, dotted, twodash, longdash
g2 = ggplot(df2, aes(x = bgrid, y = value, col = variable, linetype=variable)) + geom_line() +
  theme_bw() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("black", "red", "purple", "blue"), labels = c("(i) LB-DDP", "(ii) AR1sq-DDP","(iii) tsDDP (inf: indep. beta)","(iv) Copula-based DDP (inf: Fréchet lower bound)"), breaks = c("temp1", "temp2", "temp3", "temp4")) +
  scale_linetype_manual(values = c(1,2,3,4), labels = c("(i) LB-DDP", "(ii) AR1sq-DDP","(iii) tsDDP (inf: indep. beta)","(iv) Copula-based DDP (inf: Fréchet lower bound)"), breaks = c("temp1", "temp2", "temp3", "temp4"))+
  scale_x_continuous(trans = "log2", breaks =c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16), labels = c("1/8", "1/4","1/2","1","2","4","8","16")) + ylim(0,1) + xlab("b") + ylab(TeX("Corr$(G_{x}(B), G_{x'}(B))$"))+
  labs(color = "", linetype="") + theme_bw() + theme(text = element_text(family = "serif"))
#+ theme(legend.position = "top")
g2
# 4 linetypes: solid, dotted, twodash, longdash
# g2 = ggplot(df2, aes(x = bgrid, y = value, col = variable, linetype=variable)) + geom_line() +
#   theme_bw() + theme(legend.position = "bottom") +
#   scale_color_manual(values = c("black", "red", "purple", "blue"), labels = mylabel, breaks = c("temp1", "temp2", "temp3", "temp4")) +
#   scale_linetype_manual(values = c(1,2,3,4), labels = mylabel, breaks = c("temp1", "temp2", "temp3", "temp4"))+
#   scale_x_continuous(trans = "log2", breaks =c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16), labels = c("1/8", "1/4","1/2","1","2","4","8","16")) + ylim(0,1) + xlab("b") + ylab(TeX("inf Corr$(G_{x}(B), G_{x'}(B))$"))+
#   labs(color = "", linetype="") + theme_bw() + theme(text = element_text(family = "serif"))
#+ theme(legend.position = "top")
g2 = g2 + geom_line(aes(x = bgrid, y = 1), linetype = "aa", col = "grey60")
g2


library(ggpubr)
g012 = ggarrange(g, g1, g2, ncol = 3, common.legend = T)
g012
ggsave(g012, filename = "corr_beta_ddp_newnew_012.pdf", width = 8, height = 2.5)



#############################################

g = ggplot(df2, aes(x = bgrid, y = value, col = variable, linetype=variable)) + geom_line() +
  theme_bw() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("black", "red", "purple", "blue"), labels = c("(i) LB-DDP", "(ii) AR1sq-DDP","(iii) tsDDP","(iv) Copula-based DDP"), breaks = c("temp1", "temp2", "temp3", "temp4")) +
  scale_linetype_manual(values = c(1,2,3,4), labels = c("(i) LB-DDP", "(ii) AR1sq-DDP","(iii) tsDDP","(iv) Copula-based DDP"), breaks = c("temp1", "temp2", "temp3", "temp4"))+
  scale_x_continuous(trans = "log2", breaks =c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16), labels = c("1/8", "1/4","1/2","1","2","4","8","16")) + ylim(0,1) + xlab("b") + ylab(TeX("inf Corr$(G_{x_i}(B), G_{x_j}(B))$"))+
  labs(color = "", linetype="") + theme_bw() + theme(legend.position = "top", text = element_text(family = "serif"))


g = ggplot(df2, aes(x = bgrid, y = value, col = variable, linetype=variable)) + geom_line() +
  theme_bw() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("black", "red", "purple", "blue"), labels = mylabel, breaks = c("temp1", "temp2", "temp3", "temp4")) +
  scale_linetype_manual(values = c(1,2,3,4), labels = mylabel, breaks = c("temp1", "temp2", "temp3", "temp4"))+
  scale_x_continuous(trans = "log2", breaks =c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16), labels = c("1/8", "1/4","1/2","1","2","4","8","16")) + ylim(0,1) + xlab("b") + ylab(TeX("inf Corr$(G_{x_i}(B), G_{x_j}(B))$"))+
  labs(color = "", linetype="") + theme_bw() + theme(legend.position = "top", text = element_text(family = "serif"))


#which(bgrid==2)
g = g + geom_point(data=data.frame(bgrid = 2, variable = "temp1", value = mean(df2[c(51,52),3])) , colour="black", shape=1, size=2, stroke=1.5)
g = g + geom_point(data=data.frame(bgrid = 2, variable = "temp1", value = mean(df2[c(51,52) + 89,3])) , colour="red", shape=1, size=2, stroke=1.5)
g = g + geom_point(data=data.frame(bgrid = 2, variable = "temp1", value = mean(df2[c(51,52) + 89*2,3])) , colour="purple", shape=1, size=2, stroke=1.5)
g = g + geom_point(data=data.frame(bgrid = 2, variable = "temp1", value = mean(df2[c(51,52) + 89*3,3])) , colour="blue", shape=1, size=2, stroke=1.5)

g
ggsave(g, filename = "corr_beta_ddp_circle_newlegend.pdf", width = 6, height = 3)
#ggsave(g, filename = "corr_beta_ddp_circle_nolegend.eps", width = 6, height = 3)







b = 2
set.seed(123)
mlbdraw = 1/(1+exp(-rmlb(1000, 1, b, matrix(c(1,-0.9999,-0.9999,1),2,2))))

#
z = matrix(0, 1000, 2)
z[,1] = rnorm(1000)
z[,2] = -z[,1]
copuladraw = 1-(1-pnorm(z))^(1/b)

# ar1sq-ddp
zetasq = rnorm(1000)^2
eta1sq = rnorm(1000)^2
eta2sq = rnorm(1000)^2

ar1sqdraw = cbind(1 - exp(-(zetasq + eta1sq)/(2*b)),  1 - exp(-(zetasq + eta2sq )/(2*b)))
plot(ar1sqdraw)

# ts ddp
tsddpdraw = cbind(rbeta(1000, 1, b), rbeta(1000, 1, b))
plot(tsddpdraw)
# change above bivarate scatterplot to ggplot
# 2 plots
# draw scatterplot of mlbdraw
df11 = data.frame(x = mlbdraw[,1], y = mlbdraw[,2])
df22 = data.frame(x = ar1sqdraw[,1], y = ar1sqdraw[,2])
df33 = data.frame(x = tsddpdraw[,1], y = tsddpdraw[,2])
df44 = data.frame(x = copuladraw[,1], y = copuladraw[,2])



p1 = ggplot(df11, aes(x = x, y = y)) + geom_point(alpha = 0.3) + theme_bw() +  labs(x="",y="") #+ ggtitle("LB-DDP, bivariate beta") +  theme(plot.title = element_text(hjust = 0.5))
p2 = ggplot(df22, aes(x = x, y = y)) + geom_point( color = "red", alpha = 0.3) + theme_bw() + labs(x="",y="") #+ ggtitle("AR1sq-DDP") +  theme(plot.title = element_text(hjust = 0.5))
p3 = ggplot(df33, aes(x = x, y = y)) + geom_point( color = "purple", alpha = 0.3) + theme_bw() + labs(x="",y="") #+ ggtitle("tsDDP") +  theme(plot.title = element_text(hjust = 0.5))
p4 = ggplot(df44, aes(x = x, y = y)) + geom_point( color = "blue", alpha = 0.3) + theme_bw() + labs(x="",y="")#+ ggtitle("Fréchet lower bound") +  theme(plot.title = element_text(hjust = 0.5))

library(patchwork)
# # 2 by 2
p1234 = p1 + p2 + p3 + p4 + plot_layout(ncol = 2, nrow = 2) #+ plot_annotation(tag_levels = "i") #+ plot_annotation(subtitle = c("Corresponding bivariate betas with Beta(1,2) marginals"), tag_levels = "i")
p1234

#ggsave(p1234, filename = "corr_beta_scatter_newlegend.pdf", width = 4, height = 4)


