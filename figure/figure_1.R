# Reproduce figure 1
rm(list = ls())

# Load libraries
library(forcats)
library(tidyverse)

# Source functions
source("R/polya.R")
source("R/lb.R")

xgrid = seq(-5, 5, length = 200)
plot(xgrid, dulb(xgrid, 3, 3), type="l", col = 1)
lines(xgrid, dulb(xgrid, 2, 4), col = 2)
lines(xgrid, dulb(xgrid, 5, 1), col = 3)


lambdagrid = seq(0.01, 2.5, length = 200)
plot(lambdagrid, dpolya(lambdagrid, 3, 3), type="l", col = 1)
lines(lambdagrid, dpolya(lambdagrid, 2, 4), col= 2)
lines(lambdagrid, dpolya(lambdagrid, 5, 1), col = 3)



df = data.frame(x = rep(xgrid,3),
                lambda = rep(lambdagrid,3),
                group = rep(c("33","24","51"), each = 200),
                #linetype = rep(c("soild","longdash","twodash"), each = 200),
                lb = c(dulb(xgrid, 3, 3),dulb(xgrid, 2, 4),dulb(xgrid, 5, 1)),
                po = c(dpolya(lambdagrid, 3, 3), dpolya(lambdagrid, 2, 4),dpolya(lambdagrid, 5, 1)))
mainplot1 = ggplot(df, aes(x=x, y = lb)) +
  geom_line(aes(col = group, linetype = group)) +
  xlab("eta") + ylab("density") + theme_light() +
  facet_wrap(~"Logistic-beta distribution") +
  scale_color_manual(name = expression(eta ~" ~ LB(a,b)"),
                     values = c("33" = "black", "24" = "red","51"="blue"),
                     labels = c("33" = "(a,b) = (3,3)", "24" = "(a,b) = (2,4)", "51" = "(a,b) = (5,1)"))+
  scale_linetype_manual(name = expression(eta ~" ~ LB(a,b)"),
                        values = c("33" = 1, "24" = 5,"51"=4),
                        labels = c("33" = "(a,b) = (3,3)", "24" = "(a,b) = (2,4)", "51" = "(a,b) = (5,1)"))+
  theme(legend.justification = c(0, 1),
        legend.position = c(0.05,0.95))+
  xlab(expression(eta))+
  ylab("Density")
mainplot1

mainplot2 = ggplot(df, aes(x=lambda, y = po)) +
  geom_line(aes(col = group, linetype = group)) +
  xlab("lambda") + ylab("density") + theme_light() +
  facet_wrap(~"Polya distribution") +
  scale_color_manual(name = expression(lambda ~" ~ Polya(a,b)"),
                     values = c("33" = "black", "24" = "red","51"="blue"),
                     labels = c("33" = "(a,b) = (3,3)", "24" = "(a,b) = (2,4) or (4,2)", "51" = "(a,b) = (5,1) or (1,5)"))+
  scale_linetype_manual(name =expression(lambda ~" ~ Polya(a,b)"),
                        values = c("33" = 1, "24" = 5,"51"=4),
                        labels = c("33" = "(a,b) = (3,3)", "24" = "(a,b) = (2,4) or (4,2)", "51" = "(a,b) = (5,1) or (1,5)"))+
  theme(legend.justification = c(1, 1), legend.position = c(0.95,0.95))+
  xlab(expression(lambda))+
  ylab("Density")
mainplot2


ggpubr::ggarrange(mainplot1, mainplot2, ncol = 2)
ggsave("~/Documents/ResearchProjects/logisticbeta-reproduce/LB_Polya_density_plot.pdf", width =  10.26, height = 2.92)




