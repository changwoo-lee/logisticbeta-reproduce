# Reproduce Figure 2

# Load libraries
library(forcats)
library(tidyverse)

# Source functions
source("R/polya.R")
source("R/lb.R")

make_density_logit <- function(a = 2, b = 4, rho_values = c(-0.8, 0, 0.8)){
  logit = function(x) log(x/(1-x))
  # for density plot in (0,1)x(0,1)
  xygrid_0101 = expand.grid(seq(0.005, 0.995, length = 50), seq(0.005, 0.995, length=50))
  xygrid_trans = logit(xygrid_0101)
  # density in R2
  z_trans = numeric(nrow(xygrid_0101))
  df <- data.frame()
  for(j in 1:length(rho_values)){
    for(i in 1:nrow(xygrid_0101)) z_trans[i] = dmlb(xygrid_trans[i,], a=a, b=b, R = matrix(c(1, rho_values[j], rho_values[j], 1), 2, 2))
    # convert density in (0,1)x(0,1)
    jacobian = (exp(-xygrid_trans[,1])/(exp(-xygrid_trans[,1])+1)^2 * exp(-xygrid_trans[,2])/(exp(-xygrid_trans[,2])+1)^2)^(-1)
    df = rbind(df,
               data.frame(x = xygrid_0101[,1],
                          y = xygrid_0101[,2],
                          dens = z_trans*jacobian,
                          logdens = log(z_trans) + log(jacobian),
                          dist = paste0("Logistic transformed LB(", a, ", ", b,", ", rho_values[j] ,")")))
  }
  df$dist <- fct_relevel(df$dist, unique(df$dist))
  return(df)
}

make_density_eta <- function(a = 2, b = 4, rho_values = c(-0.8, 0, 0.8)){
  logit = function(x) log(x/(1-x))
  # for density plot in (0,1)x(0,1)
  xygrid = expand.grid(seq(-3, 3, length=70), seq(-3, 3, length=70))
  # density in R2
  z = numeric(nrow(xygrid))
  df <- data.frame()
  for(j in 1:length(rho_values)){
    for(i in 1:nrow(xygrid)) z[i] = dmlb(xygrid[i,], a=a, b=b, R = matrix(c(1, rho_values[j], rho_values[j], 1), 2, 2))
    # convert density in (0,1)x(0,1)
    df = rbind(df,
               data.frame(x = xygrid[,1],
                          y = xygrid[,2],
                          dens = z,
                          dist = paste0("LB(", a, ", ", b,", ", rho_values[j] ,")")))
  }
  df$dist <- fct_relevel(df$dist, unique(df$dist))
  return(df)
}



make_plot <- function(a = 2, b = 4, rho_values = c(-0.8, 0, 0.8)){
  df_logit <- make_density_logit(a = a, b = b, rho_values = rho_values)
  df_eta <- make_density_eta(a = a, b = b, rho_values = rho_values)

  scaleFUN <- function(x) sprintf("%.1f", x)
  pal1_plasma = viridisLite::plasma(n = 100, begin = 0, end = 0.5) # or opCon = “viridis”
  pal2_plasma = viridisLite::plasma(n = 200, begin = 0.5, end = 1) # or opCon = “viridis”
  pal12_plasma = c(pal1_plasma, pal2_plasma)

  pal1 = viridisLite::viridis(n = 100, begin = 0, end = 0.5) # or opCon = “viridis”
  pal2 = viridisLite::viridis(n = 200, begin = 0.5, end = 1) # or opCon = “viridis”
  pal12 = c(pal1, pal2)

  p_logit <- ggplot(df_logit) +
    geom_raster(mapping = aes(x=x, y=y, fill = dens), interpolate = T) +
    scale_fill_gradientn("Density", colours = pal12) +
    xlab("") + ylab("") +
    expand_limits(fill = 0) +
    theme_light() +
    theme(aspect.ratio = 1, panel.grid = element_blank())+
    facet_wrap(~dist, scales = "free")+
    xlab(expression(sigma(eta[1])))+
    ylab(expression(sigma(eta[2])))

  p_eta <- ggplot(df_eta) +
    geom_raster(mapping = aes(x=x, y=y, fill = dens), interpolate = T) +
    scale_fill_gradientn("Density", colours = pal12_plasma) +
    xlab("") + ylab("") +
    expand_limits(fill = 0) +
    theme_light() +
    theme(aspect.ratio = 1, panel.grid = element_blank())+
    facet_wrap(~dist, scales = "free")+
    xlab(expression(eta[1]))+
    ylab(expression(eta[2]))+
    scale_y_continuous(labels=scaleFUN)+
    scale_x_continuous(labels=scaleFUN)
  ggpubr::ggarrange(p_eta, p_logit, ncol = 1)

}

make_plot()
ggsave("Bivariate_density_plot.pdf", width =  9.44, height = 5.75)



