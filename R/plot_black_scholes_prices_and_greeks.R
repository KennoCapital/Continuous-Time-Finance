# Author: Nicklas Kenno Hansen
# Date:   14-03-2023

rm(list=ls())

require(dplyr)
require(gridExtra)
require(ggplot2)
require(colorspace)

source('black_scholes.R')


# Parameter values
S <- seq(50.0, 150.0, by=1)
K <- 100.0
r <- c(-0.01, 0.0, 0.01, 0.03, 0.05, 0.1)
sigma <- seq(from=0.1, to=0.5, by=0.1)
t <- 0.0
T_ <- c(0.01, 0.1, 0.25, 0.5, 1.0, 2.0)
type <- c('C', 'P')


# Set "default values"
default_r <- 0.03
default_sigma <- 0.2
default_T <- 1.0


# Setup grid
df <- expand.grid(S, K, r, sigma, t, T_, type)
colnames(df) <- c('S', 'K', 'r', 'sigma', 't', 'T_', 'type')


# Calculate option price and greeks
df <- df %>% rowwise() %>% mutate(
  moneyness = if(type=='C'){S/K} else{K/S},
  price = bs_price(S, K, r, sigma, t, T_, type),
  delta = bs_delta(S, K, r, sigma, t, T_, type),
  gamma = bs_gamma(S, K, r, sigma, t, T_, type),
  vega = bs_vega(S, K, r, sigma, t, T_, type),
  theta = bs_theta(S, K, r, sigma, t, T_, type),
  rho = bs_rho(S, K, r, sigma, t, T_, type)
)

# Function for making each subplot
subplot <- function(df, x_ax, y_ax, color_by){
  p <- ggplot(df, aes(x=!!sym(x_ax))) +
    geom_point(aes(y=!!sym(y_ax), color=!!sym(color_by))) +
    facet_wrap(vars(type)) +
    scale_color_continuous_sequential("ag_Sunset") + 
    theme_bw()
  return(p)
}

# Vary `T`, `sigma`, and `r` one at the time
df_vary_T <- df %>% filter(
  (r == default_r & sigma == default_sigma)
)

df_vary_sigma <- df %>% filter(
  (r == default_r & T_ == default_T)
)

df_vary_r <- df %>% filter(
  (T_ == default_T & sigma == default_sigma)
)

# Make lists for storing plots
p_vary_T <- list()
p_vary_sigma <- list()
p_vary_r <- list()

# Make each subplot
greeks <- c("price", "delta", "gamma", "vega", "theta", "rho")
for (greek in greeks){
  p_vary_T[[greek]] = subplot(df=df_vary_T, x_ax='S', y_ax=paste0(greek), color_by='T_')
  p_vary_sigma[[greek]] = subplot(df=df_vary_sigma, x_ax='S', y_ax=paste0(greek), color_by='sigma')
  p_vary_r[[greek]] = subplot(df=df_vary_r, x_ax='S', y_ax=paste0(greek), color_by='r')
}

# Plot the results
do.call("grid.arrange", c(p_vary_T, ncol=2))
do.call("grid.arrange", c(p_vary_sigma, ncol=2))
do.call("grid.arrange", c(p_vary_r, ncol=2))

