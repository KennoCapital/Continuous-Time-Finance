# Author: Nicklas Kenno Hansen
# Date:   08-Mar-2023

require(ggplot2)
require(reshape2)
require(gridExtra)
require(gtable)

# ---------------------------------------------------------------------------- #
# Function for measuring the time it takes to execute a function call
# ---------------------------------------------------------------------------- #
exec_time <- function(func_call){
  time_start <- Sys.time()
  func_call
  return(Sys.time() - time_start)
}

# ---------------------------------------------------------------------------- #
# Function for plotting Monte Carlo simulations
# ---------------------------------------------------------------------------- #

plot_mc <- function(A, T_, M, title='', legend=F){
  # A:      Matrix, each column is a realization
  # T_:     Terminal time
  # M:      Number of steps
  # Title:  Title of plot
  
  df <- as.data.frame(A)
  df['t'] <- T_ * 1:(M+1)/(M+1)
  df <- melt(df, id.vars='t', variable.name = 'Simulation')
  p <- ggplot(df, aes(t, value)) +
    geom_line(aes(colour = Simulation), size=0.75) +
    ggtitle(title) +
    theme_bw()
  
  if (!(legend)){
    p <- p + theme(legend.position="none")
  }
  
  return(p)
}


# ---------------------------------------------------------------------------- #
# Algorithm 1.11
# Euler discretization of an SDE
# ---------------------------------------------------------------------------- #

sde_euler <- function(t0, x0, T_, M, N, a, b, seed=0){
  # Approximate the Ito process X given by the SDE 
  # dX(t) = a(X(t),t) * dt + b(X(t),t) * dW(t)
  
  # t0: Initial time
  # x0: Initial value of X
  # T_: Expiry
  # M:  Number of steps
  # N:  Number of simulations
  # a:  Drift function
  # b:  Diffusion function
  
  if (seed > 0){set.seed(seed)}
  
  # Initilaztion
  dt = T_ / M
  t <- vector(); t[1] <- t0
  y <- matrix(data=NA, nrow=M+1, ncol=N); y[1, ] <- x0
  
  # Simulate
  for (j in 1:M){
    t[j+1] <- t[j] + dt
    Z <- rnorm(N, mean=0, sd=1)
    dW <- Z * sqrt(dt)
    y[j+1,] <- y[j,] + a(y[j,], t[j]) * dt + b(y[j,], t[j]) * dW
  }
  return(y)
}

# ---------------------------------------------------------------------------- #
# Simulation SDE using Milstein (algorithm 3.5, page 134) 
#                  and Runge-Kutta (eq. 3.15, page 136)
# ---------------------------------------------------------------------------- #

sde_milstein <- function(t0, x0, T_, M, N, a, b, seed=0){
  # Milstein simulation of SDE using RUnge-Kutta
  
  # t0:       Initial time
  # x0:       Initial value of X
  # T_:       Expiry
  # M:        Number of steps
  # N:        Number of simulations
  # a:        Drift function
  # b:        Diffusion function
  
  if (seed > 0){set.seed(seed)}
  
  # Initilaztion
  dt = T_ / M
  t <- vector(); t[1] <- t0
  y <- matrix(data=NA, nrow=M+1, ncol=N); y[1, ] <- x0
  
  # Simulate
  for (j in 1:M){
    t[j+1] <- t[j] + dt
    Z <- rnorm(N, mean=0, sd=1)
    dW <- Z * sqrt(dt)
    y_hat <- y[j, ] + a(y[j, ], t[j]) * dt + b(y[j, ], t[j]) * sqrt(dt)
    y[j+1, ] <- y[j, ] + a(y[j, ], t[j]) * dt + b(y[j, ], t[j]) * dW +
      1/(2*sqrt(dt)) * (dW^2 - dt)*(b(y_hat, t) - b(y[j, ], t))
  }
  return(y)
}


# ---------------------------------------------------------------------------- #
# Simulation of "exact" GBM
# page 147
# ---------------------------------------------------------------------------- #

sim_gbm <- function(t0, x0, T_, M, N, mu, sigma, seed=0){
  # "Exact" simulation of GBM
  # Page 147
  # X(t) = x(0) * exp{(mu-0.5*sigma^2)*dt + sigma*dW(t)}
  
  # t0: Initial time
  # x0: Initial value of X
  # T_: Expiry
  # M:  Number of steps
  # N:  Number of simulations
  # a:  Drift coefficient
  # b:  Diffusion coefficient
  
  if (seed > 0){set.seed(seed)}
  
  # Initilaztion
  dt = T_ / M
  t <- vector(); t[1] <- t0
  x <- matrix(data=NA, nrow=M+1, ncol=N); x[1, ] <- x0
  
  # Simulate
  for (j in 1:M){
    t[j+1] <- t[j] + dt
    Z <- rnorm(N, mean=0, sd=1)
    dW <- Z * sqrt(dt)
    x[j+1, ] <- x[j, ] * exp((mu-0.5*sigma^2)*dt + sigma*dW)
  }
  return(x)
}

# ---------------------------------------------------------------------------- #
# Mean-reverting volatility tandem
# Note: Only 1 realization is being sampled
# ---------------------------------------------------------------------------- #

sde_mrvt <- function(t0, s0, T_, M, sigma0, zeta0, alpha, beta, seed=0){
  # SDE simulation of mean-reverting volatility tandem
  # Example 1.15
  if (seed > 0){set.seed(seed)}
  
  dt = T_ / M
  t <- vector(); t[1] <- t0
  S <- vector(); S[1] <- s0
  sigma <- vector(); sigma[1] <- sigma0
  zeta <- vector(); zeta[1] <- zeta0
  
  for (j in 1:M){
    t[j+1] <- t[j] + dt
    Z <- rnorm(2, mean=0, sd=1)
    dW <- Z * sqrt(dt)
    
    S[j+1] <- S[j] + sigma[j] * S[j] * dW[1]
    sigma[j+1] <- sigma[j] - (sigma[j]-zeta[j]) * dt + alpha * sigma[j] * dW[2]
    zeta[j+1] <- zeta[j] + beta * (sigma[j] - zeta[j])* dt
  }
  
  output <- list(t=t, S=S, sigma=sigma, zeta=zeta)
  return(output)
}



# ---------------------------------------------------------------------------- #
# Problem 3.1
# ---------------------------------------------------------------------------- #

if (sys.nframe() == 0) {

  # Part A (with Rolf's bells & whistles)
  
  # Set parameters
  t0 <- 0.0
  x0 <- 50.0
  T_ <- 5.0
  M <- 250
  N <- 10
  a = function(x, t){0.07*x}
  b = function(x, t){0.2*x}
  mu <- a(1, 0)
  sigma <- b(1, 0)
  seed <- 1
  
  # SDE for log of GBM:
  # d[log(X(t))] = (mu-0.5*sigma^2)*dt + sigma*dW(t)
  a_log <- function(x,t){mu-0.5*sigma^2}
  b_log <- function(x,t){sigma}
  
  
  mc_euler_sde <- sde_euler(t0=t0, x0=x0, T_=T_, M=M, N=N, a=a, b=b, seed=seed)
  p_euler <- plot_mc(A=mc_euler_sde, T_=T_, M=M, title='Euler Discretization of SDE')
  
  mc_milstein <- sde_milstein(t0=t0, x0=x0, T_=T_, M=M, N=N, a=a, b=b, seed=seed)
  p_milstein <- plot_mc(A=mc_milstein, T_=T_, M=M, title='Milstein simulation of SDE with Runge-Kutta')
  
  mc_log_sde <- exp(sde_euler(t0=t0, x0=log(x0), T_=T_, M=M, N=N, a=a_log, b=b_log, seed=seed))
  p_log <- plot_mc(A=mc_log_sde, T_=T_, M=M, title='Exp{Log dynamics for Euler SDE}')
    
  mc_exact_gbm <- sim_gbm(t0=t0, x0=x0, T_=T_, M=M, N=N, mu=mu, sigma=sigma, seed=seed)
  p_exact_gbm <- plot_mc(A=mc_exact_gbm, T_=T_, M=M, title='Simulation using "Exact" GBM')
  
  grid.arrange(p_euler, p_milstein, p_log, p_exact_gbm, nrow=2)
  
  # Mean Absolute Error
  MAE <- list()
  MAE$euler <- mean(abs(mc_euler_sde - mc_exact_gbm))
  MAE$milstein <- mean(abs(mc_milstein - mc_exact_gbm))
  MAE$exp_of_log_sde <- mean(abs(mc_log_sde - mc_exact_gbm))
  MAE
  
  
  # Part B 
  
  # Set additional parameters
  s0 <- 50
  alpha <- 0.3
  beta <- 10
  sigma0 <- 0.2
  zeta0 <- 0.2
  data <- as.data.frame(
    sde_mrvt(t0=t0, s0=s0, T_=T_, M=M, sigma0=sigma0, zeta0=zeta0, alpha=alpha, beta=beta, seed=seed)
  )
  
  # Plot
  p1 <- ggplot(data, aes(x=t)) +
    geom_line(aes(y=S, color='S'), size=0.75) +
    scale_color_manual(name = "Underlying asset (Stock)",
                       values = c("S" = "black")) +
    ylab('') +
    xlab('') +
    ggtitle('Mean-Reverting Volatility Tandem') +
    theme_bw()
  
  p2 <- ggplot(data, aes(x=t)) +
    geom_line(aes(y=sigma, color='sigma'), size=0.75) +
    geom_line(aes(y=zeta, color='zeta'), size=0.75) +
    scale_color_manual(name = "Volatility processes",
                       values = c("sigma" = "blue", "zeta" = "red")) +
    ylab('') +
    theme_bw() 
  
  grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "first"))
   
  
  # Measure execution time of 50K simulations
  N <- 50000
  
  exec_time(sde_euler(t0=t0, x0=x0, T_=T_, M=M, N=N, a=a, b=b, seed=seed))
  exec_time(sde_milstein(t0=t0, x0=x0, T_=T_, M=M, N=N, a=a, b=b, seed=seed))
  exec_time(exp(sde_euler(t0=t0, x0=log(x0), T_=T_, M=M, N=N, a=a_log, b=b_log, seed=seed)))
  exec_time(sim_gbm(t0=t0, x0=x0, T_=T_, M=M, N=N, mu=mu, sigma=sigma, seed=seed))

}
