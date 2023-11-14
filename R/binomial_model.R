# AUTHOR: Nicklas Kenno Hansen
# 
# DESCRIPTION:
#   Implementation of the Binomial (tree / lattice) options pricing model
#   with a helper-function for the Black-Scholes specification.
#   This
#
# NOTATION:
#   pr:     risk neutral probability
#   S0:     initial spot price of the underlying stock
#   r:      risk free rate, i.e drift
#   q:      continuously compounded dividends
#   sigma:  volatility, i.e. diffusion coefficient
#   T_:     expiry
#   M:      number of time steps
#   u:      factor for up move
#   d:      factor for down move
# 
# NOTE:
#   The function works for both European (EUR) and American (AMR) options.
#   The delta calculated is the "forward delta".
#
# REFERENCES:
#   Hull, J.    Options, Futures, and Other Derivatives
#   Seydel, R.  Tools for computational finance
#   Poulsen, R. Finance 1 and Beyond
#   https://en.wikipedia.org/wiki/Binomial_options_pricing_model


source('payoff_func.R')


binomial_tree <- function(K, T_, S0, r, q, M, u, d, payoff_func, type = "PUT", eur_amr = "EUR") {
  # Auxiliary variables
  dt <- T_ / M
  a <- exp((r - q) * dt)
  pr <- (a - d) / (u - d)
  df <- exp(-r * dt)
  
  if (dt >= sigma^2 / (r-q)^2){
    stop('Condition for positive probabilites is not satisfied: pr =', pr)
  }
  
  # Initialize stock prices and option payoff at expiry; 
  S <- S0 * d^seq(M, 0, -1) * u^seq(0, M, 1)
  V <- payoff_func(S, K, type = type)
  
  # Delta and early exercise boundary
  delta <- NaN
  B <- rep(NaN, length.out = M + 1)
  B[M+1] <- K
  
  # Backward recursion through the tree
  for (i in M:1) {
    S <- S0 * d^seq(i-1, 0, -1) * u^seq(0, i-1, 1)
    V[1:i] <- df * (pr * V[2:(i + 1)] + (1 - pr) * V[1:i])
    V <- V[1:length(V)-1]
    
    if (eur_amr == "AMR") {
      payoff <- payoff_func(S, K, type = type)
      ex <- V < payoff
      if (sum(ex) > 0) {
        B[i] <- max(S[ex])
      }
      V <- pmax(V, payoff)
    }
    
    if (i == 2) {
      delta <- (V[1] - V[2]) / (S[1] - S[2])
    }
  }
  
  res <- list(price=V[1], delta=delta, eeb=B)
  return (res)
}

binomial_tree_bs <- function(K, T_, S0, r, q, sigma, M, payoff_func, type = "P", eur_amr = "EUR") {
  dt <- T_ / M
  u <- exp(sigma * sqrt(dt))
  d <- 1 / u
  
  res <- binomial_tree(K=K, T_=T_, S0=S0, r=r, q=q, M=M, u=u, d=d, 
                       payoff_func=payoff_func, type=type, eur_amr=eur_amr)
  return(res)
}

# Main function
if (sys.nframe() == 0){
  M <- 5000
  T_ <- 0.25
  r <- 0.06
  q <- 0.0
  sigma <- 0.2
  S0 <- 40.0
  K <- 40.0
  type <- "P"
  eur_amr <- "AMR"
  
  # Calculate and print results
  res <- binomial_tree_bs(
      K=K, T_=T_, S0=S0, r=r, q=q, sigma=sigma, M=M, 
      payoff_func=payoff_european, type=type, eur_amr=eur_amr
    )
  print(c(res$price, res$delta))
  
  # Plot 
  plot(1:(M+1), res$eeb, type = "l", xlab = 'Time step', ylab = 'Price')
  abline(h = K, lwd = 1, lty=4)
  title("Early Exercise Boundary")
}



