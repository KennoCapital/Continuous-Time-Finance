# AUTHOR: Nicklas Kenno Hansen
# 
# DESCRIPTION:
#   Implementation of standard Black-Scholes formula, Greeks, and implicit vol.
#
# NOTATION:
#   S:      spot price of the underlying stock
#   r:      risk free rate, i.e drift
#   sigma:  volatility, i.e. diffusion coefficient
#   T_:     expiry
#   t:      current time
#   q:      continuously compounded dividends (optional argument)
#   fwd:    forward price
# 
# NOTE:
#   The function `bs_analytical` provides the most common values of interest,
#   i.e. value / price and (spot) Greeks. 
#   
#   The function `bs_imp_vol` finds the implicit volatility
#   using the Bisection-method
#
# REFERENCE:
#   https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model#Instruments_paying_continuous_yield_dividends

source('payoff_func.R')


bs_fwd = function(S, r, T_, t, q=0.0){
  fwd = S * exp((r-q) * (T_-t))
  return(fwd)
}

bs_d <- function(S, K, r, sigma, t, T_, q=0.0){
  # Using `ifelse` allows for easy vectorization
  d1 <- ifelse(t >= T_, NA, (log(S/K) + (r - q + 0.5*sigma^2)*(T_-t)) / (sigma*sqrt(T_-t)))
  d2 <- ifelse(t >= T_, NA, d1 - (sigma*sqrt(T_-t)))
  
  d <- list(d1=d1, d2=d2)
  return(d)
}
 

bs_price <- function(S, K, r, sigma, t, T_, q=0.0, type='C'){
  d <- bs_d(S, K, r, sigma, t, T_, q)
  d1 <- d$d1
  d2 <- d$d2
  
  fwd = bs_fwd(S, r, T_, t, q)
  
  if (type=='C'){
    return(ifelse(t==T_, 
                  payoff_european(S=S, K=K, type=type),
                  pnorm(d1) * fwd - pnorm(d2) * K * exp(-r*(T_-t))
                  ))
    }
  
  if (type=='P'){
    return(ifelse(t==T_, 
                  payoff_european(S=S, K=K, type=type),
                  pnorm(-d2) * K * exp(-r*(T_-t)) - fwd * pnorm(-d1)
                  ))
  }
  
  return(NA)
}


bs_delta <- function(S, K, r, sigma, t, T_, q=0.0, type='C'){
  d <- bs_d(S, K, r, sigma, t, T_, q)
  d1 <- d$d1
  d2 <- d$d2
  if (type=='C'){return(pnorm(d1))}
  if (type=='P'){return(pnorm(d1) - 1)}
  
  return(NA)
}


bs_gamma <- function(S, K, r, sigma, t, T_, q=0.0, type='C'){
  d <- bs_d(S, K, r, sigma, t, T_, q)
  d1 <- d$d1
  d2 <- d$d2
  if (type=='C' | type=='P'){return(dnorm(d1) / (S*sigma*sqrt(T_-t)))}
  
  return(NA)
}


bs_vega <- function(S, K, r, sigma, t, T_, q=0.0, type='C'){
  d <- bs_d(S, K, r, sigma, t, T_, q)
  d1 <- d$d1
  d2 <- d$d2
  if (type=='C' | type=='P'){return(S*dnorm(d1)*sqrt(T_-t))}

  return(NA)
}


bs_theta <- function(S, K, r, sigma, t, T_, q=0.0, type='C'){
  d <- bs_d(S, K, r, sigma, t, T_, q)
  d1 <- d$d1
  d2 <- d$d2
  if (type=='C'){
    return(-S*dnorm(d1)*sigma/(2*sqrt(T_-t)) - r*K*exp(-r*(T_-t))*pnorm(d2))
  }
  if (type=='P'){
    return(-S*dnorm(d1)*sigma/(2*sqrt(T_-t)) + r*K*exp(-r*(T_-t))*pnorm(-d2))
  }
  
  return(NA)
}


bs_rho <- function(S, K, r, sigma, t, T_, q=0.0, type='C'){
  d <- bs_d(S, K, r, sigma, t, T_, q)
  d1 <- d$d1
  d2 <- d$d2
  if (type=='C'){return(K*(T_-t)*exp(-r*(T_-t))*pnorm(d2))}
  if (type=='P'){return(-K*(T_-t)*exp(-r*(T_-t))*pnorm(-d2))}
  
  return(NA)
}


bs_analytical <- function(S, K, r, sigma, t, T_, q=0.0, type='C'){
  param <- c(S, K, r, sigma, t, T_, q, type)
  return(
      list(
        price=bs_price(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, q=q, type=type),
        delta=bs_delta(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, q=q, type=type),
        gamma=bs_gamma(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, q=q, type=type),
        vega=bs_vega(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, q=q, type=type),
        theta=bs_theta(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, q=q, type=type),
        rho=bs_rho(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, q=q, type=type),
        param=param
        )
    )
}


bs_imp_vol <- function(price, S, K, r, t, T_, q=0.0, type='C', ...){
  obj_func <- function(x){
    return(price - bs_price(S=S, K=K, r=r, sigma=x, t=t, T_=T_, q=q, type='C'))
  }
  
  return(
    uniroot(function(x) obj_func(x), ...)$root
  )
}



if (sys.nframe() == 0) {
  S <- 100.0
  K <- 100.0
  r <- 0.03
  sigma <- 0.2
  t <- 0.0
  T_ <- 1.0
  q <- 0.0
  
  call <- bs_analytical(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, type='C')
  put <- bs_analytical(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, type='P')
  imp_vol <- bs_imp_vol(price=call$price, S=S, K=K, r=r, t=t, T_=T_, lower=0.05, upper=1.0)
}  