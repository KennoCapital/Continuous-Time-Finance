# Author: Nicklas Kenno Hansen
# Date:   09-03-2023

payoff_european <- function(S, K, type='C'){
  if (type=='C'){
    return(pmax(S-K, 0))
  }
  
  if (type=='P'){
    return(pmax(K-S, 0))
  }
  
  return(NA)
}


bs_d <- function(S, K, r, sigma, t, T_){
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*(T_-t)) / (sigma*sqrt(T_-t))
  d2 <- d1 - (sigma*sqrt(T_-t))
  return(c(d1, d2))
}
 

bs_price <- function(S, K, r, sigma, t, T_, type='C'){
  
  d1 <- ifelse(t==T_, NA, (log(S/K) + (r + 0.5*sigma^2)*(T_-t)) / (sigma*sqrt(T_-t)))
  d2 <- ifelse(t==T_, NA, d1 - (sigma*sqrt(T_-t)))
  
  if (type=='C'){
    return(ifelse(t==T_, 
                  payoff_european(S=S, K=K, type=type),
                  pnorm(d1) * S - pnorm(d2) * K * exp(-r*(T_-t))
                  ))
    }
  
  if (type=='P'){
    return(ifelse(t==T_, 
                  payoff_european(S=S, K=K, type=type),
                  pnorm(-d2) * K * exp(-r*(T_-t)) - pnorm(-d1) * S
                  ))
  }
  
  return(NA)
}


bs_delta <- function(S, K, r, sigma, t, T_, type='C'){
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*(T_-t)) / (sigma*sqrt(T_-t))
  d2 <- d1 - (sigma*sqrt(T_-t))
  if (type=='C'){return(pnorm(d1))}
  if (type=='P'){return(pnorm(d1) - 1)}
  
  return(NA)
}


bs_gamma <- function(S, K, r, sigma, t, T_, type='C'){
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*(T_-t)) / (sigma*sqrt(T_-t))
  d2 <- d1 - (sigma*sqrt(T_-t))
  if (type=='C' | type=='P'){return(dnorm(d1) / (S*sigma*sqrt(T_-t)))}
  
  return(NA)
}


bs_vega <- function(S, K, r, sigma, t, T_, type='C'){
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*(T_-t)) / (sigma*sqrt(T_-t))
  d2 <- d1 - (sigma*sqrt(T_-t))
  if (type=='C' | type=='P'){return(S*dnorm(d1)*sqrt(T_-t))}

  return(NA)
}


bs_theta <- function(S, K, r, sigma, t, T_, type='C'){
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*(T_-t)) / (sigma*sqrt(T_-t))
  d2 <- d1 - (sigma*sqrt(T_-t))
  if (type=='C'){
    return(-S*dnorm(d1)*sigma/(2*sqrt(T_-t)) - r*K*exp(-r*(T_-t))*pnorm(d2))
  }
  if (type=='P'){
    return(-S*dnorm(d1)*sigma/(2*sqrt(T_-t)) + r*K*exp(-r*(T_-t))*pnorm(-d2))
  }
  
  return(NA)
}


bs_rho <- function(S, K, r, sigma, t, T_, type='C'){
  d1 <- (log(S/K) + (r + 0.5*sigma^2)*(T_-t)) / (sigma*sqrt(T_-t))
  d2 <- d1 - (sigma*sqrt(T_-t))
  if (type=='C'){return(K*(T_-t)*exp(-r*(T_-t))*pnorm(d2))}
  if (type=='P'){return(-K*(T_-t)*exp(-r*(T_-t))*pnorm(-d2))}
  
  return(NA)
}


bs_analytical <- function(S, K, r, sigma, t, T_, type='C'){
  param <- c(S, K, r, sigma, t, T_, type)
  return(
      list(
        price=bs_price(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, type=type),
        delta=bs_delta(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, type=type),
        gamma=bs_gamma(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, type=type),
        vega=bs_vega(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, type=type),
        theta=bs_theta(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, type=type),
        rho=bs_rho(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, type=type),
        param=param
        )
    )
}


bs_imp_vol <- function(price, S, K, r, t, T_, type='C', ...){
  return(
    uniroot(function(x) price-bs_price(S=S, K=K, r=r, sigma=x, t=t, T_=T_, type='C'), ...)$root
  )
}



if (sys.nframe() == 0) {
  S <- 100.0
  K <- 100.0
  r <- 0.03
  sigma <- 0.2
  t <- 0.0
  T_ <- 1.0
  
  call <- bs_analytical(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, type='C')
  put <- bs_analytical(S=S, K=K, r=r, sigma=sigma, t=t, T_=T_, type='P')
  imp_vol <- bs_imp_vol(price=call$price, S=S0, K=K, r=r, t=t, T_=T_, lower=0.05, upper=1)
}  