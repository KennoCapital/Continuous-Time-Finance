# Author: Nicklas Kenno Hansen
# Date:   26-03-2023

rm(list=ls())

source('black_scholes.R')
source('SDE_and_GBM_simulations.R')

require(latex2exp)
require(ggplot2)
require(gridExtra)
require(grid)

delta_hedge <- function(t0, S0, K, T_, M, N, mu, sigma, r, type='C', seed=0, epsilon=0.0, sigma_imp=NA, sigma_hedge=NA){
  # Perform delta-hedge experiment for replicating the payoff of an European option:
  #   We sell (short) the option for `price` and hedge with the portfolio `V`
  #   Hence,  PnL = V - price 
  
  # PARAMETERS
  # t0:             Initial time
  # S0:             Initial value of S (the stock)
  # K:              Strike
  # T_:             Expiry
  # M:              Number of steps
  # N:              Number of simulations
  # mu:             Drift coefficient
  # sigma:          Diffusion coefficient (volatility)
  # r:              Risk-free rate
  # type:           Option type
  # epsilon:        Arbitrage made immediately on the option (at time t=0).
  # sigma_imp:      Volatility implied by the market (standard is `sigma`, but we can choose otherwise)
  # sigma_hedge:    Volatility used to hedge (standard is `sigma`, but we can choose otherwise)
  
  # RETURN (list containing):
  #   t       Vector of time steps
  #   S       Matrix of simulated stocks
  #   V       Matrix of the value of the hedging portfolio
  #   price   Matrix of the BS-price of the option for all t and simulated S(t)
  #   PnL     Matrix of the profit and loss (PnL = V-price)
  #   h_s     Matrix of the number of stocks hold (i.e. the delta)
  #   h_b     Matrix of the amount held in the risk-free asset (the back account)
  
  
  if (seed > 0){set.seed(seed)}
  if (is.na(sigma_hedge)){sigma_hedge <- sigma}
  if (is.na(sigma_imp)){sigma_imp <- sigma}
  
  dt <- T_/M
  
  # Simulate stock paths, calculate the analytical prices and deltas,
  t <- matrix(seq(from=t0, to=T_, by=dt), M+1, N)
  S <- sim_gbm(t0=t0, x0=S0, T_=T_, M=M, N=N, mu=mu, sigma=sigma, seed=seed)
  price <- bs_price(S=S, K=K, r=r, sigma=sigma_imp, t=t, T_=T_, type=type)
  delta <- bs_delta(S=S, K=K, r=r, sigma=sigma_hedge, t=t, T_=T_, type=type)

  # Initialize objects for calculating the
  #    value of the hedge portfolio (V)
  #    no. of stocks in the portfolio (h_s), amount of cash in the bank (h_b)
  V <- matrix(NA, M+1, N); V[1, ] <- price[1, ] + epsilon
  h_s <- matrix(NA, M+1, N); h_s[1, ] <- delta[1, ]
  h_b <- matrix(NA, M+1, N); h_b[1, ] <- price[1, ] - h_s[1, ] * S[1, ] + epsilon
  
  # Perform hedging experiment
  for (j in 2:(M+1)){
    V[j, ] <- h_s[j-1, ] * S[j, ] + h_b[j-1, ] * exp(dt*r)
    h_s[j, ] <- delta[j, ]
    h_b[j, ] <- V[j, ] - h_s[j, ]*S[j, ]
  }
  
  PnL = V - price
  
  return(list(t=t[, 1], S=S, V=V, price=price, PnL=PnL, h_s=h_s, h_b=h_b))

}


plot_compare_hedge_payoff <- function(V, S, K, sigma, sigma_hedge, M, type='C'){
  df <- data.frame(V=V, S=S, Option=payoff_european(S=S, K=K, type=type))
  
  if (is.na(sigma_hedge)){sigma_hedge <- sigma}
  
  p <- ggplot(df, aes(x=S)) + 
    geom_point(aes(y=V, color='Delta_hedge'), alpha=0.75, size=2) +
    geom_line(aes(y=Option, color='Payoff'), linewidth=1) +
    ggtitle(paste0('type=', type, ', M=',M, ', sigma=', sigma, ', sigma_hedge=', sigma_hedge)) +
    scale_color_manual(values=c(Delta_hedge='blue', Payoff='black'),
                       name = "") +
    theme_bw()
  
  return(p) 
}


if (sys.nframe() == 0) {
  # Set parameters
  t0 <- 0.0
  S0 <- 100
  K <- 100 
  T_ <- 0.25
  M <- 1000 
  N <- 5000 
  mu <- 0.07 
  r <- 0.03
  epsilon <- 0.0
  seed <- 1
  epsilon <- 0.0
  type <- 'C'
  
  sigma <- 0.2 
  
  
  # Use correct sigma to hedge
  sigma_hedge <- 0.2
  hedge_1 <- delta_hedge(t0=t0, S0=S0, K=K, T_=T_, M=M, N=N, mu=mu, sigma=sigma, r=r, type=type, seed=seed, epsilon=epsilon, sigma_hedge=sigma_hedge)
  p1 <- plot_compare_hedge_payoff(V=hedge_1$V[M+1, ], S=hedge_1$S[M+1, ], K=K, sigma=sigma, sigma_hedge=sigma_hedge, M=M, type=type)
  
  
  # Use wrong sigma to hedge
  sigma_hedge <- 0.3
  hedge_2 <- delta_hedge(t0=t0, S0=S0, K=K, T_=T_, M=M, N=N, mu=mu, sigma=sigma, r=r, type=type, seed=seed, epsilon=epsilon, sigma_hedge=sigma_hedge)
  p2 <- plot_compare_hedge_payoff(V=hedge_2$V[M+1, ], S=hedge_2$S[M+1, ], K=K, sigma=sigma, sigma_hedge=sigma_hedge, M=M, type=type)
  
  grid.arrange(p1, p2, nrow=2)

  # ----------------------------------------------- #
  #                                                 #
  # What if we made an immediately arbitrage        #
  #                                                 #
  # ----------------------------------------------- #
  
  epsilon <- 1.982
  sigma <- 0.2
  hedge_freq <- c(5, 10, 25, 50, 75, 100, 250, 500, 750, 1000)
  
  # Perform hedging experiment for different hedge frequencies
  PnL <- matrix(NA, nrow=length(hedge_freq), ncol=N)
  for (i in 1:length(hedge_freq)){
    res <- delta_hedge(t0=t0, S0=S0, K=K, T_=T_, M=hedge_freq[i], N=N, mu=mu, sigma=sigma, r=r, type=type, seed=seed, epsilon=epsilon)
    PnL[i, ] <- tail(res$PnL, 1)
  }
  
  # Plot some analytics of the PnL
  df = data.frame(M=hedge_freq,
                  avg=apply(PnL, 1 ,mean),
                  std=apply(PnL, 1 ,sd),
                  freq_loss=apply(PnL, 1, function(x) sum(x < 0) / N)
  )
  
  p3 <- ggplot(df, aes(M, avg)) + geom_point(size=2.5) +
    theme_bw() + 
    theme(text = element_text(size=14),
          axis.text = element_text(size=14)) +
    geom_hline(yintercept=epsilon*exp(T_*r), color='blue', linetype='dashed', linewidth=1, alpha=0.8) +
    labs(x='No. of hedge points (M)', y=TeX('$E(V_n(T))$'), 
         title='Avergage PnL for different hedge frequencies',
         subtitle = paste0('(',N, ' simulations)'))
  
  p4 <- ggplot(df, aes(M, std)) + geom_point(size=2.5) +
    theme_bw() +
    geom_smooth(method = lm, se=FALSE, size=1, alpha=0.8) +
    theme_bw() +
    xlab("No. hedgepoints") +
    ylab("Error") +
    scale_y_log10() +
    scale_x_log10() +
    theme(text = element_text(size=14),
          axis.text = element_text(size=14)) +
    labs(x='No. of hedge points (M)', y=TeX('$\\sqrt{var^p(V_n(T))}$'), 
         title='Standard deviation of PnL for different hedge frequencies',
         subtitle = paste0('(',N, ' simulations)'))
  
  p5 <- ggplot(df, aes(M, freq_loss)) + geom_point(size=2.5) +
    geom_hline(yintercept=0.0, color='blue', linetype='dashed', linewidth=1, alpha=0.8) +
    theme_bw() +
    theme(text = element_text(size=14),
          axis.text = element_text(size=14)) +
    labs(x='No. of hedge points (M)', y=TeX('$P(V_n(T) < 0)$'), 
         title='Frequency of loss for different hedge frequencies',
         subtitle = paste0('(',N, ' simulations)'))
  
  grid.arrange(p3, p4, p5, ncol = 1)
  
  
  # ----------------------------------------------- #
  #                                                 #
  # Willmott's Hedgeing Expirement                  #
  #                                                 #
  # ----------------------------------------------- #
  
  # We sell the option (at a premium, measured by implied volatility) and delta hedge

  epsilon <- 0.0
  N <- 10

  # Hedge with "model" volatility
  sigma_model <- 0.2
  sigma_imp <- 0.3
  sigma_hedge <- 0.2
  arb <- bs_price(S=S0, K=K, r=r, sigma=sigma_imp, t=t0, T_=T_, type=type) - bs_price(S=S0, K=K, r=r, sigma=sigma_model, t=t0, T_=T_, type=type)
  
  willmott <- delta_hedge(t0=t0, S0=S0, K=K, T_=T_, M=M, N=N, mu=mu, sigma=sigma_model, r=r, type=type, seed=seed, epsilon=epsilon, sigma_imp=sigma_imp, sigma_hedge=sigma_hedge)
  df <- data.frame(cbind(t=willmott$t, willmott$PnL))
  df <- melt(df, id.vars='t', value.name='PnL')
  p6 <- ggplot(df) + geom_line(aes(x=t, y=PnL, color=variable)) + theme_bw() + theme(legend.position="none") + ggtitle('Hedge with "model" volatility') +
    geom_point(aes(x=T_, y=arb), color='black', size=3) + ylim(-0.2, 3.3)
  
  
  # Hedge with "market implied" volatility
  sigma_model <- 0.2
  sigma_imp <- 0.3
  sigma_hedge <- 0.3
  
  willmott <- delta_hedge(t0=t0, S0=S0, K=K, T_=T_, M=M, N=5, mu=mu, sigma=sigma_model, r=r, type=type, seed=seed, epsilon=epsilon, sigma_imp=sigma_imp, sigma_hedge=sigma_hedge)
  df <- data.frame(cbind(t=willmott$t, willmott$PnL))
  df <- melt(df, id.vars='t', value.name='PnL')
  p7 <- ggplot(df) + geom_line(aes(x=t, y=PnL, color=variable)) + theme_bw() + theme(legend.position="none") + ggtitle('Hedge with "implied" volatility') + ylim(-0.2, 3.3)
 
  grid.arrange(p6, p7, ncol=2) 
  
  
}

