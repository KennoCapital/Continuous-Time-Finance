set.seed(2)

require(ggplot2)

n <- 100
T_ <- 5

x <- rnorm(n)
t = seq(from=0, to=T_, length.out=n+1)
w <- vector(mode='numeric', length=n+1)
w2 <- vector(mode='numeric', length=n+1)
dt = T_ / n

w[1] <- 0


for (i in 1:n){
  w[i+1] <- w[i] + sqrt(dt) * x[i] # Simulating THE CORRECT way
  w2[i+1] <- sqrt(t[i]) * x[i]     # Simulating THE WORNG way
}

# PLOT
data <- data.frame(t, w, w2)
ggplot(data, aes(x=t)) + 
  geom_line(aes(x=t, y=w), color='blue') + 
  geom_line(aes(x=t, y=w2), color='red') +
  theme_bw()

