#Inverse Gamma
install.packages("invgamma")
library(invgamma)
X1 <- rinvgamma(100, 0.5,10)

#Laplace
library(rmutil)
X2 <- rlaplace(100, 5, 10)

#Exponential
X3 <- rexp(100, 0.1)


#Inverse-Gaussian Distribution
install.packages("statmod")
library(statmod)
X4 <- rinvgauss(100, 10,1 )



#Accept-reject sampling from truncated normal
# density of target function: normal density truncated to region theta > k
p.density <- function(theta, k){
  ifelse(theta>k, exp(-theta^2 / 2), 0)
}

# majorizing function: translated exponential density
m.density <- function(theta, k, alpha){
  alpha * exp( -alpha * (theta - k))
}

# sampling from a standard normal density truncated to lie above k=1

T <- 100000 # number of random draws
theta <- rep(NA, T) # vector to be filled with random draws
alpha <- 1.62 # rescaling coefficient 
k <- 1 # truncation point

for (t in 1:T){
  z <- rexp(1, rate=alpha) + k # majorizing function (translated exponential)
  u <- runif(1, 0, 1) # will be used later to accept with probability r
  r <- p.density(z, k) / m.density(z, k, alpha) # accept probability
  if (u <= r){
    theta[t] <- z
  }
  else {
    next
  }
}


