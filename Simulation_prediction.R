library(mvtnorm)
library(MASS)
library(truncnorm) ## rtruncnorm(N0,mean=mu.z[y.data==0],sd=1,a=-Inf,b=0)
library(statmod) ## rinvgauss(10,1,1)
library(VGAM) ## rinv.gaussian(10,1,1)
library(LearnBayes)
## Prepare Simulation Data ##
set.seed(100)
n=1024
d=15
x1<-rnorm(n,3,1)
x2<-rnorm(n,x1,1)
x3<-rnorm(n,x2,2)
x4<-runif(n,5,10)
x5<-runif(n,x4,x4+3)
x6<-rnorm(n,3.5,1)
x7<-rnorm(n,x6,1)
x8<-rnorm(n,5.2,2)
x9<-runif(n,x8,x8+3)
x10<-runif(n,x9,x9+1)
x11<-rnorm(n,5,1)
x12<-rnorm(n,x11,1)
x13<-rnorm(n,x12,2)
x14<-runif(n,5,10)
x15<-runif(n,x14,x14+3)
X<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15)
X[1:5,]
head(X)
summary(x8)
summary(x15)
apply(X,2,mean)
beta<-c(4,-4,5,7,-6,rep(0,10));beta
length(beta)
Xb<-X%*%beta                               
Xb
# Obtain the vector with probabilities of success p using the probit link
p<-pnorm(Xb)
p
y.data<-rbinom(1024,1,p)
y.data
length(y.data)
simulation.data <- data.frame(X,y.data)
head(simulation.data)

split.sim <- 0.8
sample.sim <- floor(split.sim*nrow(simulation.data)) #training data size = 819

train_indices.sim <- sample(seq_len(nrow(simulation.data)), size = sample.sim)
simulation.train <- simulation.data[train_indices.sim, ]    # These are 80% of data
simulation.test <- simulation.data[-train_indices.sim, ]    # These are rest of them, 20%.  
simulation.train <- as.data.frame(simulation.train)

#X matrix of training data
x.train <- cbind(simulation.train$x1,simulation.train$x2,simulation.train$x3,simulation.train$x4,simulation.train$x5,simulation.train$x6,
                 simulation.train$x7,simulation.train$x8,simulation.train$x9,simulation.train$x10,simulation.train$x11,simulation.train$x12,
                 simulation.train$x13,simulation.train$x14,simulation.train$x15)
y.train <- simulation.train$y.data

fit.sim <- glm(y.train~0+x.train,family=binomial(link=probit),data = simulation.train)
summary(fit.sim)
fit.sim$coefficients
n1 <- nrow(simulation.train)
## Z
z<-rep(NA,n1)
N.sim<-270000
N1<-sum(y.train);N1
N0<-n1-N1;N0
## beta
beta.sim<-matrix(NA,nrow=N.sim+1,ncol=d)
beta.sim[1,] <- c(fit.sim$coefficients)
beta.sim
## tau
tau.2.sim<-matrix(NA,nrow=N.sim+1,ncol=d)
tau.2.sim[1,]<-c(rep(1,d))
tau.2.sim
## lamda
lamda.sim<-rep(NA,n1)
lamda.sim[1]<-2
lamda.sim
## sampling
for (j in 1:N.sim){
  mu.z<-x.train%*%beta.sim[j,]
  z[y.train==0]<-rtruncnorm(N0,mean=mu.z[y.train==0],sd=1,a=-Inf,b=0)
  z[y.train==1]<-rtruncnorm(N1,mean=mu.z[y.train==1],sd=1,a=0,b=Inf)
  tau<-diag(tau.2.sim[j,])
  E<-solve(solve(tau)+t(x.train)%*%x.train)
  beta.sim[j+1,]<-rmvnorm(1,E%*%t(x.train)%*%z,E)
  tau.2.sim[j+1,]<-1/rinv.gaussian(d,(lamda.sim[j]^2/beta.sim[j+1,]^2)^0.5,lamda.sim[j]^2)
  lamda.sim[j+1]<-(rgamma(1,d+1,0.5*sum(tau.2.sim[j+1,])))^0.5
}
## estimation
summary(lamda.sim)
plot(lamda.sim)## very stable
z
beta.sim[,1]
tau.2.sim[,1]
ind<-seq(from=20000,to=270000,by=125)
summary(lamda.sim[ind])##estimate lambda
sd(lamda.sim[ind])
ts.plot(lamda.sim[ind]);acf(lamda.sim[ind])
beta.est.sim <- colMeans(beta.sim[ind,])##estimate beta
tau.est.sim <- colMeans(tau.2.sim[ind,])##estimate tau
newb.sim<-beta.sim[ind,]
sd.beta.sim<-rep(NA,15)
for (k in 1:15){
  sd.beta.sim[k]<-sd(newb.sim[,k])
}
sd.beta.sim##sd of beta estimates

simulation.test
#X matrix of testing data
x.test.sim <- cbind(simulation.test$x1,simulation.test$x2,simulation.test$x3,simulation.test$x4,simulation.test$x5,simulation.test$x6,
                 simulation.test$x7,simulation.test$x8,simulation.test$x9,simulation.test$x10,simulation.test$x11,simulation.test$x12,
                 simulation.test$x13,simulation.test$x14,simulation.test$x15)
y.test.sim <- simulation.test$y.data
length(y.test.sim)
beta.est.sim
XB.sim <- x.test.sim%*%beta.est.sim
XB.sim
y.sim.prob <- pnorm(XB.sim)
y.sim.prob
y.sim.pred <- ifelse(y.sim.prob > 0.5, "1","0")

y.sim.pred

accuracy.sim <- (sum(y.sim.pred==y.test.sim)) / length(y.test.sim); accuracy.sim



quantile(beta.sim[ind,1],c(0.025,0.5,0.975))
quantile(beta.sim[ind,2],c(0.025,0.5,0.975))
quantile(beta.sim[ind,3],c(0.025,0.5,0.975))
quantile(beta.sim[ind,4],c(0.025,0.5,0.975))
quantile(beta.sim[ind,5],c(0.025,0.5,0.975))
quantile(beta.sim[ind,6],c(0.025,0.5,0.975))
quantile(beta.sim[ind,7],c(0.025,0.5,0.975))
quantile(beta.sim[ind,8],c(0.025,0.5,0.975))
quantile(beta.sim[ind,9],c(0.025,0.5,0.975))
quantile(beta.sim[ind,10],c(0.025,0.5,0.975))
quantile(beta.sim[ind,11],c(0.025,0.5,0.975))
quantile(beta.sim[ind,12],c(0.025,0.5,0.975))
quantile(beta.sim[ind,13],c(0.025,0.5,0.975))
quantile(beta.sim[ind,14],c(0.025,0.5,0.975))
quantile(beta.sim[ind,15],c(0.025,0.5,0.975))
## plots check
par(mfrow=c(2,2))
ts.plot(lamda.sim);acf(lamda.sim,lag=5000)
ts.plot(beta.sim[,1]);acf(beta.sim[,1],lag=500)
ts.plot(beta.sim[,5]);acf(beta.sim[,5])
ts.plot(tau.2.sim[,1]);acf(tau.2.sim[,1])
ts.plot(tau.2.sim[,5]);acf(tau.2.sim[,5])
## after burn in and slicing plots
ts.plot(beta.sim[ind,1]);acf(beta.sim[ind,1],lag=100)

