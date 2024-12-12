##################################
## Title: Bayesian Lasso Approch #
## Name: Vindyani Herath #########
## Start Date: August 17, 2020 ###
## Last Updated: March 27, 2021 ##
##################################
library(mvtnorm)
library(MASS)
library(truncnorm) 
library(statmod) 
library(VGAM) 
library(LearnBayes)
## Simulation Study ##
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
## find some initial betas
fit<-glm(y.data~0+X,family=binomial(link=probit))
summary(fit) 
z<-rep(NA,n)
N.sim<-270000
N1<-sum(y.data);N1
N0<-n-N1;N0
## beta
beta<-matrix(NA,nrow=N.sim+1,ncol=d)
beta[1,]<-c(fit$coefficients)
beta
## tau
tau.2<-matrix(NA,nrow=N.sim+1,ncol=d)
tau.2[1,]<-c(rep(1,d))
tau.2
## lamda
lamda<-rep(NA,n)
lamda[1]<-2
lamda
## sampling
for (j in 1:N.sim){
  mu.z<-X%*%beta[j,]
  z[y.data==0]<-rtruncnorm(N0,mean=mu.z[y.data==0],sd=1,a=-Inf,b=0)
  z[y.data==1]<-rtruncnorm(N1,mean=mu.z[y.data==1],sd=1,a=0,b=Inf)
  tau<-diag(tau.2[j,])
  E<-solve(solve(tau)+t(X)%*%X)
  beta[j+1,]<-rmvnorm(1,E%*%t(X)%*%z,E)
  tau.2[j+1,]<-1/rinv.gaussian(d,(lamda[j]^2/beta[j+1,]^2)^0.5,lamda[j]^2)
  lamda[j+1]<-(rgamma(1,d+1,0.5*sum(tau.2[j+1,])))^0.5
}
## estimation
summary(lamda)
plot(lamda)
beta[,1]
tau.2[,1]
ind<-seq(from=20000,to=270000,by=125)
summary(lamda[ind])##estimate lambda
sd(lamda[ind])
ts.plot(lamda[ind]);acf(lamda[ind])
beta.x8Changed <- colMeans(beta[ind,])##estimate beta
colMeans(tau.2[ind,])##estimate tau
newb<-beta[ind,]
sd.beta<-rep(NA,15)
for (k in 1:15){
  sd.beta[k]<-sd(newb[,k])
}
sd.beta##sd of beta estimates
quantile(beta[ind,1],c(0.025,0.5,0.975))
quantile(beta[ind,2],c(0.025,0.5,0.975))
quantile(beta[ind,3],c(0.025,0.5,0.975))
quantile(beta[ind,4],c(0.025,0.5,0.975))
quantile(beta[ind,5],c(0.025,0.5,0.975))
quantile(beta[ind,6],c(0.025,0.5,0.975))
quantile(beta[ind,7],c(0.025,0.5,0.975))
quantile(beta[ind,8],c(0.025,0.5,0.975))
quantile(beta[ind,9],c(0.025,0.5,0.975))
quantile(beta[ind,10],c(0.025,0.5,0.975))
quantile(beta[ind,11],c(0.025,0.5,0.975))
quantile(beta[ind,12],c(0.025,0.5,0.975))
quantile(beta[ind,13],c(0.025,0.5,0.975))
quantile(beta[ind,14],c(0.025,0.5,0.975))
quantile(beta[ind,15],c(0.025,0.5,0.975))
## plots check
par(mfrow=c(2,2))
ts.plot(lamda);acf(lamda,lag=5000)
ts.plot(beta[,1]);acf(beta[,1],lag=500)
ts.plot(beta[,5]);acf(beta[,5])
ts.plot(tau.2[,1]);acf(tau.2[,1])
ts.plot(tau.2[,5]);acf(tau.2[,5])
## after burn in and slicing plots
ts.plot(beta[ind,1]);acf(beta[ind,1],lag=100)

## Boxplot of distribution of beta values ##
set.seed(0204)
A.X8 <- replicate(40,{
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
  #X[1:5,]
  #summary(x8)
  #summary(x15)
  apply(X,2,mean)
  beta<-c(4,-4,5,7,-6,rep(0,10));beta
  #length(beta)
  Xb<-X%*%beta                               
  #Xb
  
  # Obtain the vector with probabilities of success p using the probit link
  p<-pnorm(Xb)

  y.data<-rbinom(1024,1,p)
  
  fit<-glm(y.data~0+X,family=binomial(link=logit))
  
  z<-rep(NA,n)
  N.sim<-270000
  N1<-sum(y.data);N1
  N0<-n-N1;N0
  ## beta
  beta<-matrix(NA,nrow=N.sim+1,ncol=d)
  (beta[1,]<-c(fit$coefficients))
  #beta
  ## tau
  tau.2<-matrix(NA,nrow=N.sim+1,ncol=d)
  tau.2[1,]<-c(rep(1,d))
  lamda<-rep(NA,n)
  lamda[1]<-2
  ## sampling
  for (j in 1:N.sim){
    mu.z<-X%*%beta[j,]
    z[y.data==0]<-rtruncnorm(N0,mean=mu.z[y.data==0],sd=1,a=-Inf,b=0)
    z[y.data==1]<-rtruncnorm(N1,mean=mu.z[y.data==1],sd=1,a=0,b=Inf)
    tau<-diag(tau.2[j,])
    E<-solve(solve(tau)+t(X)%*%X)
    beta[j+1,]<-rmvnorm(1,E%*%t(X)%*%z,E)
    tau.2[j+1,]<-1/rinv.gaussian(d,(lamda[j]^2/beta[j+1,]^2)^0.5,lamda[j]^2)
    lamda[j+1]<-(rgamma(1,d+1,0.5*sum(tau.2[j+1,])))^0.5
  }
 
ind<-seq(from=20000,to=270000,by=125)
bb <- colMeans(beta[ind,])##estimate beta
})


b1.X8 <- c(A.X8[1,26],A.X8[1,2],A.X8[1,30],A.X8[1,4],A.X8[1,5],A.X8[1,6],A.X8[1,7],A.X8[1,31],A.X8[1,9],A.X8[1,10],A.X8[1,11],A.X8[1,37],A.X8[1,38],A.X8[1,15],
           A.X8[1,16],A.X8[1,17],A.X8[1,18],A.X8[1,33],A.X8[1,20],A.X8[1,21],A.X8[1,22],A.X8[1,23],A.X8[1,24],A.X8[1,25],A.X8[1,27],A.X8[1,28],A.X8[1,29],
           A.X8[1,34],A.X8[1,35],A.X8[1,36])


b2.X8 <- c(A.X8[2,26],A.X8[2,2],A.X8[2,30],A.X8[2,4],A.X8[2,5],A.X8[2,6],A.X8[2,7],A.X8[2,31],A.X8[2,9],A.X8[2,10],A.X8[2,11],A.X8[2,37],A.X8[2,38],A.X8[2,15],
           A.X8[2,16],A.X8[2,17],A.X8[2,18],A.X8[2,33],A.X8[2,20],A.X8[2,21],A.X8[2,22],A.X8[2,23],A.X8[2,24],A.X8[2,25],A.X8[2,27],A.X8[2,28],A.X8[2,29],
           A.X8[2,34],A.X8[2,35],A.X8[2,36])

b3.X8 <- c(A.X8[3,26],A.X8[3,2],A.X8[3,30],A.X8[3,4],A.X8[3,5],A.X8[3,6],A.X8[3,7],A.X8[3,31],A.X8[3,9],A.X8[3,10],A.X8[3,11],A.X8[3,37],A.X8[3,38],A.X8[3,15],
           A.X8[3,16],A.X8[3,17],A.X8[3,18],A.X8[3,33],A.X8[3,20],A.X8[3,21],A.X8[3,22],A.X8[3,23],A.X8[3,24],A.X8[3,25],A.X8[3,27],A.X8[3,28],A.X8[3,29],
           A.X8[3,34],A.X8[3,35],A.X8[3,36])

b4.X8 <- c(A.X8[4,26],A.X8[4,2],A.X8[4,30],A.X8[4,4],A.X8[4,5],A.X8[4,6],A.X8[4,7],A.X8[4,31],A.X8[4,9],A.X8[4,10],A.X8[4,11],A.X8[4,37],A.X8[4,38],A.X8[4,15],
           A.X8[4,16],A.X8[4,17],A.X8[4,18],A.X8[4,33],A.X8[4,20],A.X8[4,21],A.X8[4,22],A.X8[4,23],A.X8[4,24],A.X8[4,25],A.X8[4,27],A.X8[4,28],A.X8[4,29],
           A.X8[4,34],A.X8[4,35],A.X8[4,36])

b5.X8 <- c(A.X8[5,26],A.X8[5,2],A.X8[5,30],A.X8[5,4],A.X8[5,5],A.X8[5,6],A.X8[5,7],A.X8[5,31],A.X8[5,9],A.X8[5,10],A.X8[5,11],A.X8[5,37],A.X8[5,38],A.X8[5,15],
           A.X8[5,16],A.X8[5,17],A.X8[5,18],A.X8[5,33],A.X8[5,20],A.X8[5,21],A.X8[5,22],A.X8[5,23],A.X8[5,24],A.X8[5,25],A.X8[5,27],A.X8[5,28],A.X8[5,29],
           A.X8[5,34],A.X8[5,35],A.X8[5,36])

b6.X8 <- c(A.X8[6,26],A.X8[6,2],A.X8[6,30],A.X8[6,4],A.X8[6,5],A.X8[6,6],A.X8[6,7],A.X8[6,31],A.X8[6,9],A.X8[6,10],A.X8[6,11],A.X8[6,37],A.X8[6,38],A.X8[6,15],
           A.X8[6,16],A.X8[6,17],A.X8[6,18],A.X8[6,33],A.X8[6,20],A.X8[6,21],A.X8[6,22],A.X8[6,23],A.X8[6,24],A.X8[6,25],A.X8[6,27],A.X8[6,28],A.X8[6,29],
           A.X8[6,34],A.X8[6,35],A.X8[6,36])

b7.X8 <- c(A.X8[7,26],A.X8[7,2],A.X8[7,30],A.X8[7,4],A.X8[7,5],A.X8[7,6],A.X8[7,7],A.X8[7,31],A.X8[7,9],A.X8[7,10],A.X8[7,11],A.X8[7,37],A.X8[7,38],A.X8[7,15],
           A.X8[7,16],A.X8[7,17],A.X8[7,18],A.X8[7,33],A.X8[7,20],A.X8[7,21],A.X8[7,22],A.X8[7,23],A.X8[7,24],A.X8[7,25],A.X8[7,27],A.X8[7,28],A.X8[7,29],
           A.X8[7,34],A.X8[7,35],A.X8[7,36])

b8.X8 <- c(A.X8[8,26],A.X8[8,2],A.X8[8,30],A.X8[8,4],A.X8[8,5],A.X8[8,6],A.X8[8,7],A.X8[8,31],A.X8[8,9],A.X8[8,10],A.X8[8,11],A.X8[8,37],A.X8[8,38],A.X8[8,15],
           A.X8[8,16],A.X8[8,17],A.X8[8,18],A.X8[8,33],A.X8[8,20],A.X8[8,21],A.X8[8,22],A.X8[8,23],A.X8[8,24],A.X8[8,25],A.X8[8,27],A.X8[8,28],A.X8[8,29],
           A.X8[8,34],A.X8[8,35],A.X8[8,36])

b9.X8 <- c(A.X8[9,26],A.X8[9,2],A.X8[9,30],A.X8[9,4],A.X8[9,5],A.X8[9,6],A.X8[9,7],A.X8[9,31],A.X8[9,9],A.X8[9,10],A.X8[9,11],A.X8[9,37],A.X8[9,38],A.X8[9,15],
           A.X8[9,16],A.X8[9,17],A.X8[9,18],A.X8[9,33],A.X8[9,20],A.X8[9,21],A.X8[9,22],A.X8[9,23],A.X8[9,24],A.X8[9,25],A.X8[9,27],A.X8[9,28],A.X8[9,29],
           A.X8[9,34],A.X8[9,35],A.X8[9,36])

b10.X8 <- c(A.X8[10,26],A.X8[10,2],A.X8[10,30],A.X8[10,4],A.X8[10,5],A.X8[10,6],A.X8[10,7],A.X8[10,31],A.X8[10,9],A.X8[10,10],A.X8[10,11],A.X8[10,37],A.X8[10,38],A.X8[10,15],
            A.X8[10,16],A.X8[10,17],A.X8[10,18],A.X8[10,33],A.X8[10,20],A.X8[10,21],A.X8[10,22],A.X8[10,23],A.X8[10,24],A.X8[10,25],A.X8[10,27],A.X8[10,28],A.X8[10,29],
            A.X8[10,34],A.X8[10,35],A.X8[10,36])

b11.X8 <- c(A.X8[11,26],A.X8[11,2],A.X8[11,30],A.X8[11,4],A.X8[11,5],A.X8[11,6],A.X8[11,7],A.X8[11,31],A.X8[11,9],A.X8[11,10],A.X8[11,11],A.X8[11,37],A.X8[11,38],A.X8[11,15],
            A.X8[11,16],A.X8[11,17],A.X8[11,18],A.X8[11,33],A.X8[11,20],A.X8[11,21],A.X8[11,22],A.X8[11,23],A.X8[11,24],A.X8[11,25],A.X8[11,27],A.X8[11,28],A.X8[11,29],
            A.X8[11,34],A.X8[11,35],A.X8[11,36])

b12.X8 <- c(A.X8[12,26],A.X8[12,2],A.X8[12,30],A.X8[12,4],A.X8[12,5],A.X8[12,6],A.X8[12,7],A.X8[12,31],A.X8[12,9],A.X8[12,10],A.X8[12,11],A.X8[12,37],A.X8[12,38],A.X8[12,15],
            A.X8[12,16],A.X8[12,17],A.X8[12,18],A.X8[12,33],A.X8[12,20],A.X8[12,21],A.X8[12,22],A.X8[12,23],A.X8[12,24],A.X8[12,25],A.X8[12,27],A.X8[12,28],A.X8[12,29],
            A.X8[12,34],A.X8[12,35],A.X8[12,36])

b13.X8 <- c(A.X8[13,26],A.X8[13,2],A.X8[13,30],A.X8[13,4],A.X8[13,5],A.X8[13,6],A.X8[13,7],A.X8[13,31],A.X8[13,9],A.X8[13,10],A.X8[13,11],A.X8[13,37],A.X8[13,38],A.X8[13,15],
            A.X8[13,16],A.X8[13,17],A.X8[13,18],A.X8[13,33],A.X8[13,20],A.X8[13,21],A.X8[13,22],A.X8[13,23],A.X8[13,24],A.X8[13,25],A.X8[13,27],A.X8[13,28],A.X8[13,29],
            A.X8[13,34],A.X8[13,35],A.X8[13,36])

b14.X8 <- c(A.X8[14,26],A.X8[14,2],A.X8[14,30],A.X8[14,4],A.X8[14,5],A.X8[14,6],A.X8[14,7],A.X8[14,31],A.X8[14,9],A.X8[14,10],A.X8[14,11],A.X8[14,37],A.X8[14,38],A.X8[14,15],
            A.X8[14,16],A.X8[14,17],A.X8[14,18],A.X8[14,33],A.X8[14,20],A.X8[14,21],A.X8[14,22],A.X8[14,23],A.X8[14,24],A.X8[14,25],A.X8[14,27],A.X8[14,28],A.X8[14,29],
            A.X8[14,34],A.X8[14,35],A.X8[14,36])

b15.X8 <- c(A.X8[15,26],A.X8[15,2],A.X8[15,30],A.X8[15,4],A.X8[15,5],A.X8[15,6],A.X8[15,7],A.X8[15,31],A.X8[15,9],A.X8[15,10],A.X8[15,11],A.X8[15,37],A.X8[15,38],A.X8[15,15],
            A.X8[15,16],A.X8[15,17],A.X8[15,18],A.X8[15,33],A.X8[15,20],A.X8[15,21],A.X8[15,22],A.X8[15,23],A.X8[15,24],A.X8[15,25],A.X8[15,27],A.X8[15,28],A.X8[15,29],
            A.X8[15,34],A.X8[15,35],A.X8[15,36])

length(b1.X8)
max(b1.X8)
beta.df.X8 <- data.frame(b1.X8,b2.X8,b3.X8,b4.X8,b5.X8,b6.X8,b7.X8,b8.X8,b9.X8,b10.X8,b11.X8,b12.X8,b13.X8,b14.X8,b15.X8)
boxplot(beta.df.X8)
library(ggplot2)
library(reshape)
beta_long.X8 <- melt(beta.df.X8)

ggplot(beta_long.X8, aes(x = variable, y = value)) +            # Applying ggplot function
  geom_boxplot(fill="steelblue") + 
  labs(title = "Distribution of Beta values", x = "Beta")




## Application using real world data: Cardiovascular Disease Data ##
attach(cardio_train)
length(cardio_train)
nrow(cardio_train)
cardio <- cardio_train
nrow(cardio)
library(mvtnorm)
library(MASS)
library(truncnorm) 
library(statmod) 
library(VGAM)
library(LearnBayes)
n = 70000
str(cardio_train)
d = 11
X1 <- scale(cardio$age)
X2 <- cardio$gender
X3 <- scale(cardio$height)
X4 <- scale(cardio$weight)
X5 <- scale(cardio$ap_hi)
X6 <- scale(cardio$ap_lo)
X7 <- cardio$cholesterol
X8 <- cardio$gluc
X9 <- cardio$smoke
X10 <- cardio$alco
X11 <- cardio$active
head(X1)
X.data.scaled <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11) #Combine predictor variables by columns
colnames(X.data.scaled) <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11")
X.data.scaled[1:5,] #Display first 5 observations
summary(X3) #Summary of height variable
summary(X4) #summary of weight variable
summary(X7) #summary of cholesterol level
apply(X.data,2,mean) #Mean of predictor variables
y.data <- cardio$cardio
length(y.data)
cardio <- data.frame(X.data.scaled,y.data)
head(cardio)
split.size <- 0.70 
n1 = n*split.size ; n1
n2 = n*(1- split.size) ;n2
sample.size <- floor(split.size*nrow(cardio)) #training data size = 49000
set.seed(1234)
train_indices <- sample(seq_len(nrow(cardio)), size = sample.size)
cardio.train <- cardio[train_indices, ]    # 70% of data
cardio.test <- cardio[-train_indices, ]    # 30% of data.  
#X matrix of training data
X.train.data <- cbind(cardio.train$X1,cardio.train$X2,cardio.train$X3,cardio.train$X4,cardio.train$X5,
                      cardio.train$X6,cardio.train$X7,cardio.train$X8,cardio.train$X9,cardio.train$X10,
                      cardio.train$X11)
head(X.train.data)
Y.train.data <- cardio.train$y.data
length(Y.train.data)
model.01 <- glm(Y.train.data~0+X.train.data, family=binomial(link=probit), data = cardio.train)
summary(model.01) #summary of the model
model.01$coefficients
z<-rep(NA,n1)
N.sim<-270000
N1<-sum(Y.train.data);N1
N0<-n1-N1;N0
beta.train<-matrix(NA,nrow=N.sim+1,ncol=d) #beta matrix 
beta.train[1,]<-c(model.01$coefficients)
beta.train
## tau
tau.2<-matrix(NA,nrow=N.sim+1,ncol=d)
tau.2[1,]<-c(rep(1,d))
tau.2
lamda<-rep(NA,n1)
lamda[1]<-2
lamda
## sampling
for (j in 1:N.sim){
  mu.z<-X.train.data%*%beta.train[j,]
  z[Y.train.data==0]<-rtruncnorm(N0,mean=mu.z[Y.train.data==0],sd=1,a=-Inf,b=0)
  z[Y.train.data==1]<-rtruncnorm(N1,mean=mu.z[Y.train.data==1],sd=1,a=0,b=Inf)
  tau<-diag(tau.2[j,])
  E<-solve(solve(tau)+t(X.train.data)%*%X.train.data)
  beta.train[j+1,]<-rmvnorm(1,E%*%t(X.train.data)%*%z,E)
  tau.2[j+1,]<-1/rinv.gaussian(d,(lamda[j]^2/beta.train[j+1,]^2)^0.5,lamda[j]^2)
  lamda[j+1]<-(rgamma(1,d+1,0.5*sum(tau.2[j+1,])))^0.5
}

summary(lamda)
plot(lamda)
z
beta.train[,1]
tau.2[,1]
ind<-seq(from=20000,to=270000,by=125)
summary(lamda[ind])##estimate lambda
sd(lamda[ind])
ts.plot(lamda[ind]);acf(lamda[ind])
beta.means <-colMeans(beta.train[ind,])##estimate beta.train
tau.means <- colMeans(tau.2[ind,])##estimate tau
newb<-beta.train[ind,]
sd.beta.train<-rep(NA,11)
for (k in 1:11){
  sd.beta.train[k]<-sd(newb[,k])
}
sd.beta.train##sd of beta.train estimates
quantile(beta.train[ind,1],c(0.025,0.5,0.975))
quantile(beta.train[ind,2],c(0.025,0.5,0.975))
quantile(beta.train[ind,3],c(0.025,0.5,0.975))
quantile(beta.train[ind,4],c(0.025,0.5,0.975))
quantile(beta.train[ind,5],c(0.025,0.5,0.975))
quantile(beta.train[ind,6],c(0.025,0.5,0.975))
quantile(beta.train[ind,7],c(0.025,0.5,0.975))
quantile(beta.train[ind,8],c(0.025,0.5,0.975))
quantile(beta.train[ind,9],c(0.025,0.5,0.975))
quantile(beta.train[ind,10],c(0.025,0.5,0.975))
quantile(beta.train[ind,11],c(0.025,0.5,0.975))

par(mfrow=c(2,2))
ts.plot(lamda);acf(lamda,lag=50)
ts.plot(beta.train[,1]);acf(beta.train[,1],lag=50)
ts.plot(beta.train[,5]);acf(beta.train[,5])
ts.plot(tau.2[,1]);acf(tau.2[,1])
ts.plot(tau.2[,5]);acf(tau.2[,5])
## after burn in and slicing plots
ts.plot(beta.train[ind,1]);acf(beta.train[ind,1],lag=100)
cardio.test
X1.test <- cardio.test$X1
X2.test <- cardio.test$X2
X3.test <- cardio.test$X3
X4.test <- cardio.test$X4
X5.test <- cardio.test$X5
X6.test <- cardio.test$X6
X7.test <- cardio.test$X7
X8.test <- cardio.test$X8
X9.test <- cardio.test$X9
X10.test <- cardio.test$X10
X11.test <- cardio.test$X11
y.test.data <- cardio.test$y.data
#X matrix of testing data
X.test.data <- cbind(X1.test,X2.test,X3.test,X4.test,X5.test,X6.test,X7.test,X8.test,X9.test,X10.test,X11.test)
head(X.test.data)
#y data of testing data
y.test.data <- cardio.test$y.data
beta.means <- colMeans(beta.train[ind,])
xb.train <- X.test.data%*%beta.means
y.prob <- pnorm(xb.train)
y.pred <- ifelse(y.prob > 0.5, 1,0)
head(cardio.test$y.data)
head(y.pred)
accuracy_01 <- (sum(y.pred==y.test.data)) / length(y.test.data) ;accuracy_01
mean(y.pred == y.test.data)


