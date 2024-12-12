attach(cardio_train)
length(cardio_train)
nrow(cardio_train)
cardio <- cardio_train
nrow(cardio)

library(mvtnorm)
library(MASS)
library(truncnorm) ## rtruncnorm(N0,mean=mu.z[y.data==0],sd=1,a=-Inf,b=0)
library(statmod) ## rinvgauss(10,1,1)
library(VGAM) ## rinv.gaussian(10,1,1)
library(LearnBayes)
n = 70000

d = 11
X1 <- cardio$age
X2 <- cardio$gender
X3 <- cardio$height
X4 <- cardio$weight
X5 <- cardio$ap_hi
X6 <- cardio$ap_lo
X7 <- cardio$cholesterol
X8 <- cardio$gluc
X9 <- cardio$smoke
X10 <- cardio$alco
X11 <- cardio$active

X.data <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11) #Combine predictor variables by columns

X.data[1:5,] #Display first 5 observations
summary(X3) #Summary of height variable
summary(X4) #summary of weight variable
summary(X7) #summary of cholesterol level

apply(X.data,2,mean) #Mean of predictor variables

y.data <- cardio$cardio
length(y.data)
cardio <- data.frame(X.data,y.data)
head(cardio)


split.size <- 0.70 
n1 = n*split.size ; n1
n2 = n*(1- split.size) ;n2
sample.size <- floor(split.size*nrow(cardio)) #training data size = 66500
set.seed(0211)
train_indices <- sample(seq_len(nrow(cardio)), size = sample.size)
cardio.train <- cardio[train_indices, ]    # These are 70% of data
cardio.test <- cardio[-train_indices, ]    # These are rest of them, 30%.  

#X matrix of training data
X.train.data <- cbind(cardio.train$X1,cardio.train$X2,cardio.train$X3,cardio.train$X4,cardio.train$X5,
            cardio.train$X6,cardio.train$X7,cardio.train$X8,cardio.train$X9,cardio.train$X10,
            cardio.train$X11)
head(X.train.data)
Y.train.data <- cardio.train$y.data
length(Y.train.data)
model.01 <- glm(Y.train.data~0+X.train.data, family=binomial(link=logit), data = cardio.train)
summary(model.01) #summary of the model
model.01$coefficients
## Z
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

## lamda
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

#y data of testing data
y.test.data <- cardio.test$y.data

beta.means <- colMeans(beta.train[ind,])

library(Rfast)
beta.medians <- colMedians(beta.train[ind,])

xb.train <- X.test.data%*%beta.means

y.prob <- pnorm(xb.train)
y.pred <- ifelse(y.prob > 0.5, 1,0)

head(cardio.test$y.data)
head(y.pred)
accuracy_01 <- (sum(y.pred==y.test.data)) / length(y.test.data) ;accuracy_01
mean(y.pred == y.test.data)


set.seed(0301)
random.pick <- rbinom(21000, 1,0.5)
head(random.pick)
mean(random.pick == y.test.data)



for (k in 1:11){
  sd.beta.train[k]<-var(newb[,k])
}

sort(sd.beta.train)
ab10 <- beta.means[10]
ab9 <- beta.means[9]
ab11 <- beta.means[11]
ab2 <- beta.means[2]
ab8 <- beta.means[8]
ab7 <- beta.means[7]
beta.red <- cbind(ab2,ab7,ab8,ab9,ab10,ab11)

ay2 <- cbind(X2.test,X7.test,X8.test,X9.test,X10.test,X11.test)

ay.train <- ay2%*%as.vector(beta.red)

ay.prob <- pnorm(ay.train)
ay.pred <- ifelse(ay.prob > 0.5, 1,0)

head(cardio.test$y.data)
head(ay.pred)
accuracy_01 <- (sum(y.pred==y.test.data)) / length(y.test.data) ;accuracy_01
mean(ay.pred == y.test.data)







