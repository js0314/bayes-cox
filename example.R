rm(list=ls())
library(MASS); library(caret); library(dplyr); 
library(glmnet); library(ncvreg) 
library(survival); library(survcomp); library(Rcpp); library(RcppArmadillo); library(mvtnorm) 

source("GdPrior.R")
sourceCpp("LogParLik.cpp")


#=========================================================#
#   DATA for sc 1-3        
#=========================================================#
n=100; p = 20; sc = 2; tau0 <- 0.0001; iter <-100; KAP.j = seq(from=0.1, to=3, length.out = 20);

if(sc==1){
  Beta <- c(rep(1, 2), rep(0, 2), rep(1, 2))
  # Beta <- rep(1, 5)
  Cor <- diag(p)
}else  if(sc==2){
  Beta <- c(rep(1, 5), rep(0, 5), rep(1, 5)) 
  rho <- 0.5
  Cor <- rho^(abs(matrix(1:p,p,p)-t(matrix(1:p,p,p))))
}else if(sc==3){
  Beta <- c(4, 4, 4, -6*sqrt(2), 4/3)/4 
  Cor <- matrix(1/2, p, p);  Cor[,4] <- Cor[4,] <- 1/sqrt(2.2); Cor[,5] <- Cor[5,] <- 0; diag(Cor) <- 1
}

set.seed(1234)
cholmat = chol(Cor)
X = matrix(rnorm(n*p, mean=0, sd=1), n, p)
X_train = X%*%cholmat
beta <- numeric(p); beta[c(1:length(Beta))] <- Beta
myrates = as.vector(exp(X_train %*% beta))
Sur = rexp(n, myrates); CT = rexp(n, 0.1)
time_train = pmin(Sur,CT); event_train = as.numeric(Sur<=CT)

X_train <- X_train[order(time_train),]
event_train <- event_train[order(time_train)]
time_train <- time_train[order(time_train)]
y_train <- cbind(time=time_train, status=event_train)
y_train2 <- Surv(time=time_train, event=event_train)
tXX.j <- apply(as.matrix(X_train)^2, 2, sum)


#=========================================================#
#   GD prior 
#=========================================================#
fit <- GdPrior(X_train, time_train, event_train, KAP.j=KAP.j, xi_j=0.01, mc.size=5000)

# est beta.map
fit$beta.map

# est dmap
fit$d.map

# selected kappa
fit$opt.k

# selected bic
fit$opt.bic
