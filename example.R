rm(list=ls())
library(MASS); library(modelr); library(numDeriv); library(caret); library(dplyr); library(tidyverse); library(tidyr); 
library(ggplot2); library(reshape2); library(cowplot); library(gridExtra); library(flextable); library(forcats) # figure
library(glmnet); library(ncvreg) # var sel
library(survival); library(survcomp); library(survAUC) # survival
library(HDInterval); library(LearnBayes); library(BhGLM) # bayesian
library(Rcpp); library(RcppArmadillo); library(mvtnorm) # cpp


source("GdPrior.R")
sourceCpp("LogParLik.cpp")


#=========================================================#
#   DATA         
#=========================================================#
n=100; p = 20; sc = 1; tau0 <- 0.0001; iter <-100

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


set.seed(4342)
cholmat = chol(Cor)
X = matrix(rnorm(n*p, mean=0, sd=1), n, p)
X_train = X%*%cholmat
beta <- numeric(p); beta[c(1:length(Beta))] <- Beta
myrates = as.vector(exp(X_train %*% beta))
Sur = rexp(n, myrates); CT = rexp(n, 0.1)
time_train = pmin(Sur,CT); event_train = as.numeric(Sur<=CT)
time_train = pmin(Sur,CT); event_train = as.numeric(Sur<=CT)
y_train <- cbind(time=time_train, status=event_train)
y_train2 <- Surv(time=time_train, event=event_train)
tXX.j <- apply(as.matrix(X_train)^2, 2, sum)


#=========================================================#
#   GD prior 
#=========================================================#
f4 = GdPrior(X_train, time_train, event_train, KAP.j = seq(0.1, 3, length.out=30), xi_j=0.01, mc.size=3000)
# est coef
f4$beta.map
# est d map
f4$d.map
# selected kap
f4$opt.k


