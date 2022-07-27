rm(list=ls())
library(MASS); library(modelr); library(numDeriv); library(caret); library(dplyr); library(tidyverse); library(tidyr); 
library(ggplot2); library(reshape2); library(cowplot); library(gridExtra); library(flextable); library(forcats) # figure
library(glmnet); library(ncvreg) # var sel
library(survival); library(survcomp); library(survAUC) # survival
library(HDInterval); library(LearnBayes); library(BhGLM) # bayesian
library(Rcpp); library(RcppArmadillo); library(mvtnorm) # cpp


source("GdPrior.R")
source("GdPrior_post.R")
source("MyCriteriaList.R")
sourceCpp("LogParLik.cpp")


#=========================================================#
#   DATA SETUP         
#=========================================================#
n=100; p = 20; sc = 1; KAP.j = seq(from=0.1, to=3, length.out = 20);
xi_j = seq(from = 0.1, to = 3, length = 20); 
B = 100; split = 0.5; mc.size=3000; iter=100; tau0 <- 0.0001
lambda <- 0.8; cens_rate <- 0.27
image.title <- paste0(Sys.Date(),"_final_sc", sc, "_n", n, "_p", p, "_kap_.1_3",".RData")


if(sc==1){
  Beta <- c(rep(1, 2), rep(0, 2), rep(1, 2))
  Cor <- diag(p)
}else  if(sc==2){
  Beta <- c(rep(1, 5), rep(0, 5), rep(1, 5)) 
  rho <- 0.5
  Cor <- rho^(abs(matrix(1:p,p,p)-t(matrix(1:p,p,p))))
}else if(sc==3){
  Beta <- c(1, 1, 1, -1.5*sqrt(2), 1/3)
  Cor <- matrix(1/2, p, p);  Cor[,4] <- Cor[4,] <- 1/sqrt(2.2); Cor[,5] <- Cor[5,] <- 0; diag(Cor) <- 1
}


#=========================================================#
#   OUPTUT LIST
#=========================================================#
SEL.BETA.oracle <- SEL.BETA.GD <- SEL.D.GD <- matrix(0, B, p)
MSE.TABLE <- FPR.TABLE <- FNR.TABLE <- FDR.TABLE <- CINDEX.TABLE <-  matrix(0, B, 2)
SEL.ID <- SEL.KAP <- rep(-Inf, B); MC.COV <- GD.ACCPT.RATE <- matrix(-Inf, B, p); GD.POST <- list()


#=========================================================#
#   Algorithm
#=========================================================#


for(song in 1:B){

  print(paste("sim=",song))
  set.seed(song + 1234)
  cholmat = chol(Cor)
  X = matrix(rnorm(n*p, mean=0, sd=1), n, p)
  X_train = X%*%cholmat
  beta <- numeric(p); beta[c(1:length(Beta))] <- Beta
  myrates = as.vector(exp(X_train %*% beta))
  Sur = rexp(n,myrates); CT = rexp(n, 0.1)
  time_train = pmin(Sur,CT); event_train = as.numeric(Sur<=CT)
  y_train <- cbind(time=time_train, status=event_train)
  y_train2 <- Surv(time=time_train, event=event_train)
  tXX.j <- apply(as.matrix(X_train)^2, 2, sum)
  
  
  # test
  set.seed(song + 3333)
  X_test = matrix(rnorm(n*p, mean=0, sd=1), n, p)
  X_test = X_test%*%cholmat
  myrates = as.vector(exp(X_test %*% beta))
  Sur_Test = rexp(n,myrates); CT_Test = rexp(n, 0.1)
  time_test = pmin(Sur_Test,CT_Test); event_test = as.numeric(Sur<=CT)
  y_test <- cbind(time=time_test, status=event_test)
  y_test2 <- Surv(time=time_test, event=event_test)
  
  #----------------------------------------------------------------------------------------------------------------
  
  beta.oracle <- rep(0, p)
  SEL.BETA.oracle[song,][which(beta!=0)] <- as.numeric(coxph(Surv(time_train, event_train) ~ X_train[, which(beta!=0)])$coefficient) #a-lasso estimate
  
  f4 = GdPrior(X_train, time_train, event_train, KAP.j = seq(0.1, 3, length.out=20), xi_j=0.01, mc.size=3000)
  SEL.BETA.GD[song,] <- f4$beta.map
  SEL.D.GD[song,] <- f4$d.map
  SEL.KAP[song] <- f4$opt.k
  
  # # # # posterior generator
  # post.gd <- GdPrior_post(X = X_train, time = time_train, event = event_train, init = f4$beta.map, beta.map = f4$beta.map, d.map = f4$d.map, variance.cont = 1.2)
  # GD.POST[[song]] <- post.gd$BETA.mc
  
  #----------------------------------------------------------------------------------------------------------------

  # evaluation
  MSE.TABLE[song,] <- c(mse.cal(beta.hat = SEL.BETA.oracle[song,], true.beta = beta),
                        mse.cal(beta.hat = SEL.BETA.GD[song,], true.beta = beta))
  
  FPR.TABLE[song,] <- c(fpr.cal(beta.hat = SEL.BETA.oracle[song,], true.beta = beta),
                        fpr.cal(beta.hat = SEL.BETA.GD[song,], true.beta = beta)) 
  
  FNR.TABLE[song,] <- c(fnr.cal(beta.hat = SEL.BETA.oracle[song,], true.beta = beta),
                        fnr.cal(beta.hat = SEL.BETA.GD[song,], true.beta = beta))
  
  FDR.TABLE[song,] <- c(fdr.cal(beta.hat = SEL.BETA.oracle[song,], true.beta = beta),
                        fdr.cal(beta.hat = SEL.BETA.GD[song,], true.beta = beta))
  
  
  CINDEX.TABLE[song,] <- c(
    concordance.index(x=as.numeric(X_test%*%SEL.BETA.oracle[song,]), surv.time= time_test, surv.event=event_test, method="noether")$c.index,
    concordance.index(x=as.numeric(X_test%*%SEL.BETA.GD[song,]), surv.time= time_test, surv.event=event_test, method="noether")$c.index)

}


#=====================================================================
# Summary
#=====================================================================
SEL.BETA.oracle <- as.data.frame(SEL.BETA.oracle); SEL.BETA.GD <- as.data.frame(SEL.BETA.GD);
colnames(SEL.BETA.oracle) <- colnames(SEL.BETA.GD) <- c(paste0("Beta0", 1:9), paste0("Beta", 10:p)) 


# Summary 1 -------------------------------
NUM.SEL <- c( round(mean(apply(SEL.BETA.oracle!=0, 1, sum)), 0),
            round(mean(apply(SEL.BETA.GD!=0, 1, sum)), 0))

final.result <- cbind(MSE = apply(MSE.TABLE, 2, mean, na.rm = TRUE),
                       CINDEX = apply(CINDEX.TABLE, 2, mean, na.rm = TRUE),
                       FPR = apply(FPR.TABLE, 2, mean, na.rm = TRUE),
                       FNR = apply(FNR.TABLE, 2, mean, na.rm = TRUE),
                       FDR = apply(FDR.TABLE, 2, mean, na.rm = TRUE), 
                      NUM.SEL = NUM.SEL)
rownames(final.result) <- c("oracle",  "gd")

data.frame(final.result) %>% rownames_to_column() %>% flextable() %>% colformat_double(big.mark=",", digits = 3, na_str = "N/A")



# Summary 2 -------------------------------
SEL.BETA.GD.coll <- data.frame(SEL.BETA.GD[, which(beta!=0) ],
                               noise = apply(SEL.BETA.GD[, which(beta==0)], 1, mean))

beta.coll <- c(beta[which(beta!=0)],0)

final.result_v2 <- cbind(
  MEAN = apply(SEL.BETA.GD.coll, 2, mean),
  SE = apply(SEL.BETA.GD.coll, 2, sd),
  MSE = apply((SEL.BETA.GD.coll - matrix(beta.coll, by=T, nrow=B, ncol=length(beta.coll)))^2, 2, mean))

flextable(data.frame(final.result_v2)) %>% theme_box() %>% 
  colformat_double(big.mark=",", digits = 5, na_str = "N/A")

