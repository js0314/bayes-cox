GdPrior <- function(X_train, time_train, event_train, KAP.j, xi_j, mc.size=3000){
  
  y_train <- cbind(time=time_train, status=event_train)
  y_train2 <- Surv(time=time_train, event=event_train)
  
  hat.beta <- as.numeric(glmnet(x = X_train, y = y_train2, family="cox", alpha=0, lambda = KAP.j[1])$beta)
  beta.GD0 <- hat.beta
  
  
  # Calculate beta for each kappa
  LL <-  length(KAP.j)
  BETA.GD <- matrix(0,LL,p); D.GD <- matrix(0,LL,p)
  NEG.MAR.GD <- rep(NA, LL)
  POST.SAM.GD <- POST.SAM.GD2 <- list()
  
  # GD.rep = 1
  for(GD.rep in 1:LL){
    print(GD.rep)
    kappa <- KAP.j[GD.rep]
    hat.beta <- hat.beta2 <- as.numeric(glmnet(x = X_train, y = y_train, family="cox", alpha=0, lambda=kappa)$beta)
    out.threshold <- which(hat.beta!=0)
    hat.d <- rep(0, p)
    
    # ICM Algorithm
    for(i in 1:iter){
      hat.d <- 2*kappa / ((1/tXX.j)*tau0^2+(hat.beta)^2)
      for(j in out.threshold){
        d_j <- hat.d[j]
        Likehood <- function(b.j){ 
          hat.beta[j] <- b.j
          Xbeta <- as.numeric(X_train%*%hat.beta)
          (-1)*(lplik_fun(X = X_train, time_train, event_train, hat.beta))+(-1)*dnorm(b.j, 0,1/sqrt(d_j),log=TRUE)
        }
        ini.b.j <- hat.beta[j]
        hat.beta[j] <- hat.beta2[j] <- optim(ini.b.j, Likehood, method="BFGS")$par
        hat.beta[j] <- as.numeric(hat.beta[j])*as.numeric(abs(hat.beta[j])>xi_j)
      }
      out.threshold <- which(hat.beta!=0)
      if(max((hat.beta-beta.GD0)^2)<0.0001){break}
      beta.GD0 <- hat.beta
    }
    
    Post <- function(bb){
      dd_j <- length(out.threshold)
      (-1)*(lplik_fun(X = as.matrix(X_train[,out.threshold]), time_train, event_train, bb))+ (-1)*sum(dnorm(bb, 0, 1/sqrt(dd_j), log=TRUE))
    }
    
    if(sum(abs(hat.beta))==0){
      next
    } else{
      # MCMC for selected variable only
      Post.opt <- optim(hat.beta[out.threshold], Post, method="BFGS", hessian=TRUE)
      Post.mean <- Post.opt$par
      Post.var <- solve(Post.opt$hessian)
      Post.sam <-  mvrnorm(mc.size, Post.mean, Post.var) #mvrnorm
      log.g <- dmnorm(Post.sam, Post.mean, Post.var, log=TRUE)
      num_cal <- function(X, time, event, out.threshold, bb){
        lplik_fun(X = as.matrix(X[,out.threshold]), time, event, bb) + sum(dnorm(bb, 0, 1/sqrt(length(out.threshold)), log=TRUE))
        }
      log.like <- sapply(1:mc.size, function(i) num_cal(X = as.matrix(X_train), time = time_train, event = event_train, out.threshold = out.threshold, bb = Post.sam[i,]))
      log.like.fin.tmp <- (log.like-log.g) 
      max.log <- max(log.like.fin.tmp)
      NEG.MAR.GD[GD.rep] <- (-1)*( log(sum(exp(log.like.fin.tmp - max.log) )) + max.log - log(nrow(X_train)))
      BETA.GD[GD.rep,] <- hat.beta; 
      D.GD[GD.rep,] <- hat.d
      POST.SAM.GD[[GD.rep]] <- Post.sam
    }
  }
  
  sel.kap.id <- which.min(NEG.MAR.GD)
  
  return(list(
    beta.map = BETA.GD[sel.kap.id, ],
    d.map = D.GD[sel.kap.id, ],
    opt.k = KAP.j[sel.kap.id],
    opt.bic =  NEG.MAR.GD[sel.kap.id],
    NEG.MAR.GD = NEG.MAR.GD,
    BETA.GD = BETA.GD
  ))
}



