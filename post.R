"GdPrior_post" <- function(X, time, event, init, beta.map, d.map, mc.ite = 3000, variance.cont = 1){
  
  X <- as.matrix(X)
  beta.map <- as.numeric(beta.map)
  d.map <- as.numeric(d.map)
  init <- as.numeric(init)
  
  n <- nrow(X); p <- ncol(X)
  BETA.mc <- matrix(0, mc.ite, p); 
  BETA.mc[1,] <- beta.map

 
  neg_ddloglik <- matrix(0, nrow=p, ncol=p)
  neg_ddloglik_cpp_ext2(neg_ddloglik,  X = X , time = time , delta = event, beta = beta.map) 
  neg_ddloglik_diag <- diag(neg_ddloglik)
  
  mh.accept.rate <- mh.u <- matrix(0, mc.ite, p)
  for(i in 2:mc.ite){
    for(loc in 1:p){
      beta.old <- beta.new <- BETA.mc[i-1, ]
      beta_j_old <- beta.old[loc]
      beta_j_new <- rnorm(1, mean = beta.old[loc], sd = sqrt(variance.cont/(neg_ddloglik_diag[loc] + d.map[loc])) )
      beta.new[loc] <- beta_j_new
      log.ratio <- loglik_cpp_ext2(X = X, time = time, delta = event, beta = beta.new) + sum(dnorm(beta.new, 0, 1/sqrt(d.map[loc]), log=TRUE)) - loglik_cpp_ext2( X = X, time = time, delta = event, beta = beta.old)  - sum(dnorm(beta.old, 0, 1/sqrt(d.map[loc]), log=TRUE))
      mh.accept.rate[i, loc] <- min(1, exp(log.ratio))
      mh.u[i, loc] <- runif(1)
      if(mh.u[i, loc] < mh.accept.rate[i, loc]){
        BETA.mc[i, loc] <- beta_j_new
      }else{
        BETA.mc[i, loc] <- beta_j_old
      }
    }
  }
  
  accept.rate <- apply(mh.accept.rate > mh.u, 2, mean)
  
  return(list(BETA.mc = BETA.mc,
              accept.rate = accept.rate))
  
}
