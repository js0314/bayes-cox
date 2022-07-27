# Precision correction

# regular CI
"creint" = function(x, alpha, digit=3){
  tmp.cre <- paste("(",round(quantile(x, alpha/2), digit),",", round(quantile(x, 1- alpha/2), digit),")", sep = "")
  tmp.cre
}


# HPD
"hdint" = function(x, alpha, digit=3){
  hpd <- round(as.numeric(hdi(x, credMass=1-alpha)), digit)
  tmp.cre <- paste("(", hpd[1],",", hpd[2],")", sep = "")  # print(tmp.cre, quote = FALSE)
  tmp.cre
}

"hdint.extend" = function(x, alpha, digit=3){
  hpd <- round(as.numeric(hdi(x, credMass=1-alpha)), digit)
  hpd
}

"mc.cov.check" = function(HDINT, beta){
  p <- length(beta)
  checker <- as.numeric(sapply(1:p, function(i) between(beta[i], HDINT[i, 1], HDINT[i, 2])))
  return(checker)
}

# regular CI considering nonzero cases
"creint2" = function(x, alpha){
  x.new <- x[x!=0]
  tmp.cre <- paste("(",round(quantile(x.new, alpha/2), 4),",", round(quantile(x.new, 1- alpha/2), 4),")", sep = "")  # print(tmp.cre, quote = FALSE)
  tmp.cre
}


# HPD considering nonzero cases
"hdint2" = function(x, alpha){
  x.new <- x[x!=0]
  hpd <- round(as.numeric(hdi(x.new, credMass=1-alpha)), 4)
  tmp.cre <- paste("(", hpd[1],",", hpd[2],")", sep = "")  # print(tmp.cre, quote = FALSE)
  tmp.cre
}


## Credible interval difference
"crediff" = function(x){
  crediff <- round(quantile(x, c(0.975)), 4) - round(quantile(x, c(0.025)), 4)
  crediff
}


# MSE Calculation
"mse.cal" <- function(beta.hat, true.beta){
  mean((beta.hat-true.beta)^2)
}



#fpr
"fpr.cal" <- function(beta.hat, true.beta){
  TN.which <- which(true.beta==0)
  sum(beta.hat[TN.which]!=0)/length(TN.which)*100
}


#fnr
"fnr.cal" <- function(beta.hat, true.beta){
  TP.which <- which(true.beta!=0)
  sum(beta.hat[TP.which]==0)/length(TP.which)*100
}


#fdr
"fdr.cal" <- function(beta.hat, true.beta){
  TN.which <- which(true.beta==0)
  pred.pos <- which(beta.hat!=0)
  sum(beta.hat[TN.which]!=0)/length(pred.pos)*100
}