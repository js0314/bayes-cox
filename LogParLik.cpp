#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double sign(double a) {
  double mysign;
  if(a == 0) {mysign=0;} else {mysign = (a>0 ? 1 : -1);}
  return mysign;
}

double absf(double a) {
  double value = (a>0 ? a : ((-1)*a));
  return value;
}

// [[Rcpp::export]]
double loglik_num(arma::mat& X, arma::colvec& time, arma::colvec& delta, arma::colvec& beta) {
  
  int nrow = X.n_rows;
  
  // initialization
  double loglik = 0;
  arma::mat at_risk(nrow,1);
  arma::mat xbeta = X*beta;
  arma::mat exp_xbeta = exp(xbeta);
  double time_i;
  
  // calculation
  loglik += arma::as_scalar(xbeta.t()*delta);
  for(int i=0; i<nrow; i++) {
    time_i = time(i);
    if(delta(i) != 0) {
      for(int k=0; k<nrow; k++) {
        at_risk(k,0) = (time(k) >= time_i ? 1 : 0);
      }
      loglik += (-1.0)*log(arma::as_scalar(exp_xbeta.t()*at_risk));
    }
  }
  
  return loglik;
}


//[[Rcpp::export]]
void neg_ddloglik_num(arma::mat& neg_ddloglik,arma::mat& X, arma::colvec& time, arma::colvec& delta, arma::colvec& beta) {
  
  int ncol = X.n_cols;
  int nrow = X.n_rows;
  
  // initialization
  arma::mat xbeta = X*beta;
  arma::mat exp_xbeta = exp(xbeta);
  double mu0_i;
  arma::colvec mu1_i(ncol);
  arma::mat mu2_i(ncol, ncol);
  arma::mat at_risk(nrow, 1);
  double time_i;
  neg_ddloglik.fill(0);
  
  // calculation
  for(int i=0; i<nrow; i++) {
    time_i = time(i);
    if(delta(i) != 0) {
      for(int k=0; k<nrow; k++) {
        at_risk(k,0) = (time(k) >= time_i ? 1 : 0);
      }
      mu0_i = arma::as_scalar(exp_xbeta.t()*at_risk);
      mu1_i = trans(X)*(at_risk%exp_xbeta);
      mu2_i = trans(X)*diagmat(vectorise(at_risk%exp_xbeta))*X;
      neg_ddloglik += mu2_i/mu0_i - (mu1_i/mu0_i)*trans(mu1_i/mu0_i);
    }
  }
  neg_ddloglik = (+1.0)*neg_ddloglik;
}

