/* 

Calculate kernel according to the following formula:
K = exp(-dist(X_j, X_i) / sigma)

*/


#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat fastKgipMat(NumericMatrix AA, double sigma0) {
  // AA, adjacent matrix {0, 1}
  // sigma0, initial kernel width
	  
  NumericMatrix Ar = clone(AA);
  int m = Ar.nrow(); 
  int k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  double sigma = sigma0 * (accu(A) / m);
  arma::colvec An =  sum(square(A), 1);
  arma::mat C = -2 * (A * A.t());
  C.each_col() += An;
  C.each_row() += An.t();
  // kernel matrix
  arma::mat Kmat = exp(-(C / sigma));
  return Kmat;
}