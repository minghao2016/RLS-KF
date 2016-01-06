/*

Closed-form solution: RLS

*/


#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace arma;

//[[Rcpp::export]]
vec fastSolve(mat A, vec y, uword numTr, uword lambda) {
  // closed-form solution
  return solve(A + lambda * eye<mat>(numTr, numTr), y);
}