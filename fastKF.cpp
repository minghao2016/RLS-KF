/* 

Fast version of Kernel Fusion (KF) algorithm based on the SNFtool package
by Wang et al.[1] using Rcpp and RcppArmadillo packages developed by 
Dirk Eddelbuettel and co-workers.

[1] B Wang, et. al. Nat. Methods 11 (2014) 333-337.

*/

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// sort: from large to small
struct greater {
  template<class T>
  bool operator()(T const &a, T const &b) const {return a > b;}
};


//[[Rcpp::export]]
mat fastKF(NumericMatrix sim1, NumericMatrix sim2, int nNeig, int nIter) {
	// sim1: first similarity matrix
	// sim2: second similarity matrix
	// nNeig: number of neighbours
	// nIter: number of iterations
	
	int NM = 2;
	int nr = sim1.nrow();
	int nc = nr;
	NumericMatrix S_1 = clone(sim1);
	NumericMatrix S_2 = clone(sim2);
  
	mat S1 = mat(S_1.begin(), nr, nc, false);
	mat S2 = mat(S_2.begin(), nr, nc, false);
	
	// Normalization
	S1.each_col() /= sum(S1, 1);
	S1 = (S1 + S1.t()) / 2;
	S2.each_col() /= sum(S2, 1);
	S2 = (S2 + S2.t()) / 2;
	
	// Backup
	mat S1_o = S1;
	mat S2_o = S2;
	
	// do K nearest neighbour operation
	rowvec rv_s1(nc);
	rowvec rv2_s1(nc);
	rowvec rv_s2(nc);
	rowvec rv2_s2(nc);
	
	// kth largest value
	double kth_s1;
	double kth_s2;
	// K nearest neighbours
	for (int ii = 0; ii < nr; ii++) {
		// current row
		rv_s1 = S1.row(ii);
		rv2_s1 = rv_s1;
		rv_s2 = S2.row(ii);
		rv2_s2 = rv_s2;
		// Kth largest value
		std::nth_element(rv_s1.begin(), rv_s1.begin() + nNeig - 1, rv_s1.end(), greater());
		std::nth_element(rv_s2.begin(), rv_s2.begin() + nNeig - 1, rv_s2.end(), greater());
		kth_s1 = rv_s1[nNeig - 1];
		kth_s2 = rv_s2[nNeig - 1];
		// Keep both 'equal' and 'greater' values
		rv2_s1.elem(find(rv2_s1 < kth_s1)).zeros();
		rv2_s2.elem(find(rv2_s2 < kth_s2)).zeros();
		// put row values back to original matrix
		S1.row(ii) = rv2_s1;
		S2.row(ii) = rv2_s2;
	}
	
	S1.each_col() /= sum(S1, 1);
	S2.each_col() /= sum(S2, 1);

	// Fusion
	for (int ii = 0; ii < nIter; ii++) {
		mat S1_next = S1 * S2_o * S1.t();
		mat S2_next = S2 * S1_o * S2.t();
		// update S1_o and S2_o values
		S1_o = S1_next + eye<mat>(nr, nr);
		S2_o = S2_next + eye<mat>(nr, nr);
		// make them symmetric
		S1_o = (S1_o + S1_o.t()) / 2;
		S2_o = (S2_o + S2_o.t()) / 2;
	}
  // Fused matrix
	mat K = (S1_o + S2_o) / NM;
	K.each_col() /= sum(K, 1);
	K = (K + K.t() + eye<mat>(nr, nr)) / 2;
	return K;
}
