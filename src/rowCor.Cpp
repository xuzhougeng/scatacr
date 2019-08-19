#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
// [[Rcpp::export]]
NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, NumericMatrix X, NumericMatrix Y) {

  if(X.ncol() != Y.ncol()){
    stop("Columns of Matrix X and Y must be equal length!");
  }
  if(max(idxX) > X.nrow()){
    stop("Idx X greater than nrow of Matrix X");
  }
  if(max(idxY) > Y.nrow()){
    stop("Idx Y greater than nrow of Matrix Y");
  }
  // Transpose Matrices
  X = transpose(X);
  Y = transpose(Y);

  const int nx = X.ncol();
  const int ny = Y.ncol();
  // Centering the matrices
  for (int j = 0; j < nx; ++j) {
    X(_, j) = X(_, j) - mean(X(_, j));
  }
  for (int j = 0; j < ny; ++j) {
    Y(_, j) = Y(_, j) - mean(Y(_, j));
  }
  // Compute 1 over the sample standard deviation
  NumericVector inv_sqrt_ss_X(nx);
  for (int i = 0; i < nx; ++i) {
    inv_sqrt_ss_X(i) = 1/sqrt(sum( X(_, i) * X(_, i) ));
  }
  NumericVector inv_sqrt_ss_Y(ny);
  for (int i = 0; i < ny; ++i) {
    inv_sqrt_ss_Y(i) = 1/sqrt(sum( Y(_, i) * Y(_, i) ));
  }
  //Calculate Correlations
  const int n = idxX.size();
  NumericVector cor(n);
  for(int k = 0; k < n; k++){
    cor[k] = sum( X(_, idxX[k] - 1) * Y(_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);
  }
  return(cor);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
