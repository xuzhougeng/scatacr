#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){

  int k = z.size() ;

  IntegerMatrix  mat(nrows, ncols);

  for (int i = 0; i < k; i++){
    mat(rp[i],cp[i]) = z[i];
  }

  return mat;
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
mat <- Matrix::Matrix(1:9, nrow=3, sparse = TRUE)
row_pos <- mat@i
col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                nrows =  mat@Dim[1], ncols = mat@Dim[2])
*/



