#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerMatrix tabulate2dCpp(IntegerVector x1,
                            int xmin,
                            int xmax,
                            IntegerVector y1,
                            int ymin,
                            int ymax){
  if(x1.size() != y1.size()){
    stop("width must equal size!");
  }
  IntegerVector x = clone(x1);
  IntegerVector y = clone(y1);
  int n = x.size();
  IntegerVector rx = seq(xmin,xmax);
  IntegerVector ry = seq(ymin,ymax);
  IntegerMatrix mat( ry.size() , rx.size() );
  int xi,yi;
  for(int i = 0; i < n; i++){
    xi = (x[i] - xmin);
    yi = (y[i] - ymin);
    if(yi >= 0 && yi < ry.size()){
      if(xi >= 0 && xi < rx.size()){
        mat( yi , xi ) = mat( yi , xi ) + 1;
      }
    }
  }
  return mat;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

  */
