#include <Rcpp.h>
#include "Filter1_TV.h"

using namespace Rcpp;

// [[Rcpp::interfaces(r,cpp)]]

// time-varying 1st-order linear filter.
void filter1_TV(double* y, double y0, double* x, double* rho, int N) {
  y[0] = y0;
  for(int ii=0; ii<N-1; ii++) {
    y[ii+1] = rho[ii]*y[ii] + x[ii];
  }
  return;
}

// [[Rcpp::export("filter1.tv")]]
NumericVector Filter1_TV(NumericVector rho, NumericVector x,
			 double y0 = 0) {
  int N = x.length()+1;
  NumericVector y(N);
  filter1_TV(&(y[0]), y0, &(x[0]), &(rho[0]), N);
  return y;
}
