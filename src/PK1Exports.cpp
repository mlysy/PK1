#include <Rcpp.h>
#include "filterTV1.h"
using namespace Rcpp;

// [[Rcpp::export("filter1.tv")]]
NumericVector FilterTV1(NumericVector rho, NumericVector x,
			double y0 = 0) {
  int N = x.length()+1;
  NumericVector y(N);
  filterTV1(&(y[0]), y0, &(x[0]), &(rho[0]), N);
  return y;
}

//[[Rcpp::export(".PK1_Sim")]]
NumericMatrix PK1_Sim(int nReps,
		      NumericVector X0, NumericVector Dose,
		      NumericVector tObs,
		      NumericVector Cl, NumericVector Ka, NumericVector Ke,
		      NumericVector sigmaP) {
  int nObs = tObs.length();
  NumericMatrix Xt(nObs, nReps);
  NumericVector dt(nObs-1);
  NumericVector tObs2(nObs-1);
  double R;
  NumericVector rho(nObs-1);
  NumericVector lambda(nObs-1);
  NumericVector eps(nObs-1);
  int ii, jj;

  for(ii=0; ii<nObs-1; ii++) {
    dt[ii] = tObs[ii+1]-tObs[ii];
    tObs2[ii] = tObs[ii];
  }

  for(ii=0; ii<nReps; ii++) {
    R = Ka[ii] * Ke[ii] * Dose[ii] / (Cl[ii] * (Ke[ii]-Ka[ii]));
    rho = exp(-Ke[ii] * dt);
    lambda = R * exp(-Ka[ii]*tObs2) * (exp(-Ka[ii] * dt) - rho);
    eps = rnorm(nObs-1);
    eps = sigmaP[ii] * sqrt((1-rho*rho)/(2*Ke[ii])) * eps;
    lambda = lambda + eps;
    filterTV1(&(Xt(0,ii)), X0[ii], &(lambda[0]), &(rho[0]), nObs);
  }
  return Xt;
}

//[[Rcpp::export(".PK1_Cmax_Sim")]]
NumericMatrix PK1_Cmax_Sim(int nReps,
		      NumericVector X0, NumericVector Dose,
		      NumericVector tObs,
		      NumericVector Cl, NumericVector Ka, NumericVector Ke,
		      NumericVector sigmaP) {
  int nObs = tObs.length();
  NumericVector Xt(nObs);
  NumericVector dt(nObs-1);
  NumericVector tObs2(nObs-1);
  double R;
  NumericVector rho(nObs-1);
  NumericVector lambda(nObs-1);
  NumericVector eps(nObs-1);
  NumericMatrix Cmax(nReps,2);
  int tCmax;
  int ii;

  for(ii=0; ii<nObs-1; ii++) {
    dt[ii] = tObs[ii+1]-tObs[ii];
    tObs2[ii] = tObs[ii];
  }

  for(ii=0; ii<nReps; ii++) {
    R = Ka[ii] * Ke[ii] * Dose[ii] / (Cl[ii] * (Ke[ii]-Ka[ii]));
    rho = exp(-Ke[ii] * dt);
    lambda = R * exp(-Ka[ii]*tObs2) * (exp(-Ka[ii] * dt) - rho);
    eps = rnorm(nObs-1);
    eps = sigmaP[ii] * sqrt((1-rho*rho)/(2*Ke[ii])) * eps;
    lambda = lambda + eps;
    filterTV1(&(Xt[0]), X0[ii], &(lambda[0]), &(rho[0]), nObs);
    tCmax = which_max(Xt);
    Cmax(ii,0) = Xt[tCmax];
    Cmax(ii,1) = tObs[tCmax];
  }
  return Cmax;
}
