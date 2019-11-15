#ifndef filterTV1_h
#define filterTV1_h 1

// time-varying 1st-order linear filter.
// filter equation is:
// y[t+1] = rho[t] * y[t] + x[t]
inline void filterTV1(double* y, double y0, double* x, double* rho, int N) {
  y[0] = y0;
  for(int ii=0; ii<N-1; ii++) {
    y[ii+1] = rho[ii]*y[ii] + x[ii];
  }
  return;
}

#endif
