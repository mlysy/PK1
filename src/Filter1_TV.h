#ifndef Filter1_TV_h
#define Filter1_TV_h 1

// time-varying 1st-order linear filter.
// filter equation is:
// y[t+1] = rho[t] * y[t] + x[t]
void filter1_TV(double* y, double y0, double* x, double* rho, int N);

#endif
