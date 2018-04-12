// ODE solution
vector PK1_ODE(real D, real Cl, real Ka, real Ke, real[] t) {
  real R;
  R = Ke*Ka*D/Cl/(Ke-Ka);
  return R * (exp(-Ka * to_vector(t)) - exp(-Ke * to_vector(t)));
}

// SDE loglikelihood
// ignores first observation for now
real PK1_SDE_lpdf(real[] Xt, real Ka, real Ke, real Cl, real sigmaP,
		 real Dose, real[] t, real[] dt) {
  real R;
  real ll;
  real rho;
  real lambda;
  real tau;
  int nObs;
  nObs = num_elements(t);
  ll = 0.0;
  R = Ke*Ka*Dose/Cl/(Ke-Ka);
  for(ii in 1:(nObs-1)) {
    rho = exp(-Ke * dt[ii]);
    lambda = R * exp(-Ka * t[ii]) * (exp(-Ka * dt[ii]) - rho);
    tau = sigmaP * sqrt((1-rho^2)/(2.0*Ke));
    ll = ll + normal_lpdf(Xt[ii+1] | rho * Xt[ii] + lambda, tau);
  }
  return ll;
}
