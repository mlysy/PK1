////////////////////////////////////////////////////////////
//
// Inference for the 1-compartment PK model
// krai, jdubin, mlysy 2015
//
// Model is:
//
// Ke, Ka, Cl: one per patient
// sigmaP: one per study
// X_i(t): concentration at time t for patient i
// X_ij: concentration at time t_j for patient i
// Y_ij: measurement at time t_j for patient i
//
// (Ke, Ka, Cl) ~ prior
// dXt = -Ke(Xt - D*Ka/Cl * exp(-Ka*t)) dt + sigmaP dBt
// Y_ij = X_ij
//
////////////////////////////////////////////////////////////

data {
  int<lower=1> nObs; // number of observations
  int<lower=1> nSub; // number of subjects

  real Yt[nSub, nObs]; // noisy observations
  real t[nSub, nObs]; // observation times
  real<lower=0> D[nSub]; // doses

  real<lower=0> sdDef; // standard deviation for default prior
}

transformed data {
  real dt[nSub, nObs-1]; // time intervals
  real Xt[nSub, nObs]; // concentrations

  for (jj in 1:(nObs-1)) {
    for(ii in 1:nSub) {
      dt[ii,jj] <- t[ii,jj+1]-t[ii,jj];
      Xt[ii,jj] <- Yt[ii,jj];
    }
  }
  for(ii in 1:nSub) {
    Xt[ii,nObs] <- Yt[ii,nObs];
  }
}

parameters {
  // parameters of interest
  real<lower=0> Cl; // clearance rates (associated with the old mu)
  real<lower=0> Ka; // absorption rates (associated with the old mu)
  real<lower=0> Ke; // elimination rates (the old gamma)

  // nuisance parameters
  real<lower=0> sigmaP; // SDE volatility

}

// inference
model {
  // Simplifying expressions
  real R;
  real rho;
  real lambda;
  real tau;

  // Priors
  Cl ~ lognormal(0, sdDef);
  Ka ~ lognormal(0, sdDef);
  Ke ~ lognormal(0, sdDef);
  sigmaP ~ lognormal(0, sdDef);

  // Likelihood
  // concentrations
  for(ii in 1:nSub) {
    R <- Ke*Ka*D[ii]/Cl/(Ke-Ka);
    for(jj in 1:(nObs-1)) {
      rho <- exp(-Ke * dt[ii,jj]);
      lambda <- R * exp(-Ka * t[ii,jj]) * (exp(-Ka * dt[ii,jj]) - rho);
      tau <- sigmaP*sqrt((1-rho^2)/(2.0*Ke));
      Xt[ii,jj+1] ~ normal(rho * Xt[ii,jj] + lambda, tau);
    }
  }
}
