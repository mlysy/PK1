////////////////////////////////////////////////////////////
//
// Inference for the 1-compartment PK model
// krai, jdubin, mlysy 2015
//
// Model is:
//
// Ke, Ka, Cl: one per study
// sigmaM: one per study
// X_i(t): concentration at time t for patient i
// X_ij: concentration at time t_j for patient i
// Y_ij: measurement at time t_j for patient i
//
// (Ke, Ka, Cl) ~ prior
// dXt/dt = -Ke(Xt - D*Ka/Cl * exp(-Ka*t))
// Y_ij ~iid N(X_ij, sigmaM^2)
//
////////////////////////////////////////////////////////////

functions {
  #include "PK1_DE_Functions.stan"
}

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

  for (jj in 1:(nObs-1)) {
    for(ii in 1:nSub) {
      dt[ii,jj] <- t[ii,jj+1]-t[ii,jj];
    }
  }
}

parameters {
  // parameters of interest
  real<lower=0> Cl; // clearance rates (associated with the old mu)
  real<lower=0> Ka; // absorption rates (associated with the old mu)
  real<lower=0> Ke; // elimination rates (the old gamma)

  // nuisance parameters
  real<lower=0> sigmaM; // measurement error

}

transformed parameters {
  real Xt[nSub, nObs]; // concentrations
  for(ii in 1:nSub) {
    Xt[ii] <- to_array_1d(PK1_ODE(D[ii], Cl, Ka, Ke, t[ii]));
  }
}

// inference
model {
  // Priors
  Cl ~ lognormal(0, sdDef);
  Ka ~ lognormal(0, sdDef);
  Ke ~ lognormal(0, sdDef);
  sigmaM ~ lognormal(0, sdDef);

  // Likelihood
  // measurements
  to_array_1d(Yt) ~ normal(to_array_1d(Xt), sigmaM);
}
