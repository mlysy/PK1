////////////////////////////////////////////////////////////
//
// Inference for the 1-compartment PK model
// mlysy 2018
//
// Model is:
//
// Ke, Ka, Cl: one per patient
// sigmaM: one per study
// X_i(t): concentration at time t for patient i
// X_ij: concentration at time t_j for patient i
// Y_ij: measurement at time t_j for patient i
//
// (Ke_i, Ka_i, Cl_i) ~iid prior
// dXt/dt = -Ke(Xt - D*Ka/Cl * exp(-Ka*t))
// Y_ij ~iid N(X_ij, sigmaM^2)
//
////////////////////////////////////////////////////////////

functions {
#include /include/PK1_DE_Functions.stan
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
      dt[ii,jj] = t[ii,jj+1]-t[ii,jj];
    }
  }
}

parameters {
  // parameters of interest
  real<lower=0> Cl[nSub]; // clearance rates (associated with the old mu)
  real<lower=0> Ka[nSub]; // absorption rates (associated with the old mu)
  real<lower=0> Ke[nSub]; // elimination rates (the old gamma)

  // nuisance parameters
  real<lower=0> sigmaM; // measurement error

  // hyperparameters
  // (means and standard deviations of params on log scale)
  real muCl;
  real muKa;
  real muKe;
  real<lower=0> sigCl;
  real<lower=0> sigKa;
  real<lower=0> sigKe;

}

transformed parameters {
  real Xt[nSub, nObs]; // concentrations
  for(ii in 1:nSub) {
    Xt[ii] = to_array_1d(PK1_ODE(D[ii], Cl[ii], Ka[ii], Ke[ii], t[ii]));
  }
}

// inference
model {
  // Hyper-priors
  muCl ~ normal(0, sdDef);
  muKa ~ normal(0, sdDef);
  muKe ~ normal(0, sdDef);
  sigCl ~ lognormal(0, sdDef);
  sigKa ~ lognormal(0, sdDef);
  sigKe ~ lognormal(0, sdDef);

  // Priors
  Cl ~ lognormal(muCl, sigCl);
  Ka ~ lognormal(muKa, sigKa);
  Ke ~ lognormal(muKe, sigKe);
  sigmaM ~ lognormal(0, sdDef);

  // Likelihood
  // measurements
  to_array_1d(Yt) ~ normal(to_array_1d(Xt), sigmaM);
}
