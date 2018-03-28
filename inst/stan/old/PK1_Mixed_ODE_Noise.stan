////////////////////////////////////////////////////////////
//
// Inference for the 1-compartment PK model
// krai, jdubin, mlysy 2015
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
  real[,] muODE(real[] D, real[] Cl, real[] Ka, real[] Ke,
		real[,] t, real[,] dt, int nObs, int nSub) {
    real R;
    real rho;
    real lambda;
    real Xt[nSub, nObs];

    // concentrations
    for(ii in 1:nSub) {
      Xt[ii,1] <- 0.0;
      R <- Ke[ii]*Ka[ii]*D[ii]/Cl[ii]/(Ke[ii]-Ka[ii]);
      for(jj in 1:(nObs-1)) {
	rho <- exp(-Ke[ii] * dt[ii,jj]);
	lambda <- R * exp(-Ka[ii] * t[ii,jj]) * (exp(-Ka[ii] * dt[ii,jj]) - rho);
	Xt[ii,jj+1] <- rho * Xt[ii,jj] + lambda;
      }
    }
    return Xt;
  }
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
  real<lower=0> sdCl;
  real<lower=0> sdKa;
  real<lower=0> sdKe;

}

transformed parameters {
  real Xt[nSub, nObs]; // concentrations
  Xt <- muODE(D, Cl, Ka, Ke, t, dt, nObs, nSub);
}

// inference
model {
  // Simplifying expressions
  //real R;
  //real rho;
  //real lambda;
  //real Xt[nSub, nObs]; // concentrations

  // Hyper-priors
  muCl ~ normal(0, sdDef);
  muKa ~ normal(0, sdDef);
  muKe ~ normal(0, sdDef);
  sdCl ~ lognormal(0, sdDef);
  sdKa ~ lognormal(0, sdDef);
  sdKe ~ lognormal(0, sdDef);

  // Priors
  Cl ~ lognormal(muCl, sdCl);
  Ka ~ lognormal(muKa, sdKa);
  Ke ~ lognormal(muKe, sdKe);
  sigmaM ~ lognormal(0, sdDef);

  // Likelihood
  // concentrations
  //for(ii in 1:nSub) {
  //  Xt[ii,1] <- 0.0;
  //  R <- Ke[ii]*Ka[ii]*D[ii]/Cl[ii]/(Ke[ii]-Ka[ii]);
  //  for(jj in 1:(nObs-1)) {
  //    rho <- exp(-Ke[ii] * dt[ii,jj]);
  //    lambda <- R * exp(-Ka[ii] * t[ii,jj]) * (exp(-Ka[ii] * dt[ii,jj]) - rho);
  //    Xt[ii,jj+1] <- rho * Xt[ii,jj] + lambda;
  //  }
  //}
  // measurements
  for(ii in 1:nSub) {
    for(jj in 1:nObs) {
      Yt[ii,jj] ~ normal(Xt[ii,jj], sigmaM);
    }
  }
}
