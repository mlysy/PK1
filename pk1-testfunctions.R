# extract a single iteration from stanfit object.
# returns a list which is formatted for log_prob and log_prob_grad
extract.iter <- function(stanfit, ii) {
  lapply(extract(stanfit), function(arr) {
    lst <- apply(arr, 1, list)
    lst[[ii]][[1]]
  })
}

# log-density of PK1 sde
# no density for first observation
pk1.sde.ldens <- function(Xt, tObs, Dose, Ke, Ka, Cl, sigmaP) {
  nObs <- length(tObs)
  dt <- diff(tObs)
  R <- Ke*Ka*Dose/Cl/(Ke-Ka)
  rho <- exp(-Ke * dt)
  lambda <- R*exp(-Ka*tObs[-nObs])*(exp(-Ka * dt) - rho)
  sigma <- sigmaP*sqrt((1-rho^2)/(2*Ke))
  sum(dnorm(x = Xt[2:nObs], mean = rho * Xt[2:nObs-1] + lambda,
            sd = sigma, log = TRUE))
}

# flimsy log-posterior on the Stan scale
# just for debugging purposes
pk1.logpost <- function(Cl, Ka, Ke, sigmaP, sigmaM,
                        muCl, muKa, muKe,
                        sdCl, sdKa, sdKe,
                        Xt, tObs, yObs, Dose, sdDef = 5, descr,
                        effect = c("Fixed", "Mixed", "Indep"),
                        DE = c("ODE", "SDE"), meas = c("Pure", "Noise"),
                        ..., debug = FALSE) {
  nSub <- nrow(yObs)
  nObs <- ncol(yObs)
  if(!missing(descr)) {
    effect <- descr["effect"]
    DE <- descr["DE"]
    meas <- descr["meas"]
  }
  effect <- match.arg(effect)
  DE <- match.arg(DE)
  meas <- match.arg(meas)
  lp <- 0
  if(debug) browser()
  if(effect == "Mixed") {
    # hyperparameters
    lp <- lp + dnorm(muCl, 0, sdDef, log = TRUE)
    lp <- lp + dnorm(muKa, 0, sdDef, log = TRUE)
    lp <- lp + dnorm(muKe, 0, sdDef, log = TRUE)
    lp <- lp + dlnorm(sdCl, 0, sdDef, log = TRUE) + log(sdCl)
    lp <- lp + dlnorm(sdKa, 0, sdDef, log = TRUE) + log(sdKa)
    lp <- lp + dlnorm(sdKe, 0, sdDef, log = TRUE) + log(sdKe)
    # parameters
    lp <- lp + sum(dlnorm(Cl, muCl, sdCl, log = TRUE)) + sum(log(Cl))
    lp <- lp + sum(dlnorm(Ka, muKa, sdKa, log = TRUE)) + sum(log(Ka))
    lp <- lp + sum(dlnorm(Ke, muKe, sdKe, log = TRUE)) + sum(log(Ke))
  } else if(effect == "Fixed") {
    # parameters
    lp <- lp + dlnorm(Cl[1], 0, sdDef, log = TRUE) + log(Cl[1])
    lp <- lp + dlnorm(Ka[1], 0, sdDef, log = TRUE) + log(Ka[1])
    lp <- lp + dlnorm(Ke[1], 0, sdDef, log = TRUE) + log(Ke[1])
    Cl <- rep(Cl, len = nSub)
    Ka <- rep(Ka, len = nSub)
    Ke <- rep(Ke, len = nSub)
  }
  if(DE == "SDE") {
    if(meas == "Pure") Xt <- yObs
    lp <- lp + dlnorm(sigmaP, 0, sdDef, log = TRUE) + log(sigmaP)
    for(ii in 1:nSub) {
      lp <- lp + pk1.sde.ldens(Xt = Xt[ii,], tObs = tObs[ii,],
                               Dose = Dose[ii],
                               Ke = Ke[ii], Ka = Ka[ii], Cl = Cl[ii],
                               sigmaP = sigmaP)
    }
  } else if(DE == "ODE") {
    Xt <- matrix(0, nSub, nObs)
    for(ii in 1:nSub) {
     Xt[ii,] <- pk1.sim(tObs = tObs[ii,], X0 = 0, Dose = Dose[ii],
                        Ke = Ke[ii], Ka = Ka[ii], Cl = Cl[ii], sigmaP = 0)
    }
  }
  if(meas == "Noise") {
    lp <- lp + dlnorm(sigmaM, 0, sdDef, log = TRUE) + log(sigmaM)
    lp <- lp + sum(dnorm(yObs, mean = Xt, sd = sigmaM, log = TRUE))
  }
  lp
}

pk1.format.pars <- function(pars, descr,
                            effect = c("Fixed", "Mixed", "Indep"),
                            DE = c("ODE", "SDE"),
                            meas = c("Pure", "Noise")) {
  if(!missing(descr)) {
    effect <- descr["effect"]
    DE <- descr["DE"]
    meas <- descr["meas"]
  }
  effect <- match.arg(effect)
  DE <- match.arg(DE)
  meas <- match.arg(meas)
  par.names <- c("Cl", "Ka", "Ke")
  if(effect == "Mixed") {
    par.names <- c(par.names,
                   "muCl", "sdCl", "muKa", "sdKa", "muKe", "sdKe")
  }
  if(DE == "SDE") {
    par.names <- c(par.names, "sigmaP")
  }
  if(meas == "Noise") {
    par.names <- c(par.names, "sigmaM")
  }
  if(DE == "SDE" && meas == "Noise") {
    par.names <- c(par.names, "Xt")
  }
  pars <- pars[par.names]
  if(effect == "Fixed") {
    pars$Cl <- pars$Cl[1]
    pars$Ka <- pars$Ka[1]
    pars$Ke <- pars$Ke[1]
  }
  pars
}
