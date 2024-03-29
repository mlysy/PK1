#' Simulation of the PK1 SDE model with unequal sampling times.
#'
#' @template param-tObs
#' @template param-X0
#' @template param-Dose
#' @template param-Ke
#' @template param-Ka
#' @template param-Cl
#' @template param-sigmaP
#' @return Concentration matrix of size \code{nReps x nObs}, where \code{nObs = length(tObs)} and \code{nReps} is determined by the length of the other inputs.
#' @details Each of \code{X0}, \code{Dose}, \code{Ke}, \code{Ka}, \code{Cl}, and \code{sigmaP} can be scalars or vectors.  In case of the latter, a concentration curve is returned for each value.
#' @export
pk1.sim <- function(tObs, X0, Dose, Ke, Ka, Cl, sigmaP) {
  nObs <- length(tObs)
  # determine nReps
  ln <- c(length(X0), length(Dose), length(Ke),
          length(Ka), length(Cl), length(sigmaP))
  nReps <- max(ln)
  if(any(ln != 1 && ln != nReps))
    stop("All arguments except tObs must have common length or length 1.")
  X0 <- rep(X0, len = nReps)
  Dose <- rep(Dose, len = nReps)
  Cl <- rep(Cl, len = nReps)
  Ke <- rep(Ke, len = nReps)
  Ka <- rep(Ka, len = nReps)
  sigmaP <- rep(sigmaP, len = nReps)
  Xt <- .PK1_Sim(nReps = nReps, X0 = X0, Dose = Dose, tObs = tObs,
                 Cl = Cl, Ka = Ka, Ke = Ke, sigmaP = sigmaP)
  t(Xt)
}

#' Solution of the PK1 ODE model.
#'
#' @inheritParams pk1.sim
#' @details The initial concentration (at time \code{t = 0}) is assumed to be \code{X0 = 0}.
#' @return Solution vector at times \code{tObs}.
#' @export
pk1.ode <- function(tObs, Dose, Cl, Ka, Ke) {
  R <- Dose * Ka*Ke / Cl / (Ke - Ka)
  R * (exp(-Ka * tObs) - exp(-Ke * tObs))
}

#' Simulation of maximal concentration for the PK1 model.
#'
#' @details Generates a concentration curve at specified timepoints, and records the maximal concentration value along with the time when it happens.
#' @inheritParams pk1.sim
#' @return A matrix with columns \code{Cmax} and \code{tCmax}.
#' @export
pk1.Cmax.sim <- function(tObs, X0, Dose, Ke, Ka, Cl, sigmaP) {
  nObs <- length(tObs)
  # determine nReps
  ln <- c(length(X0), length(Dose), length(Ke),
          length(Ka), length(Cl), length(sigmaP))
  nReps <- max(ln)
  if(any(ln != 1 && ln != nReps))
    stop("All arguments except tObs must have common length or length 1.")
  X0 <- rep(X0, len = nReps)
  Dose <- rep(Dose, len = nReps)
  Cl <- rep(Cl, len = nReps)
  Ke <- rep(Ke, len = nReps)
  Ka <- rep(Ka, len = nReps)
  sigmaP <- rep(sigmaP, len = nReps)
  Cmax <- .PK1_Cmax_Sim(nReps = nReps, X0 = X0, Dose = Dose, tObs = tObs,
                        Cl = Cl, Ka = Ka, Ke = Ke, sigmaP = sigmaP)
  colnames(Cmax) <- c("Cmax", "tCmax")
  Cmax
}

#' Loglikelihood of the PK1 Model
#'
#' Calculates the loglikelihood for different formulations of the fixed-effects PK1 model, i.e., as an SDE or ODE, and with or without measurement error.
#'
#' @return If \code{yObs} is not missing, loglikelihood.  Otherwise, list with elements \code{mu} and \code{V} of \code{yObs}.
#' @param yObs Vector of measurements.
#' @param tObs Vector of observation times.
#' @param sigmaM Standard deviation of measurment error.
#' @inheritParams pk1.sim
#' @export
pk1.loglik <- function(Cl, Ka, Ke, sigmaP, sigmaM, yObs, tObs, Dose) {
  if(missing(sigmaP)) sigmaP <- 0
  if(missing(sigmaM)) sigmaM <- 0
  if((sigmaP == 0) && (sigmaM == 0)) {
    stop("sigmaP and sigmaM can't both be zero.")
  }
  nObs <- length(tObs)
  dt <- diff(tObs)
  R <- Ke*Ka*Dose/Cl/(Ke-Ka)
  rho <- exp(-Ke * dt)
  lambda <- R*exp(-Ka*tObs[-nObs])*(exp(-Ka * dt) - rho)
  if(sigmaP > 0) {
    # SDE
    sigma <- sigmaP*sqrt((1-rho^2)/(2*Ke))
    MV <- gAR1.MV(aa = rho, bb = lambda, cc = sigma)
    # add initial value
    mu <- c(0, MV$b)
    tau <- sigmaP^2/(2*Ke)
    B <- tau * MV$A
    V <- cbind(B, MV$V + B %o% MV$A)
    V <- rbind(c(tau, B), V)
    if(sigmaM > 0) {
      diag(V) <- diag(V) + sigmaM^2
    }
  } else {
    # ODE
    MV <- gAR1.MV(aa = rho, bb = lambda, cc = 0, x0 = 0)
    mu <- c(0, MV$mu)
    V <- diag(rep(sigmaM^2, len = nObs))
  }
  if(!missing(yObs)) {
    out <- .lmvn(x = yObs- mu, Sigma = V)
  } else {
    out <- list(mu = mu, V = V)
  }
  out
}

# log-pdf for a single observation x ~ N(0, Sigma)
.lmvn <- function(x, Sigma) {
  if(length(x) == 0) return(0)
  l2pi <- 1.837877066409345483560659472811 # log(2*pi)
  iV <- .solveV(Sigma, x, ldV = TRUE)
  -.5 * (crossprod(x, iV$y)[1] + iV$ldV + length(x) * l2pi)
}

# solve method for variance matrices
# optionally computes log determinant as well
.solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V)
  if(missing(x)) x <- diag(nrow(V))
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}


#' Initialize the stan sampler
#'
#' @export
pk1.init <- function(descr = NULL, effect = c("Fixed", "Mixed", "Indep"),
                     DE = c("SDE", "ODE"), meas = c("Pure", "Noise"),
                     Cl, Ka, Ke, sigmaP, sigmaM, muCl, muKa, muKe,
                     sigCl, sigKa, sigKe, Xt) {
  if(!is.null(descr)) {
    effect <- descr["effect"]
    DE <- descr["DE"]
    meas <- descr["meas"]
  }
  effect <- match.arg(effect)
  DE <- match.arg(DE)
  meas <- match.arg(meas)
  init <- NULL
  if(meas == "Noise") init <- list(Xt = Xt, sigmaM = sigmaM)
  theta <- list(Cl = Cl, Ka = Ka, Ke = Ke)
  if(effect == "Fixed") {
    #init <- c(lapply(init, function(x) rowMeans(as.matrix(x))), lapply(theta, mean))
    init <- c(init, lapply(theta, mean))
  } else {
    init <- c(init, theta)
    if(effect == "Mixed") {
      init <- c(init,
                muCl = muCl, muKa = muKa, muKe = muKe,
                sigCl = sigCl, sigKa = sigKa, sigKe = sigKe)
    }
  }
  if(DE == "SDE") {
    init <- c(init, sigmaP = sigmaP)
  }
  init
}

#--- depreciated ----------------------------------------------------------------

## if(FALSE) {
## pk1.sim <- function(tObs, X0, Dose, Ke, Ka, Cl, sigmaP, z) {
##   nObs <- length(tObs)
##   dt <- diff(tObs)
##   R <- Ka*Dose/Cl
##   rho <- exp(-Ke * dt)
##   lambda <- Ke*R*exp(-Ka*tObs[-nObs])*(exp(-Ka * dt) - rho)/(Ke-Ka)
##   sigma <- sigmaP*sqrt((1-rho^2)/(2*Ke))
##   # for-loop in C++
##   #Xt <- rep(NA, nObs)
##   #Xt[1] <- X0
##   #eps <- rnorm(nObs-1)
##   #for(jj in 1:(nObs-1)) {
##   #  Xt[jj+1] <- rho[jj] * Xt[jj] + lambda[jj] + sigma[jj] * eps[jj]
##   #}
##   #Xt
##   if(missing(z)) z <- rnorm(nObs-1)
##   filter1.tv(rho = rho, x = lambda + sigma*z, y0 = X0)
## }
## if(!exists("pk1.compile") || pk1.compile) {
##   message("--- COMPILING STAN CODE ---")
##   pk1.mixed.sde.noise <- stan_model(file = file.path("inst", "stan",
##                                       "PK1_Mixed_SDE_Noise.stan"))
##   pk1.mixed.sde.pure <- stan_model(file = file.path("inst", "stan",
##                                      "PK1_Mixed_SDE_Pure.stan"))
##   pk1.mixed.ode.noise <- stan_model(file = file.path("inst", "stan",
##                                       "PK1_Mixed_ODE_Noise.stan"))
##   pk1.fixed.ode.noise <- stan_model(file = file.path("inst", "stan",
##                                       "PK1_Fixed_ODE_Noise.stan"))
##   pk1.fixed.sde.pure <- stan_model(file = file.path("inst", "stan",
##                                      "PK1_Fixed_SDE_Pure.stan"))
##   pk1.fixed.sde.noise <- stan_model(file = file.path("inst", "stan",
##                                       "PK1_Fixed_SDE_Noise.stan"))
##   message("---------- DONE -----------")
## }
## }
