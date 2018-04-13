# check stan vs R posteriors

## require(PK1)
require(rstan)
source("PK1-testfunctions.R")

context("stanmodels")

#--- generate data -------------------------------------------------------------

# observation times
nSub <- 12
nObs <- 7
tObs <- t(replicate(nSub, c(0, sort(runif(nObs-2)*12), 12)))
# doses chosen arbitrarily between 3 and 6
Dose <- runif(n = nSub, min = 3, max = 6)
# Define parameters for (lognormal) priors
muCl <- -3.22; sigCl <- 0.1
muKa <- 0.40; sigKa <- 0.1
muKe <- -2.52; sigKe <- 0.1
muSigmaP <- sqrt(0.2); sigSigmaP <- 0.01
muSigmaM <- sqrt(0.1); sigSigmaM <- 0.01
# generate parameters
Cl <- rlnorm(n = nSub, meanlog = muCl, sdlog = sigCl)
Ka <- rlnorm(n = nSub, meanlog = muKa, sdlog = sigKa)
Ke <- rlnorm(n = nSub, meanlog = muKe, sdlog = sigKe)
sigmaP <- muSigmaP
sigmaM <- muSigmaM
# generate SDE
Xt <- t(sapply(1:nSub, function(ii) {
  pk1.sim(tObs = tObs[ii,], X0 = 0, sigmaP = sigmaP,
          Dose = Dose[ii], Ke = Ke[ii], Ka = Ka[ii], Cl = Cl[ii])
}))
# generate data
yObs <- matrix(rnorm(nObs*nSub, mean = Xt, sd = sigmaM), nSub, nObs)

#--- generate test parameters --------------------------------------------------

# probably easiest to do this by sampling from richest model

# data + initialization
pk1.data <- list(nObs = nObs, nSub = nSub, Yt = yObs, t = tObs,
                 D = Dose, sdDef = 5)

# fitting
nsamples <- 100
pk1.init.ls <- list(pk1.init(effect = "Mixed", DE = "SDE",
                             meas = "Noise",
                             Cl = Cl, Ka = Ka, Ke = Ke,
                             sigmaP = sigmaP, sigmaM = sigmaM,
                             muCl = muCl, muKa = muKa, muKe = muKe,
                             sigCl = sigCl, sigKa = sigKa, sigKe = sigKe,
                             Xt = Xt))
suppressWarnings({
  PK1.Full.Fit <-  sampling(object = pk1.mixed.sde.noise,
                            data = pk1.data, chains = 1,
                            init = pk1.init.ls, iter = nsamples)
})

#--- compare stan and R likelihoods --------------------------------------------

ntest <- nrow(as.data.frame(PK1.Full.Fit)) # number of posterior samples

# identifiers for all models
PK1.Models <- expand.grid(meas = c("Pure", "Noise"),
                          DE = c("ODE", "SDE"),
                          effect = c("Fixed", "Mixed"))
PK1.Models <- PK1.Models[!with(PK1.Models, {
  iout <- (DE == "ODE") & (meas == "Pure")
  iout
}),]
nMod <- nrow(PK1.Models)
rownames(PK1.Models) <- 1:nMod
PK1.Models <- as.matrix(PK1.Models[,3:1])

# run through each specific stan model (imod)
for(imod in 1:nMod) {
  model.name <- paste(c("pk1", tolower(PK1.Models[imod,])),
                      collapse = ".")
  test_that(paste0("logpost_R == logpost_Stan for model: ", model.name), {
    pk1.model <- get(x = model.name)
    # fit the stan model to get log_post
    pk1.init.ls <- list(pk1.init(descr = PK1.Models[imod,],
                                 Cl = Cl, Ka = Ka, Ke = Ke,
                                 sigmaP = sigmaP, sigmaM = sigmaM,
                                 muCl = muCl, muKa = muKa, muKe = muKe,
                                 sigCl = sigCl, sigKa = sigKa, sigKe = sigKe,
                                 Xt = Xt))
    suppressWarnings({
      pk1.fit <- sampling(object = pk1.model, data = pk1.data,
                          init = pk1.init.ls,
                          iter = 1, warmup = 0, chains = 1,
                          algorithm = "Fixed_param")
    })
    # stan logposterior (regular scale)
    lp.stan <- sapply(1:ntest, function(ii) {
      pars <- extract.iter(PK1.Full.Fit, ii)
      # keep only relevant pars
      pars <- pk1.format.pars(pars = pars, descr = PK1.Models[imod,])
      upars <- unconstrain_pars(object = pk1.fit, pars = pars)
      log_prob(object = pk1.fit, upars = upars, adjust_transform = TRUE)
    })
    # R logposterior
    lp.R <- sapply(1:ntest, function(ii) {
      pars <- extract.iter(PK1.Full.Fit, ii)
      # keep only relevant pars
      pars <- pk1.format.pars(pars = pars, descr = PK1.Models[imod,])
      do.call(pk1.logpost, args = c(pars,
                                    list(yObs = yObs, tObs = tObs, Dose = Dose,
                                         descr = PK1.Models[imod,])))
    })
    # maximum difference (excluding additive constant)
    expect_equal(max(abs(diff(lp.stan-lp.R))), 0)
  })
}


## max.diff <- sapply(1:nrow(PK1.Models), function(imod) {
##   pk1.model <- paste(c("pk1", tolower(PK1.Models[imod,])),
##                      collapse = ".")
##   pk1.model <- get(x = pk1.model)
##   # fit the stan model to get log_post
##   pk1.init.ls <- list(pk1.init(descr = PK1.Models[imod,],
##                                Cl = Cl, Ka = Ka, Ke = Ke,
##                                sigmaP = sigmaP, sigmaM = sigmaM,
##                                muCl = muCl, muKa = muKa, muKe = muKe,
##                                sigCl = sigCl, sigKa = sigKa, sigKe = sigKe,
##                                Xt = Xt))
##   pk1.fit <- sampling(object = pk1.model, data = pk1.data,
##                       init = pk1.init.ls,
##                       iter = 1, warmup = 0, chains = 1,
##                       algorithm = "Fixed_param")
##   # stan logposterior (regular scale)
##   lp.stan <- sapply(1:nsamples, function(ii) {
##     pars <- extract.iter(PK1.Full.Fit, ii)
##     # keep only relevant pars
##     pars <- pk1.format.pars(pars = pars, descr = PK1.Models[imod,])
##     upars <- unconstrain_pars(object = pk1.fit, pars = pars)
##     log_prob(object = pk1.fit, upars = upars, adjust_transform = TRUE)
##   })
##   # R logposterior
##   lp.R <- sapply(1:nsamples, function(ii) {
##     pars <- extract.iter(PK1.Full.Fit, ii)
##     # keep only relevant pars
##     pars <- pk1.format.pars(pars = pars, descr = PK1.Models[imod,])
##     do.call(pk1.logpost, args = c(pars,
##                            list(yObs = yObs, tObs = tObs, Dose = Dose,
##                                 descr = PK1.Models[imod,])))
##   })
##   # maximum difference (excluding additive constant)
##   max(abs(diff(lp.stan-lp.R)))
## })

## max.diff
