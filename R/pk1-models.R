#' @name pk1-models
#' @title Variations on the One-Compartment Pharmacokinetic (PK1) Model.
#' @description Stan implementations of the PK1 model, with ODE or SDE formulation, fixed or mixed effects, and with or without additional measurement error.
#' @details Each of these objects are of class \code{stanmodel}, such that parameter inference via MCMC can be done with the \code{rstan} package, e.g., with the call
#' \preformatted{
#' sampling(object = stanmodel, iter = 1e4, chains = 1)
#' }
#' @include stanmodels.R

#' @rdname pk1-models
#' @usage sampling(object = pk1.mixed.sde.noise, ...)
#' @export
pk1.mixed.sde.noise <- stanmodels$PK1_Mixed_SDE_Noise

#' @rdname pk1-models
#' @usage sampling(object = pk1.mixed.sde.pure, ...)
#' @export
pk1.mixed.sde.pure <- stanmodels$PK1_Mixed_SDE_Pure

#' @rdname pk1-models
#' @usage sampling(object = pk1.mixed.ode.noise, ...)
#' @export
pk1.mixed.ode.noise <- stanmodels$PK1_Mixed_ODE_Noise

#' @rdname pk1-models
#' @usage sampling(object = pk1.fixed.ode.noise, ...)
#' @export
pk1.fixed.ode.noise <- stanmodels$PK1_Fixed_ODE_Noise

#' @rdname pk1-models
#' @usage sampling(object = pk1.fixed.sde.pure, ...)
#' @export
pk1.fixed.sde.pure <- stanmodels$PK1_Fixed_SDE_Pure

#' @rdname pk1-models
#' @usage sampling(object = pk1.fixed.sde.noise, ...)
#' @export
pk1.fixed.sde.noise <- stanmodels$PK1_Fixed_SDE_Noise
