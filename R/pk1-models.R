#' @name pk1-models
#' @title Variations on the One-Compartment Pharmacokinetic (PK1) Model.
#' @description Stan implementations of the PK1 model, with ODE or SDE formulation, fixed or mixed effects, and with or without additional measurement error.
#' @details Each of these objects are of class \code{stanmodel}, such that parameter inference via MCMC can be done with the \code{rstan} package, e.g., with the call
#'
#' \code{sampling(object = stanmodel, iter = 1e4, chains = 1)}

#' @rdname pk1-models
#' @usage sampling(object = pk1.mixed.sde.noise, ...)
#' @docType data
#' @export
"pk1.mixed.sde.noise"

#' @rdname pk1-models
#' @usage sampling(object = pk1.mixed.sde.pure, ...)
#' @docType data
#' @export
"pk1.mixed.sde.pure"

#' @rdname pk1-models
#' @usage sampling(object = pk1.mixed.ode.noise, ...)
#' @docType data
#' @export
"pk1.mixed.ode.noise"

#' @rdname pk1-models
#' @usage sampling(object = pk1.fixed.ode.noise, ...)
#' @docType data
#' @export
"pk1.fixed.ode.noise"

#' @rdname pk1-models
#' @usage sampling(object = pk1.fixed.sde.pure, ...)
#' @docType data
#' @export
"pk1.fixed.sde.pure"

#' @rdname pk1-models
#' @usage sampling(object = pk1.fixed.sde.noise, ...)
#' @docType data
#' @export
"pk1.fixed.sde.noise"
