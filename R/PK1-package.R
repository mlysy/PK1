#' Fixed/Mixed Effect ODE/SDE Implementations of the PK1 Model
#'
#' @details The Stochastic Differential Equation (SDE) representing the One-Compartment Pharmacokinetic (PK1) model is
#' \preformatted{
#' d Xt = - Ke (Xt - Dose/Cl * Ka exp(-Ka t)) dt + sigmaP dBt,
#' }
#' where \code{Bt} is Brownian motion.  When \code{sigmaP = 0} the SDE above reduces to an Ordinary Differential Equation (ODE).
#' @docType package
#' @name PK1
#' @import Rcpp
#' @import methods
#' @import rstan
#' @useDynLib PK1, .registration = TRUE
NULL
