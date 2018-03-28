#' Fixed/Mixed Effect ODE/SDE Implementations of the PK1 Model
#'
#' @details The Stochastic Differential Equation (SDE) representing the One-Compartment Pharmacokinetic (PK1) model is
#'
#' \deqn{d Xt = - Ke (Xt - Dose/Cl * Ka exp(-Ka t)) dt + gamma dBt,}
#'
#' where \code{Bt} is Brownian motion.  When \code{gamma = 0} the SDE above reduces to an Ordinary Differential Equation (ODE).
#' @docType package
#' @name PK1
#' @importFrom Rcpp evalCpp
#' @useDynLib PK1
NULL
