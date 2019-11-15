#' Mean and Variance of a Generalized AR(1) Process
#'
#' @param aa Scalar or vector mean-reversion parameter (see Details).
#' @param bb Scalar or vector offset parameter (see Details).
#' @param cc Scalar or vector noise standard deviation parameter (see Details).
#' @param x0 Scalar initial value.
#' @return A list with elements \code{mu} and \code{V} if \code{x0} is given.  Otherwise, a list with elements \code{A}, \code{b} and \code{V}, where \code{mu = A x0 + b}.  The sizes of \code{mu} and \code{V} are \code{N = max(length(aa), length(bb), length(cc))}.
#' @details The GAR(1) process is given by
#' \preformatted{
#' x[n] = aa[n] * x[n-1] + bb[n] + cc[n] * e[n],   e[n] ~iid N(0,1).
#' }
#' @export
gAR1.MV <- function(aa, bb, cc, x0) {
  # determine length
  N <- c(length(aa), length(bb), length(cc))
  if(length(unique(c(1,N))) > 2) {
    stop("aa, bb, cc have incompatible lengths.")
  }
  N <- max(N)
  aa <- rep(aa, len = N)
  bb <- rep(bb, len = N)
  cc <- rep(cc, len = N)
  # cholesky
  L <- diag(N)
  for(ii in 1:(N-1)) {
    L[(ii+1):N,ii] <- cumprod(aa[(ii+1):N])
  }
  # variance
  V <- tcrossprod(sweep(L, 2, cc, FUN = "*"))
  # mean
  alpha <- aa[1] * L[,1]
  beta <- c(L %*% bb)
  if(!missing(x0)) {
    mu <- alpha * x0 + beta
    out <- list(mu = mu, V = V)
  } else {
    out <- list(A = alpha, b = beta, V = V)
  }
}
