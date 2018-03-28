#' Mean and Variance of a Generalized AR(1) Process
#'
#' @return A list with elements \code{mu} and \code{V} if \code{x0} is given.  Otherwise, a list with elements \code{A}, \code{b} and \code{V}, where \code{mu = A x0 + b}.
#' @details The GAR(1) process is given by
#' \deqn{x[n] = a[n] * x[n-1] + b[n] + c[n] * e[n], \qquad e[n] ~iid N(0,1)}.
#' @export
gAR1.MV <- function(aa, bb, cc, x0, debug = FALSE) {
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
  if(debug) browser()
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
