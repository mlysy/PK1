#--- check mean and variance ----------------------------------------------------

require(PK1)
require(mniw)

N <- 5
aa <- rnorm(N)
bb <- rnorm(N)
cc <- abs(rnorm(N))
x0 <- rnorm(1)

# simulate and calculate log-likelihood
nreps <- 10
X <- matrix(NA, nreps, N)
ll1 <- 0
X0 <- x0
for(ii in 1:N) {
  X1 <- aa[ii] * X0 + bb[ii] + cc[ii] * rnorm(nreps)
  ll1 <- ll1 + dnorm(X1, aa[ii] * X0 + bb[ii], cc[ii], log = TRUE)
  X[,ii] <- X1
  X0 <- X1
}

# now using MV
MV <- gAR1.MV(aa = aa, bb = bb, cc = cc, x0 = x0)

ll2 <- dmNorm(x = X, mu = MV$mu, V = MV$V, log = TRUE)

ll1-ll2
