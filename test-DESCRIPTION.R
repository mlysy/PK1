# ok fresh rstan_package_skeleton

require(rstantools)
rstan_package_skeleton(name = "PK1stan",
                       code_files = "/Users/mlysy/Documents/R/test/rstantools/R/utils.R",
                       stan_files = "PK1_Fixed_ODE_Noise.stan")

source("/Users/mlysy/Documents/R/test/rstantools/R/utils.R")


descr <- read.dcf("DESCRIPTION")

format_pkgfield <- function(pkg_names) {
  pkgs <- strsplit(pkg_names, "[ ]*,[ ]*")[[1]]
  ver <- strsplit(pkgs, "[ ]*")[[1]]
}


pkg_names <- paste0(descr[1,"LinkingTo"], "   , Matrix  , Vector")
pkg_vers <- ver_split(pkg_names)

ver_comb(pkg_vers)


tmp <- "asd    (fjdakl a.  da dsa"
strsplit(tmp, "[ ]+")[[1]]

tmp <- "fea fanda,, afd,a dacdkal;me( CldaCVD    a"
strsplit(tmp, "[ ]*,[ ]*")[[1]][1]
grepl

format_pkgfield()

pkg_names
new_names <- "Matrix > 1.3 a, rstan (>= 2), msde"

append_pkgfield(pkg_names, new_names)

# ok use_rstan does the following things:

# 1. update description file
# 2. create stan cc files (i.e., make_cc)
# 3. check/add to namespace for things Rcpp.package.skeleton would add
# 4. create Makevars file...or at least do something with the stan_files code
# 5. tell you about everything it couldn't do
