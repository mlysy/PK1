#--- package setup --------------------------------------------------------------

require(devtools)
require(Rcpp)
require(rstan)
pkg.name <- "PK1"

# package update
# save.stan <- TRUE # uncomment this the first time
compileAttributes()
document()
devtools::install()

#--- create package folder (one-time only) --------------------------------------

#pkg.path <- getwd()
#remove.packages(pkg.name)
#unlink(file.path(pkg.path, pkg.name), recursive = TRUE)
#Rcpp.package.skeleton(name = pkg.name, example_code = FALSE)
## add inst folder for stan code
#dir.create(path = file.path(pkg.path, pkg.name, "inst", "stan"),
#           recursive = TRUE)
## add rstan dependency
#descr <- read.dcf(file = file.path(pkg.path, pkg.name, "DESCRIPTION"))
#if("Depends" %in% colnames(descr)) {
#  descr[,"Depends"] <- paste0(descr[,"Depends"], ", rstan")
#} else {
#  descr <- cbind(descr, Depends = "rstan")
#}
#write.dcf(descr, file = file.path(pkg.path, pkg.name, "DESCRIPTION"))
#compileAttributes(pkgdir = file.path(pkg.path, pkg.name))
