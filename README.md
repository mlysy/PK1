## PK1: Mixed Effect SDE Inference for a One-Compartment Pharmacodynamic Model

*Martin Lysy, Kamal Rai*

<!-- *April 21, 2018* -->

---

The primary purpose of this package is to test the [**rstantools**](https://github.com/stan-dev/rstantools) build mechanism for Stan-enabled R packages.  To install and run unit tests:

```r
if (!require("devtools")) {
  # requires the devtools package
  install.packages("devtools")
}

# install the package itself
devtools::install_github("mlysy/PK1")

# run unit tests
if (!require("testthat")) {
  # requires the testthat package
  install.packages("testthat")
}

testthat::test_package("PK1", reporter = "progress")

```
