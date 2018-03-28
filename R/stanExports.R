#--- compile the stan models ----------------------------------------------------

# * always goes to inst/stan
# * creates an RData file containing all stan functions.
# * then steps through each function.  compile stan code unless stanmodel already exists and has same c++ code.
# * TODO: setup .Rbuildignore

if(file.exists(file.path("inst", "stan", "PK1-stan.RData"))) {
  load(file.path("inst", "stan", "PK1-stan.RData"))
}

stan_model_full <- function(stan.name, stan.mod = NULL) {
  fl <- file.path("inst", "stan", stan.name)
  st <- stanc_builder(file = fl)
  if(!is.null(stan.mod) &&
     (stan.mod@model_cpp$model_cppcode == st$cppcode)) {
    md <- stan.mod
  } else {
    md <- stan_model(stanc_ret = st)
  }
  md
}

# compile models
stan.names <- c("PK1_Mixed_SDE_Noise", "PK1_Mixed_SDE_Pure",
                "PK1_Mixed_ODE_Noise", "PK1_Fixed_ODE_Noise",
                "PK1_Fixed_SDE_Pure", "PK1_Fixed_SDE_Noise")
R.names <- c("pk1.mixed.sde.noise", "pk1.mixed.sde.pure",
             "pk1.mixed.ode.noise", "pk1.fixed.ode.noise",
             "pk1.fixed.sde.pure", "pk1.fixed.sde.noise")
message("--- COMPILING STAN CODE ---")
for(ii in 1:length(stan.names)) {
  message(stan.names[ii], ".stan (", ii, "/", length(stan.names), ")")
  assign(R.names[ii],
         value = stan_model_full(paste0(stan.names[ii], ".stan"),
           get0(R.names[ii])))
}
message("---------- DONE -----------")

# save only if already exists, or explicitly instructed to do so.
if(file.exists(file.path("inst", "stan", "PK1-stan.RData")) ||
   (exists("save.stan") && save.stan)) {
  save(list = R.names, file = file.path("inst", "stan", "PK1-stan.RData"))
}
rm(stan.names, R.names, ii, stan_model_full) # cleanup
