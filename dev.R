
library("devtools")
load_all(recompile = TRUE)

nml_scr <- run_nml(fun= SCR, parl=10, ks=rep(3,5), Ns=c(rep(100,4),200), 
                   fits = 2, batchsize=5000, burn=10000, precision=1)

library("examplemodels")
nml_scr <- run_nml(fun= SCR2, parl=10, ks=rep(3,5), Ns=c(rep(100,4),200), 
                   packages_multicore = "examplemodels",
                   fits = 2, batchsize=5000, burn=10000, precision=1)
