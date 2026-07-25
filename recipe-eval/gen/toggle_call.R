
dyn.load("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/libtoggle.so")
res <- .Fortran("toggle_mfpt",
   seed   = as.integer(11),
   ttot   = as.double(100000.0),
   dt     = as.double(0.01),
   trelax = as.double(100.0),
   g0     = as.double(10.0),
   g1     = as.double(40.0),
   Kthr   = as.double(100.0),
   kdeg   = as.double(0.1),
   nexp   = as.integer(4),
   b      = as.double(20.0),
   tau    = as.double(0),
   rate   = as.double(0),
   ncross = as.integer(0))
cat(sprintf("R_NCROSS %d\n", res$ncross))
cat(sprintf("R_TAU %.15e\n", res$tau))
cat(sprintf("R_RATE %.15e\n", res$rate))
