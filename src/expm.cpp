// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// -------- pilfered from https://github.com/eddelbuettel/rcppkalman/blob/master/src/expmMat.cpp ----------

//#include <RcppArmadillo.h>

//copy the following to RcppExports, in place of existing R_init_*, also include expm.h in RcppExports

//arma::mat expm(arma::mat x) {
//    arma::mat z(x.n_rows, x.n_cols);
//    (*expmatrix)(x.begin(), x.n_rows, z.begin(), Ward_2);
//    return z;
//}

//extern "C" void R_init_coaldecoder(DllInfo *dll) { 
//    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
//    R_useDynamicSymbols(dll, FALSE);
//    expmat = (void (*) (double*, int, double*, precond_type)) R_GetCCallable("expm", "expm"); 
//} 
