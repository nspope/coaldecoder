#include <RcppArmadillo.h>

#ifndef _EXPM_H
#define _EXPM_H

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp1z)]]

typedef enum { Ward_2, Ward_1, Ward_buggy_octave } precond_type;
inline void (*expmatrix)(double *x, int n, double *z, precond_type precond_kind); //req c++17

inline arma::mat expm(arma::mat x) {
    arma::mat z(x.n_rows, x.n_cols);
    (*expmatrix)(x.begin(), x.n_rows, z.begin(), Ward_2);
    return z;
}

#endif
