# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

test_matrix_exponential_multiply <- function(A, B, t, g) {
    .Call(`_coaldecoder_test_matrix_exponential_multiply`, A, B, t, g)
}

test_coalescent_epoch_gaussian <- function(states, migr_mat, gradient) {
    .Call(`_coaldecoder_test_coalescent_epoch_gaussian`, states, migr_mat, gradient)
}

test_SparseMatrixExponentialMultiply <- function(A, B, t, g) {
    .Call(`_coaldecoder_test_SparseMatrixExponentialMultiply`, A, B, t, g)
}

test_TrioTransitionRates <- function(M, G) {
    .Call(`_coaldecoder_test_TrioTransitionRates`, M, G)
}

test_TrioAdmixtureProportions <- function(A, G) {
    .Call(`_coaldecoder_test_TrioAdmixtureProportions`, A, G)
}

test_CoalescentEpoch <- function(states, M, A, gradient) {
    .Call(`_coaldecoder_test_CoalescentEpoch`, states, M, A, gradient)
}

