library(coaldecoder)
library(Matrix)
library(expm)

set.seed(1)
P <- 5
A <- matrix(0, 3^P, 3^P)
for(i in 1:nrow(A)){
  nz <- sample((1:nrow(A))[-i], 3*(P-1))
  A[i,nz] <- rexp(length(nz))
}
diag(A) <- -rowSums(A)
B <- matrix(rexp(nrow(A)*3), nrow(A), 3)
g <- matrix(rnorm(nrow(A)*3), nrow(A), 3)
spA <- as(A, "dgCMatrix")

one <- test_SparseMatrixExponentialMultiply(t(spA), B, 0.51, g)
two <- test_matrix_exponential_multiply(t(A), B, 0.51, g)
three <- expm::expm(0.51 * t(A)) %*% B
#TODO: use numDeriv::jacobian to test gradient

all(abs(one$X - two$X) < sqrt(.Machine$double.eps))
all(abs(one$dA[as.vector(one$dA[]!=0)] - two$dA[as.vector(one$dA[]!=0)]) < sqrt(.Machine$double.eps))
all(abs(one$dB - two$dB) < sqrt(.Machine$double.eps))
all(abs(one$X - three) < sqrt(.Machine$double.eps))

#microbenchmark::microbenchmark(one <- test_SparseMatrixExponentialMultiply(t(spA), B, 0.51, g))
#microbenchmark::microbenchmark(two <- test_matrix_exponential_multiply(t(A), B, 0.51, g))
